/*
   asc2eph - program to convert JPL ASCII ephemerides to binary format.

         Reads ASCII ephemeris header file and ASCII ephemeris data file,
         and writes binary ephemeris file.

         Format:

            asc2eph header-input-file data-input-file ephemeris-output-file
*/

#include <stdio.h>
#include <stdlib.h>    //exit()
#include <string.h>    //strcmp()
#include "ephcom.h"


int ephcom_readascii_header(FILE *infp, struct ephcom_Header *header);
int ephcom_readascii_block(FILE *infp, struct ephcom_Header *header, double *datablock);
int ephcom_writeascii_header(FILE *outfp, struct ephcom_Header *header);
int ephcom_writebinary_header(FILE *outfp, struct ephcom_Header *header);
int ephcom_writebinary_block(FILE *outfp, struct ephcom_Header *header, 
                             int blocknum, double *datablock);


int main(int argc, char *argv[]){

    struct ephcom_Header header1;
    double *datablock;  /* Will hold coefficients from a data block */
    static int nblocks=0; /* Read 0 data blocks so far - first time through */
    double startjd, stopjd; /* First and last desired JD */
    double laststart, laststop;
/*
   Output file parameters.
*/
    FILE *infp, *outfp;

    if (argc < 4) {
        fprintf(stderr,
           "\nFormat:\n\n         %s ascii-header ascii-data binary-output [startJD [stopJD]]\n\n",
           argv[0]);
        exit(1);
    }

    if ((infp = fopen(argv[1],"r")) == NULL) {
        fprintf(stderr,"\nERROR: Can't open %s for input.\n\n", argv[1]);
        exit(1);
    }

    if ((outfp = fopen(argv[3],"r")) == NULL) {    //先用只读方式打开，判断文件是否存在
        if ((outfp = fopen(argv[3],"wb")) == NULL) {
            fprintf(stderr,"\nERROR: Can't open ASCII header file %s for output.\n\n", argv[3]);
            exit(1);
        }   
    }
    else {
        fprintf(stderr,"\nERROR: Output ephemeris file %s already exists.\n\n", argv[3]);
        exit(1);
    }

    if (argc > 4) {        
        startjd = atof(argv[4]);
        if (argc > 5)
            stopjd = atof(argv[5]);
        else 
            stopjd = EPHCOM_MAXJD;
    }
    else {
        startjd = EPHCOM_MINJD;
        stopjd = EPHCOM_MAXJD;
    }

/*
   Make sure ephemeris is within the desired range.
*/
    ephcom_readascii_header(infp, &header1);
    if (header1.ss[0] > stopjd) {
        fprintf(stderr,"\nERROR: ephemeris begins after desired end JD.\n\n");
        exit(1);
    }
    if (header1.ss[1] < startjd) {
        fprintf(stderr, "\nERROR: ephemeris ends before desired start JD.\n\n");
        exit(1);
    }

/*
   Done reading ASCII header.  Now we'll read and write data blocks.
*/
    fclose(infp);

    if ((infp = fopen(argv[2],"r")) == NULL) {
        fprintf(stderr,"\nERROR: Can't open ASCII data file %s for input.\n\n", argv[2]);
        exit(1);
    }

    datablock = (double *)malloc(header1.ncoeff * sizeof(double));
    while (ephcom_readascii_block(infp, &header1, datablock) > 0 &&
           datablock[0] <= stopjd) {
        if (nblocks == 0) { /* First time through */
            if (datablock[1] >= startjd) {
                if (datablock[1] - datablock[0] != header1.ss[2]) { /* Days / block */
                    fprintf(stderr,
                            "ERROR: Header block span %g doesn't match block 1 (%g to %g).\n",
                            header1.ss[2], datablock[0], datablock[1]);
                    fprintf(stderr,
                            "       Wrong header file in use.  Can't continue.\n");
                    exit(1);
                }
                /* At least write one block */
                ephcom_writebinary_block(outfp, &header1, nblocks, datablock);
                header1.ss[0] = datablock[0]; /* Start Julian Day */
                header1.ss[1] = datablock[1]; /* Stop Julian Day */
                nblocks++;
            }
        }
        else { /* Not the first data block */
      /*
         Blocks overlap in some ASCII files at their common start/end;
         ignore duplicated block but don't report - this is normal.
      */
            if (laststart == datablock[0] && laststop == datablock[1]) {
            }
            else if (laststop == datablock[0] &&
                    (datablock[1] - datablock[0]) == header1.ss[2]) {
                if (datablock[1] >= startjd) {
                    ephcom_writebinary_block(outfp, &header1, nblocks, datablock);
                    header1.ss[1] = datablock[1]; /* New last day */
                    nblocks++;
                }
            }
            else if (laststop != datablock[0]) {
                fprintf(stderr,
                        "ERROR: Blocks %d (%g to %g) and %d (%g to %g) not adjacent.\n",
                        nblocks, laststart, laststop,
                        nblocks+1, datablock[0], datablock[1]);
                fprintf(stderr,"       Can't continue.\n");
                exit(1);
            }
      /* If blocks are contiguous and spans match header, continue */
            else if ((laststop - laststart) == header1.ss[2] &&
                     (datablock[1] - datablock[0]) != header1.ss[2]) {
                fprintf(stderr,
                        "ERROR: Blocks %d (%g to %g) and %d (%g to %g) have differnt spans.\n",
                        nblocks, laststart, laststop,
                        nblocks+1, datablock[0], datablock[1]);
                 fprintf(stderr,"       Can't continue.\n");
                 exit(1);
            }
            else {
                fprintf(stderr,
                        "ERROR: Unexpected condition in block %d (%g to %g).\n",
                        nblocks+1, datablock[0], datablock[1]);
                fprintf(stderr,
                        "       Previous block was (%g to %g).\n", laststart, laststop);
                fprintf(stderr,
                        "       Can't continue.\n");
                exit(1);
            }
        }
        laststart = datablock[0];
        laststop = datablock[1];
    }
/*
   Now write the header information, with updated start/stop dates.
*/
    ephcom_writeascii_header(stdout, &header1);
    ephcom_writebinary_header(outfp, &header1);

    fclose(outfp);
    fclose(infp);

    printf("Wrote 2 header blocks + %d data blocks, %d coefficients per data block.\n\n",
           nblocks, header1.ncoeff);
    printf("%s should be exactly (2+%d)*(%d)*8 = %d bytes.\n\n",
           argv[3], nblocks, header1.ncoeff, (2+nblocks)*header1.ncoeff*8);

    return 0;
}