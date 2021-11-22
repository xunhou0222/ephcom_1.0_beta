#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "ephcom.h"
/*
   Read a JPL Ephemeris ASCII header from the file pointed to by infp
   and store values in header structure.  Write any errors to stderr.
*/
int ephcom_readascii_header(FILE *infp, struct ephcom_Header *header) {

    char group[13]; /* To store the "GROUP" header line information */
    int i, j;
    int iword; /* word number we're reading in a line */
    int blockout; /* number of bytes we've written to current block/rec in file */
    char readbuf[EPHCOM_MAXLINE + 1]; //To store the characters in each line

//Declaration of used functions
    void ephcom_nxtgrp(char *, char *, FILE *);
    char *fgets(char *, int, FILE *);
    size_t fwrite(const void *ptr, size_t size, size_t  nmemb, FILE *stream);

/*
   First header line: KSIZE= # NCOEFF= #
*/
    rewind(infp);
    fgets(readbuf, EPHCOM_MAXLINE, infp);
    sscanf(readbuf, "%*6s%6d%*11s%6d", &header->ksize, &header->ncoeff);
    if (header->ksize != 2*header->ncoeff) {
        fprintf(stderr, "Badly formed header; KSIZE <> 2*NCOEFF\n\n");
        exit(1);
    }
/*
   GROUP 1010: Title of ephemeris (DE/LE number, start JD, end JD)
*/
    ephcom_nxtgrp(group, "GROUP   1010", infp);
    fgets(header->ttl[0], EPHCOM_MAXLINE, infp);  /* JPL Ephemeris title line */
    if (strncmp(header->ttl[0], "JPL ", 4) != 0) {
        fprintf(stderr,"\nERROR: file is not a JPL ASCII header file.\n\n");
        exit(1);
    }
    fgets(header->ttl[1], EPHCOM_MAXLINE, infp);  /* Start epoch */
    fgets(header->ttl[2], EPHCOM_MAXLINE, infp);  /* Finish epoch */
/*
   Convert any newlines or tabs to single spaces.
*/
    for (i=0; i<3; i++) {
        for (j=0; j<84; j++)
            if (isspace(header->ttl[i][j]))  header->ttl[i][j] = ' ';    //ctype.h
        header->ttl[i][84] = '\0';
    }
/*
   GROUP 1030: Start and End JD, timestep (in JD) per block.
*/
    ephcom_nxtgrp(group, "GROUP   1030", infp);
    fgets(readbuf, EPHCOM_MAXLINE, infp);
    sscanf(readbuf, " %lE %lE %lE", &header->ss[0], &header->ss[1], &header->ss[2]);
/*
   GROUP 1040: Constant names.
*/
    ephcom_nxtgrp(group, "GROUP   1040", infp);
    fgets(readbuf, EPHCOM_MAXLINE, infp);
    header->ncon = atoi(readbuf);//convert str to int
/*
   Now read the constant names, 10 per line, each 6 characters long
   preceded by 2 blanks.  Pad names with blanks to make 6 characters.
*/
    for (i=0; i<header->ncon;) {
        fgets(readbuf, EPHCOM_MAXLINE, infp);
        for (iword=0; iword<10 && i<header->ncon; iword++, i++) {
            strncpy(header->cnam[i], &readbuf[2 + iword*8], 6);
            header->cnam[i][6] = '\0';
        }
    }
/*
   GROUP 1041: Constant values.
*/
    ephcom_nxtgrp(group, "GROUP   1041", infp);
    fgets(readbuf, EPHCOM_MAXLINE, infp);
    header->nval = atoi(readbuf);
    if (header->nval != header->ncon) {
        fprintf(stderr,"Error: number of constants and values not equal.\n\n");
        exit(1);
    }
/*
   Now read constant values, 3 per line, 26 characters each.
*/
    for (i = 0; i < header->ncon; i += 3) {
        fgets(readbuf, EPHCOM_MAXLINE, infp);
        for (j=0; j<strlen(readbuf); j++)
            if (tolower(readbuf[j]) == 'd') 
                readbuf[j] = 'E'; /* exponent is 'E' */
        if (i+2 < header->ncon)
            sscanf(readbuf, "%lE %lE %lE",
                   &header->cval[i], &header->cval[i+1], &header->cval[i+2]);
        else if (i+1 < header->ncon)
            sscanf(readbuf, "%lE %lE", &header->cval[i], &header->cval[i+1]);
    }
/*
   GROUP 1050: Constant values.
*/
    ephcom_nxtgrp(group, "GROUP   1050", infp);
    for (i =0; i < 3; i++) {
        fgets(readbuf, EPHCOM_MAXLINE, infp); /* Read line of 13 6-digit integers */
        for (j = 0; j < 12; j++) 
            header->ipt[j][i] = atoi(&readbuf[6*j]);
        header->lpt[i] = atoi(&readbuf[6*12]);
    }
/*
   If there are no coefficients for an ipt[i][] object (i.e., ipt[i][1]==0),
   then ipt[i][0] should contain the value of the next available coefficient
   number rather than 0, as per communication of Myles Standish to Paul Hardy
   on preferred format of ephemeris headers.

   If there are no libration coefficients (i.e., lpt[1]==0), then lpt[0]
   should contain the value of the next available coefficient number rather
   than 0 as well, as per the same communication from Myles Standish.
*/
/* First set j to maximum index into ipt[] that has coefficients */
    j = 0;
    for (i=0; i<12; i++)
        if (header->ipt[i][1] > 0 && header->ipt[i][0] > j)
            j = i;
/* Now set j to next available index count. */
    if (header->lpt[1] > 0 && header->lpt[0] > j)
        j = header->lpt[1] + header->lpt[1] * header->lpt[2] * 3;
    else
        j = header->ipt[j][0] +
            header->ipt[j][1] * header->ipt[j][2] * (j==11 ? 2 : 3);  //if j == 11, *2; otherwise, *3.
    for (i=1; i<12; i++)
        if (header->ipt[i][0] == 0) header->ipt[i][0] = j;
    if (header->lpt[0] == 0) header->lpt[0] = j;
/*
   Set the maximum number of Chebyshev coefficients possible for this file,
   to initialize position and velocity Chebyshev coefficient arrays during
   Chebyshev interpolation.
*/
    header->maxcheby = 0;
    for (i=0; i<12; i++)
        if (header->ipt[i][1] > header->maxcheby)
            header->maxcheby = header->ipt[i][1];
    if (header->lpt[1] > header->maxcheby)
        header->maxcheby = header->lpt[1];

    header->au = 0.0;
    header->emrat = 0.0;
    header->numde = 0;
    header->clight = 0.0;
    header->numle = 0;
    for (i = 0; i < header->ncon; i++) {
        if (strncmp(header->cnam[i], "AU    ", 6) == 0)
            header->au = header->cval[i];
        else if (strncmp(header->cnam[i], "EMRAT ", 6) == 0)
            header->emrat = header->cval[i];
        else if (strncmp(header->cnam[i], "DENUM ", 6) == 0)
            header->numde = header->cval[i];
        else if (strncmp(header->cnam[i], "CLIGHT", 6) == 0)
            header->clight = header->cval[i];
        else if (strncmp(header->cnam[i], "LENUM ", 6) == 0)
            header->numle = header->cval[i];
    }
    if (header->numle == 0) header->numle = header->numde;
/*
   GROUP 1070: Constant values.
*/
    ephcom_nxtgrp(group, "GROUP   1070", infp);
/*
   Now we're pointing to the first block of coefficient data, after header.
   Return at the point where we can start reading coefficients.
*/
    return(0);
}




/*
   Read a block of data coefficients from a JPL ASCII ephemeris file.
   Returns number of coefficients read, 0 at EOF.
*/
int ephcom_readascii_block(FILE *infp, struct ephcom_Header *header, double *datablock) {

    int i,j;
    int datapoints; /* points of data we've read/converted/written */
    char readbuf[EPHCOM_MAXLINE + 1];
    double val1, val2, val3; /* To read text line with 3 double precision words */

/*
   First line in an ASCII block will be the block number, followed by
   the number of coefficients.
*/
    datapoints = 0;
    if (fgets(readbuf, EPHCOM_MAXLINE, infp) && !feof(infp)) {
        sscanf(readbuf, "%d %d", &i, &j);
        if (j != header->ncoeff) {
            fprintf(stderr,
                    "\nERROR: ASCII data file's %d coefficients/block\n", j);
            fprintf(stderr,
                    "       doesn't match ASCII header's %d coefficients/block.\n\n",
                    header->ncoeff);
            exit(1);
        }
        for (i=0; i < header->ncoeff && !feof(infp); i += 3) {
            fgets(readbuf, EPHCOM_MAXLINE, infp);
            for (j=0; j<strlen(readbuf); j++)
                if (tolower(readbuf[j]) == 'd')  
                    readbuf[j] = 'E';
            sscanf(readbuf, " %lE %lE %lE", &val1, &val2, &val3);
            datablock[i] = val1;
            datapoints++;
            if ((i+1) < header->ncoeff) {
                datablock[i+1] = val2;
                datapoints++;
                if ((i+2) < header->ncoeff) {
                    datablock[i+2] = val3;
                    datapoints++;
                }
            }
        }
    }

    return(datapoints);
}




/*
   Read a JPL Ephemeris header in binary format.  Store values in
   an ephcom_Header struct.
*/
int ephcom_readbinary_header(FILE *infp, struct ephcom_Header *header) {

    int i, j;

// Declaration of used functions
    char *fgets(char *, int, FILE *);
    int fseek(FILE *, long, int);
    int fgetc(FILE *);
    double ephcom_indouble(FILE *);
    int ephcom_inint(FILE *);

    rewind(infp);
/*
   Read title lines.
*/
    for (i=0; i<3; i++) {
        for (j=0; j<84; j++) 
            header->ttl[i][j] = fgetc(infp);
        if (i == 0 && strncmp(header->ttl[0], "JPL ", 4) != 0) {
            fprintf(stderr,"\nERROR: file is not a JPL ephemeris file.\n\n");
            if (strncmp(header->ttl[0], "KSIZE", 5) == 0)
                fprintf(stderr,"File is an ASCII JPL ephemeris header instead.\n\n");
            exit(1);
        }
        header->ttl[i][j] = '\0';
    }
/*
   Read constant names.
*/
    for (i=0; i<400; i++) {
        for (j=0; j<6; j++) 
            header->cnam[i][j] = fgetc(infp);
        header->cnam[i][j] = '\0';
    }
/*
   Read ephemeris start epoch, stop epoch, and step size (in Julian Days).
*/
    for (i=0; i<3; i++)
        header->ss[i] = ephcom_indouble(infp);
/*
   Read NCON, AU, EMRAT.
*/
    header->ncon  = ephcom_inint(infp);
    header->au    = ephcom_indouble(infp);
    header->emrat = ephcom_indouble(infp);
    header->nval  = header->ncon;
/*
   Read indexes for coefficients in data block.  Written in transposed
   order (Fortran and C matrices are transposed).
*/
    for (i=0; i<12; i++) {
        for (j=0; j<3; j++)
            header->ipt[i][j] = ephcom_inint(infp);
    }
    header->numde = ephcom_inint(infp);  /* Get ephemeris number */
    for (i=0; i<3; i++) 
        header->lpt[i] = ephcom_inint(infp);
/*
   If there are no coefficients for an ipt[i][] object (i.e., ipt[i][1]==0),
   then ipt[i][0] should contain the value of the next available coefficient
   number rather than 0, as per communication of Myles Standish to Paul Hardy
   on preferred format of ephemeris headers.

   If there are no libration coefficients (i.e., lpt[1]==0), then lpt[0]
   should contain the value of the next available coefficient number rather
   than 0 as well, as per the same communication from Myles Standish.
*/
/* First set j to maximum index into ipt[] that has coefficients */
    j = 0;
    for (i=1; i<12; i++)
        if (header->ipt[i][1] > 0 && header->ipt[i][0] > j)
    j = i;
/* Now set j to next available index count. */
    if (header->lpt[1] > 0 && header->lpt[0] > j)
        j = header->lpt[1] + header->lpt[1] * header->lpt[2] * 3;
    else
        j = header->ipt[j][0] + header->ipt[j][1] * header->ipt[j][2] * (j==11 ? 2 : 3);
    for (i=1; i<12; i++)
        if (header->ipt[i][0] == 0) 
            header->ipt[i][0] = j;
    if (header->lpt[0] == 0) 
        header->lpt[0] = j;
/*
   Set the maximum number of Chebyshev coefficients possible for this file,
   to initialize position and velocity Chebyshev coefficient arrays during
   Chebyshev interpolation.
*/
    header->maxcheby = 0;
    for (i=0; i<12; i++)
        if (header->ipt[i][1] > header->maxcheby)
            header->maxcheby = header->ipt[i][1];
    if (header->lpt[1] > header->maxcheby)
        header->maxcheby = header->lpt[1];

/*
   Calculate number of coefficients, starting with
   highest index into a data block (stored in j).
*/
    j = 0;
    for (i=1; i<12; i++)
        if (header->ipt[i][1] > 0 && header->ipt[i][0] > header->ipt[j][0]) 
            j = i;

/*
   Now see if the starting point we found is lower than where
   lpt[] starts.  If not, use lpt[] for largest value.
*/
    if (header->lpt[1] > 0 && header->lpt[0] > header->ipt[j][0]) {
        header->ncoeff = header->lpt[0] - 1 +     /* starting point */
                         (header->lpt[1] *        /* coefficients per coordinate */
                         header->lpt[2]) *       /* subblocks per block */
                         3;                      /* coordinates */
    }
    else {
        header->ncoeff = header->ipt[j][0] - 1 +  /* starting point */
                         (header->ipt[j][1] *     /* coefficients per coordinate */
                         header->ipt[j][2]) *    /* subblocks per block */
                         (j == 11 ? 2 : 3);       /* coordinates */
    }
    header->ksize  = 2 * header->ncoeff;
/*
   Skip to second block in file.
*/
    fseek(infp, header->ncoeff * 8, SEEK_SET);
/*
   Read ephemeris constants.
*/
    for (i=0; i<header->ncon; i++) {
        header->cval[i] = ephcom_indouble(infp);
        if (strncmp(header->cnam[i], "LENUM ", 6) == 0)
            header->numle = header->cval[i];
        else if (strncmp(header->cnam[i], "CLIGHT", 6) == 0)
            header->clight = header->cval[i];
    }
    if (header->numle == 0) 
        header->numle = header->numde;

    return(0);
}




/*
   Read a JPL Ephemeris data block in binary format.

   This is the only routine in this library that accesses a file
   as a direct access file, with a specified block number.  The
   block number ranges from 0 on up (starting at first data block,
   after the 2 header blocks).  Returns the number of coefficients
   read, or 0 at EOF.
*/
int ephcom_readbinary_block(FILE *infp, struct ephcom_Header *header,
                            int blocknum, double *datablock) {

    int i;
    long filebyte;
    double ephcom_indouble(FILE *);
    int fseek(FILE *, long, int);

    filebyte = (blocknum + 2) * header->ncoeff * 8; /* 8 bytes per coefficient */
    fseek(infp, filebyte, SEEK_SET);
    for (i=0; !feof(infp) && i<header->ncoeff; i++)
        datablock[i] = ephcom_indouble(infp);
    if (i < header->ncoeff && feof(infp)) 
        i = -1; /* 0 --> EOF */
    
    i += 1;
    return(i); /* Number of coefficients successfuly read (all or nohing). */
}




/*
   Write header information in ASCII format.
*/
int ephcom_writeascii_header(FILE *outfp, struct ephcom_Header *header) {

    char group[13];
    double val1, val2, val3; /* To read text line with 3 double precision words */
    int i, j, k, n;
    int iword; /* word number we're reading in a line */
    int blockout; /* number of bytes we've written to current block/rec in file */
    static char spaces[84]="                                                                                \r\n";
    int idate[6];
    char *month[12] = {"JAN","FEB","MAR","APR","MAY","JUN",
                   "JUL","AUG","SEP","OCT","NOV","DEC"};

    char writebuf[EPHCOM_MAXLINE + 1];
    char outhcars[EPHCOM_MAXLINE + 1];

//used functions
    size_t fwrite(const void *ptr, size_t size, size_t  nmemb, FILE *stream);
    int ephcom_jd2cal(double tjd, int idate[6], int calendar_type);
    int ephcom_doublestrc2f(char *buf);

/*
   First header line: KSIZE= # NCOEFF= #
*/
    sprintf(writebuf, "KSIZE=%5d    NCOEFF=%5d", header->ksize, header->ncoeff);
    k = strlen(writebuf);
    strcpy(&writebuf[k], &spaces[k]);
    fprintf(outfp, writebuf);
    if (header->ksize != 2*header->ncoeff) {
        fprintf(stderr, "Badly formed header; KSIZE <> 2*NCOEFF\r\n\r\n");
        exit(1);
    }
/*
   GROUP 1010: Title of ephemeris (DE/LE number, start JD, end JD)
*/
    fprintf(outfp, spaces); /* blank line */
    sprintf(writebuf, "GROUP   1010");
    k = strlen(writebuf);
    strcpy(&writebuf[k], &spaces[k]);
    fprintf(outfp, writebuf);
    fprintf(outfp, spaces); /* blank line */
/*
   Header title lines with dates, for example:

      JPL Planetary Ephemeris DE405/LE405
      Start Epoch: JED=  2305424.5 1599 DEC 09 00:00:00
      Final Epoch: JED=  2525008.5 2201 FEB 20 00:00:00
*/
    sprintf(header->ttl[0],"JPL Planetary Ephemeris DE%03d/LE%03d",
            header->numde, header->numle);
    k = strlen(header->ttl[0]);
    strcpy(&header->ttl[0][k], &spaces[k]);
    ephcom_jd2cal(header->ss[0], idate, 0);
    sprintf(header->ttl[1],"Start Epoch: JED=%11.1f%5d %3s %02d %02d:%02d:%02d",
            header->ss[0], idate[0], month[idate[1]-1], idate[2],
            idate[3], idate[4], idate[5]);
    k = strlen(header->ttl[1]);
    strcpy(&header->ttl[1][k], &spaces[k]);
    ephcom_jd2cal(header->ss[1], idate, 0);
    sprintf(header->ttl[2],"Final Epoch: JED=%11.1f%5d %3s %02d %02d:%02d:%02d",
            header->ss[1], idate[0], month[idate[1]-1], idate[2],
            idate[3], idate[4], idate[5]);
    k = strlen(header->ttl[2]);
    strcpy(&header->ttl[2][k], &spaces[k]);
/*
   Don't print trailing blanks at the end of these 3 lines.
*/
    for (i=0; i<3; i++) {
        strncpy(writebuf, header->ttl[i], 80);
        writebuf[80] = '\r';
        writebuf[81] = '\n';
        writebuf[82] = '\0';
        fprintf(outfp, writebuf);
    }
/*
   GROUP 1030: Start and End JD, timestep (in JD) per block.
*/
    fprintf(outfp, spaces); /* blank line */
    sprintf(writebuf, "GROUP   1030");
    k = strlen(writebuf);
    strcpy(&writebuf[k], &spaces[k]);
    fprintf(outfp, writebuf);
    fprintf(outfp, spaces); /* blank line */

    sprintf(writebuf, "%12.2f%12.2f%12.0f.",
            header->ss[0], header->ss[1], header->ss[2]);
    k = strlen(writebuf);
    strcpy(&writebuf[k], &spaces[k]); /* pad with spaces */
    fprintf(outfp, writebuf);
/*
   GROUP 1040: Constant names.
*/
    fprintf(outfp, spaces); /* blank line */
    sprintf(writebuf, "GROUP   1040");
    k = strlen(writebuf);
    strcpy(&writebuf[k], &spaces[k]);
    fprintf(outfp, writebuf);
    fprintf(outfp, spaces); /* blank line */

    sprintf(writebuf, "%6d", header->ncon);
    k = strlen(writebuf);
    strcpy(&writebuf[k], &spaces[k]);
    fprintf(outfp, writebuf);

/*
   Now write the constant names, 10 per line, each 6 characters long
   preceded by 2 blanks.  Pad names with blanks to make 6 characters.
*/
    for (i=0; i<header->ncon; i++) {
       fprintf(outfp, "  %-6s", header->cnam[i]);
       if (i % 10 == 9) fprintf(outfp, "\r\n");
    }
    if (i % 10 != 0) {  /* Pad last line with spaces (i is 1 more than above) */
        for ( ; i % 10 != 0; i++) 
            fprintf(outfp, "        ");
        fprintf(outfp, "\r\n");
    }
/*
   GROUP 1041: Constant values.
*/
    fprintf(outfp, spaces); /* blank line */
    sprintf(writebuf, "GROUP   1041");
    k = strlen(writebuf);
    strcpy(&writebuf[k], &spaces[k]);
    fprintf(outfp, writebuf);
    fprintf(outfp, spaces); /* blank line */

    sprintf(writebuf, "%6d", header->nval);
    k = strlen(writebuf);
    strcpy(&writebuf[k], &spaces[k]);
    fprintf(outfp, writebuf);

    if (header->nval != header->ncon) {
        fprintf(stderr,"Error: number of constants and values not equal.\n\n");
        exit(1);
    }
/*
   Now read constant values, 3 per line, 26 characters each.
*/
    for (i = 0; i < header->ncon; i += 3) {
        val1 = header->cval[i];
        val2 = (i+1 < header->ncon) ? header->cval[i+1] : 0.0;
        val3 = (i+2 < header->ncon) ? header->cval[i+2] : 0.0;
        sprintf(writebuf, "%25.17E %25.17E %25.17E   \r\n", val1, val2, val3);
   /* Note that the character holding the sign for each # is left as is. */
        ephcom_doublestrc2f(&writebuf[1]);
        ephcom_doublestrc2f(&writebuf[27]);
        ephcom_doublestrc2f(&writebuf[53]);
        fprintf(outfp, "%s", writebuf);
    }
/*
   GROUP 1050: Constant values.
*/
    fprintf(outfp, spaces); /* blank line */
    sprintf(writebuf, "GROUP   1050");
    k = strlen(writebuf);
    strcpy(&writebuf[k], &spaces[k]);
    fprintf(outfp, writebuf);
    fprintf(outfp, spaces); /* blank line */
/*
   If there are no coefficients for an ipt[i][] object (i.e., ipt[i][1]==0),
   then ipt[i][0] should contain the value of the next available coefficient
   number rather than 0, as per communication of Myles Standish to Paul Hardy
   on preferred format of ephemeris headers.

   If there are no libration coefficients (i.e., lpt[1]==0), then lpt[0]
   should contain the value of the next available coefficient number rather
   than 0 as well, as per the same communication from Myles Standish.
*/
/* First set j to maximum index into ipt[] that has coefficients */
    j = 0;
    for (i=1; i<12; i++)
        if (header->ipt[i][1] > 0 && header->ipt[i][0] > j)
            j = i;
/* Now set j to next available index count. */
    if (header->lpt[1] > 0 && header->lpt[0] > j)
        j = header->lpt[1] + header->lpt[1] * header->lpt[2] * 3;
    else
        j = header->ipt[j][0] + header->ipt[j][1] * header->ipt[j][2] * (j==11 ? 2 : 3);

    for (i=1; i<12; i++)
        if (header->ipt[i][0] == 0) 
            header->ipt[i][0] = j;

    if (header->lpt[0] == 0) 
        header->lpt[0] = j;
/*
   Write ipt array in transposed order (arrays are transposed in FORTRAN
   from their order in C).
*/
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 12; j++)
            fprintf(outfp, "%6d", header->ipt[j][i]);
        fprintf(outfp, "%6d  \r\n", header->lpt[i]);
    }
/*
   GROUP 1070: Constant values.
*/
    fprintf(outfp, spaces); /* blank line */
    sprintf(writebuf, "GROUP   1070");
    k = strlen(writebuf);
    strcpy(&writebuf[k], &spaces[k]);
    fprintf(outfp, writebuf);
    fprintf(outfp, spaces); /* blank line */
/*
   Now we're pointing to the first block of coefficient data, after header.
*/
return(0);
}




/*
   Write coefficient block information in ASCII format.
*/
int ephcom_writeascii_block(FILE *outfp, struct ephcom_Header *header,
                            int blocknum, double *datablock) {

    double val1, val2, val3; /* To write text line with 3 double precision words */
    int i;
    char writebuf[EPHCOM_MAXLINE + 1];

    size_t fwrite(const void *ptr, size_t size, size_t  nmemb, FILE *stream);
    int fputc(int, FILE *);
    double ephcom_indouble(FILE *);
    int ephcom_doublestrc2f(char *); /* Convert C formatted double to FORTRAN format */

/*
   Write first line in block, which is block number and ncoeff.
   Note that lines in the data block files end in "\r\n", while
   lines in the header files just end in "\n".
*/
    fprintf(outfp, "%6d%6d", blocknum + 1, header->ncoeff);
    for (i=0; i<68; i++) 
        fputc(' ', outfp);
    fprintf(outfp, "\r\n");
/*
   Now write the data, 3 coefficients per line, 26 characters each.
   Convert format to match what appears in JPL Ephemeris ASCII files.
*/
    for (i = 0; i < header->ncoeff; i += 3) {
        val1 = datablock[i];
        val2 = (i+1) < header->ncoeff ? datablock[i+1] : 0.0;
        val3 = (i+2) < header->ncoeff ? datablock[i+2] : 0.0;
   /*
      Write values, 3 coefficients per line, pad lines with 0.0E+00
   */
        sprintf(writebuf, "%25.17E %25.17E %25.17E   \r\n", val1, val2, val3);
   /*
      Now re-format numbers the way the JPL header file writes them:
      all with a leading "0.", so the exponent is one greater.
      Note that here we start at strlen(writebuf)-6, but in the section
      that handles header data coefficients we start at strlen(writebuf)-5.
      This is because the data blocks end ASCII lines with "\r\n" whereas
      the header ends ASCII lines with just "\n".
   */
        ephcom_doublestrc2f(&writebuf[1]);  /* Reformat first number */
        ephcom_doublestrc2f(&writebuf[27]);  /* Reformat second number */
        ephcom_doublestrc2f(&writebuf[53]);  /* Reformat third number */
// printf("[%d]%s", strlen(writebuf), writebuf);
        fprintf(outfp, "%s", writebuf);
    }

    return(0);
}




/*
   Write a JPL Ephemeris header in binary format.
*/
int ephcom_writebinary_header(FILE *outfp, struct ephcom_Header *header) {

    char group[13];  /* To hold "GROUP" header line */
    int blockout; /* number of bytes we've written to current block/rec in file */
    int blockbytes; /* number of bytes in a block, equals 8 * ncoeff */
    double val1, val2, val3; /* To read text line with 3 double precision words */
    int i, j;
    int idate[6];
    char *month[12] = {"JAN", "FEB", "MAR", "APR", "MAY", "JUN",
                       "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};

    int ephcom_outdouble(FILE *outfp, double x);
    int ephcom_outint(FILE *outfp, unsigned u);
    void ephcom_nxtgrp(char *group, char *expected, FILE *infile);
    int ephcom_jd2cal(double tjd, int idate[6], int calendar_type);
    size_t fwrite(const void *ptr, size_t size, size_t  nmemb, FILE *stream);


    rewind(outfp);    //Point to beginning of output file.
    blockbytes = sizeof(double) * header->ncoeff;    //First header line: KSIZE= # NCOEFF= #
/*
   Start writing output ephemeris, beginning with header.
*/
/*
   Header title lines with dates, for example:

      JPL Planetary Ephemeris DE405/LE405
      Start Epoch: JED=  2305424.5 1599 DEC 09 00:00:00
      Final Epoch: JED=  2525008.5 2201 FEB 20 00:00:00
*/
    sprintf(header->ttl[0],"JPL Planetary Ephemeris DE%03d/LE%03d",
            header->numde, header->numle);
    for (i=strlen(header->ttl[0]); i<84; i++) 
        header->ttl[1][i] = ' ';
    ephcom_jd2cal(header->ss[0], idate, 0);
    sprintf(header->ttl[1],"Start Epoch: JED=%11.1f%5d %3s %02d %02d:%02d:%02d",
            header->ss[0], idate[0], month[idate[1]-1], idate[2], idate[3], idate[4], idate[5]);
    for (i=strlen(header->ttl[1]); i<84; i++)
        header->ttl[1][i] = ' ';
    ephcom_jd2cal(header->ss[1], idate, 0);
    sprintf(header->ttl[2],"Final Epoch: JED=%11.1f%5d %3s %02d %02d:%02d:%02d",
            header->ss[1], idate[0], month[idate[1]-1], idate[2], idate[3], idate[4], idate[5]);
    for (i=strlen(header->ttl[2]); i<84; i++)
            header->ttl[2][i] = ' ';
    header->ttl[0][84] = header->ttl[1][84] = header->ttl[2][84] = '\0';

/*
   ephcom_Header title lines.

   Write the three title lines to the output file, padded with blanks,
   84 characters long (84 is the first even multiple of 6 that is > 80,
   so the JPL software uses that value because it reads in Fortran 'A6'
   character strings.
*/
    fprintf(outfp, "%-84s%-84s%-84s", header->ttl[0], header->ttl[1], header->ttl[2]);
    blockout = 3*84;  /* Just wrote 3 84-byte strings to start output file */
/*
   Now output 400 cnam entries to the output file.
*/
    for (i = 0; i < header->ncon; i++) {
        fprintf(outfp, "%-6s", header->cnam[i]);
        blockout += 6;
    }
    for ( ; i < 400; i++) {
        fprintf(outfp, "      ");  /* Round out to 400 entries, all blank at end */
        blockout += 6;
    }
/*
   Binary values: Make sure bytes are in big-endian (network) order for file.
*/
    for (i=0; i<3; i++) {
        ephcom_outdouble(outfp, header->ss[i]);  /* Write net-order bytes from double precision */
        blockout += 8;
    }
    ephcom_outint(outfp, header->ncon);
    blockout += 4;
    ephcom_outdouble(outfp, header->au);
    blockout += 8;
    ephcom_outdouble(outfp, header->emrat);
    blockout += 8;
/*
   If there are no coefficients for an ipt[i][] object (i.e., ipt[i][1]==0),
   then ipt[i][0] should contain the value of the next available coefficient
   number rather than 0, as per communication of Myles Standish to Paul Hardy
   on preferred format of ephemeris headers.

   If there are no libration coefficients (i.e., lpt[1]==0), then lpt[0]
   should contain the value of the next available coefficient number rather
   than 0 as well, as per the same communication from Myles Standish.
*/
/* First set j to maximum index into ipt[] that has coefficients */
    j = 0;
    for (i=1; i<12; i++)
        if (header->ipt[i][1] > 0 && header->ipt[i][0] > j)
            j = i;
/* Now set j to next available index count. */
    if (header->lpt[1] > 0 && header->lpt[0] > j)
        j = header->lpt[1] + header->lpt[1] * header->lpt[2] * 3;
    else
        j = header->ipt[j][0] + header->ipt[j][1] * header->ipt[j][2] * (j==11 ? 2 : 3);
    for (i=1; i<12; i++)
        if (header->ipt[i][0] == 0) 
            header->ipt[i][0] = j;
    if (header->lpt[0] == 0) 
        header->lpt[0] = j;

    for (j = 0; j < 12; j++) {
        for (i = 0; i < 3; i++) {
            ephcom_outint(outfp, header->ipt[j][i]);
            blockout += 4;
        }
    }
    ephcom_outint(outfp, header->numde);
    blockout += 4;
    for (i = 0; i < 3; i++) {
        ephcom_outint(outfp, header->lpt[i]);
        blockout += 4;
    }
/*
   Now pad the end of the first record with null bytes.  Note: the
   JPL Fortran software just skips to next record at this point.
*/
    for (i = blockout; i < blockbytes; i++) {
        fputc('\0', outfp);
    }
/*
   End of first block.  Now set blockout to 0 and start with next block.
*/
    blockout = 0;
    for (i=0; i<header->ncon; i++) {
        ephcom_outdouble(outfp, header->cval[i]);
        blockout += 8;
    }
/*
   Pad with double-precision zeroes for rest of array.
*/
    for ( ; i < 400; i++) {
        ephcom_outdouble(outfp, (double)0.0);
        blockout += 8;
    }
/*
   Pad with nulls for rest of block.
*/
    for (i = blockout; i < blockbytes; i++) {
        fputc('\0', outfp);
    }
/*
   Finished normally.
*/
    return(0);
}




/*
   Write a block of data coefficients in JPL binary file format.
*/
int ephcom_writebinary_block(FILE *outfp, struct ephcom_Header *header, 
                             int blocknum, double *datablock) {

    int i;
    int ephcom_outdouble(FILE *, double);
    int filebyte;
    int filepos;

/*
   Find out where we need to point in the binary file.
*/
    filebyte = (blocknum + 2) * header->ncoeff * 8; /* 8 bytes per coefficient */
/*
   If the file isn't that large, pad it with null bytes
*/
    fseek(outfp, 0L, SEEK_END);
    filepos = ftell(outfp);
    if (filepos < filebyte) {
        for (i=0; i < (filebyte - filepos); i++) {
            fputc('\0', outfp);
        }
    }
/*
   Now go to position where we want to start writing.
*/
    fseek(outfp, filebyte, SEEK_SET);
    for (i = 0; i < header->ncoeff; i++) {
        ephcom_outdouble(outfp, datablock[i]);
    }

    return(0);
}




/*
   ephcom_parse_block() - Parse a binary block of data.  Warning: verbose!
                          Writes parsed output to file pointer outfp.
*/
int ephcom_parse_block(FILE *outfp, struct ephcom_Header *header, double *datablock) {

    int i0, i1, i2, i3;
    int blockword;
/*
   Names of the objects in Chebyshev coefficient arrays.
*/
    static char *ephcom_coeffname[13] = {
        "Mercury", "Venus", "EMBary", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune",
        "Pluto", "Moon", "Sun", "Nutation", "Libration"};

    blockword=0;
    fprintf(outfp, "@%04d StartJD\t%12.2f\n", blockword++, datablock[0]);
    fprintf(outfp, "@%04d EndJD\t%12.2f\n", blockword++, datablock[1]);

    for (i0=0; i0<12; i0++) {  /* For first 12 bodies */
        fprintf(outfp, "Body\t%d (%s)\n", i0+1, ephcom_coeffname[i0]);
        for (i1=0; i1<header->ipt[i0][2]; i1++) { /* For all subintervals */
            fprintf(outfp, "  Subinterval %d of %d\n", i1+1, header->ipt[i0][2]);
            for (i2=0; i2 < (i0==11 ? 2 : 3); i2++) {  /* For all coordinates */
                fprintf(outfp, "    %cCoefficients\n", 'X' + i2);
                for (i3=0; i3<header->ipt[i0][1]; i3++) {  /* For all coefficients */
                    blockword = header->ipt[i0][0] +
                                i1*header->ipt[i0][1] * (i0==11 ? 2 : 3) +
                                i2*header->ipt[i0][1] + i3 - 1;
                    fprintf(outfp, "      @%04d [%2d of %2d] %25.17E\n",
                            blockword, i3+1, header->ipt[i0][1], datablock[blockword]);
                }
            }
        }
    }
// For libration    
    fprintf(outfp, "Body\t13 (%s)\n", ephcom_coeffname[12]);
        for (i1=0; i1<header->lpt[2]; i1++) { /* For all subintervals */
            fprintf(outfp, "  Subinterval %d of %d\n", i1+1, header->lpt[2]);
            for (i2=0; i2 < 3; i2++) {  /* For all coordinates */
                switch (i2){
                    case 0: 
                        fprintf(outfp, "    %sCoefficients\n", "phi");
                        break;
                    case 1:
                        fprintf(outfp, "    %sCoefficients\n", "theta");
                        break;
                    case 2:
                        fprintf(outfp, "    %sCoefficients\n", "psi");
                        break;
                }
                for (i3=0; i3<header->lpt[1]; i3++) {  /* For all coefficients */
                    blockword = header->lpt[0] +
                                i1*header->lpt[1] *3 +
                                i2*header->lpt[1] + i3 - 1;
                    fprintf(outfp, "      @%04d [%2d of %2d] %25.17E\n",
                            blockword, i3+1, header->lpt[1], datablock[blockword]);
                }
            }
        }

    return(0);
}




/*
   Get the next group header in the file.
   Group header lines are 12 characters long, of the form "GROUP   nnnn".

   Parameters:
      group - the GROUP header we read
      expected - the header we expected
      infile - the file pointer to read
*/
void ephcom_nxtgrp(char *group, char *expected, FILE *infile) {

    char readbuf[EPHCOM_MAXLINE + 1];
    char *fgets(char *, int, FILE *);

    fgets(readbuf, EPHCOM_MAXLINE, infile); /* Blank Line     */
    fgets(readbuf, EPHCOM_MAXLINE, infile); /* "GROUP   dddd\n" */
    strncpy(group, readbuf, 12);
    group[12] = '\0';
    if (strncmp(group, expected, 12) != 0) {
        fprintf(stderr, "Badly formed header; \"%s\" not found.\n\n", expected);
        exit(1);
    }
    fgets(readbuf, EPHCOM_MAXLINE, infile); /* Blank Line    */
}




/*
   Print a double precision value to the given file with bytes swapped
   if necessary to match network order (Big Endian).  On Intel 80x86
   the bytes will get swapped, on Motorola or SPARC they won't.
*/
int ephcom_outdouble(FILE *outfp, double x) {
   double retval;
   unsigned char ch[8];
   unsigned char *gnulliver64c(unsigned char *);
   void *memcpy(void *, const void *, size_t);

   memcpy((void *)ch, (const void *)&x, 8);
// printf("Got OD: %25.17e => %02x %02x %02x %02x %02x %02x %02x %02x, ",
//        x, ch[0], ch[1], ch[2], ch[3], ch[4], ch[5], ch[6], ch[7]);
   (void)gnulliver64c(ch);
// printf("Sending OD: %02x %02x %02x %02x %02x %02x %02x %02x\n",
//        ch[0], ch[1], ch[2], ch[3], ch[4], ch[5], ch[6], ch[7]);
   fwrite(ch, 1, 8, outfp);
   return(0);
}




/*
   Print an integer (4-byte) value to the given file with bytes swapped
   if necessary to match network order (Big Endian).  On Intel 80x86
   the bytes will get swapped, on Motorola or SPARC they won't.
*/
int ephcom_outint(FILE *outfp, unsigned u) {
   unsigned u2;
   unsigned gnulliver32(unsigned);

   u2 = gnulliver32(u);
   fwrite(&u2, 4, 1, outfp);
   return(0);
}




/*
   Read in a double precision value from the given file with bytes swapped
   if necessary to match host order (Little- or DEC- Endian).  On Intel 80x86
   the bytes will get swapped, on Motorola or SPARC they won't.
*/
double ephcom_indouble(FILE *infp) {
   double x;
   static double retval;
   size_t fread(void *ptr,  size_t size, size_t nmemb, FILE *stream);
   unsigned char ch[8];
   unsigned char *gnulliver64c(unsigned char *);
   void *memcpy(void *, const void *, size_t);
   
   /*
      Handle as character string until bytes are in correct order,
      then copy to double once they are.
   */
   fread(ch, 1, 8, infp);
// printf("Got ID: %02x %02x %02x %02x %02x %02x %02x %02x, ",
//        ch[0], ch[1], ch[2], ch[3], ch[4], ch[5], ch[6], ch[7]);
   (void)gnulliver64c(ch);
   memcpy((void *)&retval, (const void *)ch, (size_t)8);
// printf("Sending ID: %02x %02x %02x %02x %02x %02x %02x %02x [%25.18E]\n",
//        ch[0], ch[1], ch[2], ch[3], ch[4], ch[5], ch[6], ch[7], retval);
   return(retval);
}




/*
   Read in an integer (4--byte) value to the given file with bytes swapped
   if necessary to match host order (Little- or DEC- Endian).  On Intel 80x86
   the bytes will get swapped, on Motorola or SPARC they won't.
*/
int ephcom_inint(FILE *infp) {
   unsigned u;
   static int retval;
   unsigned gnulliver32(unsigned);

   fread(&u, 4, 1, infp);
   retval = (int)gnulliver32(u);
   return(retval);
   }




/*
   ephcom_doublstrc2f() - function to convert a string with a double precision
                          value written in C to a double precision value that
                          FORTRAN creates.  Conversion happens in place.
*/

int ephcom_doublestrc2f(char *buf) {

    int i, j, istart, istop, exp, edigits;
    double x;

    for (istart=0; isspace(buf[istart]); istart++);
    x = atof(&buf[istart]);
    for (istop=istart; toupper(buf[istop]) != 'E'; istop++);
    exp = atoi(&buf[istop+1]);
    exp++;
    if (exp < 0) {
        buf[istop+2] = '-';
        exp = -exp;
    }
    else {
        buf[istop+2] = '+';
    }
    if (x == 0.0) 
        exp=0;
    if (exp < 100) 
        edigits = 2;
    else if (exp < 1000) 
        edigits = 3;
    else edigits = 4;

    while (edigits > 0) {
        buf[istop + edigits + 2] = exp % 10 + '0';
        exp /= 10;
        edigits--;
    }

    buf[istop+1] = 'D';

    while (istop > istart && buf[istop-1] != '.') {
        buf[istop] = buf[istop-1];
        istop--;
    }

    buf[istop] = buf[istop-2];  /* buf[istop-1] == '.' */
    buf[istop-2] = '0';         /* leading zero */

    return(0);
}




/*
   Planetary Ephemeris.  Takes coordinates already calculated in
   coords structure an converts to vectors and vector dot in testr[].
   Bodies start at 1 for Mercury, to match the JPL PLEPH() numbering.
   Values for ntarg and ncntr correspond to locations ntarg-1 and
   ncntr-1 in coords->pv[].
*/
int ephcom_pleph(struct ephcom_Coords *coords, int ntarg, int ncntr, double *r) {

    int i,j;

    if (ntarg != 14 && ncntr != 14) { /* If not nutation, handle normally */
        if (ntarg == 15 || ncntr == 15) { /* Libration */
            for (i=0; i<6; i++)
                r[i] = coords->pv[14][i];
        }
        else {
            for (i=0; i<6; i++)
                r[i] = coords->pv[ntarg-1][i] - coords->pv[ncntr-1][i];
        }
    }
    else { /* Nutation */
        r[0] = coords->pv[13][0];
        r[1] = coords->pv[13][1];
        r[2] = coords->pv[13][2];
        r[3] = coords->pv[13][3];
        r[4] = 0.0;
        r[5] = 0.0;
   }

    return(0);
}




/*
   ephcom_get_coords() - Interpolate positions and velocities at given time.
*/
int ephcom_get_coords(FILE *infp, struct ephcom_Header *header,
                      struct ephcom_Coords *coords, double *datablock) {

    double et2[2];    /* Ephemeris time, as coarse (whole) and fine time  in JD */
    double totaltime; /* Sum of whole and fractional JD */
    double filetime;  /* JDs since start of ephemeris file */
    double blocktime; /* JDs since start of data block */
    double subtime;   /* JDs since start of subinterval in block */
    int i, j, k;
    int blocknum;
    int nsub; /* Number of subintervals in data block for this body */
    int subinterval; /* Number of subinterval for this body */
    int dataoffset; /* Offset in datablock for current body and subinterval */
    double subspan; /* Span of one subinterval in days */
    double chebytime; /* Normalized Chebyshev time, in interval [-1,1]. */
    int ncoords; /* Number of coordinates for position and velocity */
    int ncf; /* Number of Chebyshev coefficients per coordinate */
    int retval; /* Return value */

//Declaration of used functions
    int ephcom_cheby(int maxcoeffs, double x, double span, double *y, 
                     int ncoords, int ncoeffs, double *pv);

    retval = 0; /* Assume normal return */
/*
   Split time JD into whole JDs (et2[0]) and fractional JD (et2[1]).
*/
    totaltime = coords->et2[0] + coords->et2[1];
    if (totaltime < header->ss[0] || totaltime > header->ss[1]) {
        fprintf(stderr,"Time is outside ephemeris range.\n");
        retval = -1;
    }
    else {
        et2[0] = (int)totaltime;
        et2[1] = (coords->et2[0] - et2[0]) + coords->et2[1];
        filetime = totaltime - header->ss[0]; /* Days from start of file */
        blocknum = (int)(filetime / header->ss[2]); /* Data block in file, 0.. */
   /*
      Read the data block that contains coefficients for desired date
   */
        ephcom_readbinary_block(infp, header, blocknum, datablock);
   /*
      Now step through the bodies and interpolate positions and velocities.
   */
        blocktime = totaltime - datablock[0]; /* Days from block start */
        for (i=0; i<13; i++) {
            if (i == 12)
                subspan = header->ss[2] / header->lpt[2];
            else
                subspan = header->ss[2] / header->ipt[i][2]; /* Days/subinterval */
            subinterval = (int)((totaltime - datablock[0]) / subspan);

            ncoords = (i == 11 ? 2 : 3); /* 2 coords for nutation, else 3 */

            if (i == 12)
                dataoffset = header->lpt[0] - 1 +
                             ncoords * header->lpt[1] * subinterval;
            else
                dataoffset = header->ipt[i][0] - 1 +
                             ncoords * header->ipt[i][1] * subinterval;

            subtime = blocktime - subinterval * subspan;
      /*
         Divide days in this subblock by total days in subblock
         to get interval [0,1].  The right part of the expression
         will evaluate to a whole number: subinterval lengths are
         all integer multiples of days in a block (all powers of 2).
      */
            chebytime = subtime / subspan;
            chebytime = 2.0*chebytime - 1.0;
            if (chebytime < -1.0 || chebytime > 1.0) {
                fprintf(stderr, "Chebyshev time is beyond [-1,1] interval.\n");
                fprintf(stderr, "filetime=%f, blocktime=%f, subtime=%f, chebytime=%f\n",
                        filetime, blocktime, subtime, chebytime);
            }
            else {
                if (i == 12)
                    ephcom_cheby(header->maxcheby, chebytime, subspan, &datablock[dataoffset],
                             ncoords, header->lpt[1], coords->pv[i]);
                else 
                    ephcom_cheby(header->maxcheby, chebytime, subspan, &datablock[dataoffset],
                             ncoords, header->ipt[i][1], coords->pv[i]);
            }
         /*
            Everything is as expected.  Interpolate coefficients.
         */
                
        }
   /*
      With interpolations complete, calculate Earth from EMBary and
      Sun from SSBary.  Preserve other coordinates.
   */
        for (j=0; j<6; j++) {
            coords->pv[15][j] = coords->pv[ 9][j]; /* Save original lunar coords */
            coords->pv[14][j] = coords->pv[12][j]; /* Librations if on file */
            coords->pv[13][j] = coords->pv[11][j]; /* Nutations if on file */
            coords->pv[11][j] = 0.0;
      /*
         Calculate Earth and Moon from EMBary and geocentric Moon.
      */
            coords->pv[12][j] = coords->pv[2][j]; /* Move EMBary from Earth spot */
            coords->pv[2][j] -= coords->pv[9][j] / (1.0 + header->emrat); /* Earth */
            coords->pv[9][j] += coords->pv[2][j]; /* Moon (change geo->SS-centric) */
        }

        if (!coords->km) { /* Calculate AU, not kilometers */
            for (i=0; i<16; i++ ) {
                if (i == 13) i = 15; /* Skip over nutations and librations */
                for (j=0; j<6; j++)
                    coords->pv[i][j] /= header->au;
            }
        }

        if (coords->seconds){ /* for km/sec, 86400 sec/day */
            for (i=0; i<16; i++){
                if (i==13 || i==14){
                    coords->pv[i][2] /= 86400.0;
                    coords->pv[i][3] /= 86400.0;
                }
                else{
                    for (j=3; j<6; j++)
                        coords->pv[i][j] /= 86400.0;
                }
            }
        }

    }

    return(retval);
}




/*
   ephcom_cheby() - interpolate at a point using Chebyshev coefficients
*/
inline int ephcom_cheby(
    int maxcoeffs, /* Maximum number of Chebyshev components possible */
    double x,      /* Value of x over [-1,1] for Chebyshev interpolation */
    double span,   /* Span in time of subinterval, for velocity */
    double *y,     /* Chebyshev coefficients */
    int ncoords,   /* Total number of coordinates to interpolate */
    int ncoeffs,   /* Number of Chebyshev coefficients per coordinate */
    double *pv     /* Array to hold position in 1st half, velocity in 2nd */
    ) {

    int i, j;
    static double *pc, *vc; /* Position and velocity polynomial coefficients. */
    static double lastx=2.0; /* x from last call; initialize to impossible value */
    static int init=1; /* Need to initialize pc[] and vc[] */

/*
   Allocate position and velocity Chebyshev coefficients.
*/
    if (init) {
        pc = (double *)malloc(maxcoeffs * sizeof(double));
        vc = (double *)malloc(maxcoeffs * sizeof(double));
        init = 0;
    }
/*
   This need only be called once for each Julian Date,
   saving a lot of time initializing polynomial coefficients.
*/
    if (lastx != x) {
        lastx = x;
   /*
      Initialize position polynomial coefficients
   */
        pc[0] = 1.0;    /* Chebyshev T[0](x) = 1 */
        pc[1] = x;      /* Chebyshev T[1](x) = x */
        for (i=2; i<maxcoeffs; i++) {
            pc[i] = 2.0*x * pc[i-1] - pc[i-2];
       /*
         Resolve bug with gcc generating -0.0 (also makes
         the smallest represented number equal to zero).
      */
            if (pc[i]*pc[i] == 0.0) 
                pc[i] = 0.0;
        }
   /*
      Initialize derivative polynomial coefficients
   */
        vc[0] = 0.0;          /* d(1)/dx        = 0  */
        vc[1] = 1.0;          /* d(x)/dx        = 1  */
        for (i=2; i<maxcoeffs; i++) 
            vc[i] = 2.0*x * vc[i-1] + 2*pc[i-1] - vc[i-2];
    }
/*
   Interpolate to get position for each component
*/
    for (i=0; i<ncoords; i++) { /* Once each for x, y, and z */
        pv[i] = 0.0;
        for (j=ncoeffs-1; j >= 0; j--) 
            pv[i] += pc[j] * y[i*ncoeffs + j];
    }
/*
   Interpolate velocity (first derivative)
*/
    for (i=0; i<ncoords; i++) {
        pv[ncoords + i] = 0.0;
        for (j=ncoeffs-1; j >= 0; j--) 
            pv[ncoords + i] += vc[j] * y[i*ncoeffs + j];
            pv[ncoords + i] *= 2.0 / span;
    }

    return(0);
}




/*
   ephcom_jd2cal() - convert Julian Day to calendar date and time.

      tjd: double precision Julian Day
      idate: integer year, month, day, hour, minute, second of tjd
      calendar_type: -1=Julian; 0=Automatic; 1=Gregorian

   If automatic, use Julian calendar for dates before 15 October 1582.

   From pp. 604, 606 in the Explanatory Supplement to the Astronomical Almanac.
*/
int ephcom_jd2cal(double tjd, int idate[6], int calendar_type) {

    int ihour, imin, isec;
    int j;
    int I, J, K, L, N, D, M, Y;

    tjd += 0.5 + 0.5/86400.0; /* Round to nearest second */
    j = tjd;  /* Integer Julian Day */
    tjd = (tjd - j) * 24.0;
    ihour = tjd;
    tjd = (tjd - ihour) * 60.0;
    imin = tjd;
    tjd = (tjd - imin) * 60.0;
    isec = tjd;
/*
   Julian calendar.  Explanatory Supplement to Astronomical Alamanac, p. 606.
   If automatic, use Julian calendar for dates before 15 October 1582.
*/
    if (calendar_type == -1 || (calendar_type == 0 && j <= 2299160)) {
        J = j + 1402;
        K = (J - 1) / 1461;
        L = J - 1461 * K;
        N = (L - 1) / 365 - L / 1461;
        I = L - 365 * N + 30;
        J = (80 * I) / 2447;
        D = I - (2447 * J) / 80;
        I = J / 11;
        M = J + 2 - 12 * I;
        Y = 4 * K + N + I - 4716;
    }
/*
   Gregorian calendar.
*/
    else { /* Explanatory Supplement to Astronomical Almanac, p. 604 */
        L = j + 68569;
        N = (4 * L) / 146097;
        L = L - (146097 * N + 3) / 4;
        I = (4000 * (L + 1)) / 1461001;
        L = L - (1461 * I) / 4 + 31;
        J = (80 * L) / 2447;
        D = L - (2447 * J) / 80;
        L = J / 11;
        M = J + 2 - 12 * L;
        Y = 100 * (N - 49) + I + L;
    }

    idate[0] = Y;
    idate[1] = M;
    idate[2] = D;
    idate[3] = ihour;
    idate[4] = imin;
    idate[5] = isec;

    return(0);
}




/*
   ephcom_cal2jd() - convert calendar date and time to JD.

      idate: integer year, month, day, hour, minute, second
      calendar_type: -1=Julian; 0=Automatic; 1=Gregorian
      return value: double precision Julian Day of idate[]

   From pp. 604, 606 in the Explanatory Supplement to the Astronomical Almanac.
*/
double ephcom_cal2jd(int idate[6], int calendar_type) {

    double tjd;
    int jd;

/*
   Convert hours, minutes, seconds to fractional JD.
*/
    tjd = (idate[3] + (idate[4] + idate[5] / 60.0) / 60.0) / 24.0 - 0.5;
/*
   Julian calendar.  Explanatory Supplement to Astronomical Alamanac, p. 606.
   If automatic, use Julian calendar for dates before 15 October 1582.
*/
    if (calendar_type == -1 ||
        (calendar_type == 0 && 
         (idate[0] < 1582 ||                         /* Before 1582 */
          (idate[0] == 1582 &&
           (idate[1] < 10 ||                         /* Before October 1582 */
             (idate[1] == 10 && idate[2] < 15)))))) { /* Before 15 October 1582 */
        jd = 367 * idate[0] -
             (7 * (idate[0] + 5001 + (idate[1] - 9) / 7)) / 4 +
             (275 * idate[1]) / 9 +
             idate[2] + 1729777;
    }
/*
   Gregorian calendar.
*/
    else { /* Explanatory Supplement to Astronomical Almanac, p. 604 */
        jd = (1461 * (idate[0] + 4800 + (idate[1] - 14) / 12)) / 4 +
             (367 * (idate[1] - 2 - 12 * ((idate[1] - 14) / 12))) / 12 -
             (3 * ((idate[0] + 4900 + (idate[1] - 14) / 12) / 100)) / 4 +
             idate[2] - 32075;
    }
/*
   Return value is whole JD number plus fractional JD number.
*/
    tjd += (double)jd;

    return(tjd);
}