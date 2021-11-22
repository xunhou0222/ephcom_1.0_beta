                                ephcom_1.0_beta        By Jiang Tingwei


1、读取的星历不可以用fortran程序生成的二进制星历文件，而要用eph文件夹里的
     二进制文件”JPLEPH421"。这是因为C语言与Fortran在二进制数据输出方面存在
     big/little endian的差异。

2、所用到的实现二进制星历文件读取、计算位置速度的程序在src文件夹下，相关函
     数在ephcom.c里，但是三个文件都要用到。

3、函数的用法可以参考example里的程序（来自于源码包），必要时可以阅读函数
     源码，或者与我交流。
     用到的ephcom中的主要函数如下：

         (1) int ephcom_readbinary_header(FILE *infp, struct ephcom_Header *header)

              从二进制星历文件（e.g., JPLEPH421）中读取头文件里的信息

              infp:：文件指针，指向二进制星历文件；

              header：结构体指针（详见头文件ephcom.h），用于存储读取到的头文件信息；


         (2) int ephcom_get_coords(FILE *infp,
                      struct ephcom_Header *header,
                      struct ephcom_Coords *coords,
                      double *datablock);
              
              利用header里的信息，和coords中的设定值，计算某一历元时刻的所有天体的质心天球坐标。

              infp、header：含义同上
             
              coords：结构体指针（详见头文件ephcom.h），设置物理量单位、历元时刻的儒略日，存储
      该时刻所有天体的质心天球坐标。
                 
              datablock：double类型指针，需要用下列语句动态分配内存
                        datablock = (double *)malloc(header.ncoeff * sizeof(double));

              注意：该函数调用之前，务必对header、coords初始化


         （3) int ephcom_pleph(struct ephcom_Coords *coords, int ntarg, int ncntr, double *r);
              
                计算目标天体相对于中心天体的位置速度

                coords：含义同上
   
                ntarg：目标天体代号（天体与代号的对应关系详见头文件ephcom.h）

                ncntr：中心天体代号，同上

                r：即r[6]，六元素double数组，用于存储目标天体相对于中心天体的位置速度矢量