#include <iostream>
#include "sintering.hpp"
#include <stdio.h>
#include "EasyBMP.h"
#include <limits.h>

int main(int argc, char **argv){
    
    register int i, j, k;
    double endtime = 1000.0;
    double dcdt,conc,time,prefactor;
    double temp[NUMOP][XSIZE][YSIZE],chempot[NUMOP][XSIZE][YSIZE];
    double bulkderiv[NUMOP][XSIZE][YSIZE];
    double maxop;
    double blue,green,red;
    char oldtime[128],newtime[128];
    struct Vec term[NUMOP][XSIZE][YSIZE];

    //seed random number generator
    //
    // FILE *seed = fopen("/dev/urandom", "r");
    FILE *fpout;
    // unsigned int random_data;
    // size_t seedsize = fread(&random_data, sizeof(random_data), 1, seed);
    // fclose(seed);

    // srand(random_data);

    // initialize microstructure
    
    std::cout << "Initializing microstructure... ";
    std::cout.flush();
    init_mic();
    std::cout << "Done!" << std::endl;
    std::cout.flush();

    //Generate a blank image that can be updated and output as an image as the program runs
    BMP* img = new BMP();
    img -> SetSize(XSIZE,YSIZE);
    //Declare a path length for naming the image files
    char imgpath[PATH_MAX];

    std::cout << "Entering main loop... " << std::endl;
    std::cout.flush();
    //Advance through time
    sprintf(oldtime,"0");

    sprintf(newtime,"%d",(int)(time));
    /*
    sprintf(imgpath,"%s.txt", argv[1]);
    if ((fpout = fopen(imgpath,"w")) == NULL) {
        printf("\nERROR: Could not open %s",imgpath);
        exit(1);
    }
    */

    for (time = 0; time <= endtime; time += Gdt){
  //      std::cout << "Calculating chempot... ";
  //      std::cout.flush();
        for (i = 0; i < XSIZE; i++){
            for (j = 0; j < YSIZE; j++){
                for (k = 0; k < NUMOP; k++) {
                    bulkderiv[k][i][j]  = calcBulkDeriv(Gop,k,i,j);
                    chempot[k][i][j] = bulkderiv[k][i][j] - (Geps2[k] * laplac(Gop,k,i,j));
                }
            }
        }
 //       std::cout << "Done!" << std::endl;
  //      std::cout << "Calculating Gfield... ";
 //       std::cout.flush();

        // Now the inner scalar term is calculated everywhere, we next need
        // to create the inner vector field

        for (i = 0; i < XSIZE; i++){
            for (j = 0; j < YSIZE; j++){
                Gfield[0][i][j] = grad(chempot,0,i,j);
            }
        }

  //      std::cout << "Done!" << std::endl;
  //      std::cout << "Updating order parameters... ";
  //      std::cout.flush();

        // Now all we need to do is take the divergence of the vector field Gfield,
        // already calculated at each point above, and then multiply by the mobility
       
        for (i = 0; i < XSIZE; i++){
            for (j = 0; j < YSIZE; j++){
                dcdt = GMc * div(Gfield,0,i,j);
                dcdt = GMc * div(Gfield,0,i,j);
                temp[0][i][j] = Gop[0][i][j] + (dcdt * Gdt);
                for (k = 1; k < NUMOP; k++) {
                    temp[k][i][j] = Gop[k][i][j] - ((GMnc * chempot[k][i][j]) * Gdt);
                }
            }
        }

        // The new concentrations are stored in the temp array. Transfer them to the
        // concentration array.

  //      std::cout << "Done!" << std::endl;
  //      std::cout << "Updating order parameters... ";
  //      std::cout.flush();

        sprintf(newtime,"%d",(int)(time));
        if (argc > 1 && (strcmp(oldtime,newtime) || time < Gdt)) { 
            for (j = 0; j < YSIZE; ++j) {
                for (i = 0; i < XSIZE; ++i) {
                    maxop= -5000.0;
                    for (k = 1; k < NUMOP; k++) {
                        if (Gop[k][i][j] > maxop) maxop = Gop[k][i][j];
                    }
                    maxop = pow(maxop*maxop,0.5);


                    if (maxop > 1.0) maxop = 1.0;

                    blue = 1.0-Gop[0][i][j];
                    red = green = maxop;

                    //set pixel colors, can replace with commented lines for grayscale instead of blue/pink
                    img -> SetPixel(i, j, (RGBApixel){
                        .Blue = (ebmpBYTE)(255*(blue)),
                        .Green =(ebmpBYTE)(255*(green)),
                        .Red = (ebmpBYTE)(255*(red)),
                        .Alpha = 0,
                    });
                    snprintf(imgpath, sizeof(imgpath), "%s_t=%04d.bmp", argv[1], atoi(newtime));
                    img -> WriteToFile(imgpath);
                    strcpy(oldtime,newtime);

                    /*
                    fprintf(fpout,"Time %04d: ",(int)(time));
                    for (j = 0; j < YSIZE; j++){
                        for (i = 0; i < XSIZE; i++) {
                            fprintf(fpout,"%.2f ",Gop[0][i][j]);
                        }
                        fprintf(fpout,"\n");
                        fprintf(fpout,"           ");
                    }
                    strcpy(oldtime,newtime);
                    fprintf(fpout,"\n");
                    fflush(fpout);
                    */
                }
            }
        }

        for (i = 0; i < XSIZE; i++){
            for (j = 0; j < YSIZE; j++) {
                for (k = 0; k < NUMOP; k++) {
                    Gop[k][i][j] = temp[k][i][j];
                }
            }
        }
       
    }

    fclose(fpout);

    exit(0);
}
