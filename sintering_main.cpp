
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
    char oldtime[128],newtime[128];
    struct Vec term[NUMOP][XSIZE][YSIZE];

    //seed random number generator
    FILE *seed = fopen("/dev/urandom", "r");
    unsigned int random_data;
    size_t seedsize = fread(&random_data, sizeof(random_data), 1, seed);
    fclose(seed);

    srand(random_data);

    // initialize microstructure

    init_mic();


    //Generate a blank image that can be updated and output as an image as the program runs
    BMP* img = new BMP();
    img -> SetSize(XSIZE,YSIZE);
    //Declare a path length for naming the image files
    char imgpath[PATH_MAX];

    //Advance through time
    sprintf(oldtime,"0");
    for (time = 0; time <= endtime; time += Gdt){

        for (i = 0; i < XSIZE; i++){
            for (j = 0; j < YSIZE; j++){
                for (k = 1; k <= NUMOP; k++) {
                    bulkderiv[k][i][j]  = calcBulkDeriv(Gop,k,i,j);
                    chempot[k][i][j] = bulkderiv[k][i][j] - (Geps2[k] * laplac(Gop,k,i,j));
                }
            }
        }


        // Now the inner scalar term is calculated everywhere, we next need
        // to create the inner vector field, just the conserved parameter

        for (i = 0; i < XSIZE; i++){
            for (j = 0; j < YSIZE; j++){
                Gfield[0][i][j] = grad(chempot,0,i,j);
            }
        }



        // Now all we need to do is take the divergence of the vector field Gfield,
        // already calculated at each point above, and then multiply by the mobility, just the nonconserved parameters
       
        for (i = 0; i < XSIZE; i++){
            for (j = 0; j < YSIZE; j++){
                dcdt = GMc * div(Gfield,0,i,j);
                temp[0][i][j] = Gop[0][i][j] + (dcdt * Gdt);
                for (k = 1; k < NUMOP; k++) {
                    temp[k][i][j] = Gop[k][i][j] - ((GMnc * chempot[k][i][j]) * Gdt);
                }
            }
        }

        // The new concentrations are stored in the temp array. Transfer them to the
        // concentration array.


        for (i = 0; i < XSIZE; i++){
            for (j = 0; j < YSIZE; j++) {
                maxop= -5000.0;
                for (k = 1; k < NUMOP; k++) {
                    if (Gop[k][i][j] > maxop) maxop = Gop[k][i][j];
                }
                maxop = maxop*maxop;
                if (maxop > 1.0) maxop = 1.0;

                //set pixel colors
                img -> SetPixel(i, j, (RGBApixel){
                    .Blue = (ebmpBYTE)(255*(Gop[0][i][j] * Gop[0][i][j] - (maxop))),
                    .Green =(ebmpBYTE)(255*maxop),
                    .Red = 0,//(ebmpBYTE)(255*maxop),
                    .Alpha = 0,
                });

                for (k = 0; k < NUMOP; k++) {
                    Gop[k][i][j] = temp[k][i][j];
                }
            }
        }
       
        //write image to file iterating name based on timestep. Numbers are padded with leading zeroes so they sort correctly
        sprintf(newtime,"%d",(int)(time));
        if (argc > 1 && (strcmp(oldtime,newtime) || time < Gdt)) { 
            snprintf(imgpath, sizeof(imgpath), "%s_t=%04d.bmp", argv[1], atoi(newtime));
             img -> WriteToFile(imgpath);
            strcpy(oldtime,newtime);
        }
    }

    exit(0);

}
