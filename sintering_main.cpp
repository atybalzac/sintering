#include "sintering.hpp"

#include <iostream>
#include <limits.h>
#include <stdio.h>

#include "EasyBMP.h"

int main(int argc, char **argv){
    
    register int i, j, k, m;
    double endtime = 900.0;
    double dcdt,conc,time,prefactor,sum;
    double temp[NUMOP][XSIZE][YSIZE],chempot[NUMOP][XSIZE][YSIZE];
    double bulkderiv[NUMOP][XSIZE][YSIZE];
    double maxop;
    double blue,green,red;
    char oldtime[128],newtime[128];
    struct Vec term[NUMOP][XSIZE][YSIZE];
    struct Vec vecterm[NUMOP][XSIZE][YSIZE];
    

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

    sprintf(newtime,"%d",(int)(10.0*time));
    /*
    sprintf(imgpath,"%s.txt", argv[1]);
    if ((fpout = fopen(imgpath,"w")) == NULL) {
        printf("\nERROR: Could not open %s",imgpath);
        exit(1);
    }
    */

    for (time = 0; time <= endtime; time += Gdt){
    // resetting values of volume, force, particle center, torque, velocity term
        initVolume();
        initForce();
        initCenter();
        initTorque();
        for(i = 0; i < XSIZE; i++){
            for(j = 0; j < YSIZE; j++){
                for(k = 0; k < NUMOP; k++){
                    vecterm[k][i][j].x = 0;
                    vecterm[k][i][j].y = 0;
                    vecterm[k][i][j].z = 0;
                }
            }
        }
    // calculates volume of each particle
        for (i = 0; i < XSIZE; i++){
            for(j = 0; j < YSIZE; j++){
                for (k = 1; k < NUMOP; k++){
                    grain_volume(Gop, k, i, j);
                }
            }
        }
    // Calculates center of each particle
        for(i = 0; i<XSIZE; i++){
            for(j = 0; j< YSIZE; j++){
                for (k = 1; k < NUMOP; k++){
                    grain_center(Gop, k, i, j);
                }
            }
        }


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

        //df/dc - eps^2 del^2 c and df/dn - eps^2 del^2 n

        // Now the inner scalar term is calculated everywhere, we next need
        // to create the inner vector field

        for (i = 0; i < XSIZE; i++){
            for (j = 0; j < YSIZE; j++){
                Gfield[0][i][j] = grad(chempot,0,i,j);
                Gfield[0][i][j].x = mob(i,j) * Gfield[0][i][j].x;
                Gfield[0][i][j].y = mob(i,j) * Gfield[0][i][j].y;
            }
        }

  //      std::cout << "Done!" << std::endl;
  //      std::cout << "Updating order parameters... ";
  //      std::cout.flush();

  //      grad(df/dc - eps^2 del^2 c)

        // Now all we need to do is take the divergence of the vector field Gfield,
        // already calculated at each point above, and then multiply by the mobility

        for (i = 0; i < XSIZE; i++){
            for (j = 0; j < YSIZE; j++){
                for(k = 1; k < NUMOP; k++){
                    vecterm[0][i][j].x += velocity(Gop, k, i, j).x;
                    vecterm[0][i][j].y += velocity(Gop, k, i, j).y;
                }
                vecterm[0][i][j].x = Gfield[0][i][j].x - (vecterm[0][i][j].x * Gop[0][i][j]);
                vecterm[0][i][j].y = Gfield[0][i][j].y - (vecterm[0][i][j].y * Gop[0][i][j]);
                dcdt = div(vecterm,0,i,j);
                temp[0][i][j] = Gop[0][i][j] + (dcdt * Gdt);
                for (k = 1; k < NUMOP; k++) {
                    vecterm[k][i][j].x = Gop[k][i][j] * velocity(Gop, k, i, j).x; // velocity term for non conserved op
                    vecterm[k][i][j].y = Gop[k][i][j] * velocity(Gop, k, i, j).y;
                    
                    temp[k][i][j] = Gop[k][i][j] - (GMnc * chempot[k][i][j] * Gdt) - div(vecterm, k, i, j);
                }
            }
        }

        // The new concentrations are stored in the temp array. Transfer them to the
        // concentration array.

  //      std::cout << "Done!" << std::endl;
  //      std::cout << "Updating order parameters... ";
  //      std::cout.flush();

        sprintf(newtime,"%d",(int)(10.0*time));
        if (argc > 1 && (strcmp(oldtime,newtime) || time < Gdt)) { 
            sum = 0.0;
            for (j = 0; j < YSIZE; ++j) {
                for (i = 0; i < XSIZE; ++i) {
                    sum  += Gop[0][i][j];
                    maxop= -5000.0;
                    for (k = 1; k < NUMOP; k++) {
                        if (Gop[k][i][j] > maxop) maxop = Gop[k][i][j];
                    }
                    maxop = pow(maxop*maxop,0.5);


                    if (maxop > 1.0) maxop = 1.0;

                    blue = 1.0-Gop[0][i][j];
                    red = green = maxop;

                    //set pixel colors
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
            std::cout << "Tot conc = " << sum << std::endl;
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
