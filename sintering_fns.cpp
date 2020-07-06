#include <iostream>
#include <math.h>
#include "sintering.hpp"

// Global Variables

const double GMnc = 0.6;
const double GMc = 0.6;
const double Gh = 1;
const double Gdt = .001;
const double GFluct = 0.02;
const int Gnumop = (int)(NUMOP);
const int Gnumcop = 1;

double Gop[NUMOP][XSIZE][YSIZE];          /* concentration field */
struct Vec Gfield[NUMOP][XSIZE][YSIZE];  /* chemical potential field */
double Geps2[NUMOP];

// Function definitions

//Checks whether we're at the boundary and loops back to the other side if we are for periodic boundary
int checkbc(int i, int lambda){

    if (i > (lambda - 1)){
        return i-lambda;
    }
    else if (i < 0){
        return i+lambda;
    }
    else {
        return i;
    }
}

//calculates the gradient of c (a vector)
struct Vec grad(double func[NUMOP][XSIZE][YSIZE], int op, int i, int j) {
    struct Vec gradient;

    gradient.x = 0.0;
    gradient.y = 0.0;

    gradient.x = (func[op][checkbc(i+1, XSIZE)][j] - func[op][checkbc(i-1,XSIZE)][j])
        / (2.0 * Gh);
    gradient.y = (func[op][i][checkbc(j+1,YSIZE)] - func[op][i][checkbc(j-1,YSIZE)])
        / (2.0 * Gh);
    return (gradient);
}

// calculates laplacian of c

double laplac(double func[NUMOP][XSIZE][YSIZE], int op, int i, int j) {
    double xderiv = func[op][checkbc(i+1,XSIZE)][j]
        - 2.0 * func[op][i][j]
        + func[op][checkbc(i-1,XSIZE)][j];
    double yderiv = func[op][i][checkbc(j+1,YSIZE)]
        - 2.0 * func[op][i][j]
        + func[op][i][checkbc(j-1,YSIZE)];
    return (xderiv+yderiv)/(Gh*Gh);
}

//calculates the divergence of a vector field at a point

double div(struct Vec v[NUMOP][XSIZE][YSIZE], int op, int i, int j) {

    double diverg = 0.0;

    diverg = (v[op][checkbc(i+1, XSIZE)][j].x - v[op][checkbc(i-1,XSIZE)][j].x) / (2.0 * Gh);
    diverg += (v[op][i][checkbc(j+1,YSIZE)].y - v[op][i][checkbc(j-1,YSIZE)].y) / (2.0 * Gh);
    return (diverg);
}

// initializes the microstructure

void init_mic(void) {

    register int i,j,k;
    float xc,yc,radius,dist,xdist2,ydist2;
    int irad,nxc,nyc,opnum;

    // Start with a blank slate of order parameters
   
    for (k = 0; k < NUMOP; ++k) {
        for (j = 0; j < YSIZE; ++j) {
            for (i = 0; i < XSIZE; ++i) {
                if (k == 0) {
                    Gop[k][i][j] = 0.99;
                } else {
                    Gop[k][i][j] = 0.0;
                }
                Gfield[k][i][j].x = 0.0;
                Gfield[k][i][j].y = 0.0;
            }
        }
    }

    // Set interface energies
    
    Geps2[0] = 0.5;
    for (k = 1; k < NUMOP; ++k) {
        Geps2[k] = 0.2;//0.5;
    }

    // Set center of first circle
   
    xc = 50.0;
    nxc = (int)(xc);
    yc = 50.0;
    nyc = (int)(yc);


    // Set radius of first circle
    
    radius = 7.5;
    irad = (int)(radius);
    opnum = 1;

    // Place first circle
    
    for (j = nyc-irad; j <= nyc+irad; ++j) {
        ydist2 = (float)(j - yc) * (float)(j - yc);
        for (i = nxc-irad; i <= nxc+irad; ++i) {
            xdist2 = (float)(i - xc) * (float)(i - xc);
            dist = sqrt(xdist2 + ydist2);
            if ((dist - 0.5) <= radius) {
                Gop[opnum][i][j] = 1.0; // Belongs to particle opnum
                Gop[0][i][j] = 0.01;     // Is solid
            } else {
                Gop[opnum][i][j] = 0.0; // Does not belong to particle opnum
                Gop[0][i][j] = 0.99;     // Is vapor
            }
        }
    }
            
    // Set center of second circle
   
    xc = 35.0;
    nxc = (int)(xc);
    yc = 50.0;
    nyc = (int)(yc);


    // Set radius of second circle
    
    radius = 7.5;
    irad = (int)(radius);
    opnum = 2;

    // Place second circle
    
    for (j = nyc-irad; j <= nyc+irad; ++j) {
        ydist2 = (float)(j - yc) * (float)(j - yc);
        for (i = nxc-irad; i <= nxc+irad; ++i) {
            xdist2 = (float)(i - xc) * (float)(i - xc);
            dist = sqrt(xdist2 + ydist2);
            if ((dist - 0.5) <= radius) {
                if (Gop[1][i][j] > 0.0) {
                    Gop[1][i][j] = 0.5;
                    Gop[opnum][i][j] = 0.5;
                } else {
                    Gop[opnum][i][j] = 1.0; // Belongs to particle opnum
                }
                Gop[0][i][j] = 0.01;     // Is solid
            } else {
                Gop[opnum][i][j] = 0.0; // Does not belong to particle opnum
                Gop[0][i][j] = 0.99;     // Is vapor
            }
        }
    }
            
    // Set center of third circle
   
    xc = 67.0;
    nxc = (int)(xc);
    yc = 50.0;
    nyc = (int)(yc);


    // Set radius of third circle
    
    radius = 9.5;
    irad = (int)(radius);
    opnum = 3;

    // Place third circle
    
    for (j = nyc-irad; j <= nyc+irad; ++j) {
        ydist2 = (float)(j - yc) * (float)(j - yc);
        for (i = nxc-irad; i <= nxc+irad; ++i) {
            xdist2 = (float)(i - xc) * (float)(i - xc);
            dist = sqrt(xdist2 + ydist2);
            if ((dist - 0.5) <= radius) {
                if (Gop[1][i][j] > 0.0) {
                    Gop[1][i][j] = 0.5;
                    Gop[opnum][i][j] = 0.5;
                } else {
                    Gop[opnum][i][j] = 1.0; // Belongs to particle opnum
                }
                Gop[0][i][j] = 0.01;     // Is solid
            } else {
                Gop[opnum][i][j] = 0.0; // Does not belong to particle opnum
                Gop[0][i][j] = 0.99;     // Is vapor
            }
        }
    }
            
    return;
}

// calculates the derivative of the bulk free energy w/r/t an order parameter

double calcBulkDeriv(double func[NUMOP][XSIZE][YSIZE], int op, int i, int j) {

    register int m, p, q;
    double bd = 0.0;
    double c;
    double term = 0.0;
    double term1 = 0.0;
    double term2 = 0.0;
    double a = 1.0;
    double b = 1.0;
    double dalpha = 1.0;
    double dbeta = 1.0;
    double gamma = 1.0;
    double delta = 1.0;
    double calpha = 0.99;
    double cbeta = 0.01;
    double cmean = (calpha+cbeta)/2;
    double epsilon[Gnumop-Gnumcop][Gnumop-Gnumcop];
    
    //set epsilon
    for (p=0; p < (Gnumop-Gnumcop); p++){
        for (q=0; q < (Gnumop-Gnumcop); q++){
            epsilon[p][q]=1;
        }
    }

    c = func[0][i][j];  // volume fraction of vacancies

    for (m = 1; m < NUMOP; m++) {
        term += pow(func[m][i][j],2.0); //works for conserved op
    }
    for (m = 1; m < NUMOP; m++){ //for nonconserved ops
        term1 += 2*epsilon[op][m] * pow(func[m][i][j],2.0);
    }
    term2 = term1 - pow(func[op][i][j],2.0);
    switch (op) {
        //case for conserved order parameter 
        case 0:
            bd = b*pow((c-cmean),3.0) - a*(c-cmean) + dalpha * pow((c-calpha),3.0) + dbeta*pow((c-cbeta),3.0);
            bd += -((c-calpha) * gamma * term);
            break;
        //case for non conserved order parameter
        default:
            bd = delta*pow(func[op][i][j],3.0) - (pow((c-calpha),2.0)*func[op][i][j]*gamma)
                 + (func[op][i][j]*term2);
            break;
    }

    return (bd);
}

