#include <iostream>
#include <math.h>
#include "sintering.hpp"

// Global Variables

const double GA = 11.0;
double Gcm = 1.0;
double GB = 1.0;
double Ggam = 1.0; 
double GDa = 1.0;
double GDb = 1.0;

const double GMnc = 5.0;
const double GMc = 5.0;
const double Gh = 0.5;
const double Gdt = .0001;
const double GFluct = 0.02;
const int Gnumop = (int)(NUMOP);
const int Gnumcop = 1;

double Gop[NUMOP][XSIZE][YSIZE];          /* concentration field */
struct Vec Gfield[NUMOP][XSIZE][YSIZE];  /* chemical potential field */
double Geps2[NUMOP];
double Gvolume[NUMOP];
struct Vec Gforce[NUMOP];
struct Vec Gcenter[NUMOP];
struct Vec Gtorque[NUMOP];
// Function definitions

//initializes Gvolume, Gforce, Gcenter, Gtorque to 0
void initVolume(void){
    for (int i = 1; i < NUMOP; i++){
        Gvolume[i]=0;
    }
    return;
}

void initForce(void){
    for (int i = 0; i < NUMOP; i++){
        Gforce[i].x = 0;
        Gforce[i].y = 0;
        Gforce[i].z = 0;
    }
    return;
}

void initCenter(void){
    for (int i = 0; i < NUMOP; i++){
        Gcenter[i].x=0;
        Gcenter[i].y=0;
        Gcenter[i].z=0;
    }
    return;
}

void initTorque(void){
    for (int i = 0; i < NUMOP; i++){
        Gtorque[i].x=0;
        Gtorque[i].y=0;
        Gtorque[i].z=0;
    }
    return;
}

//Checks whether we're at the boundary and loops back to the other side if we have periodic boundaries
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
    gradient.z = 0.0;

    gradient.x = (func[op][checkbc(i+1, XSIZE)][j] - func[op][checkbc(i-1,XSIZE)][j])
        / (2.0 * Gh);
    gradient.y = (func[op][i][checkbc(j+1,YSIZE)] - func[op][i][checkbc(j-1,YSIZE)])
        / (2.0 * Gh);
    return (gradient);
}

// Calculate the mobility at a point to simulate surface diffusion

double mob(int i, int j) {
    double term = GMc * (1.0 + pow(Gop[0][i][j],2.0));
    return (term);
}

// calculates laplacian of c

double laplac(double func[NUMOP][XSIZE][YSIZE], int op, int i, int j) {
    double term = 0.5 * func[op][checkbc(i+1,XSIZE)][j];
    term += (0.5 * func[op][checkbc(i-1,XSIZE)][j]);
    term += (0.5 * func[op][i][checkbc(j+1,YSIZE)]);
    term += (0.5 * func[op][i][checkbc(j-1,YSIZE)]);
    term += (0.25 * func[op][checkbc(i+1,XSIZE)][checkbc(j+1,YSIZE)]);
    term += (0.25 * func[op][checkbc(i+1,XSIZE)][checkbc(j-1,YSIZE)]);
    term += (0.25 * func[op][checkbc(i-1,XSIZE)][checkbc(j+1,YSIZE)]);
    term += (0.25 * func[op][checkbc(i-1,XSIZE)][checkbc(j-1,YSIZE)]);
    term -= (3.0 * func[op][i][j]);

    return (term)/(Gh*Gh);
}

//calculates the divergence of a vector field at a point

double div(struct Vec v[NUMOP][XSIZE][YSIZE], int op, int i, int j) {

    double diverg = 0.0;

    diverg = (v[op][checkbc(i+1, XSIZE)][j].x - v[op][checkbc(i-1,XSIZE)][j].x) / (2.0 * Gh);
    diverg += (v[op][i][checkbc(j+1,YSIZE)].y - v[op][i][checkbc(j-1,YSIZE)].y) / (2.0 * Gh);
    return (diverg);
}


//calculates force density

struct Vec forcedensity(double func[NUMOP][XSIZE][YSIZE], int op1, int op2, int i, int j){
    double c = func[0][i][j];  // volume fraction of vacancies
    double p1 = func[op1][i][j];  // non-conserved order parameter
    double p2 = func[op2][i][j];
    double cgb = 0.1;
    double c0 = 1; // the equilibrium density?
    double term = (c-c0);
    double k = 1; //stiffness constant
    struct Vec density;
    density.x = 0;
    density.y = 0;
    density.z = 0;
    struct Vec gradterm;
    gradterm.x = 0;
    gradterm.y = 0;
    gradterm.z = 0;
    if(op1 > 0 && op2 > 0 && (p1*p2) >= cgb){
        gradterm.x += grad(func, op1, i, j).x-grad(func,op2,i,j).x;
        gradterm.y += grad(func, op1, i, j).y-grad(func,op2,i,j).y;
    }
    density.x += k*term*gradterm.x;
    density.y += k*term*gradterm.y;
    return density;
}

//Calculates force acting on particle
void force(double func[NUMOP][XSIZE][YSIZE], int op1, int op2, int i, int j){
    Gforce[op1].x += forcedensity(func, op1, op2, i, j).x;
    Gforce[op1].y += forcedensity(func, op1, op2, i, j).y;
    return;
}

//calculates volume of particle
void grain_volume(double func[NUMOP][XSIZE][YSIZE], int op, int i, int j){
    double p = func[op][i][j];
    Gvolume[op] += p;
    return;
}

//calculates center of mass of particles
void grain_center(double func[NUMOP][XSIZE][YSIZE], int op, int i, int j){
    double p = func[op][i][j];
    Gcenter[op].x += i*p/Gvolume[op];
    Gcenter[op].y += j*p/Gvolume[op];
    return;
}

//calculates torque on particle
void torque(double func[NUMOP][XSIZE][YSIZE], int op1, int op2, int i, int j){
    struct Vec diff;
    diff.x = i - Gcenter[op1].x;
    diff.y = j - Gcenter[op1].y;
    diff.z = 0;

    Gtorque[op1].z += i * forcedensity(func, op1, op2, i, j).y - j * forcedensity(func, op1, op2, i, j).x;
    return;
}

//calculates velocity of particle
struct Vec velocity(double func[NUMOP][XSIZE][YSIZE], int op, int i, int j){
    struct Vec vtrans;
    struct Vec vrot;
    double mtrans = 1;
    double mrot = 1;
    double p = func[op][i][j];
    double term1 = (mtrans/Gvolume[op]) * p;

    vtrans.x = term1 * Gforce[op].x;
    vtrans.y = term1 * Gforce[op].y;
    vtrans.z = 0;

    struct Vec diff;
    diff.x = i - Gcenter[op].x;
    diff.y = j - Gcenter[op].y;
    diff.z = 0;
    
    double term2 = (mrot/Gvolume[op]) * p;
    
    vrot.x = -term2*diff.y*Gtorque[op].z;
    vrot.y = term2*diff.x*Gtorque[op].z;
    vrot.z = 0;

    struct Vec velocity;
    velocity.x = vtrans.x + vrot.x;
    velocity.y = vtrans.y + vrot.y;
    velocity.z = 0;
    
    return velocity;
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
                    Gop[k][i][j] = 0.0;
                } else {
                    Gop[k][i][j] = 0.0;
                }
                Gfield[k][i][j].x = 0.0;
                Gfield[k][i][j].y = 0.0;
            }
        }
    }

    // Set all non-constant parameters
    
    GB = 1.0;

    Geps2[0] = 3.0;
    for (k = 1; k < NUMOP; ++k) {
        Geps2[k] = 1.5;
    }

    // Set center of first circle
   
    xc = 60.0;
    nxc = (int)(xc);
    yc = 30.0;
    nyc = (int)(yc);

    // Set radius of first circle
    
    radius = 11.0;
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
                Gop[0][i][j] = 1.0;     // Is solid
            } else {
                Gop[opnum][i][j] = 0.0; // Does not belong to particle opnum
                Gop[0][i][j] = 0.0;     // Is vapor
            }
        }
    }
            
    // Set center of second circle
   
    xc = 27.0;
    nxc = (int)(xc);
    yc = 30.0;
    nyc = (int)(yc);


    // Set radius of second circle
    
    radius = 22.0;
    irad = (int)(radius);
    opnum = 2;

    /*
    xc = 7.0;
    nxc = (int)(xc);
    yc = 10.0;
    nyc = (int)(yc);


    // Set radius of second circle
    
    radius = 3.0;
    irad = (int)(radius);
    opnum = 2;
    */

    // Place second circle
    
    for (j = nyc-irad; j <= nyc+irad; ++j) {
        ydist2 = (float)(j - yc) * (float)(j - yc);
        for (i = nxc-irad; i <= nxc+irad; ++i) {
            xdist2 = (float)(i - xc) * (float)(i - xc);
            dist = sqrt(xdist2 + ydist2);
            if ((dist - 0.5) <= radius) {
                Gop[opnum][i][j] = 1.0; // Belongs to particle opnum
                Gop[0][i][j] = 1.0;     // Is solid
            } else {
                Gop[opnum][i][j] = 0.0; // Does not belong to particle opnum
                Gop[0][i][j] = 0.0;     // Is vapor
            }
        }
    }
            
    // Set center of third circle
   
    xc = 93.0;
    nxc = (int)(xc);
    yc = 30.0;
    nyc = (int)(yc);


    // Set radius of third circle
    
    radius = 22.0;
    irad = (int)(radius);
    opnum = 3;

    // Place third circle
    
    for (j = nyc-irad; j <= nyc+irad; ++j) {
        ydist2 = (float)(j - yc) * (float)(j - yc);
        for (i = nxc-irad; i <= nxc+irad; ++i) {
            xdist2 = (float)(i - xc) * (float)(i - xc);
            dist = sqrt(xdist2 + ydist2);
            if ((dist - 0.5) <= radius) {
                Gop[opnum][i][j] = 1.0; // Belongs to particle opnum
                Gop[0][i][j] = 1.0;     // Is solid
            } else {
                Gop[opnum][i][j] = 0.0; // Does not belong to particle opnum
                Gop[0][i][j] = 0.0;     // Is vapor
            }
        }
    }
            
    return;
}

// calculates the derivative of the bulk free energy w/r/t an order parameter

double calcBulkDeriv(double func[NUMOP][XSIZE][YSIZE], int op, int i, int j) {

    register int m;
    double bd = 0.0;
    double c,p;
    double term2 = 0.0;
    double term3 = 0.0;

    c = func[0][i][j];  // volume fraction of vacancies
    p = func[op][i][j];  // non-conserved order parameter

    for (m = 1; m < NUMOP; m++) {
        term2 += pow(func[m][i][j],2.0);
        term3 += pow(func[m][i][j],3.0);
    }

    switch (op) {
        case 0:
            bd = 2.0 * GA * c * (1.0 - c) * (1.0 - c);
            bd -= (2.0 * GA * c * c * (1.0 - c));
            bd += (2.0 * GB * c);
            bd -= (6.0 * GB * term2);
            bd += (4.0 * GB * term3);
            break;

        default:
            bd = 12.0 * GB * (1.0 - c) * p;
            bd -= (12.0 * GB * (2.0 - c) * p * p);
            bd += (12.0 * GB * term2 * p);
            break;
    }

    return (bd);
}

