#pragma once
#define XSIZE 120
#define YSIZE 60
#define NUMOP 4

struct Vec {
    double x;   // x component of a 2D vector
    double y;   // y component of a 2D vector
    double z;
};

// Global variables

extern const double Gh;
extern const double Gdt;

extern const double GMnc;
extern const double GMc;
extern const double GFluct;
extern const int Gnumcop;
extern const int Gnumop;

extern const double GA;
extern double GB;
extern double Gop[NUMOP][XSIZE][YSIZE];
extern struct Vec Gfield[NUMOP][XSIZE][YSIZE];
extern double Geps2[NUMOP];
extern double Gvolume[NUMOP];
extern struct Vec Gforce[NUMOP];
extern struct Vec Gcenter[NUMOP];
extern struct Vec Gtorque[NUMOP];

// Function declarations
int checkbc(int i, int lambda);
void initTorque(void);
void initCenter(void);
void initForce(void);
void initVolume(void);
struct Vec grad(double func[NUMOP][XSIZE][YSIZE], int op, int i, int j);
double mob(int i, int j);
double laplac(double func[NUMOP][XSIZE][YSIZE], int op, int i, int j);
double div(struct Vec v[NUMOP][XSIZE][YSIZE], int op, int i, int j);
double calcBulkDeriv(double func[NUMOP][XSIZE][YSIZE], int op, int i, int j);
void init_mic(void);
struct Vec forcedensity(double func[NUMOP][XSIZE][YSIZE], int op1, int op2, int i, int j);
void force(double func[NUMOP][XSIZE][YSIZE], int op1, int op2, int i, int j);
void grain_volume(double func[NUMOP][XSIZE][YSIZE], int op, int i, int j);
void grain_center(double func[NUMOP][XSIZE][YSIZE], int op, int i, int j);
void torque(double func[NUMOP][XSIZE][YSIZE], int op1, int op2, int i, int j);
struct Vec velocity(double func[NUMOP][XSIZE][YSIZE], int op, int i, int j);