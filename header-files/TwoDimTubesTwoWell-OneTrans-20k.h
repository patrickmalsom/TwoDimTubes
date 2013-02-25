/*
 *   TwoDimTubesConstants.h
 *   Written Spring 2013 -- Patrick Malsom
 *   Constants for Two Dimensional Tubes HMC
 */

// ==================================================================================
// Preprocessor Definitions
// ==================================================================================
// Path Definitions
#define NUMDIM    2           // Spacial Dimensions
#define NUMBEAD   20001       // Path Points
#define DU        0.0005      // path time step size
#define PREDT     6000.0       // DT=PreDT*DU^2 (path time)

// Temperature Definition
#define TEMP      0.15

// Incrimenter Definitions
#define NUMMD     50        // Number of MD steps 
//      NUMMD     ~3/(2*sqrt(2*PreDT*DU^2)) <- Approx optimal value of NUMMD
#define NUMMC     5000    // Number of Metropolis Hastings MC steps
#define NUMTUBE   100        // Number of tube steepest descent steps

// Constants for writing to stdout and config
#define WRITESTDOUT  50       // How often to print to stdout (# of MD loops)
//current implimentation writes a file at every steepest descent step
const char PotentialString[]="2WellTubes";// Potential Description 

//=============================================================================
//                   Tubes Definitions 
//=============================================================================

/* Mathematica code to generate the definitions (because CForm is broken...)
Cdef[fun__]:=StringReplace[StringReplace[ToString[CForm[fun]],{"Sinh("->"sinh(","Cosh("->"cosh(","Tanh("->"tanh(","Csch("->"1/sinh(","Sech("->"1/cosh(","Coth("->"1/tanh(","Power("->"pow("," "->""}],"pow(E,"->"exp("]
*/

// Potential Function
#define VFunc(x,y)           1.0-2.0*x*x+1.0*x*x*x*x+1.0*y*y
//potentials for calculation of <I-Iou>
#define FxFunc(x,y)            4.0*x-4.0*x*x*x
#define ddVxFunc(x,y)          -4.0+12.0*x*x
#define FyFunc(x,y)            -2.0*y
#define ddVyFunc(x,y)          2.0

// Smooth functions for use in Tubes HMC
#define meanxFunc(t)      MU2*1/tanh(5.*MU1)*tanh(MU1*(-5.+t))
#define meanyFunc(t)      0.0
#define dmeanxFunc(t)     MU1*MU2*1/tanh(5.*MU1)*pow(1/cosh(MU1*(-5.+t)),2)
#define dmeanyFunc(t)     0.0
#define ddmeanxFunc(t)    -2*pow(MU1,2)*MU2*1/tanh(5.*MU1)*pow(1/cosh(MU1*(-5.+t)),2)*tanh(MU1*(-5.+t))
#define ddmeanyFunc(t)    0.0
//#define BxxFunc(t)        7.52*7.52
#define BxxFunc(t)        pow(SIGMA1+(-SIGMA1+SIGMA2)/exp(SIGMA3*pow(-5.+t,2)),2)
#define ByyFunc(t)        4.0
#define BxyFunc(t)        0.0

//derivatives wrt the parameters
#define meanxDParx1Func(t)   MU2*(-5.+t)*1/tanh(5.*MU1)*pow(1/cosh(MU1*(-5.+t)),2)-5.*MU2*pow(1/sinh(5.*MU1),2)*tanh(MU1*(-5.+t))
#define meanxDParx2Func(t)   1/tanh(5.*MU1)*tanh(MU1*(-5.+t))
#define BxxDParx1Func(t)     2*(1-exp(-(SIGMA3*pow(-5.+t,2))))*(SIGMA1+(-SIGMA1+SIGMA2)/exp(SIGMA3*pow(-5.+t,2)))
#define BxxDParx2Func(t)     (2*(SIGMA1+(-SIGMA1+SIGMA2)*exp(-SIGMA3*pow(-5.+t,2))))*exp(-SIGMA3*pow(-5.+t,2))
#define BxxDParx3Func(t)     (-2*(-SIGMA1+SIGMA2)*(SIGMA1+(-SIGMA1+SIGMA2)/exp(SIGMA3*pow(-5.+t,2)))*pow(-5.+t,2))/exp(SIGMA3*pow(-5.+t,2)) 

