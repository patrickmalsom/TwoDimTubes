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

//Potential Function
#define VFunc(x,y)      1.-2.*gsl_pow_int(x,2)+gsl_pow_int(x,4)+gsl_pow_int(y,2)

//potentials for calculation of<I-Iou>
#define FxFunc(x,y)      4.*x-4*gsl_pow_int(x,3)
#define ddVxFunc(x,y)     -4.+12*gsl_pow_int(x,2)
#define FyFunc(x,y)     -2*y
#define ddVyFunc(x,y)      2

//Smooth functions for use in Tubes HMC
#define meanxFunc(t)   (7.5-1.5*t+0.967*exp(MU3*gsl_pow_int(-5.+t,2))*tanh(MU1*(-5+t)))/exp(1.*MU3*gsl_pow_int(-5.+t,2))
#define meanyFunc(t)   0.
#define dmeanxFunc(t)   (-1.5+MU3*gsl_pow_int(8.660254037844386-1.7320508075688772*t,2)+0.967*exp(MU3*gsl_pow_int(-5.+t,2))*MU1*gsl_pow_int(1/cosh(MU1*(-5+t)),2))/exp(1.*MU3*gsl_pow_int(-5.+t,2))
#define dmeanyFunc(t)      0
#define ddmeanxFunc(t)     (MU3*(-45.+MU3*gsl_pow_int(9.085602964160698-1.8171205928321397*t,3)+9.*t)-1.934*exp(MU3*gsl_pow_int(-5.+t,2))*gsl_pow_int(MU1,2)*gsl_pow_int(1/cosh(MU1*(-5+t)),2)*tanh(MU1*(-5+t)))/exp(1.*MU3*gsl_pow_int(-5.+t,2))
#define ddmeanyFunc(t)     0

//B parameters
#define BxxFunc(t)   gsl_pow_int((-1+exp(SIGMA3*(-5.+t)))*SIGMA1+SIGMA2,2)/exp(2.*SIGMA3*(-5.+t))
#define ByyFunc(t)   4
#define BxyFunc(t)   0

//derivatives wrt the parameters
#define meanxDParx1Func(t)     (-4.835+0.967*t)*gsl_pow_int(1/cosh(MU1*(-5+t)),2)
#define meanxDParx2Func(t)     0
#define BxxDParx1Func(t)        (2*(-1+exp(SIGMA3*(-5.+t)))*((-1+exp(SIGMA3*(-5.+t)))*SIGMA1+SIGMA2))/exp(2.*SIGMA3*(-5.+t))
#define BxxDParx2Func(t)      (2*((-1+exp(SIGMA3*(-5.+t)))*SIGMA1+SIGMA2))/exp(2.*SIGMA3*(-5.+t))
#define BxxDParx3Func(t)      (gsl_pow_int(SIGMA2,2)*(10.-2.*t)+SIGMA1*SIGMA2*(-20.+exp(SIGMA3*(-5.+t))*(10.-2.*t)+4.*t)+gsl_pow_int(SIGMA1,2)*(10.-2.*t+exp(SIGMA3*(-5.+t))*(-10.+2.*t)))/exp(2.*SIGMA3*(-5.+t))
