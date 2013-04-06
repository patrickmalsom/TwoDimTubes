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

// Initial Parameters for the means. To be optimized using KL gradient descent.
double MU1 = 2.00; // steepness of the transition
double MU2 = 0.97; // well location at finite temp

double SIGMA1 = 8.28537; // approx well hessian
double SIGMA2 = -0.90;
double SIGMA3 = 2.72;
//Potential Function
#define VFunc(x,y)      4.+gsl_pow_int(x,4)+0.5*gsl_pow_int(y,2)+gsl_pow_int(x,2)*(-4.+0.5*gsl_pow_int(y,2))

//potentials for calculation of<I-Iou>
#define FxFunc(x,y)      x*(8.-4.*gsl_pow_int(x,2)-1.*gsl_pow_int(y,2))
#define ddVxFunc(x,y)     -8.+12*gsl_pow_int(x,2)+1.*gsl_pow_int(y,2)
#define FyFunc(x,y)     (-1.-1.*gsl_pow_int(x,2))*y
#define ddVyFunc(x,y)      1.+1.*gsl_pow_int(x,2)

//Smooth functions for use in Tubes HMC
#define meanxFunc(t)   MU2*1/tanh(5.*MU1)*tanh(MU1*(-5.+t))
#define meanyFunc(t)   0.
#define dmeanxFunc(t)   MU1*MU2*1/tanh(5.*MU1)*gsl_pow_int(1/cosh(MU1*(-5.+t)),2)
#define dmeanyFunc(t)      0
#define ddmeanxFunc(t)     -2*gsl_pow_int(MU1,2)*MU2*1/tanh(5.*MU1)*gsl_pow_int(1/cosh(MU1*(-5.+t)),2)*tanh(MU1*(-5.+t))
#define ddmeanyFunc(t)     0

//B parameters
#define BxxFunc(t)   gsl_pow_int((-1+exp(SIGMA3*gsl_pow_int(-5.+t,2)))*SIGMA1+SIGMA2,2)/exp(2.*SIGMA3*gsl_pow_int(-5.+t,2))
#define ByyFunc(t)   4.
#define BxyFunc(t)   0.

//derivatives wrt the parameters
#define meanxDParx1Func(t)     MU2*((-5.+t)*1/tanh(5.*MU1)*gsl_pow_int(1/cosh(MU1*(-5.+t)),2)-5.*gsl_pow_int(1/sinh(5.*MU1),2)*tanh(MU1*(-5.+t)))
#define meanxDParx2Func(t)     1/tanh(5.*MU1)*tanh(MU1*(-5.+t))
#define BxxDParx1Func(t)        (2*(-1+exp(SIGMA3*gsl_pow_int(-5.+t,2)))*((-1+exp(SIGMA3*gsl_pow_int(-5.+t,2)))*SIGMA1+SIGMA2))/exp(2.*SIGMA3*gsl_pow_int(-5.+t,2))
#define BxxDParx2Func(t)      (2*((-1+exp(SIGMA3*gsl_pow_int(-5.+t,2)))*SIGMA1+SIGMA2))/exp(2.*SIGMA3*gsl_pow_int(-5.+t,2))
#define BxxDParx3Func(t)      (SIGMA1*SIGMA2*(gsl_pow_int(10.-2.*t,2)-1.*exp(SIGMA3*gsl_pow_int(-5.+t,2))*gsl_pow_int(7.0710678118654755-1.4142135623730951*t,2))+(-1.+exp(SIGMA3*gsl_pow_int(-5.+t,2)))*gsl_pow_int(SIGMA1,2)*gsl_pow_int(7.0710678118654755-1.4142135623730951*t,2)-1.*gsl_pow_int(SIGMA2,2)*gsl_pow_int(7.0710678118654755-1.4142135623730951*t,2))/exp(2.*SIGMA3*gsl_pow_int(-5.+t,2))
