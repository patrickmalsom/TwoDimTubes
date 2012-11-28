// =========================================
// Preprocessor Definitions
// =========================================

// Path Definitions
#define NUMDIM    2           // Spacial Dimensions
#define NUMBEAD   20001       // Path Points
#define DU        0.0005      // du (see potentials.h for estimate of DU)
#define PREDT     200.0         // dt=PreDT*du^2 (path time)

// Temperature Definition
#define TEMP      0.15

// Incrimenter Definitions
#define NUMMD     30        // Number of MD steps 
//      NUMMD     ~3/(2*sqrt(2*PreDT*DU^2)) <- Approx optimal value of NUMMD
#define NUMMC     100      // Number of Metropolis Hastings MC steps
#define NUMTUBE   10       // Number of tube steepest descent steps

// Constants for writing to stdout and config
#define WRITESTDOUT  30       // How often to print to stdout (# of MD loops)
//current implimentation writes a file at every steepest descent step
const char PotentialString[]="2WellSymTubes";// Potential Description 

// Potential Function
#define VFunc(x,y)           1.0l-2.0l*gsl_pow_int(x,2)+1.0l*gsl_pow_int(x,4)+1.0l*gsl_pow_int(y,2)

//=============================================================================
//                   Tubes Definitions 
//=============================================================================

// Initial Parameters for the means
//these parameters will be optimized using KL gradient descent
extern double MPARX1 = 0.97;
extern double MPARX2 = 2.00;

extern double BPARX1 = 7.5;
extern double BPARX2 = 0.42;
extern double BPARX3 = 5.50;

// Smooth functions for use in Tubes HMC
#define meanxFunc(t)    MPARX2/tanh(5*MPARX1)*tanh(MPARX1*(-5 + t))
#define meanyFunc(t)    0.0
#define dmeanxFunc(t)   MPARX1*MPARX2/tanh(5*MPARX1)*pow(cosh(MPARX1*(-5 + t)),-2.0)
#define dmeanyFunc(t)   0.0
#define ddmeanxFunc(t)  -2*pow(MPARX1,2.0)*MPARX2/tanh(5*MPARX1)*pow(cosh(MPARX1*(-5+t)),-2.0)*tanh(MPARX1*(-5+t))
#define ddmeanyFunc(t)  0.0
#define BxxFunc(t)      pow(BPARX1+(BPARX2-BPARX1)*exp(-BPARX3*pow(-5.+t,2)),2.0)
#define ByyFunc(t)      4.0
#define BxyFunc(t)      0.0

//derivatives wrt the parameters
#define BxxDParx1Func(t)  -2*(BPARX2 - BPARX3) * (BPARX3 + (BPARX2 - BPARX3)*exp(-BPARX1*pow(-5. + t,2))) *pow(-5. + t,2) * exp(-BPARX1*pow(-5. + t,2))
#define BxxDParx2Func(t)  2*exp(-BPARX1*pow(-5.+t,2))*(BPARX3+(BPARX2-BPARX3)*exp(-BPARX1*pow(-5.+t,2)))
#define BxxDParx3Func(t)  2*(1-exp(-BPARX1*pow(-5.+t,2)))*(BPARX3+(BPARX2-BPARX3)*exp(-BPARX1*pow(-5.+t,2)))
