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

// Initial Parameters for the means. To be optimized using KL gradient descent.
double MU1 = 3.00; // steepness of the transition
double MU2 = 20.0; // well location at finite temp

double SIGMA1 = 8.28537; // approx well hessian
double SIGMA2 = -0.90;
double SIGMA3 = 2.72;

//=============================================================================
//                   Tubes Definitions 
//=============================================================================
