/*
 *   TwoDimTubes.c 
 *   Written Fall 2012 -- Patrick Malsom
 *   Two Dimensional Tubes HMC
 */

// ==================================================================================
// Library Definitions
// ==================================================================================
//STD Libraries
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//OpenMP libraries
#include <omp.h>
//GNU Scientific Libraries
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>

// ==================================================================================
// Preprocessor Definitions
// ==================================================================================
// Path Definitions
#define NUMDIM    2           // Spacial Dimensions
#define NUMBEAD   20001       // Path Points
#define DU        0.0005      // path time step size
#define PREDT     1000.0       // DT=PreDT*DU^2 (path time)

// Temperature Definition
#define TEMP      0.15

// Incrimenter Definitions
#define NUMMD     50        // Number of MD steps 
//      NUMMD     ~3/(2*sqrt(2*PreDT*DU^2)) <- Approx optimal value of NUMMD
#define NUMMC     1000      // Number of Metropolis Hastings MC steps
#define NUMTUBE   100        // Number of tube steepest descent steps

// Constants for writing to stdout and config
#define WRITESTDOUT  50       // How often to print to stdout (# of MD loops)
//current implimentation writes a file at every steepest descent step
const char PotentialString[]="2WellTubes";// Potential Description 

// Potential Function
#define VFunc(x,y)           1.0l-2.0l*gsl_pow_int(x,2)+1.0l*gsl_pow_int(x,4)+1.0l*gsl_pow_int(y,2)

//=============================================================================
//                   Tubes Definitions 
//=============================================================================

// Smooth functions for use in Tubes HMC
#define meanxFunc(t)      MU2/tanh(5*MU1)*tanh(MU1*(-5+t))
#define meanyFunc(t)      0.0
#define dmeanxFunc(t)     MU1*MU2/tanh(5*MU1)*pow(cosh(MU1*(-5+t)),-2.0)
#define dmeanyFunc(t)     0.0
#define ddmeanxFunc(t)    -2*pow(MU1,2.0)*MU2/tanh(5*MU1)*pow(cosh(MU1*(-5+t)),-2.0)*tanh(MU1*(-5+t))
#define ddmeanyFunc(t)    0.0
#define BxxFunc(t)        pow(SIGMA1+(SIGMA2-SIGMA1)*exp(-SIGMA3*pow(-5.+t,2)),2.0)
#define ByyFunc(t)        4.0
#define BxyFunc(t)        0.0

//derivatives wrt the parameters
#define meanxDParx1Func(t)    MU2*(-5.0+t)/tanh(5.0*MU1)*pow(1.0/cosh(MU1*(-5.0+t)),2.0)-5.0*MU2*pow(1.0/sinh(5.0*MU1),2.0)*tanh(MU1*(-5.0+t))
#define meanxDParx2Func(t)    tanh(MU1*(-5.0+t))/tanh(5.0*MU1)
#define BxxDParx1Func(t)  -2*(SIGMA2-SIGMA3)*(SIGMA3+(SIGMA2-SIGMA3)*exp(-SIGMA1*pow(-5.+t,2)))*pow(-5.+t,2)*exp(-SIGMA1*pow(-5.+t,2))
#define BxxDParx2Func(t)  2*exp(-SIGMA1*pow(-5.+t,2))*(SIGMA3+(SIGMA2-SIGMA3)*exp(-SIGMA1*pow(-5.+t,2)))
#define BxxDParx3Func(t)  2*(1-exp(-SIGMA1*pow(-5.+t,2)))*(SIGMA3+(SIGMA2-SIGMA3)*exp(-SIGMA1*pow(-5.+t,2)))

//=============================================================================
//                   Global Variables
//=============================================================================
const int NUMu=NUMBEAD-1;  //num of DU's
const int NUMl=NUMBEAD-2;  //used in Linverse function

const double DT=PREDT*DU*DU;

// Initial Parameters for the means. To be optimized using KL gradient descent.
double MU1 = 2.00; // steepness of the transition
double MU2 = 0.97; // well location at finite temp

double SIGMA1 = 7.52; // approx well hessian
double SIGMA2 = 0.40;
double SIGMA3 = 5.20;

FILE *pStdOut;

// ==================================================================================
// Structure Definitions
// ==================================================================================

//Only stores the positions
typedef struct _position
{
  double pos[NUMDIM];
} position;

//Stores the KL distance expectation values
typedef struct _averages
{
  double mean[4];
  double meanB[10];
} averages;

//Stores positions and all potentials
typedef struct _config
{
  double pos[NUMDIM];
  double posm[NUMDIM];
  double posdm[NUMDIM];
  double posddm[NUMDIM];
  double posz[NUMDIM];
  double B[(NUMDIM*NUMDIM+NUMDIM)/2];
  double Energy;
  double G;
  double gradG[NUMDIM];
  double LinvG[NUMDIM];
} config;

// ==================================================================================
// Function Prototypes
// ==================================================================================

void rotateConfig(config **a, config **b, config **c);
//rotation of pointers.
// b       -> a,    c   -> b,        a   -> c
//(current -> old,  new -> current,  old -> new)
//can now write over new to calculate a new state (discarding old state)

void GenGaussRand(double GaussRand[NUMu], gsl_rng *RanNumPointer, double StdDev);
//generate NUMBEAD of Gaussian random nums; save to GaussRandArray[]

void generateBB(double bb[NUMBEAD], double GaussRandArray[NUMu], gsl_rng *RanNumPointer);
// generates a standard browinan bridge (starts at zero ends at zero)
// does not renormalize the bridge to correct quadratic variation (see renormBB)

void renormBB(double bb[NUMBEAD]);
//Renormalize the brownian bridge (velocities) for the correct thermal fluctuation

void renorm(position *currentpos);
//Renormalize the starting path (positions) for the correct thermal fluctuation

void calcPotentials(config *currentConfig, int beadIndex);
//Takes the config 'currentConfig' and calculates
// G, gradG, LinvG, Energy for currentConfig[beadIndex].
//Note: calculates potentials for a single bead (beadIndex).

void LInverse(config* currentConfig, double vecdg[NUMl], double veci1[NUMl], double veci0[NUMl]);
//L^(-1) y = x. finds the vector y
// function solves for all currentConfig.LinvG of the current config

void preconditionSPDE(config* currentConfig, config* newConfig, double bb[NUMBEAD], double GaussRandArray[NUMu], gsl_rng *RanNumPointer);
//performs the SPDE step. uses current config and a renormalized auxillary velocities to calculate a new config. Requires random numbers from gsl.

void MolecularDynamics(config *oldConfig, config *currentConfig, config *newConfig);
//performs molecular dynamics. uses oldConfig and currentConfig to calculate newConfig

void ProbAccRatio(config *currentConfig, config *newConfig, double *ratio);
// calculate the probability accaptance ratio for the step
//used to determine the acc/rej of the Metropolis-Hastings MC test

// ============================== Utility functions =================================
void saveConfig(config *currentConfig, config *saveConfig);
//copies ALL elements of currentto saveConfig
//leaves currentConfig unmodified.

void savePosition(position *currentpos, position *savepos);
//copies ALL elements of currentpos to savepos
//leaves currentpos unmodified.

void saveConfigtoPos(config *currentConfig, position *savepos);
//copies ONLY pos elements of currentConfig to savepos
//leaves currentConfig unmodified.

void savePostoConfig(position *currentpos, config *saveConfig);
//copies ONLY pos elements of currentpos to saveConfig
//leaves currentpos unmodified.

void writeConfig(config *newConfig, averages *tubeAve, int MCloopi);
//write the config to file with name position******.dat where ****** is MCloopi
//writes the positions and the tube averages
//only writes the iterations in constants.h

void printPositionPos(position* currentpos, int beadIndex);
//Print the positions of the position struct 

void printConfigPos(config* currentConfig, int beadIndex);
//Print the positions of the config struct 

void printConfigPot(config* currentConfig, int beadIndex);
//Print all information of the config struct calculated in calcPotentials
// (Postitions, grad G, Energy, G)

void printConfigAll(config* currentConfig, int beadIndex);
//Print all information of the config struct

void printDistance(config *newConfig, position *savePos);

// ============================== Tube Average functions =================================
void zeroAverages(averages *tubeAve, int *tau);
//zero all elements in the averages struct 

void initMeans(config *currentConfig);
// initialize the means. This is only called once per run

void accumulateAverages(averages *tubeAve, config *newConfig, int *tau);
//sum the averages. This includes the mean (sum of the positions) and the 
//square of the positions. This is enough information to calculate the covariance
//matrix in post processing

void normalizeAverages(averages *tubeAve, int *tau);
//normalize the average struct before writing to file
//simply dividing by the total number of accumulate average calls

void tubesSteepestDescent(averages *tubeAve);

// ==================================================================================
//               MAIN Program
// ==================================================================================
int main(int argc, char *argv[])
{ 
  setbuf(stdout,NULL);  //allows a pipe to from stdout to file

  //file name for writing std out
  pStdOut= fopen("StdOut.dat","w");
  setbuf(pStdOut,NULL);  //write to file immediatley

  //===============================================================
  // Declare variables and print to std output for reference
  //===============================================================
  //define the Config structs. Example: configOld[n].pos[i]
  //where n->Bead i->dimension
  config *configOld = calloc(NUMBEAD,sizeof(config));
  config *configCurrent = calloc(NUMBEAD,sizeof(config));
  config *configNew = calloc(NUMBEAD,sizeof(config));

  //Used to save positions for MHMC rejection
  position *savePos = calloc(NUMBEAD,sizeof(position));

  //averages used for minimization of KL distance
  averages *tubeAve = calloc(NUMBEAD,sizeof(averages));

  //Incrimenter Declarations
  int i;
  int acc,rej; // trackers for acceptance of MHMC loop
  int MDloopi,MCloopi, tubeloopi; 
  int tau=0; //incrimenter used in the averaging for normalization

  //Vectors for doing the L Inverse
  double *vecdg = calloc(NUMl,sizeof(double));;
  double *veci0 = calloc(NUMl,sizeof(double));;
  double *veci1 = calloc(NUMl,sizeof(double));;

  // array to store rand nums in
  double *GaussRandArray = calloc(NUMu,sizeof(double));

  // error sum for the MH-MC test
  double ratio;

  // random number sorage for the MH-MC test
  double randUniform;

  // storage for brownian bridge
  double *bb = calloc(NUMu,sizeof(double));

  //Print parameters for the run in stdout
  printf("=======================================================\n");
  printf("HMC method for 2D potentials \n");
  printf("=======================================================\n");
  printf("TEMPERATURE = %f \n",TEMP);
  printf("=======================================================\n");
  printf("Number of Metropolis Hastings steps: %i\n",NUMMC);
  printf("Number of MD steps: %i \n",NUMMD);
  printf("=======================================================\n");
  printf("Number of Dimensions: %i \n",NUMDIM);
  printf("Number of Beads: %i \n",NUMBEAD);
  printf("Path grid: du = %+.8e \n", DU);
  printf("Sampling Parameters: dt=%f \n",DT);
  printf("=======================================================\n");
  printf("MD step: h=%+.8e \n",sqrt(2.0l * DT));
  printf("MD time (n*h): %+.8e \n",NUMMD*sqrt(2.0l * DT));
  printf("=======================================================\n");

  //===============================================================
  // Reading the input configuration file into savepos
  //===============================================================
  //Input file to be read as first command line argument
  if(argv[1]==NULL) { 
    printf("No input file. Exiting!\n");
    exit(-1);
  }
  else {
    printf("Input Configuratrion File: %s\n",argv[1]);
  }
  int lineNum = 0;
  FILE *fptr = fopen(argv[1],"r");
  switch(NUMDIM){
    case 3:  //For 3 Dimensions
      while( EOF != fscanf(fptr,"%lf %lf %lf",
      &(savePos[lineNum].pos[0]),
      &(savePos[lineNum].pos[1]),
      &(savePos[lineNum].pos[2])) ) {
        lineNum++;
      }
      break;
    case 2:   //For 2 Dimensions 
      while( EOF != fscanf(fptr,"%lf %lf",
      &(savePos[lineNum].pos[0]),
      &(savePos[lineNum].pos[1])) ) {
        lineNum++;
      }
      break;
    case 1:  //For 1 Dimension
      while( EOF != fscanf(fptr,"%lf",
      &(savePos[lineNum].pos[0])) ) {
        lineNum++;
      }
      break;
    default:
      printf("ERROR: NUMDIM incorrectly defined. Exiting!\n");
      exit(-1);
  }

  //===============================================================
  // GNU Scientific Library Random Number Setup
  //===============================================================
  // Example shell command$ GSL_RNG_SEED=123 ./a.out
  printf("=======================================================\n");
  const gsl_rng_type * RanNumType;
  gsl_rng *RanNumPointer; 
  gsl_rng_env_setup();
  RanNumType = gsl_rng_default;
  RanNumPointer= gsl_rng_alloc (RanNumType);
  printf("Random Number Generator Type: %s \n", gsl_rng_name(RanNumPointer));
  printf("RNG Seed: %li \n", gsl_rng_default_seed);
  printf("=======================================================\n");

  renorm(savePos); //renormalize the position to have the correct thermal fluctuation

  // ==================================================================================
  //     Steepest Descent Loop for optimizing Tube Parameters
  // ==================================================================================

  printf("MU1:%.5f, MU2:%.5f, SIGMA1:%.5f, SIGMA2:%.5f, SIGMA3:%.5f \n",MU1,MU2,SIGMA1,SIGMA2,SIGMA3);
  for(tubeloopi=1; tubeloopi<=NUMTUBE+1; tubeloopi++)
  {
    // Initialize the means and B for all configuration
    // These will not change in the HMC loop
    initMeans(configCurrent);
    initMeans(configOld);
    initMeans(configNew);
 
    // ==================================================================================
    //     Start of HMC Loop (loops over Metropolis Hastings - MC steps)
    // ==================================================================================
 
    fprintf(pStdOut,"START Hybrid Monte Carlo MAIN LOOP\n");
    fprintf(pStdOut,"=======================================================\n");
    acc=0;
    rej=0;
    zeroAverages(tubeAve,&tau);
 
    for(MCloopi=1; MCloopi<=NUMMC; MCloopi++)
    {
      //zero ratio for MH MC test
      ratio=0.0l;
      // ==================================================================================
      //     Perform one SPDE step
      // ==================================================================================
    
      //store savePos.pos values to configCurrent.pos
      // savePos.pos stores the positions in case of rejection of the MHMC
      savePostoConfig(savePos, configCurrent);
  
      //(calculates potentials in config given the positions)
      #pragma omp parallel for
      for(i=0;i<NUMBEAD;i++) {calcPotentials(configCurrent,i);}
  
      //calculate LinvG for the config
      LInverse(configCurrent, vecdg, veci1, veci0);
 
      //do the preconditioned form of the SPDE
      preconditionSPDE(configCurrent, configNew, bb, GaussRandArray, RanNumPointer);
 
      //calculates potentials in config given the positions
      #pragma omp parallel for
      for(i=0;i<NUMBEAD;i++) {calcPotentials(configNew,i);}
 
      //calculate LinvG for the config
      LInverse(configNew, vecdg, veci1, veci0);
 
      //acc ratio of newconfig
      ProbAccRatio(configCurrent, configNew, &ratio);
 
      //calculate the averages for the tubes estimator
      accumulateAverages(tubeAve,configNew,&tau);
 
      fprintf(pStdOut,"SPDE ratio: %+0.10f \n",ratio);
      // ==================================================================================
      //     Start of MD Loop: This loop needs to be focused on for parallelization
      // ==================================================================================
 
      for(MDloopi=1;MDloopi<=NUMMD; MDloopi++)
      {
        //rotate the configuration        
        rotateConfig(&configOld, &configCurrent, &configNew);
 
        //do the MD position update
        MolecularDynamics(configOld, configCurrent, configNew);
  
        //calculates potentials in config given the positions
        #pragma omp parallel for
        for(i=0;i<NUMBEAD;i++) {calcPotentials(configNew,i);}
 
        //calculate LinvG for the config
        LInverse(configNew, vecdg, veci1, veci0);
 
        //calculate the average distance moved in the step and print to std out
        if(MDloopi%WRITESTDOUT==0){
 
          fprintf(pStdOut,"MDi: %.5d | MDi*h: %0.5f | MD ratio: %+0.5f | distance: ",MDloopi,MDloopi*sqrt(2*DT),ratio); //newline is in printDistance function
          printDistance(configNew, savePos);
        }
 
        //acc ratio of newconfig
        ProbAccRatio(configCurrent, configNew, &ratio);
 
        //calculate the averages for the tubes estimator
        accumulateAverages(tubeAve,configNew,&tau);
      }
      // ==================================================================================
      //Metropolis Hastings Monte-Carlo test
      // ==================================================================================
      randUniform = gsl_rng_uniform(RanNumPointer);
      if( exp(ratio/(2.0*TEMP)) > randUniform ){
        acc++;
        saveConfigtoPos(configNew, savePos);
      }
      else{
        rej++;
      }  //end MD loop
      fprintf(pStdOut,"rand=%+0.6f  Exp[ratio]=%+0.6f   dt= %+0.5e     acc= %i      rej= %i  \n",randUniform,exp(ratio/(2.0*TEMP)),DT,acc,rej);

    }  //end MC loop


    normalizeAverages(tubeAve,&tau);
    //call the steepest descent function and reinitialize MPAR and BPAR for the next loop
    tubesSteepestDescent(tubeAve);

    writeConfig(configNew,tubeAve,MCloopi);
    zeroAverages(tubeAve,&tau);
 
    printf("MU1:%.5f, MU2:%.5f, SIGMA1:%.5f, SIGMA2:%.5f, SIGMA3:%.5f \n",MU1,MU2,SIGMA1,SIGMA2,SIGMA3);
  }  //end tube gradient descent loop

  // GSL random number generator release memory and close files
  gsl_rng_free (RanNumPointer);
  fclose(pStdOut);

  return(0);
}

// ==================================================================================
//     Function Declarations
// ==================================================================================

//============================================
void rotateConfig(config **a, config **b, config **c)
// rotates the pointers
//configOld     ->  configNew (a->c)
//configCurrent ->  configOld (b->a)
//configNew     ->  configCurrent (c->b)
//the new configNew is then ready to be overwritten with new positions
{
  config *temp;

  temp=*c;
  *c=*a;
  *a=*b;
  *b=temp;
}

//============================================
void GenGaussRand(double GaussRand[NUMu], gsl_rng *RanNumPointer, double StdDev)
//generate NUMu of Gaussian random nums; save to GaussRandArray[]
{
  int i;
  for(i=0;i<NUMu; i++)
  {
    GaussRand[i]= gsl_ran_gaussian(RanNumPointer,StdDev);
  } 
}

//============================================
void generateBB(double bb[NUMBEAD], double GaussRandArray[NUMu], gsl_rng *RanNumPointer)
//generate a Brownian bridge and store to bb
{
  int n;
  double sqdu, xn;

  //generate NUMBEAD of Gaussian random nums; save to GaussRandArray[]
  GenGaussRand(GaussRandArray, RanNumPointer, 1.0l);

  //****************************
  // read in random numbers for testing precondition function
  // remove for real random numbers to be used
  //int linenum=0;
  //FILE *fptr = fopen("Random_Numbers.dat","r");
  //while( EOF != fscanf(fptr,"%lf",
  //&(GaussRandArray[linenum])) ) {
  //  linenum++;}
  //****************************

  sqdu=sqrt((2.0*TEMP)*DU);
  bb[0]=0.0l;
  for(n=1;n<NUMBEAD;n++)
  {
    bb[n]=bb[n-1]+sqdu*GaussRandArray[n-1];
  }

  xn=bb[NUMu]/((double)(NUMu));
  #pragma omp parallel for
  for(n=1;n<NUMu;n++)
  {
    bb[n]-=((double)(n))*xn;
  }
  bb[NUMBEAD-1]=0.0l;
}

// ============================================
void renormBB(double bb[NUMBEAD])
{
  int n;
  double endPtCorr;
  double sum, term, term0;
  double alpha;

  endPtCorr=gsl_pow_int(bb[0]-bb[NUMBEAD-1],2)/((double)(NUMu));

  sum=0.0l;
  #pragma omp parallel for reduction(+:sum)
  for(n=0;n<NUMu;n++)
  {
    sum+=gsl_pow_int(bb[n]-bb[n+1],2);
  }

  alpha=sqrt((((double)(NUMu))*DU*(2.0*TEMP)-endPtCorr)/(sum-endPtCorr));
  term=(1.0l-alpha)*(bb[NUMBEAD-1]-bb[0])/((double)(NUMu));
  term0=(1.0l-alpha)*bb[0];
  //term and term0 have a subtraction of roughly equal numbers and thus is not very accurate
  // alpha is ~1 with an error of 10^-4 or 5 for sample configs. This makes the routine 
  //nondeterministic between Fortran and C

  #pragma omp parallel for
  for(n=1;n<NUMu;n++)
  {
    bb[n]=alpha*bb[n]+term0+((double)(n-1))*term;
  }
}

// ============================================
void renorm(position *currentpos)
{
  long double alpha;
  double sum, term, term0; 
  double endPtCorr;
  int i,n;

  for(i=0;i<NUMDIM;i++)
  {
    endPtCorr=gsl_pow_int(currentpos[0].pos[i]-currentpos[NUMBEAD-1].pos[i],2)/((double)(NUMu));
    sum=0.0l;
    for(n=0;n<NUMu;n++)
    {
      sum+=gsl_pow_int(currentpos[n].pos[i]-currentpos[n+1].pos[i],2);
    }
    alpha=sqrt((((double)(NUMu))*DU*(2.0*TEMP)-endPtCorr)/(sum-endPtCorr));
    term=(1.0L-alpha)*(currentpos[NUMBEAD-1].pos[i]-currentpos[0].pos[i])/((double)(NUMu));
    term0=(1.0L-alpha)*currentpos[0].pos[i];
    //term and term0 have a subtraction of roughly equal numbers and thus is not very accurate
    // alpha is ~1 with an error of 10^-4 or 5 for sample configs. This makes the routine 
    //nondeterministic between Fortran and C
    for(n=1;n<NUMu;n++)
    {
      currentpos[n].pos[i]=alpha*currentpos[n].pos[i]+term0+(((double)(n))-1.0l)*term;
    }
  }
}

//============================================
void calcPotentials(config *currentConfig, int beadIndex)
{
  currentConfig[beadIndex].posz[0]=currentConfig[beadIndex].pos[0]-currentConfig[beadIndex].posm[0];
  currentConfig[beadIndex].posz[1]=currentConfig[beadIndex].pos[1]-currentConfig[beadIndex].posm[1];

  currentConfig[beadIndex].Energy = VFunc(currentConfig[beadIndex].posm[0],currentConfig[beadIndex].posm[1])+0.5*(currentConfig[beadIndex].B[0]*currentConfig[beadIndex].posz[0]*currentConfig[beadIndex].posz[0]+currentConfig[beadIndex].B[1]*currentConfig[beadIndex].posz[1]*currentConfig[beadIndex].posz[1]+2.0*currentConfig[beadIndex].B[2]*currentConfig[beadIndex].posz[0]*currentConfig[beadIndex].posz[1])-currentConfig[beadIndex].pos[0]*currentConfig[beadIndex].posdm[0]-currentConfig[beadIndex].pos[1]*currentConfig[beadIndex].posdm[1];
  currentConfig[beadIndex].G = 0.5*(currentConfig[beadIndex].B[0]*currentConfig[beadIndex].posz[0]*currentConfig[beadIndex].posz[0]+currentConfig[beadIndex].B[1]*currentConfig[beadIndex].posz[1]*currentConfig[beadIndex].posz[1]+2.0*currentConfig[beadIndex].B[2]*currentConfig[beadIndex].posz[0]*currentConfig[beadIndex].posz[1])+currentConfig[beadIndex].pos[0]*currentConfig[beadIndex].posddm[0]+currentConfig[beadIndex].pos[1]*currentConfig[beadIndex].posddm[1];
  currentConfig[beadIndex].gradG[0] = currentConfig[beadIndex].B[0]*currentConfig[beadIndex].posz[0]+currentConfig[beadIndex].B[2]*currentConfig[beadIndex].posz[1]+currentConfig[beadIndex].posddm[0];
  currentConfig[beadIndex].gradG[1] = currentConfig[beadIndex].B[1]*currentConfig[beadIndex].posz[1]+currentConfig[beadIndex].B[2]*currentConfig[beadIndex].posz[0]+currentConfig[beadIndex].posddm[1];
}

//============================================
void LInverse(config* currentConfig, double vecdg[NUMl], double veci1[NUMl], double veci0[NUMl])
{
  double lasti;
  int i,n;
  double du2 = DU*DU;


  for(i=0;i<NUMDIM;i++)
  {
    #pragma omp parallel for
    for(n=0;n<NUMl;n++){
      vecdg[n]=currentConfig[n+1].gradG[i]*du2;
    }

    veci0[0]=vecdg[0];
    veci1[0]=vecdg[0];
    //this for loop is recursive!!!!
    for(n=1;n<NUMl;n++){
      veci0[n]=veci0[n-1]+vecdg[n];
      veci1[n]=veci1[n-1]+((double)(n+1))*vecdg[n];
    }

    lasti=veci0[NUMl-1]-(currentConfig[NUMBEAD-1].pos[i]-currentConfig[0].pos[i]+veci1[NUMl-1])/((double)(NUMl+1));
    #pragma omp parallel for
    for(n=0;n<NUMl;n++){
      currentConfig[n+1].LinvG[i]=currentConfig[0].pos[i]+((double)(n+1))*(veci0[n]-lasti) - veci1[n];
    }
    currentConfig[0].LinvG[i]=currentConfig[0].pos[i];
    currentConfig[NUMBEAD-1].LinvG[i]=currentConfig[NUMBEAD-1].pos[i];

  }
}

//============================================
void preconditionSPDE(config* currentConfig, config* newConfig, double bb[NUMBEAD], double GaussRandArray[NUMu], gsl_rng *RanNumPointer)
//generates a preconditioned step using the Stochastic Partial Differential Equation
//currentConfig is the incoming configuration that has LinvG and GradG calculated
//newConfig is a temp array that is used to make all of the calsulations without touching currentConfig.pos[*][*]
//newConfig.pos is saved to currentConfig.pos before exiting the function
{

  double qvvel = 0.0l;
  double qvpos = 0.0l;
  int i,n;

  double h=sqrt(2.0l * DT);
  double co=(4.0l-h*h)/(4.0l+h*h);
  double si=(4.0l*h)/(4.0l+h*h); //(h/2)*si
  double hOverTwoSi=(2.0l*h*h)/(4.0l+h*h); //(h/2)*si

  //for grahm shmidt
  double alpha, alphaNum, alphaDenom;

  for(i=0;i<NUMDIM;i++)
  {
    generateBB(bb, GaussRandArray, RanNumPointer);
    //need to make the bb orthogonal to pos without the linear term
    //store pos w/o linear term in newconfig.pos temporarily
    #pragma omp parallel for
    for(n=0;n<NUMBEAD;n++){
      newConfig[n].pos[i]=currentConfig[n].pos[i]-currentConfig[0].pos[i]-(((double)(n))*(currentConfig[NUMBEAD-1].pos[i]-currentConfig[0].pos[i]))/((double)(NUMBEAD -1));
    }
  
    //Gram Schmidt orthogonalization
    alphaNum = 0.0l;
    alphaDenom = 0.0l;
    #pragma omp parallel for reduction(+:alphaNum,alphaDenom)
    for(n=1;n<NUMBEAD;n++)
    {
      alphaNum+=(bb[n]-bb[n-1])*(newConfig[n].pos[i]-newConfig[n-1].pos[i]);
      alphaDenom+=(newConfig[n].pos[i]-newConfig[n-1].pos[i])*(newConfig[n].pos[i]-newConfig[n-1].pos[i]);
    }
    alpha=alphaNum/alphaDenom;

    #pragma omp parallel for
    for(n=0;n<NUMBEAD;n++){ bb[n]=bb[n]-alpha*newConfig[n].pos[i];}

    renormBB(bb);

    #pragma omp parallel for
    for(n=0;n<NUMBEAD;n++){
      newConfig[n].pos[i]=hOverTwoSi*currentConfig[n].LinvG[i] + si*bb[n] + co*currentConfig[n].pos[i];
    }

    //calculate the quadratic variation 
    #pragma omp parallel for reduction(+:qvpos,qvvel)
    for(n=1;n<NUMBEAD;n++){
    qvpos+=gsl_pow_int((newConfig[n].pos[i]-newConfig[n-1].pos[i]),2);
    qvvel+=gsl_pow_int((bb[n]-bb[n-1]),2);
    }
  }

  //print the quadratic variation 
  qvvel *= 0.5/(2.0l*DU*((double)(NUMBEAD-1)));
  qvpos *= 0.5/(2.0l*DU*((double)(NUMBEAD-1)));
  fprintf(pStdOut,"qvvel=%0.10f      qvpos=%0.10f \n",qvvel,qvpos);

}

//============================================
void MolecularDynamics(config *oldConfig, config *currentConfig, config *newConfig)
// perform the MD step
{

  int i,n;
  double h=sqrt(2.0l*DT);
  double co=(4.0-h*h)/(4.0l+h*h);
  double si=sqrt(1.0l-co*co);
  double twoCosMinusOne = 2.0l*co-1.0l;
  double sinH = si*h;

  #pragma omp parallel for private(i)
  for(n=0;n<NUMBEAD;n++) {
    for(i=0;i<NUMDIM;i++)  {
      newConfig[n].pos[i]= (currentConfig[n].pos[i]-oldConfig[n].pos[i]) + twoCosMinusOne*currentConfig[n].pos[i] + sinH*currentConfig[n].LinvG[i];
    }
  }
}

//============================================
void ProbAccRatio(config *currentConfig, config *newConfig, double *ratio)
// calculate the probability accaptance ratio for the step
{

  //there are 3 terms in the error (see notes)
  //lambda1 with h/4 in front (l1h4)
  //lambda1 with h^/16 in front (l1hh16)
  //lambda2 (l2)
  double l1h4, l1hh16, l2;
  int n;

  double h=sqrt(2.0l*DT);
  //double co=(4.0l-(h*h))/(4.0l+(h*h));
  //double si=sqrt(1.0l-(co*co));
  double cot=(4.0l-h*h)/(4.0l*h);
  double csc=(4.0l+h*h)/(4.0l*h);
  
  l1h4=0.0l;
  l1hh16=0.0l;
  l2=0.0l;

  #pragma omp parallel for reduction(+:l1h4,l1hh16,l2)
  for(n=0;n<NUMBEAD;n++)
  {
    //this only works for the two dimensional case
    l1h4+=(cot*currentConfig[n].pos[0]-csc*newConfig[n].pos[0])*currentConfig[n].gradG[0] + (csc*currentConfig[n].pos[0]-cot*newConfig[n].pos[0])*newConfig[n].gradG[0]+(cot*currentConfig[n].pos[1]-csc*newConfig[n].pos[1])*currentConfig[n].gradG[1] + (csc*currentConfig[n].pos[1]-cot*newConfig[n].pos[1])*newConfig[n].gradG[1];
    l1hh16+=currentConfig[n].gradG[0]*currentConfig[n].LinvG[0] - newConfig[n].gradG[0]*newConfig[n].LinvG[0]+currentConfig[n].gradG[1]*currentConfig[n].LinvG[1] - newConfig[n].gradG[1]*newConfig[n].LinvG[1];
    l2+=newConfig[n].G-currentConfig[n].G;
  }
  l1h4*=h*0.25l;
  l1hh16*=h*h*0.0625l;
  l2*=0.5l;

  *ratio+=DU*(l1h4+l1hh16+l2);
}

// ============================== Utility functions =================================

void saveConfig(config *currentConfig, config *saveConfig)
//copies ALL elements of currentConfig to saveConfig
{
  int i,n;
  for(n=0;n<NUMBEAD;n++){
    saveConfig[n].Energy=currentConfig[n].Energy;
    saveConfig[n].G=currentConfig[n].G;
      for(i=0;i<NUMDIM;i++){
        saveConfig[n].pos[i]=currentConfig[n].pos[i];
        saveConfig[n].gradG[i]=currentConfig[n].gradG[i];
        saveConfig[n].LinvG[i]=currentConfig[n].LinvG[i];
} } }

// ============================================
void savePosition(position *currentpos, position *savepos)
//copies ALL elements of currentpos to save pos
{
  int i,n;
  for(n=0;n<NUMBEAD;n++){
    for(i=0;i<NUMDIM;i++){
      savepos[n].pos[i]=currentpos[n].pos[i];
} } }

// ============================================
void saveConfigtoPos(config *currentConfig, position *savepos)
//copies ONLY pos elements of currentConfig to savepos
{
  int i,n;
  for(n=0;n<NUMBEAD;n++){
    for(i=0;i<NUMDIM;i++){
      savepos[n].pos[i]=currentConfig[n].pos[i];
} } }

// ============================================
void savePostoConfig(position *currentpos, config *saveConfig)
//copies ONLY pos elements of currentpos to saveConfig
{
  int i,n;
  for(n=0;n<NUMBEAD;n++){
    for(i=0;i<NUMDIM;i++){
      saveConfig[n].pos[i]=currentpos[n].pos[i];
} } }

//============================================
void writeConfig(config *newConfig, averages *tubeAve, int MCloopi)
{
//print the configs to file
//File order: posx   posy   meanx   meany   posx^2   posx*posy   posy^2
  char filename[80];
  int i,n;

  sprintf(filename,"%s-T%.2f-pos%07d-h0%.5f-b%.5f-a%.5f.dat",PotentialString,TEMP,MCloopi,SIGMA1,SIGMA2,SIGMA3);
  FILE * pWritePos;
  pWritePos = fopen(filename,"w");
  for(n=0;n<NUMBEAD;n++){
    for(i=0;i<NUMDIM;i++){
    fprintf(pWritePos, "%+.15e \t",newConfig[n].pos[i]);
    }
    for(i=0;i<4;i++){
    fprintf(pWritePos, "%+.15e \t",tubeAve[n].mean[i]);
    }
    for(i=0;i<10;i++){
    fprintf(pWritePos, "%+.15e \t",tubeAve[n].meanB[i]);
    }
    fprintf(pWritePos, "%+.15e \t",newConfig[n].B[0]);
    fprintf(pWritePos, "\n");
  }
  fclose(pWritePos);

  /*
  //call the python script to plot the paths
  char * path;
  char pyCall[300];
  path=getenv("PWD");
  if(path != NULL){
    strcpy(pyCall,"2dHMC-plot ");
    strcat(pyCall,path);
    strcat(pyCall,"/");
    strcat(pyCall,filename);
    strcat(pyCall," 1000");
    system(pyCall); 
  }
  */

}

//============================================
void printPositionPos(position* currentpos, int beadIndex)
//Print the positions of the position struct 
{
  printf("position x             position y\n");
  printf("%+.17e %+.17e \n",currentpos[beadIndex].pos[0],currentpos[beadIndex].pos[1]);
}

//============================================
void printConfigPos(config* currentConfig, int beadIndex)
//Print the positions of the config struct 
{
  printf("position x             position y\n");
  printf("%+.17e %+.17e \n",currentConfig[beadIndex].pos[0],currentConfig[beadIndex].pos[1]);
}

//============================================
void printConfigPot(config* currentConfig, int beadIndex)
//Print all information of the config struct calculated in calcPotentials
// (Postitions, grad G, Energy, G)
{
  printf("position x            position y               gradG x               gradG y               Energy                      G \n");
  printf("%+.15e %+.15e %+.15e %+.15e %+.15e %+.15e \n",currentConfig[beadIndex].pos[0],currentConfig[beadIndex].pos[1],currentConfig[beadIndex].gradG[0],currentConfig[beadIndex].gradG[1],currentConfig[beadIndex].Energy,currentConfig[beadIndex].G);
}

//============================================
void printConfigAll(config* currentConfig, int beadIndex)
//Print all information of the config struct
{
  printf("position x          position y           gradG x          gradG y         Energy               G              LinvG x            LinvG y\n");
  printf("%+.10e %+.10e %+.10e %+.10e %+.10e %+.10e %+.10e %+.10e\n",currentConfig[beadIndex].pos[0],currentConfig[beadIndex].pos[1],currentConfig[beadIndex].gradG[0],currentConfig[beadIndex].gradG[1],currentConfig[beadIndex].Energy,currentConfig[beadIndex].G,currentConfig[beadIndex].LinvG[0],currentConfig[beadIndex].LinvG[1]);
}

//============================================
void printDistance(config *newConfig, position *savePos)
{
  int i,n;
  double tempSum;

  tempSum=0.0l;
  #pragma omp parallel for private(i) reduction(+:tempSum)
  for(n=1;n<NUMBEAD-1;n++)
  {
    for(i=0;i<NUMDIM;i++)
    {
      tempSum+=gsl_pow_int(newConfig[n].pos[i]-savePos[n].pos[i],2);
    }
  }
  fprintf(pStdOut,"%+.10f \n",tempSum/(NUMBEAD-2));

}

// ============================== Tube Average functions =================================

void zeroAverages(averages *tubeAve, int *tau)
  {
  int n,m;
  for(n=0;n<NUMBEAD;n++)
  {
    for(m=0;m<4;m++)
    {
      tubeAve[n].mean[m]=0.0l;
    }
    for(m=0;m<10;m++)
    {
      tubeAve[n].meanB[m]=0.0l;
    }
  }

  *tau=0;
}

// ============================================
void initMeans(config *currentConfig)
//initializes the means and B 
{
  int n;
  double dun;
  for(n=0;n<NUMBEAD;n++){
    dun=DU*((double)(n));
    currentConfig[n].posm[0]=meanxFunc(dun);
    currentConfig[n].posm[1]=meanyFunc(dun);
    currentConfig[n].posdm[0]=dmeanxFunc(dun);
    currentConfig[n].posdm[1]=dmeanyFunc(dun);
    currentConfig[n].posddm[0]=ddmeanxFunc(dun);
    currentConfig[n].posddm[1]=ddmeanyFunc(dun);
    currentConfig[n].B[0]=BxxFunc(dun);
    currentConfig[n].B[1]=ByyFunc(dun);
    currentConfig[n].B[2]=BxyFunc(dun);
  }
}

//============================================
void accumulateAverages(averages *tubeAve, config *newConfig, int *tau)
{
  int n;
  double I, Iou, ImIou;
  double Fx, ddVx, Fy, ddVy;

  *tau=*tau+1;

  I=1.0;
  Iou=0.0;

  //calculate I and Iou. these are just numbers and are used below
  for(n=0;n<NUMBEAD;n++)
  {
    Fx=4.0*newConfig[n].pos[0]-4.0*gsl_pow_int(newConfig[n].pos[0],3);
    ddVx=-4.0+12.0*newConfig[n].pos[0]*newConfig[n].pos[0];
    Fy=-2.0*newConfig[n].pos[1];
    ddVy=2.0;

    I+=DT*(0.5*(Fx*Fx+Fy*Fy)-TEMP*(ddVx+ddVy));

    Iou+=DT*newConfig[n].G;
  }

  ImIou=I-Iou;

  //Calculate dm/dtau and dBij/dtau and save to an array for averaging later
  for(n=0;n<NUMBEAD;n++)
  {
    tubeAve[n].mean[0]+=newConfig[n].pos[0];
    tubeAve[n].mean[1]+=newConfig[n].pos[1];
    tubeAve[n].mean[2]+=ImIou*(newConfig[n].LinvG[0]-newConfig[n].pos[0]);
    tubeAve[n].mean[3]+=ImIou*(newConfig[n].LinvG[1]-newConfig[n].pos[1]);
    tubeAve[n].meanB[0]+=newConfig[n].pos[0]*newConfig[n].pos[0];
    tubeAve[n].meanB[1]+=newConfig[n].pos[1]*newConfig[n].pos[1];
    tubeAve[n].meanB[2]+=newConfig[n].pos[0]*newConfig[n].pos[1];
    tubeAve[n].meanB[3]+=ImIou*newConfig[n].posz[0]*newConfig[n].posz[0];
    tubeAve[n].meanB[4]+=ImIou*newConfig[n].posz[1]*newConfig[n].posz[1];
    tubeAve[n].meanB[5]+=ImIou*newConfig[n].posz[0]*newConfig[n].posz[1];
    tubeAve[n].meanB[6]+=newConfig[n].posz[0]*newConfig[n].posz[0];
    tubeAve[n].meanB[7]+=newConfig[n].posz[1]*newConfig[n].posz[1];
    tubeAve[n].meanB[8]+=newConfig[n].posz[0]*newConfig[n].posz[1];
    tubeAve[n].meanB[9]+=ImIou;
  }
}

//============================================
void normalizeAverages(averages *tubeAve, int *tau)
{
  int n,m;
  double oneOverTau=1.0l/((double)(*tau));

  #pragma omp parallel for
  for(n=0;n<NUMBEAD;n++)
  {
    for(m=0;m<4;m++)
    {
      tubeAve[n].mean[m]*=oneOverTau;
    }
    for(m=0;m<10;m++)
    {
      tubeAve[n].meanB[m]*=oneOverTau;
    }
  }
}

//============================================
void tubesSteepestDescent(averages *tubeAve)
// calculate the gradient of the KL divergence and minimize it
//passes a new set of parameters for the smooth functions defined in constants.h
{
  int n;
  double tempMU1=0.0;
  double tempMU2=0.0;
  double tempSIGMA1=0.0;
  double tempSIGMA2=0.0;
  double tempSIGMA3=0.0;
  double gammaDescent=2.0;
  double dun;

  for(n=0;n<NUMBEAD;n++)
  {
    dun=DU*((double)(n));
    //mean integration
    tempMU1+= meanxDParx1Func(dun)*(tubeAve[n].mean[2]);
    tempMU2+= meanxDParx2Func(dun)*(tubeAve[n].mean[3]);
    //B integration 
    tempSIGMA1+= BxxDParx1Func(dun)*(tubeAve[n].meanB[3]-tubeAve[n].meanB[6]*tubeAve[n].meanB[9]);
    tempSIGMA2+= BxxDParx2Func(dun)*(tubeAve[n].meanB[3]-tubeAve[n].meanB[6]*tubeAve[n].meanB[9]);
    tempSIGMA3+= BxxDParx3Func(dun)*(tubeAve[n].meanB[3]-tubeAve[n].meanB[6]*tubeAve[n].meanB[9]);
  }

  //Steepest Descent: newp = oldp - gamma * Del Function
  MU1= MU1-10000.0*gammaDescent*tempMU1/((double)(NUMBEAD));
  MU2= MU2-0.1*gammaDescent*tempMU2/((double)(NUMBEAD));
  SIGMA1= SIGMA1-gammaDescent*tempSIGMA1/((double)(NUMBEAD));
  SIGMA2= SIGMA2-gammaDescent*tempSIGMA2/((double)(NUMBEAD));
  SIGMA3= SIGMA3-gammaDescent*tempSIGMA3/((double)(NUMBEAD));
}
