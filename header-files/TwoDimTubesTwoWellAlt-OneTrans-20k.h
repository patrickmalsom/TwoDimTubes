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
const char PotentialString[]="2WellTubesAlt";// Potential Description 

//=============================================================================
//                   Tubes Definitions 
//=============================================================================

/* Mathematica code to generate the definitions (because CForm is broken...)
Cdef[fun__]:=StringReplace[StringReplace[ToString[CForm[Expand[fun]]],{"Sinh("->"sinh(","Cosh("->"cosh(","Tanh("->"tanh(","Csch("->"1/sinh(","Sech("->"1/cosh(","Coth("->"1/tanh(","Power("->"gsl_pow_int("," "->""}],"pow(E,"->"exp("]

V=(x^2-2.)^2+y^2*(x^2+1.)/2.;
mx=MU2*Tanh[MU1*(t-5.)]/Tanh[5.*MU1];
my=0.;

Print["//Potential Function"]
Print["#define VFunc(x,y)      "<>Cdef[V]]
Print[""]

Print["//potentials for calculation of<I-Iou>"]
Print["#define FxFunc(x,y)      "<>Cdef[-D[V,x]]]
Print["#define ddVxFunc(x,y)     "<>Cdef[D[D[V,x],x]]]
Print["#define FyFunc(x,y)     "<>Cdef[-D[V,y]]]
Print["#define ddVyFunc(x,y)      "<>Cdef[D[V,{y,2}]]]
Print[""]

Print["//Smooth functions for use in Tubes HMC"]
Print["#define meanxFunc(t)   "<>Cdef[mx]]
Print["#define meanyFunc(t)   "<>Cdef[my]]
Print["#define dmeanxFunc(t)   "<>Cdef[D[mx,t]]]
Print["#define dmeanyFunc(t)      "<>Cdef[D[my,t]]]
Print["#define ddmeanxFunc(t)     "<>Cdef[D[mx,{t,2}]]]
Print["#define ddmeanyFunc(t)     "<>Cdef[D[my,{t,2}]]]
Print[""]

Print["//B parameters"]
Print["//#define BxxFunc(t)    "<>Cdef[D[V,{x,2}]^2/.{x->1,y->0}]]
Print["#define BxxFunc(t) pow(SIGMA1+(-SIGMA1+SIGMA2)/exp(SIGMA3*pow(-5.+t,2)),2)     "]
Print["#define ByyFunc(t)     "<>Cdef[D[V,{y,2}]^2/.{x->1,y->0}]]
Print["#define BxyFunc(t)    "<>Cdef[D[D[V,x],y]^2/.{x->1,y->0}]]
Print[""]

Print["//derivatives wrt the parameters"]
Print["#define meanxDParx1Func(t)     "<>Cdef[D[mx,MU1]]]
Print["#define meanxDParx2Func(t)     "<>Cdef[D[mx,MU2]]]
Print["#define BxxDParx1Func(t)      2*(1-exp(-(SIGMA3*pow(-5.+t,2))))*(SIGMA1+(-SIGMA1+SIGMA2)/exp(SIGMA3*pow(-5.+t,2)))"]
Print["#define BxxDParx2Func(t)      (2*(SIGMA1+(-SIGMA1+SIGMA2)*exp(-SIGMA3*pow(-5.+t,2))))*exp(-SIGMA3*pow(-5.+t,2))"]
Print["#define BxxDParx3Func(t)      (-2*(-SIGMA1+SIGMA2)*(SIGMA1+(-SIGMA1+SIGMA2)/exp(SIGMA3*pow(-5.+t,2)))*pow(-5.+t,2))/exp(SIGMA3*pow(-5.+t,2))"]
*/

//Potential Function
#define VFunc(x,y)      4.-4.*gsl_pow_int(x,2)+gsl_pow_int(x,4)+0.5*gsl_pow_int(y,2)+0.5*gsl_pow_int(x,2)*gsl_pow_int(y,2)

//potentials for calculation of<I-Iou>
#define FxFunc(x,y)      8.*x-4*gsl_pow_int(x,3)-1.*x*gsl_pow_int(y,2)
#define ddVxFunc(x,y)     -8.+12*gsl_pow_int(x,2)+1.*gsl_pow_int(y,2)
#define FyFunc(x,y)     -1.*y-1.*gsl_pow_int(x,2)*y
#define ddVyFunc(x,y)      1.+1.*gsl_pow_int(x,2)

//Smooth functions for use in Tubes HMC
#define meanxFunc(t)   MU2*1/tanh(5.*MU1)*tanh(MU1*(-5.+t))
#define meanyFunc(t)   0.
#define dmeanxFunc(t)   MU1*MU2*1/tanh(5.*MU1)*gsl_pow_int(1/cosh(MU1*(-5.+t)),2)
#define dmeanyFunc(t)      0
#define ddmeanxFunc(t)     -2*gsl_pow_int(MU1,2)*MU2*1/tanh(5.*MU1)*gsl_pow_int(1/cosh(MU1*(-5.+t)),2)*tanh(MU1*(-5.+t))
#define ddmeanyFunc(t)     0

//B parameters
//#define BxxFunc(t)    16.
#define BxxFunc(t) pow(SIGMA1+(-SIGMA1+SIGMA2)/exp(SIGMA3*pow(-5.+t,2)),2)     
#define ByyFunc(t)     4.
#define BxyFunc(t)    0.

//derivatives wrt the parameters
#define meanxDParx1Func(t)     -5.*MU2*1/tanh(5.*MU1)*gsl_pow_int(1/cosh(MU1*(-5.+t)),2)+MU2*t*1/tanh(5.*MU1)*gsl_pow_int(1/cosh(MU1*(-5.+t)),2)-5.*MU2*gsl_pow_int(Csch(5.*MU1),2)*tanh(MU1*(-5.+t))
#define meanxDParx2Func(t)     1/tanh(5.*MU1)*tanh(MU1*(-5.+t))
#define BxxDParx1Func(t)      2*(1-exp(-(SIGMA3*pow(-5.+t,2))))*(SIGMA1+(-SIGMA1+SIGMA2)/exp(SIGMA3*pow(-5.+t,2)))
#define BxxDParx2Func(t)      (2*(SIGMA1+(-SIGMA1+SIGMA2)*exp(-SIGMA3*pow(-5.+t,2))))*exp(-SIGMA3*pow(-5.+t,2))
#define BxxDParx3Func(t)      (-2*(-SIGMA1+SIGMA2)*(SIGMA1+(-SIGMA1+SIGMA2)/exp(SIGMA3*pow(-5.+t,2)))*pow(-5.+t,2))/exp(SIGMA3*pow(-5.+t,2))
