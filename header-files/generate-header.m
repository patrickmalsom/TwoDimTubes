#!/usr/local/bin/MathematicaScript -script

GenHeader[]:=Module[{},
  (* Mathematica script to generate potentials for the c header files *)
  Cdef[fun__]:=StringReplace[StringReplace[ToString[CForm[Simplify[Expand[fun]]]],{"Sinh("->"sinh(","Cosh("->"cosh(","Tanh("->"tanh(","Csch("->"1/sinh(","Sech("->"1/cosh(","Coth("->"1/tanh(","Power("->"gsl_pow_int("," "->""}],"gsl_pow_int(E,"->"exp("];
  file=OpenWrite["TwoDimTubes.h"];

  WriteString[file,"// TwoDimTubes.h - Constants for Two Dimensional Tubes HMC \n"];
  WriteString[file,"// ============================= Preprocessor Definitions ============================= \n"];
  WriteString[file," \n"];

  WriteString[file,"// Path Definitions \n"];
  WriteString[file,"#define NUMDIM    "<>ToString[NUMDIM]<>"           // Spacial Dimensions \n"];
  WriteString[file,"#define NUMBEAD   "<>ToString[NUMBEAD]<>"      // Path Points \n"];
  WriteString[file,"#define DU        "<>ToString[DU]<>"      // path time step size \n"];
  WriteString[file,"#define PREDT     "<>ToString[PREDT]<>"      // DT=PreDT*DU^2 (path time) \n"];
  WriteString[file," \n"];

  WriteString[file,"// Temperature Definition \n"];
  WriteString[file,"#define TEMP      "<>ToString[TEMPERATURE]<>" \n"];
  WriteString[file," \n"];

  WriteString[file,"// Incrimenter Definitions \n"];
  WriteString[file,"#define NUMMD     "<>ToString[NUMMD]<>"         // Number of MD steps  \n"];
  WriteString[file,"//      NUMMD     ~3/(2*sqrt(2*PreDT*DU^2)) <- Approx optimal value of NUMMD \n"];
  WriteString[file,"#define NUMMC     "<>ToString[NUMMC]<>"       // Number of Metropolis Hastings MC steps \n"];
  WriteString[file,"#define NUMTUBE   "<>ToString[NUMTUBE]<>"        // Number of tube steepest descent steps \n"];
  WriteString[file," \n"];

  WriteString[file,"// Constants for writing to stdout and config \n"];
  WriteString[file,"#define WRITESTDOUT  "<>ToString[WRITESTDOUT]<>"       // How often to print to stdout (# of MD loops) \n"];
  WriteString[file,"const char PotentialString[]=\""<>PotentialString<>"\"; // Potential Description  \n"];
  WriteString[file," \n"];

  WriteString[file,"// ============================= Tubes Definitions ============================= \n"];
  WriteString[file,"// Initial Parameters for the means. To be optimized using KL gradient descent. \n"];
  WriteString[file,"double MU1 = "<>ToString[mu1]<>"; // steepness of the transition \n"];
  WriteString[file,"double MU2 = "<>ToString[mu2]<>"; // well location at finite temp \n"];
  WriteString[file," \n"];

  WriteString[file,"double SIGMA1 = "<>ToString[sigma1]<>"; // approx well hessian \n"];
  WriteString[file,"double SIGMA2 = "<>ToString[sigma2]<>"; \n"];
  WriteString[file,"double SIGMA3 = "<>ToString[sigma3]<>"; \n"];
  WriteString[file," \n"];

  WriteString[file,"// ============================= Parameter Function Defn's ============================= \n"];
  WriteString[file,"//Potential Function","\n"];
  WriteString[file,"#define VFunc(x,y)      "<>Cdef[V],"\n"];
  WriteString[file,"","\n"];

  WriteString[file,"//potentials for calculation of<I-Iou>","\n"];
  WriteString[file,"#define FxFunc(x,y)      "<>Cdef[-D[V,x]],"\n"];
  WriteString[file,"#define ddVxFunc(x,y)     "<>Cdef[D[D[V,x],x]],"\n"];
  WriteString[file,"#define FyFunc(x,y)     "<>Cdef[-D[V,y]],"\n"];
  WriteString[file,"#define ddVyFunc(x,y)      "<>Cdef[D[V,{y,2}]],"\n"];
  WriteString[file,"#define gradGxFunc(x,y)      "<>Cdef[Simplify[Expand[D[0.5*(D[V,x]^2+D[V,y]^2)-TEMP*(D[V,{x,2}]+D[V,{y,2}]),x]]]],"\n"];
  WriteString[file,"#define gradGyFunc(x,y)      "<>Cdef[Simplify[Expand[D[0.5*(D[V,x]^2+D[V,y]^2)-TEMP*(D[V,{x,2}]+D[V,{y,2}]),y]]]],"\n"];
  WriteString[file,"","\n"];

  WriteString[file,"//Smooth functions for use in Tubes HMC","\n"];
  WriteString[file,"#define meanxFunc(t)   "<>Cdef[mx],"\n"];
  WriteString[file,"#define meanyFunc(t)   "<>Cdef[my],"\n"];
  WriteString[file,"#define dmeanxFunc(t)   "<>Cdef[D[mx,t]],"\n"];
  WriteString[file,"#define dmeanyFunc(t)      "<>Cdef[D[my,t]],"\n"];
  WriteString[file,"#define ddmeanxFunc(t)     "<>Cdef[D[mx,{t,2}]],"\n"];
  WriteString[file,"#define ddmeanyFunc(t)     "<>Cdef[D[my,{t,2}]],"\n"];
  WriteString[file,"","\n"];

  WriteString[file,"//B parameters","\n"];
  WriteString[file,"#define BxxFunc(t)   "<>Cdef[Bxx],"\n"];
  WriteString[file,"#define ByyFunc(t)   "<>Cdef[Byy],"\n"];
  WriteString[file,"#define BxyFunc(t)   "<>Cdef[Bxy],"\n"];
  WriteString[file,"","\n"];

  WriteString[file,"//derivatives wrt the parameters","\n"];
  WriteString[file,"#define meanxDParx1Func(t)     "<>Cdef[D[mx,MU1]],"\n"];
  WriteString[file,"#define meanxDParx2Func(t)     "<>Cdef[D[mx,MU2]],"\n"];
  WriteString[file,"#define BxxDParx1Func(t)        "<>Cdef[D[Bxx,SIGMA1]],"\n"];
  WriteString[file,"#define BxxDParx2Func(t)      "<>Cdef[D[Bxx,SIGMA2]],"\n"];
  WriteString[file,"#define BxxDParx3Func(t)      "<>Cdef[D[Bxx,SIGMA3]],"\n"];

  Close[file];
];
