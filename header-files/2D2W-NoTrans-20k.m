#!/usr/local/bin/MathematicaScript -script

(* Definitions of constants *)
NUMDIM = 2;           (* Spacial Dimensions *)
NUMBEAD = 20001;      (* Path Points *)
DU = 0.0005;          (* path time step size *)
PREDT = 6000.0;       (* DT=PreDT*DU^2 (path time) *)
TEMPERATURE = 0.15;   (* temperature *)
NUMMD = 50;           (* Number of MD steps  ~3/(2*sqrt(2*PreDT*DU^2)) *)
NUMMC = 10000;         (* Number of Metropolis Hastings MC steps *)
NUMTUBE = 100;        (* Number of tube steepest descent steps *)
WRITESTDOUT= 50;      (* How often to print to stdout (# of MD loops) *)
PotentialString ="2WellTubes";  (* Potential Description  *)

(* path parameters *)
mu1 = 2.00; (* steepness of the transition *)
mu2 = 0.97; (* well location at finite temp *)
sigma1 = 8.28537; (* approx well hessian *)
sigma2 = -0.90; (*  *)
sigma3 = 2.72; (*  *)

(* Definitions of the parameterized functions *)
V=(x^2-1.)^2+y^2;
mx=0.9756
my=0.;
Bxx=8.28537*8.28537
Byy=4.
Bxy=0.

(* generate the header file *)
Get["generate-header.m"];
GenHeader[];

Quit[];
