#!/usr/local/bin/MathematicaScript -script

(* Mathematica script to generate potentials for the c header files *)

V=(x^2-2.)^2+y^2*(x^2+1.)/2.;
mx=MU2*Tanh[MU1*(t-5.)]/Tanh[5.*MU1];
my=0.;
Bxx=(SIGMA1+(-SIGMA1+SIGMA2)/Exp[SIGMA3*(t-5)])^2
Byy=D[V,{y,2}]^2/.{x->1,y->0}
Bxy=D[D[V,y],x]^2/.{x->1,y->0}

Cdef[fun__]:=StringReplace[StringReplace[ToString[CForm[Expand[fun]]],{"Sinh("->"sinh(","Cosh("->"cosh(","Tanh("->"tanh(","Csch("->"1/sinh(","Sech("->"1/cosh(","Coth("->"1/tanh(","Power("->"gsl_pow_int("," "->""}],"pow(E,"->"exp("]
(* takes a regular mathematica function and turns it into a c function (uses gsl_pow_int) *)

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
Print["#define BxxFunc(t)   "<>Cdef[Bxx]]
Print["#define ByyFunc(t)   "<>Cdef[Byy]]
Print["#define BxyFunc(t)   "<>Cdef[Bxy]]
Print[""]

Print["//derivatives wrt the parameters"]
Print["#define meanxDParx1Func(t)     "<>Cdef[D[mx,MU1]]]
Print["#define meanxDParx2Func(t)     "<>Cdef[D[mx,MU2]]]
Print["#define BxxDParx1Func(t)        "<>Cdef[D[Bxx,SIGMA1]]]
Print["#define BxxDParx2Func(t)      "<>Cdef[D[Bxx,SIGMA2]]]
Print["#define BxxDParx3Func(t)      "<>Cdef[D[Bxx,SIGMA3]]]

Quit[];
