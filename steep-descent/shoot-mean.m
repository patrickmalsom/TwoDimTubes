#!/usr/local/bin/MathematicaScript -script
(*argv=$ScriptCommandLine;*)

(* set the directory *)
SetDirectory[Directory[]];

(* import the data set and take the reverse of G *)
data=Import["2WellTubes-T0.15-pos0001000.dat","Table"];
glist=Take[(data[[All,25]]+Reverse[data[[All,25]]])/2,{10000,20001}];
(* import the mean and take half of it *)
meandata=Import["inputMean.dat","Table"];
dmeandata=meandata[[All,2]];

(* define the starting point of forming the new glist *)
(* this newglist will be monatonic *)
tempstart=1000;
glistn=Take[glist,{1,1000}];
prev=glist[[tempstart]];
next=glist[[tempstart+1]];
(* if previous is greater than next, append next, else append previous *)
Do[
If[prev>next,
AppendTo[glistn,next],
AppendTo[glistn,prev]];
prev=Last[glistn];
next=glist[[tempstart+i+1]];
,{i,1,2000}]
(* offset the new g list by G0 *)
newglist=glistn-Last[glistn];

(* perform the shooting on the ODE dm/dt=Sqrt[2(G-G0)] *)
(* remembering to mix the new G and the old mean to form the new mdot *)
shoot[G0_,\[Alpha]_]:=Module[{},m=Table[0.0,{i,1,Length[newglist]}];
Do[m[[i+1]]=m[[i]]+0.0005 
If[newglist[[i]]-G0<0,
0,
Sqrt[2 \[Alpha](newglist[[i]]-G0)+(1-\[Alpha])*(dmeandata[[i+10000]])^2]];
,{i,1,Length[m]-1}];
m
]
(* make a small guess and iterate until the correct G0 is found *)
(* the end value of the mean should be approx the BC *)
guess=0.001;
While[Last[shoot[guess,0.2]]-Last[meandata][[1]]>0,guess=guess+0.0005;]

(* form the new mean *)
newm=Join[Drop[Reverse[-m],-1],m];
flatnum=10000-(Length[newm]-1)/2;
newm=Join[Table[First[newm],{i,1,flatnum}],newm,Table[Last[newm],{i,1,flatnum}]]*-0.9756/newm[[1]];

(* calculate dm/dt with finite differences *)
dnewm=Join[{0.0},Table[(newm[[i]]-newm[[i-1]])/(0.0005),{i,2,Length[newm]}]];

(* calculate d^2m/dt^2 with finite differences *)
newmt=MovingAverage[newm,501];
ddnewm=Join[Table[0.0,{i,1,251}],Table[(newmt[[i+1]]-2newmt[[i]]+newmt[[i-1]]),{i,2,Length[newmt]-1}]/0.0005^2,Table[0.0,{i,1,251}]];

(* Export the data *)
Export["newMean.tsv",Transpose[{newm,dnewm,ddnewm}],"TSV"];

Quit[];
