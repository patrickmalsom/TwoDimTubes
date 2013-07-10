#!/usr/local/bin/MathematicaScript -script
(*argv=$ScriptCommandLine;*)

indata=Import["2WellTubes-T0.15-pos0001000.dat","TSV"];

dDdB=-indata[[All,14]]+indata[[All,15]]*indata[[All,16]];

inpar=ToExpression[Last[Import["steep-descent/parameterHist.dat","List"]]];

func[list__]:=Table[68.6475-(list[[1]]+Sum[list[[n+1]]*(x-5)^(2n),{n,1,Length[inpar]-1,1}])*Exp[-10.(x-5)^2],{x,0,10,0.0005}];

dBdb[n_]:=If[n==0,-Table[Exp[-2(t-5)^2],{t,0,10,0.0005}],-Table[(t-5)^(2n)*Exp[-10.(t-5)^2],{t,0,10,0.0005}]];

outpar=inpar-0.01*{1,100,0,0,0,0,0}*Table[Total[dDdB*dBdb[i]],{i,0,6}];
Print[Table[Total[dDdB*dBdb[i]],{i,0,6}]]

Byy=Table[4.,{t,0,10,0.0005}];
Bxy=Table[0.,{t,0,10,0.0005}];

Export["inputB.dat",Transpose[{func[outpar],Byy,Bxy}],"TSV"];

parWrite=Import["steep-descent/parameterHist.dat","List"];
AppendTo[parWrite,outpar];
Export["steep-descent/parameterHist.dat",parWrite,"List"];


Quit[];
