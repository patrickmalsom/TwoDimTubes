#!/usr/local/bin/MathematicaScript -script
argv=$ScriptCommandLine;

(* in1 : input data file *)

data=Import[argv[[2]],"Table"];

dDdmx=Total[data[[All, 8 + 0]] - data[[All, 8 + 1]]*data[[All, 8 + 2]]];
dDdmy=Total[data[[All, 11 + 0]] - data[[All, 11 + 1]]*data[[All, 11 + 2]]];

dDdBxx=Total[-data[[All, 14 + 0]] + data[[All, 14 + 1]]*data[[All, 14 + 2]]];
dDdByy=Total[-data[[All, 17 + 0]] + data[[All, 17 + 1]]*data[[All, 17 + 2]]];

Print["sum of dD/dmx:  "<>ToString[dDdmx]]
Print["sum of dD/dmy:  "<>ToString[dDdmy]]
Print["sum of dD/dBxx:  "<>ToString[dDdBxx]]
Print["sum of dD/dByy:  "<>ToString[dDdByy]]

Quit[];
