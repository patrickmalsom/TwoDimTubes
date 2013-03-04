#!/usr/local/bin/MathematicaScript -script

(*
input: ./path-distances.m inputPos.dat
*)

argv=$ScriptCommandLine;

data=Import[argv[[2]],"Table"];
dDdBxx=Total[-data[[All, 14 + 0]] + data[[All, 14 + 1]]*data[[All, 14 + 2]]];

Print["sum of dD/dBxx:  "<>ToString[dDdBxx]]


Quit[];
