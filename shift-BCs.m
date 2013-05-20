#!/usr/local/bin/MathematicaScript -script
argv=$ScriptCommandLine;

(* script that shifts the boundary conditions on an input file *)
(* in1 : input file name *)
(* in2 : output file name *)
(* in3 : shift in the x direction *)
(* in4 : shift in the y direction *)

data=Import[argv[[2]],"Table"];
Export[argv[[3]],Transpose[{data[[All,1]]+ToExpression[argv[[4]]],data[[All,2]]+ToExpression[argv[[5]]]}],"TSV"];

Quit[];
