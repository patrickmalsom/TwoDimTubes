#!/usr/local/bin/MathematicaScript -script

argv=$ScriptCommandLine;

data=Import[argv[[2]],"Table"];
Export[argv[[3]],Transpose[{data[[All,1]]+ToExpression[argv[[4]]],data[[All,2]]}],"TSV"];

Quit[];
