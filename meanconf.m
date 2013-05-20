#!/usr/local/bin/MathematicaScript -script
argv=$ScriptCommandLine;

(* script to take an average of all output files *)
(* files MUST start at 10000 and go in increments of 10000 to be able to be averaged correctly *)
(* you must be in teh main directory where the files are stored. no absolute paths allowed *)
(* in1 : number of files to average ex:3 -> (10000+20000+30000)/3 *)

Num = ToExpression[argv[[2]]];

tempData=Import["2WellTubes-T0.15-pos0010000.dat", "Table"];

Do[
  tempData = tempData + Import["2WellTubes-T0.15-pos"<>StringDrop[ToString[NumberForm[i, 3, NumberPadding -> {"0", ""}]], 1] <> "0000.dat", "Table"];
,{i, 2, Num}];

Export["meanpath.dat",tempData/Num,"TSV"];

Quit[];
