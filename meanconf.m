#!/usr/local/bin/MathematicaScript -script

argv=$ScriptCommandLine;

importDir = argv[[2]];
Num = ToExpression[argv[[3]]];

tempData = Import[importDir <> "2WellTubes-T0.15-pos0010000.dat", "Table"];
Do[
  tempData = tempData + Import[importDir<>"2WellTubes-T0.15-pos"<>StringDrop[ToString[NumberForm[i, 3, NumberPadding -> {"0", ""}]], 1] <> "0000.dat", "Table"];
,{i, 2, Num}];

Export[importDir<>"meanpath.dat",tempData/Num,"TSV"];

Quit[];
