#!/usr/local/bin/MathematicaScript -script

argv=$ScriptCommandLine;

Num = ToExpression[argv[[2]]];

tempData=Import["2WellTubes-T0.15-pos0010000.dat", "Table"];

Do[
  tempData = tempData + Import["2WellTubes-T0.15-pos"<>StringDrop[ToString[NumberForm[i, 3, NumberPadding -> {"0", ""}]], 1] <> "0000.dat", "Table"];
,{i, 2, Num}];

Export["meanpath.dat",tempData/Num,"TSV"];

Quit[];
