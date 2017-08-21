(* ::Package:: *)

pwd="~/repos/transport-code/TransOC/";
AppendTo[$Path,pwd];
Needs["TransOC`"]





file=FileNameJoin[{Directory[],"out-Er-dw-uncoupled-Ntraj.dat"}];
Ntot00=8;
Ntraj=30;
(* Average number of Ds+Phis fraction of total sites *)
Ntotb2=IntegerPart[Ntot00/3];
outall={};
Ers=Range[0, 1, 0.05];
DW=Range[-0.5,0.5, 0.1];
Do[
outallx={};
Do[
out={};
Do[
Ds00x=Table[RandomInteger[{1,Ntot00}],{i,1,RandomInteger[{0,Ntotb2}]}];
Ds00=If[Ds00x =={},{},Transpose[Tally[Ds00x]][[1]]];
Phi00x=Complement[Table[i,{i,1,Ntot00}],Ds00];
x0=Table[RandomInteger[{1,Ntotb2}],{i,1,Ntotb2}];
x=Transpose[Tally[x0]][[1]];
Phis00=If[Length[Phi00x]< Ntotb2,Phi00x,Phi00x[[x]]];
m00=RandomInteger[{0,Ntot00}];
outtraj=Trajectory[Ntot-> Ntot00, m0-> m00, Ds0-> Ds00, Phis0->Phis00,Eb-> {0.7,0.7},Er->er ,maxiter->500, g -> 0.15,dw-> ddw];
Print[" Er, dw = ",er,", ",ddw,"  .... "];
AppendTo[out,outtraj];
,{j,1,Ntraj}];
AppendTo[outallx,out];
,{er,Ers}];
AppendTo[outall,outallx];
,{ddw,DW}];
Export[file,outall];







file=FileNameJoin[{Directory[],"out-Er-dw-coupled-Ntraj.dat"}];
Ntot00=8;
Ntraj=30;
(* Average number of Ds+Phis fraction of total sites *)
Ntotb2=IntegerPart[Ntot00/3];
outall={};
Ers=Range[0, 1, 0.05];
DW=Range[-0.5,0.5, 0.1];
Do[
outallx={};
Do[
out={};
Do[
Ds00x=Table[RandomInteger[{1,Ntot00}],{i,1,RandomInteger[{0,Ntotb2}]}];
Ds00=If[Ds00x =={},{},Transpose[Tally[Ds00x]][[1]]];
Phi00x=Complement[Table[i,{i,1,Ntot00}],Ds00];
x0=Table[RandomInteger[{1,Ntotb2}],{i,1,Ntotb2}];
x=Transpose[Tally[x0]][[1]];
Phis00=If[Length[Phi00x]< Ntotb2,Phi00x,Phi00x[[x]]];
m00=RandomInteger[{0,Ntot00}];
outtraj=Trajectory[Ntot-> Ntot00, m0-> m00, Ds0-> Ds00, Phis0->Phis00,Eb-> {0.7,0.7},Er->er ,maxiter->500, g -> 0.15,dw-> ddw];
Print[" Er, dw = ",er,", ",ddw,"  .... "];
AppendTo[out,outtraj];
,{j,1,Ntraj}];
AppendTo[outallx,out];
,{er,Ers}];
AppendTo[outall,outallx];
,{ddw,DW}];
Export[file,outall];


