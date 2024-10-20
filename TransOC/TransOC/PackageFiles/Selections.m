(* ::Package:: *)

(************************************************************************)
(* This file was generated automatically by the Mathematica front end.  *)
(* It contains Initialization cells from a Notebook file, which         *)
(* typically will have the same name as this file except ending in      *)
(* ".nb" instead of ".m".                                               *)
(*                                                                      *)
(* This file is intended to be loaded into the Mathematica kernel using *)
(* the package loading commands Get or Needs.  Doing so is equivalent   *)
(* to using the Evaluate Initialization Cells menu command in the front *)
(* end.                                                                 *)
(*                                                                      *)
(* DO NOT EDIT THIS FILE.  This entire file is regenerated              *)
(* automatically each time the parent Notebook file is saved in the     *)
(* Mathematica front end.  Any changes you make to this file will be    *)
(* overwritten.                                                         *)
(************************************************************************)



WhichSite[rlist_]:=Module[{lrz,Rlist,eta,whichsite},
lrz=Length[rlist];
Rlist=Insert[Accumulate[rlist],0,1];
eta=RandomReal[{0,Rlist[[-1]] }];
whichsite=-1;
Do[If[eta>Rlist[[i]]&&eta<=  Rlist[[i+1]],whichsite=i;Break[]],{i,1,lrz}];
If[whichsite==-1,
Print["WhichSite: not found! Aborting... "];
Print["rlist: ",rlist];
Print["Rlist: ",Rlist];
Abort[];
];
whichsite
];
WhichSiteChannel[rlist_]:=Module[{lrz,Rlist,eta,whichsite,wchannel,x},
lrz=Length[rlist];
Rlist=Insert[Accumulate[rlist],0,1];
eta=RandomReal[{0,Rlist[[-1]] }];
whichsite=-1;wchannel=-1;x=-1;
Do[If[eta>Rlist[[i]]&&eta<=  Rlist[[i+1]],x=i;Break[]],{i,1,lrz}];
If[x==-1,
Print["WhichSiteChannel: not found! Aborting... "];
Print["rlist: ",rlist];
Print["Rlist: ",Rlist];
Abort[];
];
whichsite=Quotient[3+x,4];
wchannel=Mod[x,4,1];
{whichsite,wchannel}
];

EUinds={
{1,1,2,3},
{1,1,2,3},
{1,1,2,3},
{1,1,2,3},
{4,4,5,6},
{4,4,5,6},
{7,7,8,9},
{7,7,8,9},
{11,11,11,11},
{10,10,10,10},
{11,11,11,11},
{10,10,10,10},
{11,11,11,11},
{10,10,10,10},
{11,11,11,11},
{10,10,10,10},
{13,13,13,13},
{13,13,13,13},
{12,12,12,12},
{12,12,12,12},
{13,13,13,13},
{13,13,13,13},
{12,12,12,12},
{12,12,12,12},
{2,2,2,2},
{2,2,2,2}
};
ListdNm={{0,0},{0,-1},{0,1},{2,1},{2,0},{2,2},
{-2,-1},{-2,-2},{-2,0},{1,1},{1,0},{-1,-1},{-1,0}};

DEQCsum[w0_,Ebr_,Ebl_,Er_]:=Module[{dqc,dEbulk,dEcont,signEr},
dEbulk={0,0,-w0,w0};
dEcont={
-Ebr,-Ebr+w0,-Ebl,-Ebl+w0,
Ebr-w0,Ebr,Ebl-w0,Ebl,
Ebr,-Ebr+w0,Ebr-w0,-Ebr,
Ebl,-Ebl+w0,Ebl-w0,-Ebl,
-w0,-w0};
signEr={
1,-1,-1,1,1,-1,-1,1,
1,1,-1,-1,-1,-1,1,1,
-1,1,-1,1,1,-1,1,-1,
0,0};
dqc=Table[Table[
If[wj<= 8,dEbulk[[wc]],dEcont[[wj-8]]]
-signEr[[wj]]*Er,
{wc,1,4}],{wj,1,26}];
(*Print["dqc = ",dqc];*)
dqc
];

sectors[Eigs_]:=Module[{deg,nlevels,X,Es,NElevels},
X=Tally[SetAccuracy[Eigs,8]];
Es=X[[All,1]];NElevels=X[[All,2]];nlevels=Length[NElevels];
deg=Insert[Accumulate[NElevels],0,1];
{nlevels,Es,deg}];

(*SelectFinalState[Ef_,Uf_,coeff_,wrates_]:=Module[
{X,psiff,deg,r,whichsec,Elevels,nlevels,i1,i2,Rlist,x,eta,Eff},
X=Tally[Ef];
Elevels=X[[All,2]];nlevels=Length[Elevels];
i1=0;i2=0;x=0;Rlist={0};deg={};
Do[
i2+=Elevels[[i]];
AppendTo[deg,{i1,i2}];
r=coeff[[i1+1;;i2]];
x+=Norm[r]^2;AppendTo[Rlist,x];
i1=i2;
,{i,1,nlevels}];
eta=RandomReal[{0,Rlist[[-1]] }];
whichsec=-1;
Do[If[eta>Rlist[[i]]&&eta\[LessEqual]  Rlist[[i+1]],whichsec=i;Break[]],{i,1,nlevels}];

If[whichsec\[Equal]-1,
Print["WhichSiteChannel: not found! Aborting... "];
Print["Rlist: ",Rlist];
];

If[whichsec\[Equal]-1,
{"whichsec not found! Aborting... Rlist: ",Rlist," & eta:" , eta}>>>file;
Abort[];];
{i1,i2}=deg[[whichsec]];
psiff=Sum[coeff[[i]]*Uf[[All,i]],{i,i1+1,i2}];
psiff/Norm[psiff];
Eff=XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX;
{Eff,psiff}
];*)

SelectFinalState[Ef_,Uf_,coeff_,wrates_]:=Module[
{X,psiff,deg,r,whichsec,Elevels,nlevels,i1,i2,eta,Rlist,Eff},

{nlevels,Eff,deg}=sectors[Ef];
Rlist=Insert[Accumulate[wrates],0,1];

If[Im[Rlist[[-1]]]!=0,Print["SelectFinalState: Im[Rlist[[-1]]]\[NotEqual]0, some issue.... Should I Chop?"],
Abort[];
];

eta=RandomReal[{0,Rlist[[-1]] }];
whichsec=-1;
Do[If[eta>Rlist[[i]]&&eta<=  Rlist[[i+1]],whichsec=i;Break[]],{i,1,nlevels}];

If[whichsec==-1,
Print["WhichSiteChannel: not found! Aborting... "];
Print["Rlist: ",Rlist];
Abort[];
];

i1=deg[[whichsec]];i2=deg[[whichsec+1]];
psiff=Sum[coeff[[i]]*Uf[[All,i]],{i,i1+1,i2}];
psiff/Norm[psiff];

{Eff[[whichsec]],psiff}
];




