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



(* :Context: TransOC` *)
(* :Author: Ahsan Zeb, School of Physics and Astronomy, St Andrews, UK*)
(* :Summary: The TransOC package calculates the incoherent charge transport through
			 an organic micorcavity in the strong light-matter coupling regime *)
(* :Copyright: Copyright (c) 2017 Ahsan Zeb, distributed under the XXXX license *)
(* :Package Version: 1.0 *)
(* :Mathematica Version: 10.4.0.0 *)


If[TrueQ[$VersionNumber<10],Print["Sorry, TransOC only works with version 10+ of the Wolfram Language."];Abort[]];


BeginPackage["TransOC`"];


(*auxiliary: *)
$TransOCVersion::usage="$TransOCVersion gives the version number of the TransOC package.";
$TransOCDebugQ::usage="Setting $TransOCDebugQ=True results in debugging messages being printed during evaluations.";


(* function: *)
Trajectory::usage="Trajectory runs a Montr Carlo trajectory.

\!\(\*
StyleBox[\"Optional\",\nFontSize->14]\)\!\(\*
StyleBox[\" \",\nFontSize->14]\)\!\(\*
StyleBox[\"Arguments\",\nFontSize->14]\)\!\(\*
StyleBox[\":\",\nFontSize->14]\)
Ntot: 5; Number of molecules in the system.;
m0: 1; number of excitations in the coupled light-matter sub-system;
Ds0: {}; List of sites that are D;
Phis0: {};List of sites that are \[Phi]s;					
g: 1.0; strenght of light-matter coupling for a single molecule: (\!\(\*SubscriptBox[\(\[Omega]\), \(R\)]\) = g \!\(\*SqrtBox[SubscriptBox[\(N\), \(Active\)]]\));
dw: 0.0; Detuning, \!\(\*SubscriptBox[\(\[Omega]\), \(0\)]\)-\!\(\*SubscriptBox[\(\[Omega]\), \(c\)]\);
tpar: {0.1,0.001,0.001,0.1,0.1,0.1,0.1,0.1}; 
		Hopping parameters, {th,tl,tlh,thl,JhR,JlR,JhL,JlL};
Eb: {0.0,0.0}; Barriers at the left and right ontacts: Eb={Ebl,Ebr};
kappa: 0.005; Cavity losses \[Kappa];
gamma: 0.005; Exciton non-radiative decay rate \[Gamma];
Er: 0.5; Energy E.r due to the applied electric field.
w0: 2.0; Exciton bare energy \!\(\*SubscriptBox[\(\[Omega]\), \(0\)]\).
beta: 40.0; \[Beta] = 1/\!\(\*SubscriptBox[\(k\), \(b\)]\)T;
maxiter: 20; Max number of iterations/time steps in the dynamics.
includecross: True; Include cross hops, L-H and H-L?
BlockInjection: {False,False}; Block injection of charge carriers? {BlockLeft,BlockRight}. 
				Left contact is cathode, right contact is anode. Block electron from left and holes from right.
M,wv,lambda: 3,0.2,1.0; Vibrational levels, energy, coupling to electronic states.
AlwaysLP: True; Quickly relax the quantum system to its lowest eigenstate before an other hop can take place?
VibAssis: False; Include vibrational assistance for the Zener tunneling, H-L transitions from active sites H levels?"


Begin["`Private`"];


$TransOCVersion=1.;
debug:=TrueQ[$TransOCDebugQ];


(*get code from all files in the PackageFiles directory (separate files are used to make it easier to find stuff);
they all inherit the context Crystallica`Private`, so that all files together form the actual package;
messages and option defaults can be found in the file of their respective parent function*)
Get/@FileNames["*.m",FileNameJoin[{DirectoryName[$InputFileName],"PackageFiles"}]];


End[];
EndPackage[];
