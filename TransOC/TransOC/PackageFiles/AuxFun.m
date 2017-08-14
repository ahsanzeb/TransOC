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



c[n_,m_]:=Binomial[n,m];
comb[N_,k_]:=Subsets[Table[i,{i,1,N}],{k}];
Num[x_]:=N[x];
(*Dropp[Ds_,elem_]:=Drop[Ds,{Position[Ds,elem][[1,1]]}];*)
Dropp[Ds_,elem_]:=Complement[Ds,{elem}];(*compared to old above def, complement also removes repetition, but that is not present in our case unless something is wrong! *)

