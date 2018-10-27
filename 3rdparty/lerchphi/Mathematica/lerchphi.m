(*******************************************************************

-------------------------------
Lerch's transcendent Phi(z,s,v)
-------------------------------

This program copyright by

  Sergej V. Aksenov 

(http://www.geocities.com/saksenov) and 

  Ulrich D. Jentschura 

(jentschura@physik.uni-freiburg.de), 2002.

Version 1.00 (May 1, 2002)

Calling sequence:

LerchPhiCNCT[z,s,v]    (for positive z, positive v)

The program calculates Lerch's Phi transcendent Phi(z,s,v) 
with the following parameters specified as options: relative 
accuracy (default 14 digits), maximum iterations of the convergence 
acceleration loop (default 100), transformation parameter beta 
(default 1) and n (default 0). Standard Mathematica precision 
(16 decimal figures) is used throughout the calculation.

********************************************************************)

LerchPhiCNCT::usage = "LerchPhiCNCT[z,s,v] returns Lerch's transcendent 
of positive arguments z, v, and s using the combined-nonlinear condensation 
transformation. For related formulas,
see the e-print math.NA/0202009 and references therein. 
z is restricted to the interval (0,1).  
LerchPhiCNCT can be given four options: CNCTRelativeAccuracy specifies 
the AccuracyGoal for the calculation; CNCTMaximumIterations specifies 
the maximum number of iterations of CNC transforms; and CNCTBeta and 
CNCTN specify parameters of the transformation beta (shift parameter)
and n (initial element of the transformation), respectively."

Options[LerchPhiCNCT] = {CNCTRelativeAccuracy->14,
  CNCTMaximumIterations->100, CNCTBeta->1, CNCTN->0}

LerchPhiCNCT[z_,s_,v_,opts___?OptionQ] :=
  Module[{j,omega0,storeaj,num,den,o,res,factors},

(* Subroutine that adds a term to A_j *)
ajstep[inp_List] := Module[{k,sum,bjk},
k = inp[[1]]; sum=inp[[2]];
bjk = 2^k z^(2^k (j+1)-1)/(v+2^k (j+1)-1)^s;
sum+=bjk; k++;
Return[{k,sum,bjk}]];
      
(* Subroutine that does the recurrence update 
   of a numerator or a denominator of the CNC transform S_k^n *)
recur[inp_List,pos_] := Module[{loc},
loc=inp;
loc[[pos]]=loc[[pos+1]]-loc[[pos]] factors[[pos]];
Return[loc]];
      
(* Subroutine that iterates CNC transforms S_k^n *)
sknstep[inp_List] := Module[{i,omega,sn,skn,eps},
i = inp[[1]]; omega = inp[[2]]; 
sn = inp[[3]]; skn = inp[[4]]; eps = inp[[5]];
sn += omega; AppendTo[storeaj,(-1)^i omega];
omega = (-1)^(i+1) * If[EvenQ[i],
  1/2 (storeaj[[i/2+1]]-z^(i/2)/(v+i/2)^s), j=i+1;
  NestWhile[ajstep,{1,z^j/(v+j)^s,z^j/(v+j)^s},
    (#[[3]]/#[[2]]>10^(-2-acc))&][[2]]];
AppendTo[num,sn/omega];
AppendTo[den,1/omega];
factors = Reverse[Which[i==0,{0},i==1,{0,1},True,
  Prepend[Table[(beta+n+i-1) (beta+n+i-2)/(beta+n+i+o-2)/(beta+n+i+o-3),
    {o,1,i}],0]]];
num = Fold[recur,num,Table[o,{o,Length[num]-1,1,-1}]];
den = Fold[recur,den,Table[o,{o,Length[den]-1,1,-1}]];
skn = RotateLeft[skn];
eps = RotateLeft[eps];
skn[[2]] = num[[1]]/den[[1]];
eps[[2]] = Abs[skn[[2]]-skn[[1]]];
i++;
Return[{i,omega,sn,skn,eps}]];
      
(* Read in options *)
acc = CNCTRelativeAccuracy/.{opts}/.Options[LerchPhiCNCT];
imax = CNCTMaximumIterations/.{opts}/.Options[LerchPhiCNCT];
beta = CNCTBeta/.{opts}/.Options[LerchPhiCNCT];
n = CNCTN/.{opts}/.Options[LerchPhiCNCT];
      
(* Initialize loop for CNC transforms *)
j = 0;
omega0 = NestWhile[ajstep,{1,z^j/(v+j)^s,z^j/(v+j)^s},
  (#[[3]]/#[[2]]>10^(-2-acc))&];
num={};
den={};
storeaj={};
      
(* Calculate CNC transforms S_k^n until 
   termination criteria are satisfied *)
      
res = NestWhile[sknstep,{0,omega0[[2]],0,{0,0},{0,0}},
  Not[#[[1]]>imax-1||(#[[1]]>1&&(#[[5,2]]==0||
   (#[[5,2]]<#[[5,1]]&&
   Abs[2 (#[[5,1]])^2/(#[[5,1]]-#[[5,2]])/#[[4,2]]]<10^(-acc))))]&];
      
(* Return last iterate of CNC transforms S_k^n *)
If[res[[1]]>imax-1,
  Message[LerchPhiCNCT::maxit,acc,res[[1]]]];
Return[res[[4, 2]]]] /; s \[Element] Reals && 
  (v \[Element] Reals && (v > 0 || (Not[IntegerQ[v]] && IntegerQ[s]))) && 
  (z \[Element] Reals && z > 0 && z < 1);

LerchPhiCNCT::maxit =
  "Warning: algorithm has not achieved relative accuracy of `1` digits 
   after maximum `2` iterations.";

Print[">> LerchPhiCNCT package 1.00 -- Copyright (C) 2002 Sergej Aksenov\n
and Ulrich Jentschura -- Type ?LerchPhiCNCT for information"];


