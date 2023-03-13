(* ::Package:: *)

(* ::Text:: *)
(*(*Numerical  values used in our calculation for known constants*)*)


{
(*Planck constant*)h = 6.582118989999999 10^-25, 
(*Muon mean lifetime*)\[Tau]\[Mu]=2.197 10^-6/h, 
(*tau mean lifetime*)\[Tau]\[Tau]=290.6 10^-15/h, 
(*Muon decay width*)\[CapitalGamma]\[Mu] = 1/\[Tau]\[Mu], 
(*Tau decay width*)\[CapitalGamma]\[Tau] = 1/\[Tau]\[Tau], 
(*Current bound on the br of the tau decay \[Tau]->\[Mu]\[Gamma]*)BRexp\[Tau]\[Mu]\[Gamma] = 4.4 10^-8, 
(*Fermion mass values*)me = 5.1 10^-5, mmu = 0.10565, mtau = 1.7768, mt = 173.21,mb=4.5, mc = 1.27, mZ=91.18,
(*Up and down quark charges*)Qu = 2/3,Qd=-1/3,
(*LQ charges*) QS1 = 5/3,QS0=2/3,
(*Weinberg angle sine and cosine values*) 
sW=Sqrt[0.22],
cW=Sqrt[1-0.22^2],
(*Minimum and maximum bounds on the muon anomaly at 95 % C.L. in units of 10^-11*)
\[CapitalDelta]min=251-1.96 59,
\[CapitalDelta]max=251+1.96 59,
\[Alpha]QED=1/137.,
(*Z Couplings to charged leptons, up and down quarks*)
 {gVl=g/(2cW)(-1/2 + 2 sW^2), gVu=g/(2cW)( 1/2 - 4/3 sW^2),gVd=g/(2cW)( -1/2 + 2/3 sW^2)};
 {gAl=g/(2cW)(-1/2) , gAu=g/(2cW)( 1/2) ,gAd=g/(2cW)( -1/2) };
(*Z couplings to a charge 5/3 scalar LQ*)
gSSZ1=-1/2-sW^2 QS1,
gSSZ0=1/2-sW^2 QS2,
(*Conversion factor from GeV^-1 to cm*)
ConvFac=1.97 10^-14,
mZ=91.18;
ee=Sqrt[4Pi \[Alpha]QED],g=ee/sW,v=246,mH=125,mW=80.4,\[Alpha]SmZ=0.1184,gs=Sqrt[4Pi \[Alpha]SmZ]};


(* ::Text:: *)
(*(*Functions to calculate the contribution of a scalar LQ to the muon anomaly and the decay \[Tau]->\[Mu]\[Gamma]*)*)


F1\[Gamma]p[z1_, z2_,mS_] := ( 
{mq=Sqrt[z2] mS, m=Sqrt[z1] mS};\[CapitalDelta]p=B0[0, z2 mS^2, z2 mS^2] - B0[z1 mS^2, z2 mS^2, mS^2];
\[CapitalDelta]m= B0[0, mS^2, mS^2] - B0[z1 mS^2,z2 mS^2, mS^2];
  -(1 - (2*(1 + \[CapitalDelta]m ))/z1 + (2*(1 + \[CapitalDelta]p )*z2)/z1 + 
  (2*\[CapitalDelta]m *(-1 + z1 - z2) + 4*(-1 + z1 + z2 + \[CapitalDelta]p *z2))/
   (z1^2 + (-1 + z2)^2 - 2*z1*(1 + z2)))/z1)


F2\[Gamma]p[z1_, z2_,mS_] := ( 
{mq=Sqrt[z2] mS, m=Sqrt[z1] mS};\[CapitalDelta]p=B0[0, z2 mS^2, z2 mS^2] - B0[z1 mS^2, z2 mS^2, mS^2];\[CapitalDelta]m= B0[0, mS^2, mS^2] - B0[z1 mS^2,z2 mS^2, mS^2];
-(z1^3+2*\[CapitalDelta]m *((-2+z1)*z1+(-1+z2)^2)-2*(-1+z2)^2*(-1+z2+\[CapitalDelta]p *z2)+z1*(-3+z2*(2+z2+2*\[CapitalDelta]p *(1+z2))))/(z1^2*(z1^2+(-1+z2)^2-2*z1*(1+z2))))


F\[Gamma]p[z1_, z2_,mS_, Q1_:Qu,Q2_:QS1]:=Q1 F1\[Gamma]p[z1, z2,mS]+Q2 F2\[Gamma]p[z1, z2,mS]


G1\[Gamma]p[z1_, z2_,mS_] := ( 
{mq=Sqrt[z2] mS, m=Sqrt[z1] mS};\[CapitalDelta]p=B0[0, z2 mS^2, z2 mS^2] - B0[z1 mS^2, z2 mS^2, mS^2];\[CapitalDelta]m= B0[0, mS^2, mS^2] - B0[z1 mS^2,z2 mS^2, mS^2];
-(2*((1+\[CapitalDelta]p )*z1^2+\[CapitalDelta]m *(1+z1-z2)+(-1+z2)*(-1+z2+\[CapitalDelta]p *z2)-z1*(-2+\[CapitalDelta]p +2*(1+\[CapitalDelta]p )*z2)))/(z1*(z1^2+(-1+z2)^2-2*z1*(1+z2))))
 


 G2\[Gamma]p[z1_, z2_,mS_] := ( 
{mq=Sqrt[z2] mS, m=Sqrt[z1] mS};\[CapitalDelta]p=B0[0, z2 mS^2, z2 mS^2] - B0[z1 mS^2, z2 mS^2, mS^2];\[CapitalDelta]m= B0[0, mS^2, mS^2] - B0[z1 mS^2,z2 mS^2, mS^2];
-(2*(\[CapitalDelta]m *(-1+z1+z2)+(1+z1-z2)*(-1+z1+z2+\[CapitalDelta]p *z2)))/(z1*(z1^2+(-1+z2)^2-2*z1*(1+z2))))


G\[Gamma]p[z1_, z2_,mS_, Q1_:Qu,Q2_:QS1]:=Q1 G1\[Gamma]p[z1, z2,mS]+Q2 G2\[Gamma]p[z1, z2,mS]


F\[Gamma][z1_, z2_,mS_,Q1_:Qu,Q2_:QS1]:=Q1 F1\[Gamma]p[z1, z2,mS]+Q2 F2\[Gamma]p[z1, z2,mS]


L\[Gamma]PV[mi_,mj_,mk_,mS_,Qk_,QS_,gLik_,gLkj_,gRik_,gRkj_]:=(mi(-((gLik*gLkj*mi-gRik*gRkj*mj)*(Qk-QS))/(2*(mi-mj)*(mi+mj))+(-((gRik*gRkj*mi-gLik*gLkj*mj)*mk^2*(Qk+QS))/(2*mi*(mi-mj)*mj*(mi+mj))+((gRik*gRkj*mi-gLik*gLkj*mj)*mS^2*(Qk+QS))/(2*mi*(mi-mj)*mj*(mi+mj)))*B0[0,mk^2,mS^2]+((gRik*mi*(-2*gLkj*(mi^2-mj^2)*mk*(Qk+QS)+gRkj*mj*(mi^2*(Qk-QS)+mk^2*(Qk+QS)))+gLik*gLkj*(mj^2*mk^2*(Qk+QS)-mi^2*(mj^2*(Qk-QS)+2*mk^2*(Qk+QS))))*B0[mi^2,mk^2,mS^2])/(2*mi*(mi-mj)^2*(mi+mj)^2)+((gLik*gLkj*mi*mj*(mj^2*(Qk-QS)+mk^2*(Qk+QS))+gRik*(2*gLkj*mj*(mi^2-mj^2)*mk*(Qk+QS)+gRkj*(-2*mj^2*mk^2*(Qk+QS)+mi^2*(mj^2*(-Qk+QS)+mk^2*(Qk+QS)))))*B0[mj^2,mk^2,mS^2])/(2*(mi-mj)^2*mj*(mi+mj)^2)+(mk*(gRik*gRkj*mj*mk-gLkj*(gRik*(mi^2-mj^2)+gLik*mi*mk))*Qk*C0[mi^2,mj^2,0,mk^2,mS^2,mk^2])/((mi-mj)*(mi+mj))+mS^2*(-((gRik*gRkj*mi*mj+gLik*gLkj*(-2*mi^2+mj^2))*(Qk+QS)*B0[mi^2,mk^2,mS^2])/(2*mi*(mi-mj)^2*(mi+mj)^2)-((gLik*gLkj*mi*mj+gRik*gRkj*(mi^2-2*mj^2))*(Qk+QS)*B0[mj^2,mk^2,mS^2])/(2*(mi-mj)^2*mj*(mi+mj)^2)+((gLik*gLkj*mi-gRik*gRkj*mj)*QS*C0[mi^2,mj^2,0,mS^2,mk^2,mS^2])/(mi^2-mj^2))))


R\[Gamma]PV[mi_,mj_,mk_,mS_,Qk_,QS_,gLik_,gLkj_,gRik_,gRkj_]:=(mi(-((gRik*gRkj*mi-gLik*gLkj*mj)*(Qk-QS))/(2*(mi-mj)*(mi+mj))+(-((gLik*gLkj*mi-gRik*gRkj*mj)*mk^2*(Qk+QS))/(2*mi*(mi-mj)*mj*(mi+mj))+((gLik*gLkj*mi-gRik*gRkj*mj)*mS^2*(Qk+QS))/(2*mi*(mi-mj)*mj*(mi+mj)))*B0[0,mk^2,mS^2]+((gLik*mi*(-2*gRkj*(mi^2-mj^2)*mk*(Qk+QS)+gLkj*mj*(mi^2*(Qk-QS)+mk^2*(Qk+QS)))+gRik*gRkj*(mj^2*mk^2*(Qk+QS)-mi^2*(mj^2*(Qk-QS)+2*mk^2*(Qk+QS))))*B0[mi^2,mk^2,mS^2])/(2*mi*(mi-mj)^2*(mi+mj)^2)+((gRik*gRkj*mi*mj*(mj^2*(Qk-QS)+mk^2*(Qk+QS))+gLik*(2*gRkj*mj*(mi^2-mj^2)*mk*(Qk+QS)+gLkj*(-2*mj^2*mk^2*(Qk+QS)+mi^2*(mj^2*(-Qk+QS)+mk^2*(Qk+QS)))))*B0[mj^2,mk^2,mS^2])/(2*(mi-mj)^2*mj*(mi+mj)^2)+(mk*(-(gRik*gRkj*mi*mk)+gLik*(gRkj*(-mi^2+mj^2)+gLkj*mj*mk))*Qk*C0[mi^2,mj^2,0,mk^2,mS^2,mk^2])/((mi-mj)*(mi+mj))+mS^2*(-((gLik*gLkj*mi*mj+gRik*gRkj*(-2*mi^2+mj^2))*(Qk+QS)*B0[mi^2,mk^2,mS^2])/(2*mi*(mi-mj)^2*(mi+mj)^2)-((gRik*gRkj*mi*mj+gLik*gLkj*(mi^2-2*mj^2))*(Qk+QS)*B0[mj^2,mk^2,mS^2])/(2*(mi-mj)^2*mj*(mi+mj)^2)+((gRik*gRkj*mi-gLik*gLkj*mj)*QS*C0[mi^2,mj^2,0,mS^2,mk^2,mS^2])/(mi^2-mj^2))))


\[CapitalGamma]fitofj\[Gamma][mi_,mj_,mk_,mS_,Qk_,QS_,gLik_,gLkj_,gRik_,gRkj_]:=(
L\[Gamma]=Sum[L\[Gamma]PV[mi,mj,mk[[j]],mS[[i]],Qk[[j]],QS[[i]],gLik[[i,j]],gLkj[[i,j]],gRik[[i,j]],gRkj[[i,j]]],{i,1,Length[mS]},{j,1,Length[mk]}];
R\[Gamma]=Sum[R\[Gamma]PV[mi,mj,mk[[j]],mS[[i]],Qk[[j]],QS[[i]],gLik[[i,j]],gLkj[[i,j]],gRik[[i,j]],gRkj[[i,j]]],{i,1,Length[mS]},{j,1,Length[mk]}];
LR2=Abs[L\[Gamma]]^2+Abs[R\[Gamma]]^2;(3 /(16Pi^2))^2(LR2*(mi^2-mj^2)^3)/(16*mi^5*Pi))


afPV[mi_,mk_,mS_,Qk_,QS_,gLik_,gRik_]:=(
z1=(mi/mS)^2;
z2=(mk/mS)^2;
FF\[Gamma]=Qk F1\[Gamma]p[z1,z2,mS]+QS F2\[Gamma]p[z1,z2,mS];
GG\[Gamma]=Qk G1\[Gamma]p[z1,z2,mS]+QS G2\[Gamma]p[z1,z2,mS];
-3Sqrt[z1]/(16Pi^2)(Sqrt[z1](Abs[gLik]^2+Abs[gRik]^2)FF\[Gamma]+2(gLik gRik) Sqrt[z2]GG\[Gamma]))


CF=4/3;


(* ::Text:: *)
(*(*The following functions receive as arguments the following arrays: *)
(*mS array with masses of LQs,*)
(*\[Lambda]Lik, \[Lambda]Rik arrays with the LQ couplings to the initial lepton and the internal quark,*)
(*\[Lambda]Lkj, \[Lambda]Rkj arrays with the LQ couplings to the initial lepton and the internal quark,*)
(*mk array with massed of internal quarks, *)*)


(* ::Text:: *)
(*(*Br for the contribution of a scalar LQ to the \[Tau]->\[Mu]\[Gamma] decay for a set of internal LQs and quarks*)*)


BR\[Tau]to\[Mu]\[Gamma][mk_,mS_,QS_,\[Lambda]Lik_,\[Lambda]Lkj_,\[Lambda]Rik_,\[Lambda]Rkj_]:=\[Tau]\[Tau] \[CapitalGamma]fitofj\[Gamma][mtau,mmu,mk,mS ,{2/3,2/3}ee,QS ee,\[Lambda]Lik,\[Lambda]Lkj,\[Lambda]Rik,\[Lambda]Rkj]


(* ::Text:: *)
(*(*LQ contribution to the muon anomaly for a set of internal LQs and quarks*)*)


MDMPV[mi_,mk_,mS_,Qk_,QS_,\[Lambda]Lik_,\[Lambda]Rik_]:=Sum[afPV[mi,mk[[j]],mS[[i]],Qk[[j]],QS[[i]],\[Lambda]Lik[[i,j]],\[Lambda]Rik[[i,j]]],{i,1,Length[mS]},{j,1,Length[mk]}];
