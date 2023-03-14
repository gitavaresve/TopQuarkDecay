(* ::Package:: *)

(* ::Text:: *)
(* (*We consider an effective Lagrangian for a charge 5/3 scalar leptoquark coupling with a lepton and a quark, where \[Lambda]Ltauc,\[Lambda]Rtauc,\[Lambda]Ltaut,and \[Lambda]Rtaut stand for the left-and right-LQ couplings with the tau and the c and t quarks, whereas  \[Lambda]Lmuc,\[Lambda]Rmuc,\[Lambda]Lmut,and \[Lambda]Rmut correspond to the couplings to the muon. The program calculates the branching ratio for the top quark rare decay t->c\[Gamma]\[Gamma] from the one-loop contribution from a scalar leptoquark with electric charge of 5/3 and a set of internal  leptons as described in manuscript https://arxiv.org/abs/2208.05064 by A. Bola\[NTilde]os, R. S\[AAcute]nchez V\[EAcute]lez and G. Tavares Velasco.  The program needs the LoopTools package (https://feynarts.de/looptools/) to evaluate the Passarino-Veltman scalar *)
(*functions. *)*)


(* ::Text:: *)
(*(*Scalar Passarino-Veltman functions involved in the calculation of the t->c\[Gamma]\[Gamma] decay due to a scalar leptoquark of charge 5/3 and a set of internal leptons as described in manuscript.*)*)


funesc={B1[0,mi^2*yk^2,mi^2*yS^2],B1[0,mi^2*yS^2,mi^2*yk^2],B0[0,mi^2*yk^2,mi^2*yk^2],B0[0,mi^2*yS^2,mi^2*yS^2],
B0[mi^2,mi^2*yk^2,mi^2*yS^2],B0[mi^2*xs,mi^2*yk^2,mi^2*yS^2],B0[mi^2*xt,mi^2*yk^2,mi^2*yS^2],
B0[mi^2-mi^2*xs-mi^2*xt,mi^2*yk^2,mi^2*yk^2],B0[mi^2-mi^2*xs-mi^2*xt,mi^2*yS^2,mi^2*yS^2],
C0[0,0,mi^2*xs,mi^2*yk^2,mi^2*yk^2,mi^2*yS^2],C0[0,0,mi^2*xs,mi^2*yk^2,mi^2*yS^2,mi^2*yS^2],
C0[0,0,mi^2*xt,mi^2*yk^2,mi^2*yk^2,mi^2*yS^2],C0[0,0,mi^2*xt,mi^2*yk^2,mi^2*yS^2,mi^2*yS^2],
C0[0,0,mi^2-mi^2*xs-mi^2*xt,mi^2*yk^2,mi^2*yk^2,mi^2*yk^2],
C0[0,0,mi^2-mi^2*xs-mi^2*xt,mi^2*yS^2,mi^2*yS^2,mi^2*yS^2],C0[mi^2,0,mi^2*xs,mi^2*yk^2,mi^2*yS^2,mi^2*yS^2],
C0[mi^2,0,mi^2*xs,mi^2*yS^2,mi^2*yk^2,mi^2*yk^2],C0[mi^2,0,mi^2*xt,mi^2*yk^2,mi^2*yS^2,mi^2*yS^2],
C0[mi^2,0,mi^2*xt,mi^2*yS^2,mi^2*yk^2,mi^2*yk^2],C0[mi^2,0,mi^2-mi^2*xs-mi^2*xt,mi^2*yk^2,mi^2*yS^2,mi^2*yk^2],
C0[mi^2,0,mi^2-mi^2*xs-mi^2*xt,mi^2*yS^2,mi^2*yk^2,mi^2*yS^2],
D0[mi^2,0,0,0,mi^2*xs,mi^2*xt,mi^2*yk^2,mi^2*yS^2,mi^2*yS^2,mi^2*yk^2],
D0[mi^2,0,0,0,mi^2*xs,mi^2*xt,mi^2*yS^2,mi^2*yk^2,mi^2*yk^2,mi^2*yS^2],
D0[mi^2,0,0,0,mi^2*xs,mi^2-mi^2*xs-mi^2*xt,mi^2*yk^2,mi^2*yS^2,mi^2*yS^2,mi^2*yS^2],
D0[mi^2,0,0,0,mi^2*xs,mi^2-mi^2*xs-mi^2*xt,mi^2*yS^2,mi^2*yk^2,mi^2*yk^2,mi^2*yk^2],
D0[mi^2,0,0,0,mi^2*xt,mi^2-mi^2*xs-mi^2*xt,mi^2*yk^2,mi^2*yS^2,mi^2*yS^2,mi^2*yS^2],
D0[mi^2,0,0,0,mi^2*xt,mi^2-mi^2*xs-mi^2*xt,mi^2*yS^2,mi^2*yk^2,mi^2*yk^2,mi^2*yk^2]};


(* ::Text:: *)
(*(*Squared amplitude for the decay \[Tau]->c\[Gamma]\[Gamma] in terms of the Fi and Hi form factors  defined in Eq. (B1) of our manuscript. We also use the kinematic variables xt=t hat and xs=s hat defined in Sec. III of manuscript. *)
(* *)*)


 M2=(xs*xt*(-4*xs*xt*(-1 + xu)*Abs[F1]^2 - 2*xt^2*(-1 + xu)*xu*Abs[F2]^2 + 
   xs*xt*xu^2*Abs[F3]^2 - xs*xt*xu^3*Abs[F3]^2 + 4*xs^3*xt*Abs[F4]^2 - 
   8*xs^2*xt^2*Abs[F4]^2 + 4*xs*xt^3*Abs[F4]^2 + 4*xs*xt*xu*Abs[F4]^2 - 
   4*xs*xt*xu^2*Abs[F4]^2 + 2*xs^2*xt^2*xu*Abs[F5]^2 - 
   4*xs*xt^3*xu*Abs[F5]^2 + 2*xt^4*xu*Abs[F5]^2 + 2*xt^2*xu^2*Abs[F5]^2 - 
   2*xt^2*xu^3*Abs[F5]^2 + xs^3*xt*xu^2*Abs[F6]^2 - 
   2*xs^2*xt^2*xu^2*Abs[F6]^2 + xs*xt^3*xu^2*Abs[F6]^2 + 
   xs*xt*xu^3*Abs[F6]^2 - xs*xt*xu^4*Abs[F6]^2 + 4*xs*xt*Abs[FF1]^2 + 
   4*xs^3*xt*Abs[FF1]^2 - 8*xs^2*xt^2*Abs[FF1]^2 + 4*xs*xt^3*Abs[FF1]^2 - 
   4*xs*xt*xu^2*Abs[FF1]^2 + 2*xt^2*xu*Abs[FF2]^2 + 
   2*xs^2*xt^2*xu*Abs[FF2]^2 - 4*xs*xt^3*xu*Abs[FF2]^2 + 
   2*xt^4*xu*Abs[FF2]^2 - 2*xt^2*xu^3*Abs[FF2]^2 + xs*xt*xu^2*Abs[FF3]^2 + 
   xs^3*xt*xu^2*Abs[FF3]^2 - 2*xs^2*xt^2*xu^2*Abs[FF3]^2 + 
   xs*xt^3*xu^2*Abs[FF3]^2 - xs*xt*xu^4*Abs[FF3]^2 + 4*xs*xt*Abs[H1]^2 - 
   4*xs*xt*xu*Abs[H1]^2 + 2*xs^2*xu*Abs[H2]^2 - 2*xs^2*xu^2*Abs[H2]^2 + 
   xs*xt*xu^2*Abs[H3]^2 - xs*xt*xu^3*Abs[H3]^2 + 4*xs^3*xt*Abs[H4]^2 - 
   8*xs^2*xt^2*Abs[H4]^2 + 4*xs*xt^3*Abs[H4]^2 + 4*xs*xt*xu*Abs[H4]^2 - 
   4*xs*xt*xu^2*Abs[H4]^2 + 2*xs^4*xu*Abs[H5]^2 - 4*xs^3*xt*xu*Abs[H5]^2 + 
   2*xs^2*xt^2*xu*Abs[H5]^2 + 2*xs^2*xu^2*Abs[H5]^2 - 2*xs^2*xu^3*Abs[H5]^2 + 
   xs^3*xt*xu^2*Abs[H6]^2 - 2*xs^2*xt^2*xu^2*Abs[H6]^2 + 
   xs*xt^3*xu^2*Abs[H6]^2 + xs*xt*xu^3*Abs[H6]^2 - xs*xt*xu^4*Abs[H6]^2 + 
   4*xs*xt*Abs[HH1]^2 + 4*xs^3*xt*Abs[HH1]^2 - 8*xs^2*xt^2*Abs[HH1]^2 + 
   4*xs*xt^3*Abs[HH1]^2 - 4*xs*xt*xu^2*Abs[HH1]^2 + 2*xs^2*xu*Abs[HH2]^2 + 
   2*xs^4*xu*Abs[HH2]^2 - 4*xs^3*xt*xu*Abs[HH2]^2 + 
   2*xs^2*xt^2*xu*Abs[HH2]^2 - 2*xs^2*xu^3*Abs[HH2]^2 + 
   xs*xt*xu^2*Abs[HH3]^2 + xs^3*xt*xu^2*Abs[HH3]^2 - 
   2*xs^2*xt^2*xu^2*Abs[HH3]^2 + xs*xt^3*xu^2*Abs[HH3]^2 - 
   xs*xt*xu^4*Abs[HH3]^2 + 2*xs*xt*xu*Re[F1*Conjugate[F3]] - 
   2*xs*xt*xu^2*Re[F1*Conjugate[F3]] - 2*xs^2*xt*Re[F2*Conjugate[F4]] + 
   2*xs^3*xt*Re[F2*Conjugate[F4]] + 4*xs*xt^2*Re[F2*Conjugate[F4]] - 
   2*xs^2*xt^2*Re[F2*Conjugate[F4]] - 2*xt^3*Re[F2*Conjugate[F4]] - 
   2*xs*xt^3*Re[F2*Conjugate[F4]] + 2*xt^4*Re[F2*Conjugate[F4]] + 
   2*xs*xt*xu*Re[F2*Conjugate[F4]] + 2*xt^2*xu*Re[F2*Conjugate[F4]] - 
   2*xs*xt*xu^2*Re[F2*Conjugate[F4]] - 2*xt^2*xu^2*Re[F2*Conjugate[F4]] + 
   2*xs^2*xt*Re[F1*Conjugate[F5]] - 2*xs^3*xt*Re[F1*Conjugate[F5]] - 
   4*xs*xt^2*Re[F1*Conjugate[F5]] + 2*xs^2*xt^2*Re[F1*Conjugate[F5]] + 
   2*xt^3*Re[F1*Conjugate[F5]] + 2*xs*xt^3*Re[F1*Conjugate[F5]] - 
   2*xt^4*Re[F1*Conjugate[F5]] - 2*xs*xt*xu*Re[F1*Conjugate[F5]] - 
   2*xt^2*xu*Re[F1*Conjugate[F5]] + 2*xs*xt*xu^2*Re[F1*Conjugate[F5]] + 
   2*xt^2*xu^2*Re[F1*Conjugate[F5]] + xs^2*xt*xu*Re[F3*Conjugate[F5]] - 
   xs^3*xt*xu*Re[F3*Conjugate[F5]] - 2*xs*xt^2*xu*Re[F3*Conjugate[F5]] + 
   xs^2*xt^2*xu*Re[F3*Conjugate[F5]] + xt^3*xu*Re[F3*Conjugate[F5]] + 
   xs*xt^3*xu*Re[F3*Conjugate[F5]] - xt^4*xu*Re[F3*Conjugate[F5]] - 
   xs*xt*xu^2*Re[F3*Conjugate[F5]] - xt^2*xu^2*Re[F3*Conjugate[F5]] + 
   xs*xt*xu^3*Re[F3*Conjugate[F5]] + xt^2*xu^3*Re[F3*Conjugate[F5]] - 
   xs^2*xt*xu*Re[F2*Conjugate[F6]] + xs^3*xt*xu*Re[F2*Conjugate[F6]] + 
   2*xs*xt^2*xu*Re[F2*Conjugate[F6]] - xs^2*xt^2*xu*Re[F2*Conjugate[F6]] - 
   xt^3*xu*Re[F2*Conjugate[F6]] - xs*xt^3*xu*Re[F2*Conjugate[F6]] + 
   xt^4*xu*Re[F2*Conjugate[F6]] + xs*xt*xu^2*Re[F2*Conjugate[F6]] + 
   xt^2*xu^2*Re[F2*Conjugate[F6]] - xs*xt*xu^3*Re[F2*Conjugate[F6]] - 
   xt^2*xu^3*Re[F2*Conjugate[F6]] + 2*xs^3*xt*xu*Re[F4*Conjugate[F6]] - 
   4*xs^2*xt^2*xu*Re[F4*Conjugate[F6]] + 2*xs*xt^3*xu*Re[F4*Conjugate[F6]] + 
   2*xs*xt*xu^2*Re[F4*Conjugate[F6]] - 2*xs*xt*xu^3*Re[F4*Conjugate[F6]] - 
   8*xs^2*xt*Re[F1*Conjugate[FF1]] + 8*xs*xt^2*Re[F1*Conjugate[FF1]] - 
   2*xs^2*xt*xu*Re[F3*Conjugate[FF1]] + 2*xs*xt^2*xu*Re[F3*Conjugate[FF1]] - 
   8*xs^2*xt*Re[F4*Conjugate[FF1]] + 8*xs*xt^2*Re[F4*Conjugate[FF1]] - 
   2*xs^2*xt*xu*Re[F6*Conjugate[FF1]] + 2*xs*xt^2*xu*Re[F6*Conjugate[FF1]] - 
   4*xs*xt^2*xu*Re[F2*Conjugate[FF2]] + 4*xt^3*xu*Re[F2*Conjugate[FF2]] - 
   4*xs*xt^2*xu*Re[F5*Conjugate[FF2]] + 4*xt^3*xu*Re[F5*Conjugate[FF2]] - 
   2*xs^2*xt*xu*Re[F1*Conjugate[FF3]] + 2*xs*xt^2*xu*Re[F1*Conjugate[FF3]] - 
   2*xs^2*xt*xu^2*Re[F3*Conjugate[FF3]] + 2*xs*xt^2*xu^2*
    Re[F3*Conjugate[FF3]] - 2*xs^2*xt*xu*Re[F4*Conjugate[FF3]] + 
   2*xs*xt^2*xu*Re[F4*Conjugate[FF3]] - 2*xs^2*xt*xu^2*
    Re[F6*Conjugate[FF3]] + 2*xs*xt^2*xu^2*Re[F6*Conjugate[FF3]] + 
   2*xs*xt*xu*Re[FF1*Conjugate[FF3]] + 2*xs^3*xt*xu*Re[FF1*Conjugate[FF3]] - 
   4*xs^2*xt^2*xu*Re[FF1*Conjugate[FF3]] + 
   2*xs*xt^3*xu*Re[FF1*Conjugate[FF3]] - 
   2*xs*xt*xu^3*Re[FF1*Conjugate[FF3]] + 2*xs*xt*xu*Re[F3*Conjugate[H1]] - 
   2*xs*xt*xu^2*Re[F3*Conjugate[H1]] - 2*xs^2*xt*xu*Re[FF3*Conjugate[H1]] + 
   2*xs*xt^2*xu*Re[FF3*Conjugate[H1]] - 2*xs*xt*xu*Re[F2*Conjugate[H2]] + 
   2*xs*xt*xu^2*Re[F2*Conjugate[H2]] + xs^3*xu*Re[F6*Conjugate[H2]] - 
   xs^4*xu*Re[F6*Conjugate[H2]] - 2*xs^2*xt*xu*Re[F6*Conjugate[H2]] + 
   xs^3*xt*xu*Re[F6*Conjugate[H2]] + xs*xt^2*xu*Re[F6*Conjugate[H2]] + 
   xs^2*xt^2*xu*Re[F6*Conjugate[H2]] - xs*xt^3*xu*Re[F6*Conjugate[H2]] - 
   xs^2*xu^2*Re[F6*Conjugate[H2]] - xs*xt*xu^2*Re[F6*Conjugate[H2]] + 
   xs^2*xu^3*Re[F6*Conjugate[H2]] + xs*xt*xu^3*Re[F6*Conjugate[H2]] + 
   2*xs^2*xt*xu*Re[FF2*Conjugate[H2]] - 2*xs*xt^2*xu*Re[FF2*Conjugate[H2]] + 
   2*xs*xt*xu*Re[F1*Conjugate[H3]] - 2*xs*xt*xu^2*Re[F1*Conjugate[H3]] + 
   2*xs*xt*xu^2*Re[F3*Conjugate[H3]] - 2*xs*xt*xu^3*Re[F3*Conjugate[H3]] + 
   xs^2*xt*xu*Re[F5*Conjugate[H3]] - xs^3*xt*xu*Re[F5*Conjugate[H3]] - 
   2*xs*xt^2*xu*Re[F5*Conjugate[H3]] + xs^2*xt^2*xu*Re[F5*Conjugate[H3]] + 
   xt^3*xu*Re[F5*Conjugate[H3]] + xs*xt^3*xu*Re[F5*Conjugate[H3]] - 
   xt^4*xu*Re[F5*Conjugate[H3]] - xs*xt*xu^2*Re[F5*Conjugate[H3]] - 
   xt^2*xu^2*Re[F5*Conjugate[H3]] + xs*xt*xu^3*Re[F5*Conjugate[H3]] + 
   xt^2*xu^3*Re[F5*Conjugate[H3]] - 2*xs^2*xt*xu*Re[FF1*Conjugate[H3]] + 
   2*xs*xt^2*xu*Re[FF1*Conjugate[H3]] - 2*xs^2*xt*xu^2*
    Re[FF3*Conjugate[H3]] + 2*xs*xt^2*xu^2*Re[FF3*Conjugate[H3]] + 
   2*xs*xt*xu*Re[H1*Conjugate[H3]] - 2*xs*xt*xu^2*Re[H1*Conjugate[H3]] - 
   2*xs^3*xt*xu*Re[F6*Conjugate[H4]] + 4*xs^2*xt^2*xu*Re[F6*Conjugate[H4]] - 
   2*xs*xt^3*xu*Re[F6*Conjugate[H4]] - 2*xs*xt*xu^2*Re[F6*Conjugate[H4]] + 
   2*xs*xt*xu^3*Re[F6*Conjugate[H4]] + 2*xs^2*xt*xu*Re[FF3*Conjugate[H4]] - 
   2*xs*xt^2*xu*Re[FF3*Conjugate[H4]] - 2*xs^3*Re[H2*Conjugate[H4]] + 
   2*xs^4*Re[H2*Conjugate[H4]] + 4*xs^2*xt*Re[H2*Conjugate[H4]] - 
   2*xs^3*xt*Re[H2*Conjugate[H4]] - 2*xs*xt^2*Re[H2*Conjugate[H4]] - 
   2*xs^2*xt^2*Re[H2*Conjugate[H4]] + 2*xs*xt^3*Re[H2*Conjugate[H4]] + 
   2*xs^2*xu*Re[H2*Conjugate[H4]] + 2*xs*xt*xu*Re[H2*Conjugate[H4]] - 
   2*xs^2*xu^2*Re[H2*Conjugate[H4]] - 2*xs*xt*xu^2*Re[H2*Conjugate[H4]] + 
   xs^3*xu*Re[F3*Conjugate[H5]] - xs^4*xu*Re[F3*Conjugate[H5]] - 
   2*xs^2*xt*xu*Re[F3*Conjugate[H5]] + xs^3*xt*xu*Re[F3*Conjugate[H5]] + 
   xs*xt^2*xu*Re[F3*Conjugate[H5]] + xs^2*xt^2*xu*Re[F3*Conjugate[H5]] - 
   xs*xt^3*xu*Re[F3*Conjugate[H5]] - xs^2*xu^2*Re[F3*Conjugate[H5]] - 
   xs*xt*xu^2*Re[F3*Conjugate[H5]] + xs^2*xu^3*Re[F3*Conjugate[H5]] + 
   xs*xt*xu^3*Re[F3*Conjugate[H5]] + 2*xs^3*xt*xu*Re[F5*Conjugate[H5]] - 
   4*xs^2*xt^2*xu*Re[F5*Conjugate[H5]] + 2*xs*xt^3*xu*Re[F5*Conjugate[H5]] + 
   2*xs*xt*xu^2*Re[F5*Conjugate[H5]] - 2*xs*xt*xu^3*Re[F5*Conjugate[H5]] - 
   2*xs^2*xt*xu*Re[FF2*Conjugate[H5]] + 2*xs*xt^2*xu*Re[FF2*Conjugate[H5]] + 
   2*xs^3*Re[H1*Conjugate[H5]] - 2*xs^4*Re[H1*Conjugate[H5]] - 
   4*xs^2*xt*Re[H1*Conjugate[H5]] + 2*xs^3*xt*Re[H1*Conjugate[H5]] + 
   2*xs*xt^2*Re[H1*Conjugate[H5]] + 2*xs^2*xt^2*Re[H1*Conjugate[H5]] - 
   2*xs*xt^3*Re[H1*Conjugate[H5]] - 2*xs^2*xu*Re[H1*Conjugate[H5]] - 
   2*xs*xt*xu*Re[H1*Conjugate[H5]] + 2*xs^2*xu^2*Re[H1*Conjugate[H5]] + 
   2*xs*xt*xu^2*Re[H1*Conjugate[H5]] + xs^3*xu*Re[H3*Conjugate[H5]] - 
   xs^4*xu*Re[H3*Conjugate[H5]] - 2*xs^2*xt*xu*Re[H3*Conjugate[H5]] + 
   xs^3*xt*xu*Re[H3*Conjugate[H5]] + xs*xt^2*xu*Re[H3*Conjugate[H5]] + 
   xs^2*xt^2*xu*Re[H3*Conjugate[H5]] - xs*xt^3*xu*Re[H3*Conjugate[H5]] - 
   xs^2*xu^2*Re[H3*Conjugate[H5]] - xs*xt*xu^2*Re[H3*Conjugate[H5]] + 
   xs^2*xu^3*Re[H3*Conjugate[H5]] + xs*xt*xu^3*Re[H3*Conjugate[H5]] + 
   xs^2*xt*xu*Re[F2*Conjugate[H6]] - xs^3*xt*xu*Re[F2*Conjugate[H6]] - 
   2*xs*xt^2*xu*Re[F2*Conjugate[H6]] + xs^2*xt^2*xu*Re[F2*Conjugate[H6]] + 
   xt^3*xu*Re[F2*Conjugate[H6]] + xs*xt^3*xu*Re[F2*Conjugate[H6]] - 
   xt^4*xu*Re[F2*Conjugate[H6]] - xs*xt*xu^2*Re[F2*Conjugate[H6]] - 
   xt^2*xu^2*Re[F2*Conjugate[H6]] + xs*xt*xu^3*Re[F2*Conjugate[H6]] + 
   xt^2*xu^3*Re[F2*Conjugate[H6]] - 2*xs^3*xt*xu*Re[F4*Conjugate[H6]] + 
   4*xs^2*xt^2*xu*Re[F4*Conjugate[H6]] - 2*xs*xt^3*xu*Re[F4*Conjugate[H6]] - 
   2*xs*xt*xu^2*Re[F4*Conjugate[H6]] + 2*xs*xt*xu^3*Re[F4*Conjugate[H6]] - 
   2*xs^3*xt*xu^2*Re[F6*Conjugate[H6]] + 4*xs^2*xt^2*xu^2*
    Re[F6*Conjugate[H6]] - 2*xs*xt^3*xu^2*Re[F6*Conjugate[H6]] - 
   2*xs*xt*xu^3*Re[F6*Conjugate[H6]] + 2*xs*xt*xu^4*Re[F6*Conjugate[H6]] + 
   2*xs^2*xt*xu*Re[FF1*Conjugate[H6]] - 2*xs*xt^2*xu*Re[FF1*Conjugate[H6]] + 
   2*xs^2*xt*xu^2*Re[FF3*Conjugate[H6]] - 2*xs*xt^2*xu^2*
    Re[FF3*Conjugate[H6]] - xs^3*xu*Re[H2*Conjugate[H6]] + 
   xs^4*xu*Re[H2*Conjugate[H6]] + 2*xs^2*xt*xu*Re[H2*Conjugate[H6]] - 
   xs^3*xt*xu*Re[H2*Conjugate[H6]] - xs*xt^2*xu*Re[H2*Conjugate[H6]] - 
   xs^2*xt^2*xu*Re[H2*Conjugate[H6]] + xs*xt^3*xu*Re[H2*Conjugate[H6]] + 
   xs^2*xu^2*Re[H2*Conjugate[H6]] + xs*xt*xu^2*Re[H2*Conjugate[H6]] - 
   xs^2*xu^3*Re[H2*Conjugate[H6]] - xs*xt*xu^3*Re[H2*Conjugate[H6]] + 
   2*xs^3*xt*xu*Re[H4*Conjugate[H6]] - 4*xs^2*xt^2*xu*Re[H4*Conjugate[H6]] + 
   2*xs*xt^3*xu*Re[H4*Conjugate[H6]] + 2*xs*xt*xu^2*Re[H4*Conjugate[H6]] - 
   2*xs*xt*xu^3*Re[H4*Conjugate[H6]] + 2*xs^2*xt*xu*Re[F3*Conjugate[HH1]] - 
   2*xs*xt^2*xu*Re[F3*Conjugate[HH1]] - 2*xs^2*xt*xu*Re[F6*Conjugate[HH1]] + 
   2*xs*xt^2*xu*Re[F6*Conjugate[HH1]] + 2*xs*xt*xu*Re[FF3*Conjugate[HH1]] - 
   2*xs^3*xt*xu*Re[FF3*Conjugate[HH1]] + 4*xs^2*xt^2*xu*
    Re[FF3*Conjugate[HH1]] - 2*xs*xt^3*xu*Re[FF3*Conjugate[HH1]] - 
   4*xs*xt*xu^2*Re[FF3*Conjugate[HH1]] + 
   2*xs*xt*xu^3*Re[FF3*Conjugate[HH1]] + 8*xs^2*xt*Re[H1*Conjugate[HH1]] - 
   8*xs*xt^2*Re[H1*Conjugate[HH1]] + 2*xs^2*xt*xu*Re[H3*Conjugate[HH1]] - 
   2*xs*xt^2*xu*Re[H3*Conjugate[HH1]] + 8*xs^2*xt*Re[H4*Conjugate[HH1]] - 
   8*xs*xt^2*Re[H4*Conjugate[HH1]] + 2*xs^2*xt*xu*Re[H6*Conjugate[HH1]] - 
   2*xs*xt^2*xu*Re[H6*Conjugate[HH1]] - 2*xs^2*xt*xu*Re[F2*Conjugate[HH2]] + 
   2*xs*xt^2*xu*Re[F2*Conjugate[HH2]] + 2*xs^2*xt*xu*Re[F5*Conjugate[HH2]] - 
   2*xs*xt^2*xu*Re[F5*Conjugate[HH2]] - 2*xs*xt*xu*Re[FF2*Conjugate[HH2]] + 
   2*xs^3*xt*xu*Re[FF2*Conjugate[HH2]] - 4*xs^2*xt^2*xu*
    Re[FF2*Conjugate[HH2]] + 2*xs*xt^3*xu*Re[FF2*Conjugate[HH2]] + 
   4*xs*xt*xu^2*Re[FF2*Conjugate[HH2]] - 
   2*xs*xt*xu^3*Re[FF2*Conjugate[HH2]] + 2*xs^3*xu*Re[FF3*Conjugate[HH2]] - 
   2*xs^4*xu*Re[FF3*Conjugate[HH2]] - 4*xs^2*xt*xu*Re[FF3*Conjugate[HH2]] + 
   2*xs^3*xt*xu*Re[FF3*Conjugate[HH2]] + 
   2*xs*xt^2*xu*Re[FF3*Conjugate[HH2]] + 2*xs^2*xt^2*xu*
    Re[FF3*Conjugate[HH2]] - 2*xs*xt^3*xu*Re[FF3*Conjugate[HH2]] - 
   2*xs^2*xu^2*Re[FF3*Conjugate[HH2]] - 2*xs*xt*xu^2*Re[FF3*Conjugate[HH2]] + 
   2*xs^2*xu^3*Re[FF3*Conjugate[HH2]] + 2*xs*xt*xu^3*Re[FF3*Conjugate[HH2]] + 
   4*xs^3*xu*Re[H2*Conjugate[HH2]] - 4*xs^2*xt*xu*Re[H2*Conjugate[HH2]] + 
   4*xs^3*xu*Re[H5*Conjugate[HH2]] - 4*xs^2*xt*xu*Re[H5*Conjugate[HH2]] + 
   2*xs^2*xt*xu*Re[F1*Conjugate[HH3]] - 2*xs*xt^2*xu*Re[F1*Conjugate[HH3]] + 
   2*xs^2*xt*xu^2*Re[F3*Conjugate[HH3]] - 2*xs*xt^2*xu^2*
    Re[F3*Conjugate[HH3]] - 2*xs^2*xt*xu*Re[F4*Conjugate[HH3]] + 
   2*xs*xt^2*xu*Re[F4*Conjugate[HH3]] - 2*xs^2*xt*xu^2*
    Re[F6*Conjugate[HH3]] + 2*xs*xt^2*xu^2*Re[F6*Conjugate[HH3]] + 
   2*xs*xt*xu*Re[FF1*Conjugate[HH3]] - 2*xs^3*xt*xu*Re[FF1*Conjugate[HH3]] + 
   4*xs^2*xt^2*xu*Re[FF1*Conjugate[HH3]] - 
   2*xs*xt^3*xu*Re[FF1*Conjugate[HH3]] - 
   4*xs*xt*xu^2*Re[FF1*Conjugate[HH3]] + 
   2*xs*xt*xu^3*Re[FF1*Conjugate[HH3]] + 
   2*xs^2*xt*xu*Re[FF2*Conjugate[HH3]] - 
   2*xs^3*xt*xu*Re[FF2*Conjugate[HH3]] - 
   4*xs*xt^2*xu*Re[FF2*Conjugate[HH3]] + 2*xs^2*xt^2*xu*
    Re[FF2*Conjugate[HH3]] + 2*xt^3*xu*Re[FF2*Conjugate[HH3]] + 
   2*xs*xt^3*xu*Re[FF2*Conjugate[HH3]] - 2*xt^4*xu*Re[FF2*Conjugate[HH3]] - 
   2*xs*xt*xu^2*Re[FF2*Conjugate[HH3]] - 2*xt^2*xu^2*Re[FF2*Conjugate[HH3]] + 
   2*xs*xt*xu^3*Re[FF2*Conjugate[HH3]] + 2*xt^2*xu^3*Re[FF2*Conjugate[HH3]] + 
   2*xs*xt*xu^2*Re[FF3*Conjugate[HH3]] - 2*xs^3*xt*xu^2*
    Re[FF3*Conjugate[HH3]] + 4*xs^2*xt^2*xu^2*Re[FF3*Conjugate[HH3]] - 
   2*xs*xt^3*xu^2*Re[FF3*Conjugate[HH3]] - 
   4*xs*xt*xu^3*Re[FF3*Conjugate[HH3]] + 
   2*xs*xt*xu^4*Re[FF3*Conjugate[HH3]] + 2*xs^2*xt*xu*Re[H1*Conjugate[HH3]] - 
   2*xs*xt^2*xu*Re[H1*Conjugate[HH3]] + 2*xs^2*xt*xu^2*
    Re[H3*Conjugate[HH3]] - 2*xs*xt^2*xu^2*Re[H3*Conjugate[HH3]] + 
   2*xs^2*xt*xu*Re[H4*Conjugate[HH3]] - 2*xs*xt^2*xu*Re[H4*Conjugate[HH3]] + 
   2*xs^2*xt*xu^2*Re[H6*Conjugate[HH3]] - 2*xs*xt^2*xu^2*
    Re[H6*Conjugate[HH3]] + 2*xs*xt*xu*Re[HH1*Conjugate[HH3]] + 
   2*xs^3*xt*xu*Re[HH1*Conjugate[HH3]] - 4*xs^2*xt^2*xu*
    Re[HH1*Conjugate[HH3]] + 2*xs*xt^3*xu*Re[HH1*Conjugate[HH3]] - 
   2*xs*xt*xu^3*Re[HH1*Conjugate[HH3]]))/32;


 coef1={(Qj*\[Lambda]Ljk*(2*Qk*\[Lambda]Lik + Qj*(-1 + (-3 + xt)*yk^2 - 
       (-3 + xt)*yS^2)*\[Lambda]Lik - 2*Qj*xt*yk*\[Lambda]Rik))/
   (xs*(-1 + xt)*xt^2) + (Qj^2*yk^2*\[Lambda]Ljk*
    ((-3 + xt)*yk^2*\[Lambda]Lik - (-3 + xt)*yS^2*\[Lambda]Lik - 
     2*xt*yk*\[Lambda]Rik)*BB[1])/(xs*(-1 + xt)*xt^2*(yk^2 - yS^2)) + 
  (Qj^2*yS^2*\[Lambda]Ljk*(-((-3 + xt)*yk^2*\[Lambda]Lik) + 
     (-3 + xt)*yS^2*\[Lambda]Lik + 2*xt*yk*\[Lambda]Rik)*BB[2])/
   (xs*(-1 + xt)*xt^2*(yk^2 - yS^2)) - 
  (Qj*\[Lambda]Ljk*(-2*Qk*(-2*xs*xt + 2*xs^3*xt + xt^2 + xs^2*(1 - 2*xt^2))*
      \[Lambda]Lik + Qj*(xt^2*(1 - yk^2 + yS^2)*\[Lambda]Lik + 
       xs*xt*((-2 + xt*(1 + yk^2 - yS^2))*\[Lambda]Lik - 
         2*(-1 + xt)*yk*\[Lambda]Rik) + 
       xs^2*((1 + yk^2 - yS^2 + xt^2*(-3 + yk^2 - yS^2) + 
           xt*(1 - yk^2 + yS^2))*\[Lambda]Lik + 2*(-1 + xt)^2*yk*
          \[Lambda]Rik) + xs^3*(xt*(1 + yk^2 - yS^2)*\[Lambda]Lik + 
         2*xt*yk*\[Lambda]Rik - 2*(yk^2*\[Lambda]Lik - yS^2*\[Lambda]Lik + 
           yk*\[Lambda]Rik))))*BB[3])/((-1 + xs)*xs^2*(-1 + xt)*xt^2*
    (xs + xt)) + (Qj*(2*Qk*xs - Qj*(xs - yk^2 + yS^2))*\[Lambda]Lik*
    \[Lambda]Ljk*BB[4])/((-1 + xs)*xs*xt^2) + 
  (Qj*\[Lambda]Ljk*(2*Qk*xt*(-2*xs + xt)*\[Lambda]Lik - 
     Qj*xt*(xt - yk^2 + yS^2)*\[Lambda]Lik + 
     Qj*xs*(yk^2*\[Lambda]Lik - yS^2*\[Lambda]Lik - 2*yk*\[Lambda]Rik + 
       xt*(\[Lambda]Lik + 4*yk*\[Lambda]Rik)))*BB[5])/(xs^2*(-1 + xt)*xt^2) - 
  (Qk^2*(xs - xt)^2*\[Lambda]Lik*\[Lambda]Ljk*BB[6])/(xs^2*xt^2*(xs + xt)) + 
  ((Qj - Qk)^2*(xs - xt)^2*\[Lambda]Lik*\[Lambda]Ljk*BB[7])/
   (xs^2*xt^2*(xs + xt)) - (mi^2*Qk*\[Lambda]Ljk*
    (Qk*(-1 + xs)*xs^2*(xs^2*\[Lambda]Lik + xt*(-1 + 2*yk^2 - 2*yS^2)*
        \[Lambda]Lik + xs*(-1 + 3*xt + 2*yk^2 - 2*yS^2)*\[Lambda]Lik + 
       2*(-yk^2 + yS^2)*\[Lambda]Lik - 4*xt*yk*\[Lambda]Rik + 
       4*xs*xt*yk*\[Lambda]Rik + 4*xt^2*yk*\[Lambda]Rik) + 
     Qj*((-1 + xt)^2*xt*(yk^2 - yS^2)^2*\[Lambda]Lik + 
       xs*(-1 + xt)*(-((yk^2 - yS^2)^2*\[Lambda]Lik) + xt*(yk^2 - yS^2)*
          ((-2 + 3*yk^2 - 3*yS^2)*\[Lambda]Lik - 4*yk*\[Lambda]Rik) + 
         xt^2*(4*yk^2*\[Lambda]Lik - 2*yS^2*\[Lambda]Lik + 
           4*yk^3*\[Lambda]Rik - 4*yk*yS^2*\[Lambda]Rik)) + 
       xs^3*((yk^2 - yS^2)^2*\[Lambda]Lik + xt^2*(3*\[Lambda]Lik + 
           4*yk*\[Lambda]Rik) + 2*xt*(yk^2*\[Lambda]Lik - 
           2*yS^2*\[Lambda]Lik + 2*yk^3*\[Lambda]Rik - 
           2*yk*yS^2*\[Lambda]Rik)) + xs^2*(-2*(yk^2 - yS^2)^2*\[Lambda]Lik + 
         xt^3*(\[Lambda]Lik + 4*yk*\[Lambda]Rik) + 
         xt*(3*yk^4*\[Lambda]Lik + 3*yS^2*(2 + yS^2)*\[Lambda]Lik - 
           2*yk^2*(2 + 3*yS^2)*\[Lambda]Lik - 8*yk^3*\[Lambda]Rik + 
           8*yk*yS^2*\[Lambda]Rik) + xt^2*((-2 + 6*yk^2 - 6*yS^2)*
            \[Lambda]Lik + 4*yk*(-1 + 2*yk^2 - 2*yS^2)*\[Lambda]Rik))))*
    CC[1])/(2*xs^2*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*(Qj - Qk)*\[Lambda]Ljk*
    (-(Qk*(-1 + xs)*xs^2*(xs^2*\[Lambda]Lik + 2*(yk^2 - yS^2)*\[Lambda]Lik + 
        xs*(-1 + xt - 2*yk^2 + 2*yS^2)*\[Lambda]Lik - 
        4*xs*xt*yk*\[Lambda]Rik - 2*xt^2*(\[Lambda]Lik + 2*yk*\[Lambda]Rik) + 
        xt*(\[Lambda]Lik - 2*yk^2*\[Lambda]Lik + 2*yS^2*\[Lambda]Lik + 
          4*yk*\[Lambda]Rik))) + Qj*(-1 + xs + xt)*(xs^4*\[Lambda]Lik + 
       (-1 + xt)*xt*(yk^2 - yS^2)^2*\[Lambda]Lik - 
       xs^3*(\[Lambda]Lik + 2*yk^2*\[Lambda]Lik - 2*yS^2*\[Lambda]Lik + 
         4*xt*yk*\[Lambda]Rik) + xs*((-1 + 2*xt)*yk^4*\[Lambda]Lik + 
         2*(1 - 2*xt)*yk^2*yS^2*\[Lambda]Lik + 
         yS^2*(2*xt^2 - yS^2 + 2*xt*yS^2)*\[Lambda]Lik + 
         4*(-1 + xt)*xt*yk^3*\[Lambda]Rik - 4*(-1 + xt)*xt*yk*yS^2*
          \[Lambda]Rik) + xs^2*((yk^4 + yS^2*(-2 + yS^2) - 
           2*yk^2*(-1 + yS^2))*\[Lambda]Lik - xt^2*(\[Lambda]Lik + 
           4*yk*\[Lambda]Rik) + xt*(\[Lambda]Lik - 2*yk^2*\[Lambda]Lik + 
           4*yk*(1 + yk^2 - yS^2)*\[Lambda]Rik))))*CC[2])/
   (2*xs^2*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*Qk*\[Lambda]Ljk*
    (-(Qk*(-1 + xt)*xt*((-1 + xt)^2*xt*(xt + 2*yk^2 - 2*yS^2)*\[Lambda]Lik + 
        2*xs^3*(yk^2 - yS^2)*\[Lambda]Lik - 2*xs^2*(-1 + xt)*
         (2*(-yk^2 + yS^2)*\[Lambda]Lik + xt*(\[Lambda]Lik - 
            2*yk*\[Lambda]Rik)) + xs*(-1 + xt)*
         (2*(-yk^2 + yS^2)*\[Lambda]Lik + xt*(\[Lambda]Lik + 
            4*yk^2*\[Lambda]Lik - 4*yS^2*\[Lambda]Lik - 4*yk*\[Lambda]Rik) + 
          xt^2*(\[Lambda]Lik + 4*yk*\[Lambda]Rik)))) + 
     Qj*(-((-1 + xt)^3*xt*(yk^2 - yS^2)^2*\[Lambda]Lik) + 
       4*xs^4*xt*yk*(yk*\[Lambda]Lik + (-1 + 2*xt)*\[Lambda]Rik) - 
       xs*(-1 + xt)^2*((-1 + 3*xt)*yk^4*\[Lambda]Lik + 2*(1 - 3*xt)*yk^2*yS^2*
          \[Lambda]Lik + yS^2*(-2*xt^2 - yS^2 + 3*xt*yS^2)*\[Lambda]Lik + 
         4*(-1 + xt)*xt*yk^3*\[Lambda]Rik - 4*(-1 + xt)*xt*yk*yS^2*
          \[Lambda]Rik) + xs^3*(-1 + xt)*(-((yk^2 - yS^2)^2*\[Lambda]Lik) + 
         xt^2*(\[Lambda]Lik + 12*yk*\[Lambda]Rik) + 
         2*xt*yk*(3*yk*\[Lambda]Lik - 2*yk^2*\[Lambda]Rik + 
           2*(-2 + yS^2)*\[Lambda]Rik)) - xs^2*(-1 + xt)*
        (-2*(yk^2 - yS^2)^2*\[Lambda]Lik + xt^3*(\[Lambda]Lik - 
           4*yk*\[Lambda]Rik) + xt^2*(-2*yk^2*\[Lambda]Lik - 
           2*yS^2*\[Lambda]Lik + 8*yk^3*\[Lambda]Rik - 8*yk*(-1 + yS^2)*
            \[Lambda]Rik) + xt*(3*yk^4*\[Lambda]Lik + 3*yS^4*\[Lambda]Lik + 
           yk^2*(2 - 6*yS^2)*\[Lambda]Lik - 8*yk^3*\[Lambda]Rik + 
           4*yk*(-1 + 2*yS^2)*\[Lambda]Rik))))*CC[3])/
   (2*xs^3*(-1 + xt)*xt^2*(-1 + xs + xt)^2) - 
  (mi^2*(Qj - Qk)*\[Lambda]Ljk*
    (Qk*(-1 + xt)*xt*(2*xs^3*(xt - yk^2 + yS^2)*\[Lambda]Lik + 
       (-1 + xt)^2*xt*(xt - 2*yk^2 + 2*yS^2)*\[Lambda]Lik - 
       4*xs^2*(-1 + xt)*(yk^2*\[Lambda]Lik - yS^2*\[Lambda]Lik + 
         xt*yk*\[Lambda]Rik) + xs*(-1 + xt)*(2*(yk^2 - yS^2)*\[Lambda]Lik + 
         xt^2*(\[Lambda]Lik - 4*yk*\[Lambda]Rik) + 
         xt*(\[Lambda]Lik - 4*yk^2*\[Lambda]Lik + 4*yS^2*\[Lambda]Lik + 
           4*yk*\[Lambda]Rik))) + Qj*(-1 + xs + xt)*
      (4*xs^3*xt*yS^2*\[Lambda]Lik - 
       xt*(xt^2 + yk^2 - yS^2 + xt*(-1 - yk^2 + yS^2))^2*\[Lambda]Lik + 
       xs^2*(-1 + xt)*(-((yk^2 - yS^2)^2*\[Lambda]Lik) - 
         4*xt*yk^3*\[Lambda]Rik + 2*xt*yS^2*(\[Lambda]Lik + 
           2*yk*\[Lambda]Rik) + xt^2*(\[Lambda]Lik + 4*yk*\[Lambda]Rik)) + 
       xs*(-1 + xt)*((yk^2 - yS^2)^2*\[Lambda]Lik + 4*xt^3*yk*\[Lambda]Rik - 
         2*xt*(yk^2 - yS^2)*(yk^2*\[Lambda]Lik - yS^2*\[Lambda]Lik - 
           2*yk*\[Lambda]Rik) + xt^2*((-1 + 2*yk^2)*\[Lambda]Lik - 
           4*yk*(1 + yk^2 - yS^2)*\[Lambda]Rik))))*CC[4])/
   (2*xs^3*(-1 + xt)*xt^2*(-1 + xs + xt)^2) + 
  (mi^2*Qk^2*\[Lambda]Ljk*(xs^4*\[Lambda]Lik + 
     xt*(xt^3 + xt^2*(-1 + 2*yk^2 - 2*yS^2) - 2*(yk^2 - yS^2)^2 + 
       2*xt*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2)))*\[Lambda]Lik + 
     xs^3*((-1 + 2*xt + 2*yk^2 - 2*yS^2)*\[Lambda]Lik + 
       4*xt*yk*\[Lambda]Rik) + 
     xs^2*(2*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2))*\[Lambda]Lik + 
       8*xt^2*yk*\[Lambda]Rik + xt*((-1 + 2*yk^2 - 6*yS^2)*\[Lambda]Lik + 
         4*yk*(-1 + 2*yk^2 - 2*yS^2)*\[Lambda]Rik)) + 
     xs*(-2*(yk^2 - yS^2)^2*\[Lambda]Lik + 4*xt^3*yk*\[Lambda]Rik + 
       4*xt*(yk^2 - yS^2)*((-1 + yk^2 - yS^2)*\[Lambda]Lik - 
         2*yk*\[Lambda]Rik) + xt^2*((1 + 2*yk^2 - 6*yS^2)*\[Lambda]Lik + 
         4*yk*(-1 + 2*yk^2 - 2*yS^2)*\[Lambda]Rik)))*CC[5])/(2*xs^3*xt^3) - 
  (mi^2*(Qj - Qk)^2*\[Lambda]Ljk*(xs^4*\[Lambda]Lik + 
     xt*(xt^3 - 2*(yk^2 - yS^2)^2 + xt^2*(-1 - 2*yk^2 + 2*yS^2) + 
       2*xt*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4))*\[Lambda]Lik - 
     xs^3*(\[Lambda]Lik + 2*yk^2*\[Lambda]Lik - 2*yS^2*\[Lambda]Lik + 
       4*xt*yk*\[Lambda]Rik) + xs^2*(xt*(1 - 2*yk^2 - 2*yS^2)*\[Lambda]Lik + 
       2*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4)*\[Lambda]Lik + 
       4*xt*yk*(1 + 2*yk^2 - 2*yS^2)*\[Lambda]Rik - 
       2*xt^2*(\[Lambda]Lik + 4*yk*\[Lambda]Rik)) + 
     xs*(-2*(yk^2 - yS^2)^2*\[Lambda]Lik - 4*xt^3*yk*\[Lambda]Rik + 
       4*xt*(yk^2 - yS^2)*(yk^2*\[Lambda]Lik - yS^2*\[Lambda]Lik - 
         2*yk*\[Lambda]Rik) + xt^2*((1 - 2*yk^2 - 2*yS^2)*\[Lambda]Lik + 
         4*yk*(1 + 2*yk^2 - 2*yS^2)*\[Lambda]Rik)))*CC[6])/(2*xs^3*xt^3) + 
  (mi^2*(Qj - Qk)*\[Lambda]Ljk*
    (-(Qk*(-1 + xs)*xs*(xs^4*\[Lambda]Lik + 2*(-1 + xt)*xt*(yk^2 - yS^2)*
         \[Lambda]Lik + xs^3*((-2 + xt - 2*yk^2 + 2*yS^2)*\[Lambda]Lik - 
          4*xt*yk*\[Lambda]Rik) + xs^2*((1 + 2*xt^2 + 4*yk^2 - 6*xt*yk^2 - 
            4*yS^2 + 6*xt*yS^2)*\[Lambda]Lik - 4*(-2 + xt)*xt*yk*
           \[Lambda]Rik) - xs*(2*(yk^2 - yS^2)*\[Lambda]Lik + 
          4*xt^2*(yk^2*\[Lambda]Lik - yS^2*\[Lambda]Lik - yk*\[Lambda]Rik) + 
          xt*(\[Lambda]Lik - 8*yk^2*\[Lambda]Lik + 8*yS^2*\[Lambda]Lik + 
            4*yk*\[Lambda]Rik)))) + Qj*(-1 + xs + xt)*
      (xs^5*\[Lambda]Lik - (-1 + xt)*xt*(yk^2 - yS^2)^2*\[Lambda]Lik - 
       2*xs^4*((1 + yk^2 - yS^2)*\[Lambda]Lik + 2*xt*yk*\[Lambda]Rik) + 
       xs*((1 - 3*xt + xt^2)*yk^4*\[Lambda]Lik - 2*(1 - 3*xt + xt^2)*yk^2*
          yS^2*\[Lambda]Lik + yS^2*(yS^2 - 3*xt*yS^2 + xt^2*(2 + yS^2))*
          \[Lambda]Lik - 4*(-1 + xt)*xt*yk^3*\[Lambda]Rik + 
         4*(-1 + xt)*xt*yk*yS^2*\[Lambda]Rik) + 
       xs^3*((1 + xt - xt^2 - 2*xt*yk^2 + yk^4 - 4*yS^2 + yS^4 - 
           2*yk^2*(-2 + yS^2))*\[Lambda]Lik - 4*xt*yk*(-2 + xt - yk^2 + yS^2)*
          \[Lambda]Rik) + xs^2*(-2*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4)*
          \[Lambda]Lik + xt*((-1 + 2*yk^4 + 2*yS^4 + yk^2*(2 - 4*yS^2))*
            \[Lambda]Lik - 4*yk*(1 + 2*yk^2 - 2*yS^2)*\[Lambda]Rik) + 
         xt^2*(\[Lambda]Lik - 6*yS^2*\[Lambda]Lik + 4*yk*(1 + yk^2 - yS^2)*
            \[Lambda]Rik))))*CC[7])/(2*xs^3*xt^3*(-1 + xs + xt)^2) - 
  (mi^2*Qk*\[Lambda]Ljk*(Qk*(-1 + xs)*xs*(xs^4*\[Lambda]Lik - 
       2*(-1 + xt)*xt*(yk^2 - yS^2)*\[Lambda]Lik + 
       xs^3*((-2 + 3*xt + 2*yk^2 - 2*yS^2)*\[Lambda]Lik + 
         4*xt*yk*\[Lambda]Rik) + 
       xs^2*((1 + 4*xt^2 - 4*yk^2 + 4*yS^2 + xt*(-4 + 6*yk^2 - 6*yS^2))*
          \[Lambda]Lik + 4*(-2 + xt)*xt*yk*\[Lambda]Rik) + 
       xs*(2*(yk^2 - yS^2)*\[Lambda]Lik + 
         xt^2*((-2 + 4*yk^2 - 4*yS^2)*\[Lambda]Lik - 4*yk*\[Lambda]Rik) + 
         xt*(\[Lambda]Lik - 8*yk^2*\[Lambda]Lik + 8*yS^2*\[Lambda]Lik + 
           4*yk*\[Lambda]Rik))) - 
     Qj*((-1 + xt)^2*xt*(yk^2 - yS^2)^2*\[Lambda]Lik - 
       xs*(-1 + xt)*((1 - 4*xt + xt^2)*yk^4*\[Lambda]Lik - 
         2*(1 - 4*xt + xt^2)*yk^2*yS^2*\[Lambda]Lik + 
         yS^2*(yS^2 - 4*xt*yS^2 + xt^2*(2 + yS^2))*\[Lambda]Lik - 
         4*(-1 + xt)*xt*yk^3*\[Lambda]Rik + 4*(-1 + xt)*xt*yk*yS^2*
          \[Lambda]Rik) + xs^4*(-((yk^2 - yS^2)^2*\[Lambda]Lik) + 
         xt^2*(\[Lambda]Lik - 4*yk*\[Lambda]Rik) + 
         2*xt*yk*(yk*\[Lambda]Lik - 2*yk^2*\[Lambda]Rik + 
           2*yS^2*\[Lambda]Rik)) - xs^3*(-3*(yk^2 - yS^2)^2*\[Lambda]Lik + 
         xt^3*(\[Lambda]Lik + 4*yk*\[Lambda]Rik) + 
         xt*(3*yk^4*\[Lambda]Lik + 3*yS^4*\[Lambda]Lik + yk^2*(4 - 6*yS^2)*
            \[Lambda]Lik - 12*yk^3*\[Lambda]Rik + 12*yk*yS^2*\[Lambda]Rik) + 
         xt^2*((1 - 6*yk^2 - 2*yS^2)*\[Lambda]Lik + 8*yk*(-1 + yk^2 - yS^2)*
            \[Lambda]Rik)) + xs^2*(-3*(yk^2 - yS^2)^2*\[Lambda]Lik + 
         xt*(7*yk^4*\[Lambda]Lik + 7*yS^4*\[Lambda]Lik + 2*yk^2*(1 - 7*yS^2)*
            \[Lambda]Lik - 12*yk^3*\[Lambda]Rik + 12*yk*yS^2*\[Lambda]Rik) + 
         xt^3*((1 + 4*yk^2 + 2*yS^2)*\[Lambda]Lik + 4*yk*(1 - yk^2 + yS^2)*
            \[Lambda]Rik) - xt^2*(3*yk^4*\[Lambda]Lik - 6*yk^2*(-1 + yS^2)*
            \[Lambda]Lik + yS^2*(4 + 3*yS^2)*\[Lambda]Lik - 
           16*yk^3*\[Lambda]Rik + 4*yk*(\[Lambda]Rik + 
             4*yS^2*\[Lambda]Rik)))))*CC[8])/(2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*(Qj - Qk)*\[Lambda]Ljk*
    (-(Qk*(-1 + xt)*xt*(-2*xs^3*(xt - yk^2 + yS^2)*\[Lambda]Lik + 
        (-1 + xt)^2*xt*(xt - 2*yk^2 + 2*yS^2)*\[Lambda]Lik + 
        xs*(-1 + xt)*xt*((1 + xt - 4*yk^2 + 4*yS^2)*\[Lambda]Lik - 
          4*(-1 + xt)*yk*\[Lambda]Rik) + 
        xs^2*(2*(-yk^2 + yS^2)*\[Lambda]Lik - 4*xt^2*yk*\[Lambda]Rik + 
          2*xt*(\[Lambda]Lik + 2*yk*\[Lambda]Rik)))) + 
     Qj*(-1 + xs + xt)*(4*xs^3*xt*yS^2*\[Lambda]Lik + 
       xt*(xt^2 + yk^2 - yS^2 + xt*(-1 - yk^2 + yS^2))^2*\[Lambda]Lik - 
       xs*(-1 + xt)*((yk^2 - yS^2)^2*\[Lambda]Lik + 4*xt^3*yk*\[Lambda]Rik - 
         2*xt*(yk^2 - yS^2)*(yk^2*\[Lambda]Lik - yS^2*\[Lambda]Lik - 
           2*yk*\[Lambda]Rik) + xt^2*((-1 + 2*yk^2)*\[Lambda]Lik - 
           4*yk*(1 + yk^2 - yS^2)*\[Lambda]Rik)) - 
       xs^2*((yk^2 - yS^2)^2*\[Lambda]Lik + xt^3*(\[Lambda]Lik + 
           4*yk*\[Lambda]Rik) - xt*(yk^4*\[Lambda]Lik - 
           2*yk^2*yS^2*\[Lambda]Lik + yS^2*(-2 + yS^2)*\[Lambda]Lik - 
           4*yk^3*\[Lambda]Rik + 4*yk*yS^2*\[Lambda]Rik) + 
         xt^2*((-1 + 2*yS^2)*\[Lambda]Lik - 4*yk*(1 + yk^2 - yS^2)*
            \[Lambda]Rik))))*CC[9])/(2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*Qk*\[Lambda]Ljk*
    (-(Qk*(-1 + xt)*xt*((-1 + xt)^2*xt*(xt + 2*yk^2 - 2*yS^2)*\[Lambda]Lik - 
        2*xs^3*(2*xt + yk^2 - yS^2)*\[Lambda]Lik + xs*(-1 + xt)*xt*
         ((1 + xt + 4*yk^2 - 4*yS^2)*\[Lambda]Lik + 4*(-1 + xt)*yk*
           \[Lambda]Rik) - 2*xs^2*((-yk^2 + yS^2)*\[Lambda]Lik + 
          xt^2*(\[Lambda]Lik - 2*yk*\[Lambda]Rik) - 
          2*xt*(\[Lambda]Lik - yk*\[Lambda]Rik)))) + 
     Qj*(-((-1 + xt)^3*xt*(yk^2 - yS^2)^2*\[Lambda]Lik) + 
       4*xs^4*xt*yk*(-(yk*\[Lambda]Lik) + (-1 + 2*xt)*\[Lambda]Rik) + 
       xs^3*((yk^2 - yS^2)^2*\[Lambda]Lik - 3*xt^3*(\[Lambda]Lik - 
           4*yk*\[Lambda]Rik) + xt^2*((3 - 6*yk^2 + 4*yS^2)*\[Lambda]Lik - 
           4*yk*(5 + yk^2 - yS^2)*\[Lambda]Rik) - 
         xt*(yk^4*\[Lambda]Lik + yS^2*(4 + yS^2)*\[Lambda]Lik - 
           2*yk^2*(5 + yS^2)*\[Lambda]Lik - 4*yk^3*\[Lambda]Rik + 
           4*yk*(-2 + yS^2)*\[Lambda]Rik)) - xs^2*(-1 + xt)*
        (-2*(yk^2 - yS^2)^2*\[Lambda]Lik + xt^3*(\[Lambda]Lik - 
           4*yk*\[Lambda]Rik) + 2*xt^2*((-1 + yk^2 - 3*yS^2)*\[Lambda]Lik + 
           4*yk*(1 + yk^2 - yS^2)*\[Lambda]Rik) + 
         xt*(3*yk^4*\[Lambda]Lik + 3*yS^2*(2 + yS^2)*\[Lambda]Lik - 
           2*yk^2*(4 + 3*yS^2)*\[Lambda]Lik - 8*yk^3*\[Lambda]Rik + 
           4*yk*(-1 + 2*yS^2)*\[Lambda]Rik)) - xs*(-1 + xt)^2*
        (-((yk^2 - yS^2)^2*\[Lambda]Lik) + xt*(yk^2 - yS^2)*
          ((-2 + 3*yk^2 - 3*yS^2)*\[Lambda]Lik - 4*yk*\[Lambda]Rik) - 
         2*xt^2*(-2*yk^3*\[Lambda]Rik + yS^2*(\[Lambda]Lik + 
             2*yk*\[Lambda]Rik)))))*CC[10])/(2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*Qk^2*\[Lambda]Ljk*(xs^5*\[Lambda]Lik + xt^4*(xt + 2*yk^2 - 2*yS^2)*
      \[Lambda]Lik + 4*xs^2*xt^2*(-(yk^2*\[Lambda]Lik) + yS^2*\[Lambda]Lik + 
       xt*yk*\[Lambda]Rik) + xs^4*(3*xt*\[Lambda]Lik + 2*yk^2*\[Lambda]Lik - 
       2*yS^2*\[Lambda]Lik + 4*xt*yk*\[Lambda]Rik) + 
     2*xs^3*xt*(2*(yk^2 - yS^2)*\[Lambda]Lik + 
       xt*(\[Lambda]Lik + 2*yk*\[Lambda]Rik)) + 
     xs*xt^3*(4*(yk^2 - yS^2)*\[Lambda]Lik + 
       xt*(\[Lambda]Lik + 4*yk*\[Lambda]Rik)))*CC[11])/
   (2*xs^3*xt^3*(xs + xt)) - (mi^2*(Qj - Qk)^2*\[Lambda]Ljk*
    (xs^5*\[Lambda]Lik + xt^4*(xt - 2*yk^2 + 2*yS^2)*\[Lambda]Lik + 
     xs^4*(2*(-yk^2 + yS^2)*\[Lambda]Lik + 
       xt*(\[Lambda]Lik - 4*yk*\[Lambda]Rik)) + 
     xs*xt^3*(4*(-yk^2 + yS^2)*\[Lambda]Lik + 
       xt*(\[Lambda]Lik - 4*yk*\[Lambda]Rik)) - 
     2*xs^3*xt*(2*(yk^2 - yS^2)*\[Lambda]Lik + 
       xt*(\[Lambda]Lik + 2*yk*\[Lambda]Rik)) - 
     2*xs^2*xt^2*(2*(-yk^2 + yS^2)*\[Lambda]Lik + 
       xt*(\[Lambda]Lik + 2*yk*\[Lambda]Rik)))*CC[12])/
   (2*xs^3*xt^3*(xs + xt)) + 
  (mi^4*(Qj - Qk)*Qk*((-1 + xt)^2*(yk^2 - yS^2)^2 + 
     xs^2*(xt^2 + (yk^2 - yS^2)^2 - 2*xt*(yk^2 + yS^2)) - 
     2*xs*((yk^2 - yS^2)^2 + xt^2*(yk^2 + yS^2) - 
       xt*(yk^2 + yk^4 + yS^2 - 2*yk^2*yS^2 + yS^4)))*\[Lambda]Ljk*
    ((-1 + xt)*xt*(yk^2 - yS^2)*\[Lambda]Lik + 
     xs^2*(3*xt*\[Lambda]Lik + yk^2*\[Lambda]Lik - yS^2*\[Lambda]Lik + 
       4*xt*yk*\[Lambda]Rik) + xs*(2*xt*(-1 + yk^2 - yS^2)*\[Lambda]Lik + 
       (-yk^2 + yS^2)*\[Lambda]Lik - 4*xt*yk*\[Lambda]Rik + 
       xt^2*(\[Lambda]Lik + 4*yk*\[Lambda]Rik)))*DD[1])/
   (2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^4*(Qj - Qk)*Qk*(xs*(xt - (yk - yS)^2) - (-1 + xt)*(yk - yS)^2)*
    (-((-1 + xt)*(yk + yS)^2) + xs*(xt - (yk + yS)^2))*\[Lambda]Ljk*
    ((-1 + xt)*xt*(yk^2 - yS^2)*\[Lambda]Lik + 
     xs^2*(-(xt*\[Lambda]Lik) + yk^2*\[Lambda]Lik - yS^2*\[Lambda]Lik + 
       4*xt*yk*\[Lambda]Rik) + xs*((-yk^2 + yS^2)*\[Lambda]Lik + 
       2*xt*(yk^2*\[Lambda]Lik - yS^2*\[Lambda]Lik - 2*yk*\[Lambda]Rik) + 
       xt^2*(\[Lambda]Lik + 4*yk*\[Lambda]Rik)))*DD[2])/
   (2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^4*(Qj - Qk)^2*(xs^3 + (-1 + xt)*(yk^2 - yS^2)^2 + 
     xs^2*(-1 + xt - 2*yk^2 + 2*yS^2) + xs*(yk^4 + yS^2*(-2 - 2*xt + yS^2) - 
       2*yk^2*(-1 + xt + yS^2)))*\[Lambda]Ljk*(xs^2*\[Lambda]Lik + 
     xt*(-yk^2 + yS^2)*\[Lambda]Lik - xs*((yk^2 - yS^2)*\[Lambda]Lik + 
       xt*(\[Lambda]Lik + 4*yk*\[Lambda]Rik)))*DD[3])/(2*xs^3*xt^3) - 
  (mi^4*Qk^2*(xs^3 + xs^2*(-1 + xt + 2*yk^2 - 2*yS^2) + 
     (-1 + xt)*(yk^2 - yS^2)^2 + xs*(yk^4 + yS^2*(2 - 2*xt + yS^2) - 
       2*yk^2*(1 + xt + yS^2)))*\[Lambda]Ljk*(xs^2*\[Lambda]Lik + 
     xs*(yk^2 - yS^2)*\[Lambda]Lik + xt*(yk^2 - yS^2)*\[Lambda]Lik + 
     xs*xt*(\[Lambda]Lik + 4*yk*\[Lambda]Rik))*DD[4])/(2*xs^3*xt^3) + 
  (mi^4*(Qj - Qk)^2*(xt^3 + (-1 + xs)*(yk^2 - yS^2)^2 + 
     xt^2*(-1 + xs - 2*yk^2 + 2*yS^2) + xt*(yk^4 + yS^2*(-2 - 2*xs + yS^2) - 
       2*yk^2*(-1 + xs + yS^2)))*\[Lambda]Ljk*
    (xt*(xt - yk^2 + yS^2)*\[Lambda]Lik - 
     xs*((yk^2 - yS^2)*\[Lambda]Lik + xt*(\[Lambda]Lik + 4*yk*\[Lambda]Rik)))*
    DD[5])/(2*xs^3*xt^3) + (mi^4*Qk^2*\[Lambda]Ljk*
    (-((-1 + xt)*xt*(xt + yk^2 - yS^2)^3*\[Lambda]Lik) - 
     xs*(xt + yk^2 - yS^2)*(-((yk^2 - yS^2)^2*\[Lambda]Lik) + 
       4*xt^3*yk*\[Lambda]Rik + 2*xt*(yk^2 - yS^2)*
        ((-1 + yk^2 - yS^2)*\[Lambda]Lik - 2*yk*\[Lambda]Rik) + 
       xt^2*(\[Lambda]Lik - 4*yS^2*\[Lambda]Lik + 4*yk*(-1 + yk^2 - yS^2)*
          \[Lambda]Rik)) + xs^2*(-((yk^2 - yS^2)^3*\[Lambda]Lik) + 
       xt^3*(\[Lambda]Lik + 4*yk*\[Lambda]Rik) - xt*(yk^2 - yS^2)*
        (-(yk^2*\[Lambda]Lik) - 3*yS^2*\[Lambda]Lik + 4*yk^3*\[Lambda]Rik - 
         4*yk*yS^2*\[Lambda]Rik) + xt^2*(3*yk^2*\[Lambda]Lik + 
         yS^2*\[Lambda]Lik + 8*yk^3*\[Lambda]Rik + 4*yk*(-1 + 2*yS^2)*
          \[Lambda]Rik)))*DD[6])/(2*xs^3*xt^3), 
 (2*Qj*(2*Qk*(-1 + 2*xt)*\[Lambda]Rik + Qj*(1 - yk^2 + yS^2)*\[Lambda]Rik - 
     2*Qj*xt*(2*yk*\[Lambda]Lik + \[Lambda]Rik + yk^2*\[Lambda]Rik - 
       yS^2*\[Lambda]Rik))*\[Lambda]Rjk)/(xs*(-1 + xt)*xt^2) - 
  (2*Qj^2*yk^2*((yk^2 - yS^2)*\[Lambda]Rik + 
     2*xt*(2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - yS^2*\[Lambda]Rik))*
    \[Lambda]Rjk*BB[1])/(xs*(-1 + xt)*xt^2*(yk^2 - yS^2)) + 
  (2*Qj^2*yS^2*((yk^2 - yS^2)*\[Lambda]Rik + 
     2*xt*(2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - yS^2*\[Lambda]Rik))*
    \[Lambda]Rjk*BB[2])/(xs*(-1 + xt)*xt^2*(yk^2 - yS^2)) + 
  (2*Qj*(2*Qk*(-((-2 + xt)*(-1 + xt)^2) + xs^2*xt + xs*(-2 + 3*xt - 2*xt^2))*
      \[Lambda]Rik + Qj*(-1 + xs + xt)*(2*(-1 + xs)*(1 - 3*xt + 2*xt^2)*yk*
        \[Lambda]Lik + (-4 + 7*xt - 4*xt^2 + xs*(2 - 3*xt + 2*xt^2))*yk^2*
        \[Lambda]Rik - (-2 + 2*(-2 + xs)*yS^2 + 
         xt*(1 + 7*yS^2 - 3*xs*(-1 + yS^2)) + 
         2*xt^2*(-2*yS^2 + xs*(-1 + yS^2)))*\[Lambda]Rik))*\[Lambda]Rjk*
    BB[3])/((-1 + xs)*xs*(-1 + xt)^2*xt^2*(-1 + xs + xt)) - 
  (4*Qj*(-(Qk*(-2*xs + 2*xs^2 + xt)) + Qj*(-1 + xs + xt)*(xs - yk^2 + yS^2))*
    \[Lambda]Rik*\[Lambda]Rjk*BB[4])/((-1 + xs)*xs*xt^2*(-1 + xs + xt)) - 
  (2*Qj*(2*Qk*(xs - 2*(-1 + xt)^2)*xt*\[Lambda]Rik + 
     Qj*(-1 + xs + xt)*(-2*(-1 + xt)*yk*\[Lambda]Lik + 
       (3 - 2*xt)*yk^2*\[Lambda]Rik + (2*xt^2 - 3*yS^2 + xt*(-3 + 2*yS^2))*
        \[Lambda]Rik))*\[Lambda]Rjk*BB[5])/(xs*(-1 + xt)^2*xt^2*
    (-1 + xs + xt)) - (4*Qk^2*\[Lambda]Rik*\[Lambda]Rjk*BB[6])/(xs*xt^2) + 
  (4*(Qj - Qk)^2*\[Lambda]Rik*\[Lambda]Rjk*BB[7])/(xs*xt^2) - 
  (2*mi^2*Qk*(Qj*((-1 + xt)^2*(yk^2 - yS^2)^2 + 
       xs^2*(xt^2 - 2*xt*yS^2 + (yk^2 - yS^2)^2) + 
       2*xs*((-1 + xt)*yk^4 - 2*(-1 + xt)*yk^2*yS^2 + 
         yS^2*(xt - xt^2 - yS^2 + xt*yS^2)))*\[Lambda]Rik + 
     Qk*xs*(-(xt^2*yk*\[Lambda]Lik) + (-1 + xs)^2*(xs + 2*yk^2 - 2*yS^2)*
        \[Lambda]Rik + (-1 + xs)*xt*(-(yk*\[Lambda]Lik) + 
         2*yk^2*\[Lambda]Rik + 2*(xs - yS^2)*\[Lambda]Rik)))*\[Lambda]Rjk*
    CC[1])/(xs*xt^3*(-1 + xs + xt)^2) + 
  (2*mi^2*(Qj - Qk)*(Qj*(xs^4 + (-1 + xt)^2*(yk^2 - yS^2)^2 + 
       2*xs^3*(-1 + xt - yk^2 + yS^2) + xs^2*(1 + xt^2 + yk^4 - 4*yS^2 + 
         yS^4 - 2*yk^2*(-2 + yS^2) + xt*(-2 - 4*yk^2 + 2*yS^2)) + 
       2*xs*((-1 + xt)*yk^4 + (-1 + xt)*yS^2*(-1 + yS^2) - 
         yk^2*(1 + xt^2 - 2*yS^2 + 2*xt*(-1 + yS^2))))*\[Lambda]Rik - 
     Qk*xs*((-1 + xs)^2*(xs - 2*yk^2 + 2*yS^2)*\[Lambda]Rik + 
       xt^2*(yk*\[Lambda]Lik + (-1 + 2*xs)*\[Lambda]Rik) + 
       (-1 + xs)*xt*(yk*\[Lambda]Lik - 2*yk^2*\[Lambda]Rik + 
         2*(xs + yS^2)*\[Lambda]Rik)))*\[Lambda]Rjk*CC[2])/
   (xs*xt^3*(-1 + xs + xt)^2) + 
  (2*mi^2*Qk*(-(Qk*(-1 + xt)*xt*(xt^3 + (-1 + xs)*xt*(-1 + 4*yk^2 - 4*yS^2) + 
        2*(-1 + xs)^2*(yk^2 - yS^2) + 2*xt^2*(-1 + xs + yk^2 - yS^2))*
       \[Lambda]Rik) + Qj*(-((-1 + xt)^3*(yk^2 - yS^2)^2*\[Lambda]Rik) + 
       4*xs^3*xt*yk*(\[Lambda]Lik + yk*\[Lambda]Rik) - 
       xs^2*(-1 + xt)*(xt^2*\[Lambda]Rik + (yk^2 - yS^2)^2*\[Lambda]Rik - 
         xt*yk*(7*\[Lambda]Lik + 6*yk*\[Lambda]Rik)) + 
       xs*(-1 + xt)*(2*(yk^2 - yS^2)^2*\[Lambda]Rik + 
         xt^2*(3*yk*\[Lambda]Lik + \[Lambda]Rik + 2*yk^2*\[Lambda]Rik) - 
         xt*(3*yk*\[Lambda]Lik + 2*yk^4*\[Lambda]Rik + 2*yS^4*\[Lambda]Rik + 
           yk^2*(2 - 4*yS^2)*\[Lambda]Rik))))*\[Lambda]Rjk*CC[3])/
   (xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)^2) + 
  (2*mi^2*(Qj - Qk)*(-(Qk*(-1 + xt)*xt*(xt^3 - 2*(-1 + xs)^2*(yk^2 - yS^2) + 
        2*xt^2*(-1 + xs - yk^2 + yS^2) + xt*(1 + 2*xs^2 + 4*yk^2 - 4*yS^2 + 
          xs*(-2 - 4*yk^2 + 4*yS^2)))) + 
     Qj*(xt^5 - (-1 + xs)^2*(yk^2 - yS^2)^2 + 
       xt^4*(-3 + 2*xs - 2*yk^2 + 2*yS^2) + xt^3*(3 + xs^2 + yk^4 - 6*yS^2 + 
         yS^4 - 2*yk^2*(-3 + yS^2) + xs*(-4 - 4*yk^2 + 2*yS^2)) + 
       xt*((3 - 4*xs + xs^2)*yk^4 - (-1 + xs)*yS^2*(-2 + 4*xs^2 + 3*yS^2 - 
           xs*yS^2) - 2*(-1 + xs)*yk^2*(1 - 3*yS^2 + xs*(-1 + yS^2))) - 
       xt^2*(1 + 3*yk^4 - 6*yS^2 + 3*yS^4 - 6*yk^2*(-1 + yS^2) + 
         xs^2*(1 + 2*yk^2 + 4*yS^2) - 2*xs*(yk^4 - 2*yk^2*(-2 + yS^2) + 
           (-1 + yS^2)^2))))*\[Lambda]Rik*\[Lambda]Rjk*CC[4])/
   (xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)^2) + 
  (2*mi^2*Qk^2*(xs^3*\[Lambda]Rik + xs^2*(-1 + xt + 2*yk^2 - 2*yS^2)*
      \[Lambda]Rik + (xt^3 + xt^2*(-1 + 2*yk^2 - 2*yS^2) - 
       2*(yk^2 - yS^2)^2 + 2*xt*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2)))*
      \[Lambda]Rik + xs*(xt^2*\[Lambda]Rik + 
       2*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2))*\[Lambda]Rik - 
       2*xt*(yk*\[Lambda]Lik + 2*yS^2*\[Lambda]Rik)))*\[Lambda]Rjk*CC[5])/
   (xs^2*xt^3) - (2*mi^2*(Qj - Qk)^2*(xs^3 + xt^3 - 2*(yk^2 - yS^2)^2 + 
     xt^2*(-1 - 2*yk^2 + 2*yS^2) + xs^2*(-1 + xt - 2*yk^2 + 2*yS^2) + 
     2*xt*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4) + 
     xs*(xt^2 - 4*xt*yk^2 + 2*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4)))*
    \[Lambda]Rik*\[Lambda]Rjk*CC[6])/(xs^2*xt^3) + 
  (2*mi^2*(Qj - Qk)*(-(Qk*(-1 + xs)*xs*(xs^3 - 2*(-1 + xt)^2*(yk^2 - yS^2) + 
        2*xs^2*(-1 + xt - yk^2 + yS^2) + xs*(1 + 2*xt^2 + 4*yk^2 - 4*yS^2 + 
          xt*(-2 - 4*yk^2 + 4*yS^2)))) + 
     Qj*(xs^5 - (-1 + xt)^2*(yk^2 - yS^2)^2 + 
       xs^4*(-3 + 2*xt - 2*yk^2 + 2*yS^2) + xs^3*(3 + xt^2 + yk^4 - 6*yS^2 + 
         yS^4 - 2*yk^2*(-3 + yS^2) + xt*(-4 - 4*yk^2 + 2*yS^2)) + 
       xs*((3 - 4*xt + xt^2)*yk^4 - (-1 + xt)*yS^2*(-2 + 4*xt^2 + 3*yS^2 - 
           xt*yS^2) - 2*(-1 + xt)*yk^2*(1 - 3*yS^2 + xt*(-1 + yS^2))) - 
       xs^2*(1 + 3*yk^4 - 6*yS^2 + 3*yS^4 - 6*yk^2*(-1 + yS^2) + 
         xt^2*(1 + 2*yk^2 + 4*yS^2) - 2*xt*(yk^4 - 2*yk^2*(-2 + yS^2) + 
           (-1 + yS^2)^2))))*\[Lambda]Rik*\[Lambda]Rjk*CC[7])/
   (xs^2*xt^3*(-1 + xs + xt)^2) - 
  (2*mi^2*Qk*(Qk*(-1 + xs)*xs*(xs^3 + xt^2*(1 + 2*yk^2 - 2*yS^2) + 
       2*(yk^2 - yS^2) - 4*xt*(yk^2 - yS^2) + 
       2*xs^2*(-1 + xt + yk^2 - yS^2) + xs*(1 - 4*yk^2 + 4*yS^2 + 
         xt*(-2 + 4*yk^2 - 4*yS^2)))*\[Lambda]Rik + 
     Qj*(-((-1 + xt)^2*(yk^2 - yS^2)^2*\[Lambda]Rik) + 
       xs^3*(xt^2*\[Lambda]Rik + (yk^2 - yS^2)^2*\[Lambda]Rik - 
         xt*yk*(\[Lambda]Lik + 2*yk*\[Lambda]Rik)) + 
       xs^2*(-3*(yk^2 - yS^2)^2*\[Lambda]Rik - xt^2*(3*yk*\[Lambda]Lik + 
           2*\[Lambda]Rik + 6*yk^2*\[Lambda]Rik) + 
         2*xt*(yk*\[Lambda]Lik + yk^4*\[Lambda]Rik + yS^4*\[Lambda]Rik - 
           2*yk^2*(-1 + yS^2)*\[Lambda]Rik)) + 
       xs*(3*(yk^2 - yS^2)^2*\[Lambda]Rik - 2*xt^3*yk*(\[Lambda]Lik + 
           2*yk*\[Lambda]Rik) - xt*(yk*\[Lambda]Lik + 4*yk^4*\[Lambda]Rik + 
           4*yS^4*\[Lambda]Rik + yk^2*(2 - 8*yS^2)*\[Lambda]Rik) + 
         xt^2*(3*yk*\[Lambda]Lik + yk^4*\[Lambda]Rik - 2*yk^2*(-3 + yS^2)*
            \[Lambda]Rik + (1 + yS^4)*\[Lambda]Rik))))*\[Lambda]Rjk*CC[8])/
   (xs^2*xt^3*(-1 + xs + xt)^2) + 
  (2*mi^2*(Qj - Qk)*(Qj*(xt^6 + (-1 + xs)^2*(yk^2 - yS^2)^2 + 
       2*xt^5*(-2 + xs - yk^2 + yS^2) + xt^4*(6 + xs^2 + yk^4 - 8*yS^2 + 
         yS^4 - 2*yk^2*(-4 + yS^2) + xs*(-6 - 4*yk^2 + 2*yS^2)) - 
       2*xt^3*(xs^2*(1 + yk^2) + 2*(1 + yk^4 - 3*yS^2 + yS^4 + 
           yk^2*(3 - 2*yS^2)) - xs*(3 + yk^4 - 2*yS^2 + yS^4 - 
           2*yk^2*(-3 + yS^2))) + xt^2*(1 + 6*yk^4 - 8*yS^2 + 6*yS^4 + 
         yk^2*(8 - 12*yS^2) + xs^2*(1 + yk^4 + 4*yS^2 + yS^4 - 
           2*yk^2*(-2 + yS^2)) - 2*xs*(1 + 3*yk^4 - yS^2 + 3*yS^4 - 
           6*yk^2*(-1 + yS^2))) + 2*xt*(-((2 - 3*xs + xs^2)*yk^4) + 
         (-1 + xs)*yS^2*(-1 + xs^2 + 2*yS^2 - xs*(1 + yS^2)) + 
         (-1 + xs)*yk^2*(1 - 4*yS^2 + xs*(-1 + 2*yS^2))))*\[Lambda]Rik - 
     Qk*(-1 + xt)^2*xt*((-1 + xt)^2*(xt - 2*yk^2 + 2*yS^2)*\[Lambda]Rik + 
       xs^2*(yk*\[Lambda]Lik + 2*xt*\[Lambda]Rik) + 
       xs*((-1 + xt)*yk*\[Lambda]Lik - 2*(-1 + xt)*yk^2*\[Lambda]Rik + 
         (-3*xt + 2*xt^2 - 2*yS^2 + 2*xt*yS^2)*\[Lambda]Rik)))*\[Lambda]Rjk*
    CC[9])/(xs^2*(-1 + xt)*xt^3*(-1 + xs + xt)^2) + 
  (2*mi^2*Qk*(Qk*(-1 + xt)^2*xt*(xs^2*yk*\[Lambda]Lik - 
       (-1 + xt)^2*(xt + 2*yk^2 - 2*yS^2)*\[Lambda]Rik - 
       xs*(-1 + xt)*(-(yk*\[Lambda]Lik) + 2*yk^2*\[Lambda]Rik + 
         2*(xt - yS^2)*\[Lambda]Rik)) + 
     Qj*(-((-1 + xt)^4*(yk^2 - yS^2)^2*\[Lambda]Rik) + 
       2*xs^3*xt*yk*((-1 + xt)*\[Lambda]Lik - yk*\[Lambda]Rik) - 
       xs^2*(-1 + xt)*(xt^3*\[Lambda]Rik - (yk^2 - yS^2)^2*\[Lambda]Rik - 
         xt^2*(4*yk*\[Lambda]Lik + \[Lambda]Rik + 2*yS^2*\[Lambda]Rik) + 
         xt*(4*yk*\[Lambda]Lik + yk^4*\[Lambda]Rik - 2*yk^2*(-2 + yS^2)*
            \[Lambda]Rik + yS^2*(2 + yS^2)*\[Lambda]Rik)) + 
       2*xs*(-1 + xt)^2*((yk^2 - yS^2)^2*\[Lambda]Rik + 
         xt^2*(yk*\[Lambda]Lik + yS^2*\[Lambda]Rik) - 
         xt*(yk*\[Lambda]Lik + yk^4*\[Lambda]Rik + yS^2*(1 + yS^2)*
            \[Lambda]Rik + yk^2*(\[Lambda]Rik - 2*yS^2*\[Lambda]Rik)))))*
    \[Lambda]Rjk*CC[10])/(xs^2*(-1 + xt)*xt^3*(-1 + xs + xt)^2) + 
  (2*mi^2*Qk^2*(xs^2 + xt^2)*(xs + xt + 2*yk^2 - 2*yS^2)*\[Lambda]Rik*
    \[Lambda]Rjk*CC[11])/(xs^2*xt^3) - 
  (2*mi^2*(Qj - Qk)^2*(xs^2 + xt^2)*(xs + xt - 2*yk^2 + 2*yS^2)*\[Lambda]Rik*
    \[Lambda]Rjk*CC[12])/(xs^2*xt^3) + 
  (2*mi^4*(Qj - Qk)*Qk*(xs*(xt - (yk - yS)^2) - (-1 + xt)*(yk - yS)^2)*
    ((-1 + xt)*(yk^2 - yS^2) + xs*(xt + yk^2 - yS^2))*
    (-((-1 + xt)*(yk + yS)^2) + xs*(xt - (yk + yS)^2))*\[Lambda]Rik*
    \[Lambda]Rjk*DD[1])/(xs^2*xt^3*(-1 + xs + xt)^2) + 
  (2*mi^4*(Qj - Qk)*Qk*((-1 + xt)^3*(yk^2 - yS^2)^3*\[Lambda]Rik + 
     xs^3*(xt^3*\[Lambda]Rik + (yk^2 - yS^2)^3*\[Lambda]Rik + 
       xt^2*(yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - yS^2*\[Lambda]Rik) - 
       xt*(yk^2 - yS^2)*(yk*\[Lambda]Lik + 3*yk^2*\[Lambda]Rik + 
         yS^2*\[Lambda]Rik)) + xs^2*(-3*(yk^2 - yS^2)^3*\[Lambda]Rik + 
       xt^3*(yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - (1 + yS^2)*\[Lambda]Rik) + 
       xt^2*(-2*yk^3*\[Lambda]Lik + yk*(-1 + 2*yS^2)*\[Lambda]Lik - 
         6*yk^4*\[Lambda]Rik + 4*yk^2*yS^2*\[Lambda]Rik + 
         2*yS^2*(1 + yS^2)*\[Lambda]Rik) + xt*(yk^2 - yS^2)*
        (2*yk*\[Lambda]Lik + 3*yk^4*\[Lambda]Rik - 6*yk^2*(-1 + yS^2)*
          \[Lambda]Rik + yS^2*(2 + 3*yS^2)*\[Lambda]Rik)) - 
     xs*(-1 + xt)*(3*(yk^2 - yS^2)^3*\[Lambda]Rik - xt*(yk^2 - yS^2)*
        (yk*\[Lambda]Lik + 3*yk^4*\[Lambda]Rik + yk^2*(3 - 6*yS^2)*
          \[Lambda]Rik + yS^2*(1 + 3*yS^2)*\[Lambda]Rik) + 
       xt^2*(yk^3*\[Lambda]Lik - yk*yS^2*\[Lambda]Lik + 3*yk^4*\[Lambda]Rik - 
         yS^2*(1 + yS^2)*\[Lambda]Rik - yk^2*(\[Lambda]Rik + 
           2*yS^2*\[Lambda]Rik))))*\[Lambda]Rjk*DD[2])/
   (xs^2*xt^3*(-1 + xs + xt)^2) + (2*mi^4*(Qj - Qk)^2*(xs - yk^2 + yS^2)*
    (xs^3 + (-1 + xt)*(yk^2 - yS^2)^2 + xs^2*(-1 + xt - 2*yk^2 + 2*yS^2) + 
     xs*(yk^4 + yS^2*(-2 - 2*xt + yS^2) - 2*yk^2*(-1 + xt + yS^2)))*
    \[Lambda]Rik*\[Lambda]Rjk*DD[3])/(xs^2*xt^3) - 
  (2*mi^4*Qk^2*(xs^5*\[Lambda]Rik + xs^4*(-2 + 2*xt + 3*yk^2 - 3*yS^2)*
      \[Lambda]Rik + (-1 + xt)^2*(yk^2 - yS^2)^3*\[Lambda]Rik - 
     xs^2*((-yk^6 + 3*yk^4*(2 + yS^2) - 3*yk^2*(1 + 4*yS^2 + yS^4) + 
         yS^2*(3 + 6*yS^2 + yS^4))*\[Lambda]Rik - xt*(-1 + yk^2 - yS^2)*
        (-(yk*\[Lambda]Lik) + 2*yk^2*\[Lambda]Rik - 6*yS^2*\[Lambda]Rik) + 
       xt^2*(yk*\[Lambda]Lik + 5*yk^2*\[Lambda]Rik + 3*yS^2*\[Lambda]Rik)) + 
     xs*(-((-3 + 2*yk^2 - 2*yS^2)*(yk^2 - yS^2)^2*\[Lambda]Rik) + 
       xt^2*(-(yk^3*\[Lambda]Lik) + yk*yS^2*\[Lambda]Lik - 
         yk^4*\[Lambda]Rik + 3*yS^4*\[Lambda]Rik - 2*yk^2*(-1 + yS^2)*
          \[Lambda]Rik) + xt*(yk^2 - yS^2)*(yk*\[Lambda]Lik + 
         2*yk^4*\[Lambda]Rik + 2*yS^2*(3 + yS^2)*\[Lambda]Rik - 
         2*yk^2*(\[Lambda]Rik + 2*yS^2*\[Lambda]Rik))) + 
     xs^3*(xt^2*\[Lambda]Rik + (1 + 3*yk^4 + 6*yS^2 + 3*yS^4 - 
         6*yk^2*(1 + yS^2))*\[Lambda]Rik + xt*(-(yk*\[Lambda]Lik) + 
         2*yk^2*\[Lambda]Rik - 2*(\[Lambda]Rik + 3*yS^2*\[Lambda]Rik))))*
    \[Lambda]Rjk*DD[4])/(xs^2*xt^3*(-1 + xs + xt)) + 
  (2*mi^4*(Qj - Qk)^2*(xt - yk^2 + yS^2)*(xt^3 + (-1 + xs)*(yk^2 - yS^2)^2 + 
     xt^2*(-1 + xs - 2*yk^2 + 2*yS^2) + xt*(yk^4 + yS^2*(-2 - 2*xs + yS^2) - 
       2*yk^2*(-1 + xs + yS^2)))*\[Lambda]Rik*\[Lambda]Rjk*DD[5])/
   (xs^2*xt^3) - (2*mi^4*Qk^2*((-1 + xt)^2*(xt + yk^2 - yS^2)^3*
      \[Lambda]Rik + xs^2*(xt^3*\[Lambda]Rik + (yk^2 - yS^2)^3*\[Lambda]Rik - 
       xt*(yk^2 - yS^2)*(yk*\[Lambda]Lik + yk^2*\[Lambda]Rik + 
         3*yS^2*\[Lambda]Rik) - xt^2*(3*yk*\[Lambda]Lik + 
         5*yk^2*\[Lambda]Rik + 3*yS^2*\[Lambda]Rik)) + 
     xs*(2*xt^4*\[Lambda]Rik - 2*(yk^2 - yS^2)^3*\[Lambda]Rik + 
       xt^2*(-(yk^3*\[Lambda]Lik) + yk*(3 + yS^2)*\[Lambda]Lik + 
         2*yk^4*\[Lambda]Rik - 8*yk^2*yS^2*\[Lambda]Rik + 
         6*yS^2*(1 + yS^2)*\[Lambda]Rik) + xt*(yk^2 - yS^2)*
        (yk*\[Lambda]Lik + 2*yk^4*\[Lambda]Rik + 2*yS^2*(3 + yS^2)*
          \[Lambda]Rik - 2*yk^2*(\[Lambda]Rik + 2*yS^2*\[Lambda]Rik)) + 
       xt^3*(-3*yk*\[Lambda]Lik + 2*yk^2*\[Lambda]Rik - 
         2*(\[Lambda]Rik + 3*yS^2*\[Lambda]Rik))))*\[Lambda]Rjk*DD[6])/
   (xs^2*xt^3*(-1 + xs + xt)), 
 (-2*Qj*(-2*Qk*xs*xt*(-2 + xs + xt) + Qj*(-((-1 + xt)*xt*(yk^2 - yS^2)) + 
       xs^2*(xt - yk^2 + 2*xt*yk^2 + yS^2 - 2*xt*yS^2) + 
       xs*(yk^2 - yS^2 + xt^2*(1 + 2*yk^2 - 2*yS^2) + 
         xt*(-2 - 4*yk^2 + 4*yS^2))))*\[Lambda]Lik*\[Lambda]Ljk)/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)) + 
  (2*Qj^2*(xs + xt - 2*xs*xt)*yk^2*\[Lambda]Lik*\[Lambda]Ljk*BB[1])/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2) + 
  (2*Qj^2*(-xt + xs*(-1 + 2*xt))*yS^2*\[Lambda]Lik*\[Lambda]Ljk*BB[2])/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2) + 
  (2*Qj*\[Lambda]Ljk*(-2*Qk*(2 - 2*xs + xs^2 - 2*xt + xt^2)*\[Lambda]Lik + 
     Qj*(-2 - 6*yk^2 + 6*yS^2 + xt*(4 + 8*yk^2 - 8*yS^2) + 
       xt^2*(-1 - 3*yk^2 + 3*yS^2) + xs^2*(-1 - 3*yk^2 + 3*yS^2 + 
         2*xt*(1 + yk^2 - yS^2)) + 2*xs*(2 + 4*yk^2 - 4*yS^2 - 
         4*xt*(1 + yk^2 - yS^2) + xt^2*(1 + yk^2 - yS^2)))*\[Lambda]Lik + 
     4*Qj*(-1 + xs)*(-1 + xt)*(-2 + xs + xt)*yk*\[Lambda]Rik)*BB[3])/
   ((-1 + xs)^2*xs*(-1 + xt)^2*xt*(-1 + xs + xt)) + 
  (2*Qj*\[Lambda]Ljk*(2*Qk*xs*\[Lambda]Lik + 
     Qj*((yk^2 - yS^2)*\[Lambda]Lik - 2*xs^2*(\[Lambda]Lik + 
         2*yk*\[Lambda]Rik) + xs*(\[Lambda]Lik + 4*yk*\[Lambda]Rik)))*BB[4])/
   ((-1 + xs)^2*xs^2*xt*(-1 + xs + xt)) + 
  (2*Qj*\[Lambda]Ljk*(2*Qk*xt*\[Lambda]Lik + 
     Qj*((yk^2 - yS^2)*\[Lambda]Lik - 2*xt^2*(\[Lambda]Lik + 
         2*yk*\[Lambda]Rik) + xt*(\[Lambda]Lik + 4*yk*\[Lambda]Rik)))*BB[5])/
   (xs*(-1 + xt)^2*xt^2*(-1 + xs + xt)) + 
  (2*mi^2*Qk*\[Lambda]Ljk*(Qk*xs*((-1 + xs)*(yk^2 - yS^2)*\[Lambda]Lik + 
       xt*(-1 + 2*xs + yk^2 - yS^2)*\[Lambda]Lik + 2*(-1 + xs)*xt*yk*
        \[Lambda]Rik + 2*xt^2*yk*\[Lambda]Rik) + 
     Qj*(-((-1 + xt)*xt*yk^2*\[Lambda]Lik) + 
       xs*(-(yS^2*\[Lambda]Lik) + xt*(1 - yk^2 + yS^2)*\[Lambda]Lik + 
         2*xt*yk*\[Lambda]Rik - 2*xt^2*yk*\[Lambda]Rik) + 
       xs^2*(yS^2*\[Lambda]Lik - 2*xt*(\[Lambda]Lik + yk*\[Lambda]Rik))))*
    CC[1])/(xs*xt^2*(-1 + xs + xt)^3) - 
  (2*mi^2*(Qj - Qk)*\[Lambda]Ljk*(Qj*(-xs + xs^2 + xt - xt^2)*yS^2*
      \[Lambda]Lik + Qk*xs*(xt*(-1 + yk^2 - yS^2)*\[Lambda]Lik + 
       (-1 + xs)*(yk^2 - yS^2)*\[Lambda]Lik + 2*(-1 + xs)*xt*yk*
        \[Lambda]Rik + 2*xt^2*(\[Lambda]Lik + yk*\[Lambda]Rik)))*CC[2])/
   (xs*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^2*Qk*\[Lambda]Ljk*(Qk*xt*((-1 + xt)*(yk^2 - yS^2)*\[Lambda]Lik + 
       xs*(-1 + 2*xt + yk^2 - yS^2)*\[Lambda]Lik + 2*xs^2*yk*\[Lambda]Rik + 
       2*xs*(-1 + xt)*yk*\[Lambda]Rik) + Qj*((-1 + xt)*xt*yS^2*\[Lambda]Lik - 
       xs^2*yk*(yk*\[Lambda]Lik + 2*xt*\[Lambda]Rik) + 
       xs*(yk^2*\[Lambda]Lik + xt*(1 - yk^2 + yS^2)*\[Lambda]Lik + 
         2*xt*yk*\[Lambda]Rik - 2*xt^2*(\[Lambda]Lik + yk*\[Lambda]Rik))))*
    CC[3])/(xs^2*xt*(-1 + xs + xt)^3) + 
  (2*mi^2*(Qj - Qk)*\[Lambda]Ljk*(Qj*(-xs + xs^2 + xt - xt^2)*yS^2*
      \[Lambda]Lik + Qk*xt*(-((-1 + xt)*(yk^2 - yS^2)*\[Lambda]Lik) + 
       xs*(1 - yk^2 + yS^2)*\[Lambda]Lik - 2*xs*(-1 + xt)*yk*\[Lambda]Rik - 
       2*xs^2*(\[Lambda]Lik + yk*\[Lambda]Rik)))*CC[4])/
   (xs^2*xt*(-1 + xs + xt)^3) - (4*mi^2*Qk^2*(xs + xt)*yk^2*\[Lambda]Lik*
    \[Lambda]Ljk*CC[5])/(xs^2*xt^2*(-1 + xs + xt)) + 
  (4*mi^2*(Qj - Qk)^2*(xs + xt)*yS^2*\[Lambda]Lik*\[Lambda]Ljk*CC[6])/
   (xs^2*xt^2*(-1 + xs + xt)) - (2*mi^2*(Qj - Qk)*\[Lambda]Ljk*
    (Qj*(xs^4 + (-1 + xt)*xt + xs^3*(-3 + 6*xt) + xs^2*(3 - 13*xt + 7*xt^2) + 
       xs*(-1 + 8*xt - 8*xt^2 + 2*xt^3))*yS^2*\[Lambda]Lik + 
     Qk*(-1 + xs)^2*xt*(xs*(-1 + yk^2 - yS^2)*\[Lambda]Lik + 
       (-1 + xt)*(yk^2 - yS^2)*\[Lambda]Lik + 2*xs*(-1 + xt)*yk*
        \[Lambda]Rik + 2*xs^2*(\[Lambda]Lik + yk*\[Lambda]Rik)))*CC[7])/
   ((-1 + xs)*xs^2*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^2*Qk*\[Lambda]Ljk*(Qk*(-1 + xs)^2*xt*
      ((-1 + xt)*(yk^2 - yS^2)*\[Lambda]Lik + xs*(-1 + 2*xt + yk^2 - yS^2)*
        \[Lambda]Lik + 2*xs^2*yk*\[Lambda]Rik + 2*xs*(-1 + xt)*yk*
        \[Lambda]Rik) + Qj*((-1 + xt)*xt*yS^2*\[Lambda]Lik + 
       xs^4*yk*(yk*\[Lambda]Lik - 2*xt*\[Lambda]Rik) + 
       xs^3*(-3*yk^2*\[Lambda]Lik + xt*(1 + 5*yk^2 + yS^2)*\[Lambda]Lik + 
         6*xt*yk*\[Lambda]Rik - 2*xt^2*(\[Lambda]Lik + yk*\[Lambda]Rik)) + 
       xs*(-(yk^2*\[Lambda]Lik) + 2*xt^3*yk^2*\[Lambda]Lik - 
         2*xt^2*((1 + 3*yk^2 + yS^2)*\[Lambda]Lik + yk*\[Lambda]Rik) + 
         xt*(\[Lambda]Lik + 5*yk^2*\[Lambda]Lik + 3*yS^2*\[Lambda]Lik + 
           2*yk*\[Lambda]Rik)) + xs^2*(3*yk^2*\[Lambda]Lik + 
         xt^2*((4 + 6*yk^2 + yS^2)*\[Lambda]Lik + 4*yk*\[Lambda]Rik) - 
         xt*((2 + 10*yk^2 + 3*yS^2)*\[Lambda]Lik + 6*yk*\[Lambda]Rik))))*
    CC[8])/((-1 + xs)*xs^2*xt^2*(-1 + xs + xt)^3) - 
  (2*mi^2*(Qj - Qk)*\[Lambda]Ljk*
    (Qj*(2*xs^3*xt + (-1 + xt)^3*xt + xs*(-1 + xt)^2*(-1 + 6*xt) + 
       xs^2*(1 - 8*xt + 7*xt^2))*yS^2*\[Lambda]Lik + 
     Qk*xs*(-1 + xt)^2*(xt*(-1 + yk^2 - yS^2)*\[Lambda]Lik + 
       (-1 + xs)*(yk^2 - yS^2)*\[Lambda]Lik + 2*(-1 + xs)*xt*yk*
        \[Lambda]Rik + 2*xt^2*(\[Lambda]Lik + yk*\[Lambda]Rik)))*CC[9])/
   (xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^2*Qk*\[Lambda]Ljk*(Qk*xs*(-1 + xt)^2*
      ((-1 + xs)*(yk^2 - yS^2)*\[Lambda]Lik + xt*(-1 + 2*xs + yk^2 - yS^2)*
        \[Lambda]Lik + 2*(-1 + xs)*xt*yk*\[Lambda]Rik + 
       2*xt^2*yk*\[Lambda]Rik) + Qj*(2*xs^3*xt*yk^2*\[Lambda]Lik + 
       (-1 + xt)^3*xt*yk^2*\[Lambda]Lik - xs*(-1 + xt)^2*
        (yS^2*\[Lambda]Lik + 2*xt^2*yk*\[Lambda]Rik - 
         xt*((1 + 5*yk^2 + yS^2)*\[Lambda]Lik + 2*yk*\[Lambda]Rik)) - 
       xs^2*(-1 + xt)*(yS^2*\[Lambda]Lik + 2*xt^2*(\[Lambda]Lik + 
           yk*\[Lambda]Rik) - xt*((2 + 6*yk^2 + yS^2)*\[Lambda]Lik + 
           2*yk*\[Lambda]Rik))))*CC[10])/(xs^2*(-1 + xt)*xt^2*
    (-1 + xs + xt)^3) + (2*mi^4*(Qj - Qk)*Qk*\[Lambda]Ljk*
    (-((-1 + xt)^2*xt*yk^2*(yk^2 - yS^2)*\[Lambda]Lik) - 
     xs*(-1 + xt)*(yS^2*(-yk^2 + yS^2)*\[Lambda]Lik + 
       xt^2*yk*(-(yk*\[Lambda]Lik) + 2*yk^2*\[Lambda]Rik + 
         2*yS^2*\[Lambda]Rik) + xt*(2*yk^4*\[Lambda]Lik - 
         yk^2*(1 + yS^2)*\[Lambda]Lik - yS^2*(1 + yS^2)*\[Lambda]Lik - 
         2*yk^3*\[Lambda]Rik - 2*yk*yS^2*\[Lambda]Rik)) + 
     xs^3*(yS^2*(-yk^2 + yS^2)*\[Lambda]Lik + 
       2*xt^2*(\[Lambda]Lik + yk*\[Lambda]Rik) - 
       xt*(2*yk^2*\[Lambda]Lik + 3*yS^2*\[Lambda]Lik + 2*yk^3*\[Lambda]Rik + 
         2*yk*yS^2*\[Lambda]Rik)) + xs^2*(2*yS^2*(yk^2 - yS^2)*\[Lambda]Lik + 
       2*xt^3*yk*\[Lambda]Rik + xt*(-(yk^4*\[Lambda]Lik) - 
         yk^2*(-3 + yS^2)*\[Lambda]Lik + 2*yS^2*(2 + yS^2)*\[Lambda]Lik + 
         4*yk^3*\[Lambda]Rik + 4*yk*yS^2*\[Lambda]Rik) - 
       xt^2*((1 + yk^2 + 3*yS^2)*\[Lambda]Lik + 2*yk*(1 + 2*yk^2 + 2*yS^2)*
          \[Lambda]Rik)))*DD[1])/(xs^2*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^4*(Qj - Qk)*Qk*\[Lambda]Ljk*((-1 + xt)^2*xt*yS^2*(-yk^2 + yS^2)*
      \[Lambda]Lik + xs^3*yk*(yk*(-yk^2 + yS^2)*\[Lambda]Lik + 
       2*xt^2*\[Lambda]Rik + xt*(yk*\[Lambda]Lik - 2*yk^2*\[Lambda]Rik - 
         2*yS^2*\[Lambda]Rik)) - xs*(-1 + xt)*
      (yk^2*(-yk^2 + yS^2)*\[Lambda]Lik + xt*(yk^4*\[Lambda]Lik + 
         yk^2*(-1 + yS^2)*\[Lambda]Lik - yS^2*(1 + 2*yS^2)*\[Lambda]Lik - 
         2*yk^3*\[Lambda]Rik - 2*yk*yS^2*\[Lambda]Rik) + 
       xt^2*(2*yk^2*\[Lambda]Lik + 3*yS^2*\[Lambda]Lik + 
         2*yk^3*\[Lambda]Rik + 2*yk*yS^2*\[Lambda]Rik)) + 
     xs^2*(2*yk^2*(yk^2 - yS^2)*\[Lambda]Lik + 
       2*xt^3*(\[Lambda]Lik + yk*\[Lambda]Rik) + 
       xt*(-2*yk^4*\[Lambda]Lik + yk^2*yS^2*\[Lambda]Lik + 
         yS^2*(1 + yS^2)*\[Lambda]Lik + 4*yk^3*\[Lambda]Rik + 
         4*yk*yS^2*\[Lambda]Rik) - xt^2*((1 + yk^2 + 3*yS^2)*\[Lambda]Lik + 
         2*yk*(1 + 2*yk^2 + 2*yS^2)*\[Lambda]Rik)))*DD[2])/
   (xs^2*xt^2*(-1 + xs + xt)^3) - (2*mi^4*(Qj - Qk)^2*yS^2*\[Lambda]Ljk*
    (xs^2*\[Lambda]Lik + xt*(-yk^2 + yS^2)*\[Lambda]Lik - 
     xs*((yk^2 - yS^2)*\[Lambda]Lik + xt*(\[Lambda]Lik + 4*yk*\[Lambda]Rik)))*
    DD[3])/(xs^2*xt^2*(-1 + xs + xt)) + 
  (2*mi^4*Qk^2*yk^2*\[Lambda]Ljk*(xs^2*\[Lambda]Lik + 
     xt*(yk^2 - yS^2)*\[Lambda]Lik + xs*(3*xt*\[Lambda]Lik + 
       yk^2*\[Lambda]Lik - yS^2*\[Lambda]Lik + 4*xt*yk*\[Lambda]Rik))*DD[4])/
   (xs^2*xt^2*(-1 + xs + xt)) + (2*mi^4*(Qj - Qk)^2*yS^2*\[Lambda]Ljk*
    (xs*(yk^2 - yS^2)*\[Lambda]Lik - xt*(xt - yk^2 + yS^2)*\[Lambda]Lik + 
     xs*xt*(\[Lambda]Lik + 4*yk*\[Lambda]Rik))*DD[5])/
   (xs^2*xt^2*(-1 + xs + xt)) + (2*mi^4*Qk^2*yk^2*\[Lambda]Ljk*
    (xt*(xt + yk^2 - yS^2)*\[Lambda]Lik + 
     xs*(3*xt*\[Lambda]Lik + yk^2*\[Lambda]Lik - yS^2*\[Lambda]Lik + 
       4*xt*yk*\[Lambda]Rik))*DD[6])/(xs^2*xt^2*(-1 + xs + xt)), 
 (Qj*(2*Qk*(1 - 2*xt)*\[Lambda]Rik + Qj*(-1 + yk^2 - yS^2)*\[Lambda]Rik + 
     Qj*xt*(2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - 
       (-2 + yS^2)*\[Lambda]Rik))*\[Lambda]Rjk)/(xs*(-1 + xt)*xt^2) + 
  (Qj^2*yk^2*((yk^2 - yS^2)*\[Lambda]Rik + 
     xt*(2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - yS^2*\[Lambda]Rik))*
    \[Lambda]Rjk*BB[1])/(xs*(-1 + xt)*xt^2*(yk^2 - yS^2)) + 
  (Qj^2*yS^2*((-yk^2 + yS^2)*\[Lambda]Rik + 
     xt*(-2*yk*\[Lambda]Lik - yk^2*\[Lambda]Rik + yS^2*\[Lambda]Rik))*
    \[Lambda]Rjk*BB[2])/(xs*(-1 + xt)*xt^2*(yk^2 - yS^2)) + 
  (Qj*(2*Qk*(xs - xt)*\[Lambda]Rik + Qj*(xt*(1 - yk^2 + yS^2)*\[Lambda]Rik + 
       xs*(-1 + xt)*(2*yk*\[Lambda]Lik + \[Lambda]Rik + yk^2*\[Lambda]Rik - 
         yS^2*\[Lambda]Rik) - xs^2*(2*(-1 + xt)*yk*\[Lambda]Lik + 
         (-2 + xt)*yk^2*\[Lambda]Rik + (xt + 2*yS^2 - xt*yS^2)*
          \[Lambda]Rik)))*\[Lambda]Rjk*BB[3])/((-1 + xs)*xs^2*(-1 + xt)*
    xt^2) + (Qj*(2*Qk*xs - Qj*(xs - yk^2 + yS^2))*\[Lambda]Rik*\[Lambda]Rjk*
    BB[4])/((-1 + xs)*xs*xt^2) + 
  (Qj*(-2*Qk*xt^2*\[Lambda]Rik + Qj*xt*(xt - yk^2 + yS^2)*\[Lambda]Rik + 
     Qj*xs*(-2*yk*\[Lambda]Lik - 3*yk^2*\[Lambda]Rik + 
       (xt + 3*yS^2)*\[Lambda]Rik))*\[Lambda]Rjk*BB[5])/
   (xs^2*(-1 + xt)*xt^2) + (Qk^2*(-xs + xt)*\[Lambda]Rik*\[Lambda]Rjk*BB[6])/
   (xs^2*xt^2) + ((Qj - Qk)^2*(xs - xt)*\[Lambda]Rik*\[Lambda]Rjk*BB[7])/
   (xs^2*xt^2) - 
  (mi^2*Qk*(Qk*(-1 + xs)*xs^2*(xs^2 + xt - 2*yk^2 + 2*xt*yk^2 + 2*yS^2 - 
       2*xt*yS^2 + xs*(-1 + xt + 2*yk^2 - 2*yS^2)) - 
     Qj*((-1 + xt)^2*xt*(yk^2 - yS^2)^2 + xs^3*(xt^2 + 2*xt*yk^2 - 
         (yk^2 - yS^2)^2) + xs^2*(xt^3 + 2*(yk^2 - yS^2)^2 - 
         xt*(yk^4 - 2*yk^2*(-2 + yS^2) + yS^2*(-2 + yS^2)) + 
         xt^2*(6*yk^2 - 2*(1 + yS^2))) + xs*(xt^3*(4*yk^2 - 2*yS^2) + 
         2*xt*(yk^2 - yS^2) - (yk^2 - yS^2)^2 + 
         xt^2*(yk^4 - 2*yk^2*(3 + yS^2) + yS^2*(4 + yS^2)))))*\[Lambda]Rik*
    \[Lambda]Rjk*CC[1])/(2*xs^2*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*(Qk^2*(-1 + xs)*xs^2*(xs^2 + 2*xt^2 + 2*(yk^2 - yS^2) + 
       xt*(-1 - 2*yk^2 + 2*yS^2) + xs*(-1 + 3*xt - 2*yk^2 + 2*yS^2)) - 
     Qj*Qk*(2*xs^5 - (-1 + xt)^2*xt*(yk^2 - yS^2)^2 + 
       xs^4*(-4 + 6*xt - 4*yk^2 + 4*yS^2) + xs^3*(2 + 5*xt^2 + yk^4 - 
         8*yS^2 + yS^4 - 2*yk^2*(-4 + yS^2) + xt*(-8 - 6*yk^2 + 4*yS^2)) + 
       xs*(-((-1 + xt^2)*yk^4) + 2*(-1 + xt^2)*yk^2*yS^2 - 
         (-1 + xt)*yS^2*(2*xt^2 + yS^2 + xt*yS^2)) + 
       xs^2*(xt^3 - 2*xt^2*(2 + yk^2 + yS^2) + xt*(2 + yk^4 - 4*yS^2 + yS^4 - 
           2*yk^2*(-3 + yS^2)) - 2*(yk^4 + yS^2*(-2 + yS^2) - 
           2*yk^2*(-1 + yS^2)))) + 
     Qj^2*(xs^5 - (-1 + xt)^2*xt*(yk^2 - yS^2)^2 + 
       xs^4*(-2 + 3*xt - 2*yk^2 + 2*yS^2) + xs^3*(1 + 3*xt^2 + yk^4 - 
         4*yS^2 + yS^4 - 2*yk^2*(-2 + yS^2) + xt*(-4 - 4*yk^2 + 2*yS^2)) + 
       xs*(-((-1 + xt^2)*yk^4) + 2*(-1 + xt^2)*yk^2*yS^2 - 
         (-1 + xt)*yS^2*(2*xt^2 + yS^2 + xt*yS^2)) + 
       xs^2*(xt^3 - 2*xt^2*(1 + yk^2 + yS^2) - 2*(yk^2 + yk^4 - yS^2 - 
           2*yk^2*yS^2 + yS^4) + xt*(yk^4 - 2*yk^2*(-2 + yS^2) + 
           (-1 + yS^2)^2))))*\[Lambda]Rik*\[Lambda]Rjk*CC[2])/
   (2*xs^2*xt^3*(-1 + xs + xt)^2) - 
  (mi^2*Qk*(-(Qk*(-1 + xt)*xt*((-1 + xt)^2*xt*(xt + 2*yk^2 - 2*yS^2) + 
        2*xs^3*(yk^2 - yS^2) + xs*(-1 + xt)*(3*xt^2 - 2*yk^2 + 2*yS^2 + 
          xt*(-1 + 4*yk^2 - 4*yS^2)) + 2*xs^2*(xt^2 - 2*yk^2 + 2*yS^2 + 
          xt*(-1 + 2*yk^2 - 2*yS^2)))*\[Lambda]Rik) + 
     Qj*(-((-1 + xt)^3*xt*(yk^2 - yS^2)^2*\[Lambda]Rik) - 
       xs*(-1 + xt)^2*((1 + xt)*yk^4 - 2*xt^2*yS^2 - 2*(1 + xt)*yk^2*yS^2 + 
         yS^4 + xt*yS^4)*\[Lambda]Rik + 4*xs^4*xt*yk*(\[Lambda]Lik + 
         yk*\[Lambda]Rik) - xs^3*(-1 + xt)*(xt^2*\[Lambda]Rik - 
         (yk^2 - yS^2)^2*\[Lambda]Rik - 2*xt*yk*(4*\[Lambda]Lik + 
           5*yk*\[Lambda]Rik)) - xs^2*(-1 + xt)*(xt^3*\[Lambda]Rik + 
         2*(yk^2 - yS^2)^2*\[Lambda]Rik - 2*xt^2*(2*yk*\[Lambda]Lik + 
           3*yk^2*\[Lambda]Rik + yS^2*\[Lambda]Rik) + 
         xt*(4*yk*\[Lambda]Lik - yk^4*\[Lambda]Rik - yS^4*\[Lambda]Rik + 
           2*yk^2*(3 + yS^2)*\[Lambda]Rik))))*\[Lambda]Rjk*CC[3])/
   (2*xs^3*(-1 + xt)*xt^2*(-1 + xs + xt)^2) + 
  (mi^2*(-(Qk^2*(-1 + xt)*xt*(2*xs^3*(xt - yk^2 + yS^2) + 
        (-1 + xt)^2*xt*(xt - 2*yk^2 + 2*yS^2) + 
        4*xs^2*(xt^2 + yk^2 - yS^2 + xt*(-1 - yk^2 + yS^2)) + 
        xs*(-1 + xt)*(3*xt^2 + 2*(yk^2 - yS^2) + 
          xt*(-1 - 4*yk^2 + 4*yS^2)))) + 
     Qj^2*(4*xs^4*xt*yS^2 - (-1 + xt)^3*xt*(xt - yk^2 + yS^2)^2 - 
       xs*(-1 + xt)^2*(3*xt^3 + (yk^2 - yS^2)^2 + xt*(yk^2 - yS^2)^2 + 
         xt^2*(-1 - 4*yk^2 + 2*yS^2)) - xs^3*(xt^3 + (yk^2 - yS^2)^2 - 
         xt^2*(1 + 10*yS^2) - xt*(yk^4 - 10*yS^2 - 2*yk^2*yS^2 + yS^4)) + 
       xs^2*(-3*xt^4 + 2*(yk^2 - yS^2)^2 + xt^3*(5 + 2*yk^2 + 6*yS^2) - 
         3*xt*(yk^4 - 2*yS^2 - 2*yk^2*yS^2 + yS^4) + 
         xt^2*(-2 + yk^4 - 12*yS^2 + yS^4 - 2*yk^2*(1 + yS^2)))) - 
     Qj*Qk*(4*xs^4*xt*yS^2 - (-1 + xt)^3*xt*(2*xt^2 - 4*xt*(yk^2 - yS^2) + 
         (yk^2 - yS^2)^2) - xs*(-1 + xt)^2*(6*xt^3 + (yk^2 - yS^2)^2 + 
         xt^2*(-2 - 8*yk^2 + 6*yS^2) + xt*(yk^4 + yS^2*(-2 + yS^2) - 
           2*yk^2*(-1 + yS^2))) - xs^3*(3*xt^3 + (yk^2 - yS^2)^2 - 
         xt^2*(3 + 2*yk^2 + 8*yS^2) - xt*(yk^4 + yS^2*(-8 + yS^2) - 
           2*yk^2*(1 + yS^2))) + xs^2*(-7*xt^4 + 2*(yk^2 - yS^2)^2 + 
         xt^3*(13 + 6*yk^2 + 2*yS^2) + xt^2*(-6 + yk^4 - 4*yS^2 + yS^4 - 
           2*yk^2*(5 + yS^2)) + xt*(-3*yk^4 + 2*yS^2 - 3*yS^4 + 
           yk^2*(4 + 6*yS^2)))))*\[Lambda]Rik*\[Lambda]Rjk*CC[4])/
   (2*xs^3*(-1 + xt)*xt^2*(-1 + xs + xt)^2) + 
  (mi^2*Qk^2*(xs^4 + xs^3*(-1 + 2*yk^2 - 2*yS^2) + 
     xs*(-2*xt^3 + 4*xt*(yk^2 - yS^2) - 2*(yk^2 - yS^2)^2 + 
       xt^2*(1 - 2*yk^2 + 6*yS^2)) + 
     xs^2*(-2*xt^2 + xt*(1 - 6*yk^2 + 2*yS^2) + 
       2*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2))) + 
     xt*(-xt^3 + 2*(yk^2 - yS^2)^2 + xt^2*(1 - 2*yk^2 + 2*yS^2) - 
       2*xt*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2))))*\[Lambda]Rik*
    \[Lambda]Rjk*CC[5])/(2*xs^3*xt^3) - 
  (mi^2*(Qj - Qk)^2*(xs^4 + xs^3*(-1 + 2*xt - 2*yk^2 + 2*yS^2) + 
     xs*(-2*xt^3 - 2*(yk^2 - yS^2)^2 + xt^2*(1 + 2*yk^2 + 2*yS^2)) - 
     xs^2*(xt*(1 + 2*yk^2 + 2*yS^2) - 2*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + 
         yS^4)) + xt*(-xt^3 + xt^2*(1 + 2*yk^2 - 2*yS^2) + 
       2*(yk^2 - yS^2)^2 - 2*xt*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4)))*
    \[Lambda]Rik*\[Lambda]Rjk*CC[6])/(2*xs^3*xt^3) + 
  (mi^2*(Qk^2*(-1 + xs)*xs*(xs^4 - 2*(-1 + xt)*xt*(yk^2 - yS^2) + 
       xs*(xt - 2*yk^2 + 2*yS^2) + xs^3*(-2 + 3*xt - 2*yk^2 + 2*yS^2) + 
       xs^2*(1 + 2*xt^2 + 4*yk^2 - 4*yS^2 - 2*xt*(2 + yk^2 - yS^2))) + 
     Qj*Qk*(-2*xs^6 + xs^5*(6 - 6*xt + 4*yk^2 - 4*yS^2) - 
       (-1 + xt)^2*xt*(yk^2 - yS^2)^2 - xs^4*(6 + 5*xt^2 + yk^4 - 12*yS^2 + 
         yS^4 - 2*xt*(7 + 3*yk^2 - 2*yS^2) - 2*yk^2*(-6 + yS^2)) - 
       xs^3*(-2 + xt^3 - 3*yk^4 + 12*yS^2 - 3*yS^4 + 6*yk^2*(-2 + yS^2) - 
         xt^2*(7 + 2*yk^2 + 2*yS^2) + xt*(10 + yk^4 - 6*yS^2 + yS^4 - 
           2*yk^2*(-5 + yS^2))) + xs^2*(-3*yk^4 + 4*yS^2 - 3*yS^4 + 
         xt^3*(1 + 2*yS^2) + yk^2*(-4 + 6*yS^2) + 
         xt^2*(-2 + yk^4 - 2*yS^2 - 2*yk^2*yS^2 + yS^4) + 
         xt*(2 + yk^4 + yS^4 - 2*yk^2*(-1 + yS^2))) + 
       xs*((yk^2 - yS^2)^2 + xt^3*(yk^4 + 2*yS^2 - 2*yk^2*yS^2 + yS^4) - 
         xt^2*(3*yk^4 + 3*yS^4 + yk^2*(2 - 6*yS^2)) + 
         xt*(yk^4 + yS^2*(-2 + yS^2) - 2*yk^2*(-1 + yS^2)))) + 
     Qj^2*(xs^6 + (-1 + xt)^2*xt*(yk^2 - yS^2)^2 + 
       xs^5*(-3 + 3*xt - 2*yk^2 + 2*yS^2) + xs^4*(3 + 3*xt^2 + yk^4 - 
         6*yS^2 + yS^4 - 2*yk^2*(-3 + yS^2) + xt*(-7 - 4*yk^2 + 2*yS^2)) + 
       xs^3*(-1 + xt^3 - 3*yk^4 + 6*yS^2 - 3*yS^4 + 6*yk^2*(-1 + yS^2) - 
         xt^2*(5 + 2*yk^2 + 2*yS^2) + xt*(5 + yk^4 - 4*yS^2 + yS^4 - 
           2*yk^2*(-4 + yS^2))) - xs^2*(-3*yk^4 + 2*yS^2 - 3*yS^4 + 
         xt^3*(1 + 2*yS^2) + yk^2*(-2 + 6*yS^2) + 
         xt*(yk^4 - 2*yk^2*(-2 + yS^2) + (-1 + yS^2)^2) + 
         xt^2*(-2 + yk^4 + yS^4 - 2*yk^2*(1 + yS^2))) - 
       xs*((1 + xt - 3*xt^2 + xt^3)*yk^4 - 2*(1 + xt - 3*xt^2 + xt^3)*yk^2*
          yS^2 + (-1 + xt)*yS^2*(-yS^2 - 2*xt*yS^2 + xt^2*(2 + yS^2)))))*
    \[Lambda]Rik*\[Lambda]Rjk*CC[7])/(2*xs^3*xt^3*(-1 + xs + xt)^2) - 
  (mi^2*Qk*(Qk*(-1 + xs)*xs*(xs^4 + xs^3*(-2 + xt + 2*yk^2 - 2*yS^2) + 
       xs*(-xt + 2*xt^2 + 2*yk^2 - 2*yS^2) + 2*(-1 + xt)*xt*(yk^2 - yS^2) + 
       xs^2*(1 + 2*(-2 + xt)*yk^2 - 2*(-2 + xt)*yS^2)) - 
     Qj*(-((-1 + xt)^2*xt*(yk^2 - yS^2)^2) + 
       xs^4*(xt^2 + 2*xt*yk^2 - (yk^2 - yS^2)^2) + 
       xs^3*(xt^3 + xt^2*(-1 + 6*yk^2 - 2*yS^2) + 3*(yk^2 - yS^2)^2 - 
         xt*(yk^4 + yS^4 - 2*yk^2*(-2 + yS^2))) + 
       xs*((1 + xt - 3*xt^2 + xt^3)*yk^4 - 2*(1 + xt - 3*xt^2 + xt^3)*yk^2*
          yS^2 + (-1 + xt)*yS^2*(-yS^2 - 2*xt*yS^2 + xt^2*(2 + yS^2))) + 
       xs^2*(xt^3*(-1 + 4*yk^2 - 2*yS^2) - 3*(yk^2 - yS^2)^2 + 
         xt*(yk^4 + yS^4 - 2*yk^2*(-1 + yS^2)) + 
         xt^2*(yk^4 - 2*yk^2*(3 + yS^2) + yS^2*(4 + yS^2)))))*\[Lambda]Rik*
    \[Lambda]Rjk*CC[8])/(2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*(-(Qk^2*(-1 + xt)*xt*(2*xs^3*(xt - yk^2 + yS^2) + 
        (-1 + xt)^2*xt*(xt - 2*yk^2 + 2*yS^2) + 
        xs*xt*(1 + 3*xt^2 + 4*yk^2 - 4*yS^2 - 4*xt*(1 + yk^2 - yS^2)) + 
        2*xs^2*(2*xt^2 + yk^2 - yS^2 + xt*(-1 - 2*yk^2 + 2*yS^2)))) + 
     Qj^2*(4*xs^4*xt*yS^2 - (-1 + xt)^3*xt*(xt - yk^2 + yS^2)^2 - 
       xs*(-1 + xt)^2*(3*xt^3 + (yk^2 - yS^2)^2 + xt*(yk^2 - yS^2)^2 + 
         xt^2*(-1 - 4*yk^2 + 2*yS^2)) - xs^3*(xt^3 + (yk^2 - yS^2)^2 - 
         xt^2*(1 + 10*yS^2) - xt*(yk^4 - 6*yS^2 - 2*yk^2*yS^2 + yS^4)) + 
       xs^2*(-3*xt^4 + 2*(yk^2 - yS^2)^2 + xt^3*(5 + 2*yk^2 + 6*yS^2) + 
         xt*(-3*yk^4 + 2*yS^2 + 6*yk^2*yS^2 - 3*yS^4) + 
         xt^2*(-2 + yk^4 - 8*yS^2 + yS^4 - 2*yk^2*(1 + yS^2)))) - 
     Qj*Qk*(4*xs^4*xt*yS^2 - (-1 + xt)^3*xt*(2*xt^2 - 4*xt*(yk^2 - yS^2) + 
         (yk^2 - yS^2)^2) - xs*(-1 + xt)^2*(6*xt^3 + (yk^2 - yS^2)^2 + 
         xt*(yk^2 - yS^2)^2 + xt^2*(-2 - 8*yk^2 + 6*yS^2)) - 
       xs^3*(3*xt^3 + (yk^2 - yS^2)^2 - xt^2*(3 + 2*yk^2 + 8*yS^2) - 
         xt*(yk^4 + yS^2*(-4 + yS^2) - 2*yk^2*(1 + yS^2))) + 
       xs^2*(-7*xt^4 + 2*(yk^2 - yS^2)^2 + xt^3*(11 + 6*yk^2 + 2*yS^2) + 
         xt^2*(-4 + yk^4 - 2*yS^2 + yS^4 - 2*yk^2*(4 + yS^2)) + 
         xt*(-3*yk^4 - 3*yS^4 + yk^2*(2 + 6*yS^2)))))*\[Lambda]Rik*
    \[Lambda]Rjk*CC[9])/(2*xs^3*xt^3*(-1 + xs + xt)^2) - 
  (mi^2*Qk*(-(Qk*(-1 + xt)*xt*((-1 + xt)^2*xt*(xt + 2*yk^2 - 2*yS^2) + 
        2*xs^3*(yk^2 - yS^2) + xs*xt*(1 + 3*xt^2 - 4*yk^2 + 4*yS^2 + 
          4*xt*(-1 + yk^2 - yS^2)) + 2*xs^2*(xt^2 - yk^2 + yS^2 + 
          2*xt*(yk^2 - yS^2)))*\[Lambda]Rik) + 
     Qj*(-((-1 + xt)^3*xt*(yk^2 - yS^2)^2*\[Lambda]Rik) + 
       xs*(-1 + xt)^2*(2*xt^2*yS^2 - (yk^2 - yS^2)^2 - 
         xt*(yk^4 - 2*yk^2*(1 + yS^2) + yS^2*(2 + yS^2)))*\[Lambda]Rik + 
       4*xs^4*xt*yk*(\[Lambda]Lik + yk*\[Lambda]Rik) - 
       xs^2*(-1 + xt)*(xt^3*\[Lambda]Rik + 2*(yk^2 - yS^2)^2*\[Lambda]Rik - 
         2*xt^2*(2*yk*\[Lambda]Lik + \[Lambda]Rik + 3*yk^2*\[Lambda]Rik + 
           yS^2*\[Lambda]Rik) + xt*(4*yk*\[Lambda]Lik - yk^4*\[Lambda]Rik + 
           2*yk^2*yS^2*\[Lambda]Rik - yS^2*(-2 + yS^2)*\[Lambda]Rik)) + 
       xs^3*(-(xt^3*\[Lambda]Rik) - (yk^2 - yS^2)^2*\[Lambda]Rik + 
         xt^2*(8*yk*\[Lambda]Lik + \[Lambda]Rik + 10*yk^2*\[Lambda]Rik) + 
         xt*(-8*yk*\[Lambda]Lik + yk^4*\[Lambda]Rik + yS^4*\[Lambda]Rik - 
           2*yk^2*(3 + yS^2)*\[Lambda]Rik))))*\[Lambda]Rjk*CC[10])/
   (2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*Qk^2*(xs + xt)*(xs^3 - xs*xt^2 - xt^2*(xt + 2*yk^2 - 2*yS^2) - 
     xs^2*(xt - 2*yk^2 + 2*yS^2))*\[Lambda]Rik*\[Lambda]Rjk*CC[11])/
   (2*xs^3*xt^3) - (mi^2*(Qj - Qk)^2*(xs - xt)*(xs + xt)^2*
    (xs + xt - 2*yk^2 + 2*yS^2)*\[Lambda]Rik*\[Lambda]Rjk*CC[12])/
   (2*xs^3*xt^3) - (mi^4*(Qj - Qk)*Qk*(xs*(xt - (yk - yS)^2) - 
     (-1 + xt)*(yk - yS)^2)*((-1 + xt)*xt*(yk^2 - yS^2) + 
     xs*(-2*xt + xt^2 + yk^2 - yS^2) + xs^2*(xt - yk^2 + yS^2))*
    (-((-1 + xt)*(yk + yS)^2) + xs*(xt - (yk + yS)^2))*\[Lambda]Rik*
    \[Lambda]Rjk*DD[1])/(2*xs^3*xt^3*(-1 + xs + xt)^2) - 
  (mi^4*(Qj - Qk)*Qk*(xs*(xt - (yk - yS)^2) - (-1 + xt)*(yk - yS)^2)*
    ((-1 + xt)*xt*(yk^2 - yS^2) + xs*(xt^2 + yk^2 - yS^2) + 
     xs^2*(xt - yk^2 + yS^2))*(-((-1 + xt)*(yk + yS)^2) + 
     xs*(xt - (yk + yS)^2))*\[Lambda]Rik*\[Lambda]Rjk*DD[2])/
   (2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^4*(Qj - Qk)^2*(xs^2 + xt*(yk^2 - yS^2) + xs*(xt - yk^2 + yS^2))*
    (xs^3 + (-1 + xt)*(yk^2 - yS^2)^2 + xs^2*(-1 + xt - 2*yk^2 + 2*yS^2) + 
     xs*(yk^4 + yS^2*(-2 - 2*xt + yS^2) - 2*yk^2*(-1 + xt + yS^2)))*
    \[Lambda]Rik*\[Lambda]Rjk*DD[3])/(2*xs^3*xt^3) - 
  (mi^4*Qk^2*(xs - xt)*(xs + yk^2 - yS^2)*
    (xs^3 + xs^2*(-1 + xt + 2*yk^2 - 2*yS^2) + (-1 + xt)*(yk^2 - yS^2)^2 + 
     xs*(yk^4 + yS^2*(2 - 2*xt + yS^2) - 2*yk^2*(1 + xt + yS^2)))*
    \[Lambda]Rik*\[Lambda]Rjk*DD[4])/(2*xs^3*xt^3) - 
  (mi^4*(Qj - Qk)^2*(xs*(xt + yk^2 - yS^2) + xt*(xt - yk^2 + yS^2))*
    (xt^3 + (-1 + xs)*(yk^2 - yS^2)^2 + xt^2*(-1 + xs - 2*yk^2 + 2*yS^2) + 
     xt*(yk^4 + yS^2*(-2 - 2*xs + yS^2) - 2*yk^2*(-1 + xs + yS^2)))*
    \[Lambda]Rik*\[Lambda]Rjk*DD[5])/(2*xs^3*xt^3) + 
  (mi^4*Qk^2*((-1 + xt)*xt*(xt + yk^2 - yS^2)^3*\[Lambda]Rik + 
     xs*(2*xt^4 + xt^3*(-1 + 2*yk^2 - 6*yS^2) - xt*(yk^2 - yS^2)^2 + 
       (yk^2 - yS^2)^3 - xt^2*(yk^2 - yS^2)*(3 + 4*yS^2))*\[Lambda]Rik + 
     xs^2*(xt^3*\[Lambda]Rik - (yk^2 - yS^2)^3*\[Lambda]Rik + 
       xt*(5*yk^4 - 6*yk^2*yS^2 + yS^4)*\[Lambda]Rik - 
       xt^2*(4*yk*\[Lambda]Lik + 9*yk^2*\[Lambda]Rik + 3*yS^2*\[Lambda]Rik)))*
    \[Lambda]Rjk*DD[6])/(2*xs^3*xt^3), 
 (-2*Qj*(-2*Qk + Qj*(1 + yk^2 - yS^2))*\[Lambda]Lik*\[Lambda]Ljk)/
   (xs*(-1 + xt)*xt^2) + (2*Qj^2*yk^2*\[Lambda]Lik*\[Lambda]Ljk*BB[1])/
   (xs*xt^2 - xs*xt^3) + (2*Qj^2*yS^2*\[Lambda]Lik*\[Lambda]Ljk*BB[2])/
   (xs*(-1 + xt)*xt^2) + 
  (2*Qj*\[Lambda]Ljk*(2*Qk*(xs - xs^2 + (-1 + xt)^2)*xt*\[Lambda]Lik + 
     Qj*(-1 + xs)*(-1 + xs + xt)*(xt*(1 + yk^2 - yS^2)*\[Lambda]Lik + 
       2*xt*yk*\[Lambda]Rik - 2*(yk^2*\[Lambda]Lik - yS^2*\[Lambda]Lik + 
         yk*\[Lambda]Rik)))*BB[3])/((-1 + xs)*xs*(-1 + xt)^2*xt^2*
    (-1 + xs + xt)) - (4*Qj*Qk*\[Lambda]Lik*\[Lambda]Ljk*BB[4])/
   ((-1 + xs)*xs*xt*(-1 + xs + xt)) - 
  (2*Qj*\[Lambda]Ljk*(-2*Qk*xs*xt*\[Lambda]Lik + Qj*(-1 + xs + xt)*
      (-(yk^2*\[Lambda]Lik) + yS^2*\[Lambda]Lik - 2*yk*\[Lambda]Rik + 
       xt*(\[Lambda]Lik + 2*yk*\[Lambda]Rik)))*BB[5])/
   (xs*(-1 + xt)^2*xt^2*(-1 + xs + xt)) + 
  (2*mi^2*Qk^2*yk*\[Lambda]Ljk*\[Lambda]Rik*CC[1])/(xt^2*(-1 + xs + xt)) - 
  (2*mi^2*(Qj - Qk)*Qk*\[Lambda]Ljk*((-1 + xs)*yk*\[Lambda]Rik + 
     xt*(\[Lambda]Lik + yk*\[Lambda]Rik))*CC[2])/(xt^2*(-1 + xs + xt)^2) - 
  (2*mi^2*Qk*\[Lambda]Ljk*(-(Qk*xt*\[Lambda]Lik) + 
     Qj*(-1 + xs)*yk*\[Lambda]Rik + Qj*xt*(\[Lambda]Lik + yk*\[Lambda]Rik))*
    CC[3])/(xs*xt*(-1 + xs + xt)^2) - 
  (4*mi^2*Qk^2*yk*\[Lambda]Ljk*\[Lambda]Rik*CC[5])/(xs*xt^2) + 
  (2*mi^2*Qk*\[Lambda]Ljk*(Qk*(-1 + xs)*xt*\[Lambda]Lik + 
     Qj*((-1 + xs)^2*yk*\[Lambda]Rik + 2*xt^2*yk*\[Lambda]Rik - 
       (-1 + xs)*xt*(\[Lambda]Lik - 3*yk*\[Lambda]Rik)))*CC[8])/
   (xs*xt^2*(-1 + xs + xt)^2) - (2*mi^2*(Qj - Qk)*\[Lambda]Ljk*
    (2*Qj*(-1 + xs + xt)^2*yS^2*\[Lambda]Lik + Qk*(-1 + xt)^2*
      ((-1 + xs)*yk*\[Lambda]Rik + xt*(\[Lambda]Lik + yk*\[Lambda]Rik)))*
    CC[9])/(xs*(-1 + xt)*xt^2*(-1 + xs + xt)^2) + 
  (2*mi^2*Qk*yk*\[Lambda]Ljk*(Qk*(-1 + xt)^2*\[Lambda]Rik + 
     2*Qj*(-1 + xs + xt)*(yk*\[Lambda]Lik + \[Lambda]Rik - xt*\[Lambda]Rik))*
    CC[10])/(xs*(-1 + xt)*xt^2*(-1 + xs + xt)) + 
  (2*mi^4*(Qj - Qk)*Qk*\[Lambda]Ljk*(xs^2*yk*(xt - yk^2 + yS^2)*
      \[Lambda]Rik + xs*(2*yk*(yk^2 - yS^2)*\[Lambda]Rik + 
       xt^2*(\[Lambda]Lik + yk*\[Lambda]Rik) - 
       xt*(yk^2*\[Lambda]Lik + yS^2*\[Lambda]Lik + yk*\[Lambda]Rik + 
         2*yk^3*\[Lambda]Rik - 2*yk*yS^2*\[Lambda]Rik)) - 
     (-1 + xt)*(yk*(-yk^2 + yS^2)*\[Lambda]Rik + 
       xt*(yk^2*\[Lambda]Lik + yS^2*\[Lambda]Lik + yk^3*\[Lambda]Rik - 
         yk*yS^2*\[Lambda]Rik)))*DD[2])/(xs*xt^2*(-1 + xs + xt)^2) + 
  (2*mi^4*Qk^2*yk*\[Lambda]Ljk*((xs^2 - yk^2 + yS^2 + xs*(-1 + yk^2 - yS^2))*
      \[Lambda]Rik + xt*(2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik + 
       (xs - yS^2)*\[Lambda]Rik))*DD[4])/(xs*xt^2*(-1 + xs + xt)) + 
  (2*mi^4*Qk^2*yk*\[Lambda]Ljk*(-(xt^2*\[Lambda]Rik) + 
     (-1 + xs)*(yk^2 - yS^2)*\[Lambda]Rik + 
     xt*(2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - (-1 + xs + yS^2)*
        \[Lambda]Rik))*DD[6])/(xs*xt^2*(-1 + xs + xt)), 
 (2*Qj*(xs - xt)*(-2*Qk*xs*xt + Qj*(-1 + xt)*(yk^2 - yS^2) + 
     Qj*xs*(xt + yk^2 - yS^2))*\[Lambda]Rik*\[Lambda]Rjk)/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)) + 
  (2*Qj^2*(xs - xt)*yk^2*\[Lambda]Rik*\[Lambda]Rjk*BB[1])/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2) + 
  (2*Qj^2*(-xs + xt)*yS^2*\[Lambda]Rik*\[Lambda]Rjk*BB[2])/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2) + 
  (2*Qj*(-2*xs + xs^2 - (-2 + xt)*xt)*(2*Qk + Qj*(-1 + yk^2 - yS^2))*
    \[Lambda]Rik*\[Lambda]Rjk*BB[3])/((-1 + xs)^2*xs*(-1 + xt)^2*xt*
    (-1 + xs + xt)) + 
  (2*Qj*(2*Qk*xs + Qj*(-yk^2 + yS^2 + xs*(-1 + 2*yk^2 - 2*yS^2)))*
    \[Lambda]Rik*\[Lambda]Rjk*BB[4])/((-1 + xs)^2*xs^2*xt*(-1 + xs + xt)) + 
  (2*Qj*(-2*Qk*xt + Qj*(xt + yk^2 - 2*xt*yk^2 - yS^2 + 2*xt*yS^2))*
    \[Lambda]Rik*\[Lambda]Rjk*BB[5])/(xs*(-1 + xt)^2*xt^2*(-1 + xs + xt)) + 
  (2*mi^2*Qk*(Qk*xs*((-1 + xs)*(yk^2 - yS^2) + xt*(1 + yk^2 - yS^2)) + 
     Qj*((-1 + xt)*xt*yk^2 - xs*yS^2 + xs^2*yS^2 + xs*xt*(-1 + yk^2 + yS^2)))*
    \[Lambda]Rik*\[Lambda]Rjk*CC[1])/(xs*xt^2*(-1 + xs + xt)^3) - 
  (2*mi^2*(Qj^2*(xs^2 + (-1 + xt)*xt + xs*(-1 + 2*xt))*yS^2 + 
     Qk^2*xs*(-((-1 + xs)*(yk^2 - yS^2)) + xt*(1 - yk^2 + yS^2)) + 
     Qj*Qk*(-((-1 + xt)*xt*yS^2) + xs^2*(yk^2 - 2*yS^2) + 
       xs*(-yk^2 + 2*yS^2 + xt*(-1 + yk^2 - 3*yS^2))))*\[Lambda]Rik*
    \[Lambda]Rjk*CC[2])/(xs*xt^2*(-1 + xs + xt)^3) - 
  (2*mi^2*Qk*(Qk*xt*((-1 + xt)*(yk^2 - yS^2) + xs*(1 + yk^2 - yS^2)) + 
     Qj*(-(xs*yk^2) + xs^2*yk^2 + (-1 + xt)*xt*yS^2 + 
       xs*xt*(-1 + yk^2 + yS^2)))*\[Lambda]Rik*\[Lambda]Rjk*CC[3])/
   (xs^2*xt*(-1 + xs + xt)^3) + 
  (2*mi^2*(Qj^2*(xs^2 + (-1 + xt)*xt + xs*(-1 + 2*xt))*yS^2 + 
     Qk^2*xt*(-((-1 + xt)*(yk^2 - yS^2)) + xs*(1 - yk^2 + yS^2)) + 
     Qj*Qk*(-(xs^2*yS^2) + (-1 + xt)*xt*(yk^2 - 2*yS^2) + 
       xs*(yS^2 + xt*(-1 + yk^2 - 3*yS^2))))*\[Lambda]Rik*\[Lambda]Rjk*CC[4])/
   (xs^2*xt*(-1 + xs + xt)^3) - (4*mi^2*Qk^2*(xs - xt)*yk^2*\[Lambda]Rik*
    \[Lambda]Rjk*CC[5])/(xs^2*xt^2*(-1 + xs + xt)) + 
  (4*mi^2*(Qj - Qk)^2*(xs - xt)*yS^2*\[Lambda]Rik*\[Lambda]Rjk*CC[6])/
   (xs^2*xt^2*(-1 + xs + xt)) - 
  (2*mi^2*(Qj - Qk)*(Qj*(xs^4 - (-1 + xt)*xt + xs^3*(-3 + 4*xt) + 
       xs^2*(3 - 7*xt + 5*xt^2) + xs*(-1 + 2*xt - 4*xt^2 + 2*xt^3))*yS^2 - 
     Qk*(-1 + xs)^2*xt*(xs*(-1 + yk^2 - yS^2) + (-1 + xt)*(yk^2 - yS^2)))*
    \[Lambda]Rik*\[Lambda]Rjk*CC[7])/((-1 + xs)*xs^2*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^2*Qk*(-(Qk*(-1 + xs)^2*xt*((-1 + xt)*(yk^2 - yS^2) + 
        xs*(1 + yk^2 - yS^2))) + Qj*(xs^4*yk^2 - (-1 + xt)*xt*yS^2 + 
       xs^3*(xt - 3*yk^2 + 5*xt*yk^2 - xt*yS^2) + 
       xs*(-yk^2 + 2*xt^3*yk^2 + xt*(1 + 5*yk^2 - 3*yS^2) + 
         xt^2*(-6*yk^2 + 2*yS^2)) + xs^2*(3*yk^2 + xt^2*(6*yk^2 - yS^2) + 
         xt*(-2 - 10*yk^2 + 3*yS^2))))*\[Lambda]Rik*\[Lambda]Rjk*CC[8])/
   ((-1 + xs)*xs^2*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^2*(Qj - Qk)*(Qj*(2*xs^3*xt + (-1 + xt)^3*xt + 
       xs*(-1 + xt)^2*(1 + 4*xt) + xs^2*(-1 - 4*xt + 5*xt^2))*yS^2 - 
     Qk*xs*(-1 + xt)^2*(xt*(-1 + yk^2 - yS^2) + (-1 + xs)*(yk^2 - yS^2)))*
    \[Lambda]Rik*\[Lambda]Rjk*CC[9])/(xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)^3) - 
  (2*mi^2*Qk*(-(Qk*xs*(-1 + xt)^2*((-1 + xs)*(yk^2 - yS^2) + 
        xt*(1 + yk^2 - yS^2))) + Qj*(2*xs^3*xt*yk^2 + (-1 + xt)^3*xt*yk^2 + 
       xs*(-1 + xt)^2*(xt + 5*xt*yk^2 + yS^2 - xt*yS^2) + 
       xs^2*(-1 + xt)*(6*xt*yk^2 + yS^2 - xt*yS^2)))*\[Lambda]Rik*
    \[Lambda]Rjk*CC[10])/(xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)^3) - 
  (2*mi^4*(Qj - Qk)*Qk*(-((-1 + xt)^2*xt*yk^2*(yk^2 - yS^2)) + 
     xs^3*yS^2*(xt + yk^2 - yS^2) + xs^2*(-2*yk^2*yS^2 + 2*yS^4 + 
       xt^2*(-1 + yk^2 + yS^2) + xt*(yk^2 - yk^4 + 3*yk^2*yS^2 - 2*yS^4)) + 
     xs*(xt^3*yk^2 + yS^2*(yk^2 - yS^2) + xt^2*(-2*yk^4 + yS^2 + 
         3*yk^2*yS^2 - yS^4) + xt*(2*yk^4 + yS^2*(-1 + 2*yS^2) - 
         yk^2*(1 + 4*yS^2))))*\[Lambda]Rik*\[Lambda]Rjk*DD[1])/
   (xs^2*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^4*(Qj - Qk)*Qk*((-1 + xt)^2*xt*yS^2*(yk^2 - yS^2) + 
     xs^3*yk^2*(xt - yk^2 + yS^2) + xs^2*(2*yk^2*(yk^2 - yS^2) + 
       xt^2*(-1 + yk^2 + yS^2) + xt*(-2*yk^4 + yS^2 + 3*yk^2*yS^2 - yS^4)) + 
     xs*(-yk^4 + xt^3*yS^2 + yk^2*yS^2 + xt^2*(-yk^4 - 2*yS^4 + 
         yk^2*(1 + 3*yS^2)) + xt*(2*yk^4 + yS^2*(-1 + 2*yS^2) - 
         yk^2*(1 + 4*yS^2))))*\[Lambda]Rik*\[Lambda]Rjk*DD[2])/
   (xs^2*xt^2*(-1 + xs + xt)^3) - 
  (2*mi^4*(Qj - Qk)^2*yS^2*(xs^2 + xt*(yk^2 - yS^2) + xs*(xt - yk^2 + yS^2))*
    \[Lambda]Rik*\[Lambda]Rjk*DD[3])/(xs^2*xt^2*(-1 + xs + xt)) + 
  (2*mi^4*Qk^2*yk^2*(xs^2 + xs*(xt + yk^2 - yS^2) + xt*(-yk^2 + yS^2))*
    \[Lambda]Rik*\[Lambda]Rjk*DD[4])/(xs^2*xt^2*(-1 + xs + xt)) + 
  (2*mi^4*(Qj - Qk)^2*yS^2*(xs*(xt + yk^2 - yS^2) + xt*(xt - yk^2 + yS^2))*
    \[Lambda]Rik*\[Lambda]Rjk*DD[5])/(xs^2*xt^2*(-1 + xs + xt)) - 
  (2*mi^4*Qk^2*yk^2*(xt*(xt + yk^2 - yS^2) + xs*(xt - yk^2 + yS^2))*
    \[Lambda]Rik*\[Lambda]Rjk*DD[6])/(xs^2*xt^2*(-1 + xs + xt)), 
 (Qj*(2*Qk*\[Lambda]Rik + Qj*(-1 - 3*yk^2 + 3*yS^2)*\[Lambda]Rik + 
     Qj*xt*(-2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - yS^2*\[Lambda]Rik))*
    \[Lambda]Rjk)/(xs*(-1 + xt)*xt^2) + 
  (Qj^2*yk^2*(3*(-yk^2 + yS^2)*\[Lambda]Rik + 
     xt*(-2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - yS^2*\[Lambda]Rik))*
    \[Lambda]Rjk*BB[1])/(xs*(-1 + xt)*xt^2*(yk^2 - yS^2)) + 
  (Qj^2*yS^2*(3*(yk^2 - yS^2)*\[Lambda]Rik + 
     xt*(2*yk*\[Lambda]Lik - yk^2*\[Lambda]Rik + yS^2*\[Lambda]Rik))*
    \[Lambda]Rjk*BB[2])/(xs*(-1 + xt)*xt^2*(yk^2 - yS^2)) - 
  (Qj*(-2*Qk*(-2*xs*xt + 2*xs^3*xt + xt^2 + xs^2*(1 - 2*xt^2))*\[Lambda]Rik + 
     Qj*(xt^2*(1 - yk^2 + yS^2)*\[Lambda]Rik + 
       xs*xt*(-2*(-1 + xt)*yk*\[Lambda]Lik + xt*yk^2*\[Lambda]Rik + 
         (-2 + xt - xt*yS^2)*\[Lambda]Rik) + 
       xs^3*(2*(-1 + xt)*yk*\[Lambda]Lik + (-2 + xt)*yk^2*\[Lambda]Rik + 
         (xt + 2*yS^2 - xt*yS^2)*\[Lambda]Rik) + 
       xs^2*(2*(-1 + xt)^2*yk*\[Lambda]Lik + (1 - xt + xt^2)*yk^2*
          \[Lambda]Rik + (1 - yS^2 + xt*(1 + yS^2) - xt^2*(3 + yS^2))*
          \[Lambda]Rik)))*\[Lambda]Rjk*BB[3])/((-1 + xs)*xs^2*(-1 + xt)*xt^2*
    (xs + xt)) + (Qj*(2*Qk*xs - Qj*(xs - yk^2 + yS^2))*\[Lambda]Rik*
    \[Lambda]Rjk*BB[4])/((-1 + xs)*xs*xt^2) + 
  (Qj*(2*Qk*xt*(-2*xs + xt)*\[Lambda]Rik - Qj*xt*(xt - yk^2 + yS^2)*
      \[Lambda]Rik + Qj*xs*(2*(-1 + 2*xt)*yk*\[Lambda]Lik + 
       yk^2*\[Lambda]Rik + (xt - yS^2)*\[Lambda]Rik))*\[Lambda]Rjk*BB[5])/
   (xs^2*(-1 + xt)*xt^2) - (Qk^2*(xs - xt)^2*\[Lambda]Rik*\[Lambda]Rjk*BB[6])/
   (xs^2*xt^2*(xs + xt)) + ((Qj - Qk)^2*(xs - xt)^2*\[Lambda]Rik*\[Lambda]Rjk*
    BB[7])/(xs^2*xt^2*(xs + xt)) - 
  (mi^2*Qk*(Qk*(-1 + xs)*xs^2*(4*xt^2*yk*\[Lambda]Lik + 
       (xs^2 - 2*yk^2 + 2*yS^2 + xs*(-1 + 2*yk^2 - 2*yS^2))*\[Lambda]Rik + 
       xt*(4*(-1 + xs)*yk*\[Lambda]Lik + 2*yk^2*\[Lambda]Rik + 
         (-1 + 3*xs - 2*yS^2)*\[Lambda]Rik)) + 
     Qj*((-1 + xt)^2*xt*(yk^2 - yS^2)^2*\[Lambda]Rik + 
       xs^3*((yk^2 - yS^2)^2*\[Lambda]Rik + xt^2*(4*yk*\[Lambda]Lik + 
           3*\[Lambda]Rik) + 2*xt*(2*yk^3*\[Lambda]Lik - 
           2*yk*yS^2*\[Lambda]Lik + yk^2*\[Lambda]Rik - 
           2*yS^2*\[Lambda]Rik)) + xs*(-1 + xt)*
        (-((yk^2 - yS^2)^2*\[Lambda]Rik) + xt*(yk^2 - yS^2)*
          (-4*yk*\[Lambda]Lik - 2*\[Lambda]Rik + 3*yk^2*\[Lambda]Rik - 
           3*yS^2*\[Lambda]Rik) + xt^2*(4*yk^3*\[Lambda]Lik - 
           4*yk*yS^2*\[Lambda]Lik + 4*yk^2*\[Lambda]Rik - 
           2*yS^2*\[Lambda]Rik)) + xs^2*(-2*(yk^2 - yS^2)^2*\[Lambda]Rik + 
         xt^3*(4*yk*\[Lambda]Lik + \[Lambda]Rik) + 
         xt*(-8*yk^3*\[Lambda]Lik + 8*yk*yS^2*\[Lambda]Lik + 
           3*yk^4*\[Lambda]Rik + 3*yS^2*(2 + yS^2)*\[Lambda]Rik - 
           2*yk^2*(2 + 3*yS^2)*\[Lambda]Rik) + xt^2*(8*yk^3*\[Lambda]Lik - 
           4*yk*(\[Lambda]Lik + 2*yS^2*\[Lambda]Lik) + 6*yk^2*\[Lambda]Rik - 
           2*(\[Lambda]Rik + 3*yS^2*\[Lambda]Rik)))))*\[Lambda]Rjk*CC[1])/
   (2*xs^2*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*(Qj - Qk)*(-(Qk*(-1 + xs)*xs^2*
       ((xs^2 + 2*(yk^2 - yS^2) + xs*(-1 - 2*yk^2 + 2*yS^2))*\[Lambda]Rik - 
        2*xt^2*(2*yk*\[Lambda]Lik + \[Lambda]Rik) + 
        xt*(-4*(-1 + xs)*yk*\[Lambda]Lik - 2*yk^2*\[Lambda]Rik + 
          (1 + xs + 2*yS^2)*\[Lambda]Rik))) + Qj*(-1 + xs + xt)*
      (xs^4*\[Lambda]Rik + (-1 + xt)*xt*(yk^2 - yS^2)^2*\[Lambda]Rik - 
       xs^3*(4*xt*yk*\[Lambda]Lik + \[Lambda]Rik + 2*yk^2*\[Lambda]Rik - 
         2*yS^2*\[Lambda]Rik) + 
       xs^2*((yk^4 + yS^2*(-2 + yS^2) - 2*yk^2*(-1 + yS^2))*\[Lambda]Rik - 
         xt^2*(4*yk*\[Lambda]Lik + \[Lambda]Rik) + 
         xt*(4*yk^3*\[Lambda]Lik - 4*yk*(-1 + yS^2)*\[Lambda]Lik + 
           \[Lambda]Rik - 2*yk^2*\[Lambda]Rik)) + 
       xs*(-((yk^2 - yS^2)^2*\[Lambda]Rik) + 2*xt*(yk^2 - yS^2)*
          (-2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - yS^2*\[Lambda]Rik) + 
         2*xt^2*(2*yk^3*\[Lambda]Lik - 2*yk*yS^2*\[Lambda]Lik + 
           yS^2*\[Lambda]Rik))))*\[Lambda]Rjk*CC[2])/
   (2*xs^2*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*Qk*(-(Qk*(-1 + xt)*xt*((-1 + xt)^2*xt*(xt + 2*yk^2 - 2*yS^2)*
         \[Lambda]Rik + 2*xs^3*(yk^2 - yS^2)*\[Lambda]Rik + 
        2*xs^2*(-1 + xt)*(2*xt*yk*\[Lambda]Lik - xt*\[Lambda]Rik + 
          2*yk^2*\[Lambda]Rik - 2*yS^2*\[Lambda]Rik) + 
        xs*(-1 + xt)*(2*(-yk^2 + yS^2)*\[Lambda]Rik + 
          xt^2*(4*yk*\[Lambda]Lik + \[Lambda]Rik) + 
          xt*(-4*yk*\[Lambda]Lik + \[Lambda]Rik + 4*yk^2*\[Lambda]Rik - 
            4*yS^2*\[Lambda]Rik)))) + 
     Qj*(-((-1 + xt)^3*xt*(yk^2 - yS^2)^2*\[Lambda]Rik) + 
       4*xs^4*xt*yk*((-1 + 2*xt)*\[Lambda]Lik + yk*\[Lambda]Rik) + 
       xs^3*(-1 + xt)*(-((yk^2 - yS^2)^2*\[Lambda]Rik) + 
         xt^2*(12*yk*\[Lambda]Lik + \[Lambda]Rik) + 
         2*xt*yk*(-2*(2 + yk^2 - yS^2)*\[Lambda]Lik + 3*yk*\[Lambda]Rik)) - 
       xs*(-1 + xt)^2*(-((yk^2 - yS^2)^2*\[Lambda]Rik) + 
         xt*(yk^2 - yS^2)*(-4*yk*\[Lambda]Lik + 3*yk^2*\[Lambda]Rik - 
           3*yS^2*\[Lambda]Rik) + xt^2*(4*yk^3*\[Lambda]Lik - 
           4*yk*yS^2*\[Lambda]Lik - 2*yS^2*\[Lambda]Rik)) + 
       xs^2*(-1 + xt)*(xt^3*(4*yk*\[Lambda]Lik - \[Lambda]Rik) + 
         2*(yk^2 - yS^2)^2*\[Lambda]Rik + xt^2*(-8*yk^3*\[Lambda]Lik + 
           8*yk*(-1 + yS^2)*\[Lambda]Lik + 2*yk^2*\[Lambda]Rik + 
           2*yS^2*\[Lambda]Rik) + xt*(8*yk^3*\[Lambda]Lik + 
           yk*(4 - 8*yS^2)*\[Lambda]Lik - 3*yk^4*\[Lambda]Rik - 
           3*yS^4*\[Lambda]Rik + 2*yk^2*(-1 + 3*yS^2)*\[Lambda]Rik))))*
    \[Lambda]Rjk*CC[3])/(2*xs^3*(-1 + xt)*xt^2*(-1 + xs + xt)^2) - 
  (mi^2*(Qj - Qk)*(Qk*(-1 + xt)*xt*(2*xs^3*(xt - yk^2 + yS^2)*\[Lambda]Rik + 
       (-1 + xt)^2*xt*(xt - 2*yk^2 + 2*yS^2)*\[Lambda]Rik - 
       4*xs^2*(-1 + xt)*(xt*yk*\[Lambda]Lik + (yk^2 - yS^2)*\[Lambda]Rik) - 
       xs*(-1 + xt)*(xt^2*(4*yk*\[Lambda]Lik - \[Lambda]Rik) + 
         2*(-yk^2 + yS^2)*\[Lambda]Rik + xt*(-4*yk*\[Lambda]Lik - 
           \[Lambda]Rik + 4*yk^2*\[Lambda]Rik - 4*yS^2*\[Lambda]Rik))) + 
     Qj*(-1 + xs + xt)*(4*xs^3*xt*yS^2*\[Lambda]Rik - 
       xt*(xt^2 + yk^2 - yS^2 + xt*(-1 - yk^2 + yS^2))^2*\[Lambda]Rik + 
       xs*(-1 + xt)*(4*xt^3*yk*\[Lambda]Lik + (yk^2 - yS^2)^2*\[Lambda]Rik - 
         xt^2*(4*yk^3*\[Lambda]Lik - 4*yk*(-1 + yS^2)*\[Lambda]Lik + 
           \[Lambda]Rik - 2*yk^2*\[Lambda]Rik) - 2*xt*(yk^2 - yS^2)*
          (-2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - yS^2*\[Lambda]Rik)) + 
       xs^2*(-1 + xt)*(-((yk^2 - yS^2)^2*\[Lambda]Rik) + 
         xt^2*(4*yk*\[Lambda]Lik + \[Lambda]Rik) + 
         xt*(-4*yk^3*\[Lambda]Lik + 4*yk*yS^2*\[Lambda]Lik + 
           2*yS^2*\[Lambda]Rik))))*\[Lambda]Rjk*CC[4])/
   (2*xs^3*(-1 + xt)*xt^2*(-1 + xs + xt)^2) + 
  (mi^2*Qk^2*(xs^4*\[Lambda]Rik + xt*(xt^3 + xt^2*(-1 + 2*yk^2 - 2*yS^2) - 
       2*(yk^2 - yS^2)^2 + 2*xt*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2)))*
      \[Lambda]Rik + xs^3*((-1 + 2*yk^2 - 2*yS^2)*\[Lambda]Rik + 
       2*xt*(2*yk*\[Lambda]Lik + \[Lambda]Rik)) + 
     xs*(4*xt^3*yk*\[Lambda]Lik - 2*(yk^2 - yS^2)^2*\[Lambda]Rik + 
       xt^2*(8*yk^3*\[Lambda]Lik - 4*yk*(\[Lambda]Lik + 
           2*yS^2*\[Lambda]Lik) + \[Lambda]Rik + 2*yk^2*\[Lambda]Rik - 
         6*yS^2*\[Lambda]Rik) + 4*xt*(yk^2 - yS^2)*(-2*yk*\[Lambda]Lik + 
         yk^2*\[Lambda]Rik - (1 + yS^2)*\[Lambda]Rik)) + 
     xs^2*(8*xt^2*yk*\[Lambda]Lik + 2*(yk^4 + yS^2 + yS^4 - 
         yk^2*(1 + 2*yS^2))*\[Lambda]Rik + xt*(8*yk^3*\[Lambda]Lik - 
         4*yk*(\[Lambda]Lik + 2*yS^2*\[Lambda]Lik) + 2*yk^2*\[Lambda]Rik - 
         (1 + 6*yS^2)*\[Lambda]Rik)))*\[Lambda]Rjk*CC[5])/(2*xs^3*xt^3) - 
  (mi^2*(Qj - Qk)^2*(xs^4*\[Lambda]Rik + xt*(xt^3 - 2*(yk^2 - yS^2)^2 + 
       xt^2*(-1 - 2*yk^2 + 2*yS^2) + 2*xt*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + 
         yS^4))*\[Lambda]Rik - xs^3*(4*xt*yk*\[Lambda]Lik + \[Lambda]Rik + 
       2*yk^2*\[Lambda]Rik - 2*yS^2*\[Lambda]Rik) + 
     xs^2*(2*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4)*\[Lambda]Rik - 
       2*xt^2*(4*yk*\[Lambda]Lik + \[Lambda]Rik) + 
       xt*(8*yk^3*\[Lambda]Lik + yk*(4 - 8*yS^2)*\[Lambda]Lik + 
         \[Lambda]Rik - 2*yk^2*\[Lambda]Rik - 2*yS^2*\[Lambda]Rik)) + 
     xs*(-4*xt^3*yk*\[Lambda]Lik - 2*(yk^2 - yS^2)^2*\[Lambda]Rik + 
       xt^2*(8*yk^3*\[Lambda]Lik + yk*(4 - 8*yS^2)*\[Lambda]Lik + 
         \[Lambda]Rik - 2*yk^2*\[Lambda]Rik - 2*yS^2*\[Lambda]Rik) + 
       4*xt*(yk^2 - yS^2)*(-2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - 
         yS^2*\[Lambda]Rik)))*\[Lambda]Rjk*CC[6])/(2*xs^3*xt^3) + 
  (mi^2*(Qj - Qk)*(-(Qk*(-1 + xs)*xs*(xs^4*\[Lambda]Rik + 
        2*(-1 + xt)*xt*(yk^2 - yS^2)*\[Lambda]Rik + 
        xs^3*(-2*(1 + yk^2 - yS^2)*\[Lambda]Rik + 
          xt*(-4*yk*\[Lambda]Lik + \[Lambda]Rik)) + 
        xs^2*((1 + 4*yk^2 - 4*yS^2)*\[Lambda]Rik + 
          xt^2*(-4*yk*\[Lambda]Lik + 2*\[Lambda]Rik) + 
          xt*(8*yk*\[Lambda]Lik - 6*yk^2*\[Lambda]Rik + 
            6*yS^2*\[Lambda]Rik)) - xs*(2*(yk^2 - yS^2)*\[Lambda]Rik + 
          4*xt^2*(-(yk*\[Lambda]Lik) + yk^2*\[Lambda]Rik - 
            yS^2*\[Lambda]Rik) + xt*(4*yk*\[Lambda]Lik + \[Lambda]Rik - 
            8*yk^2*\[Lambda]Rik + 8*yS^2*\[Lambda]Rik)))) + 
     Qj*(-1 + xs + xt)*(xs^5*\[Lambda]Rik - (-1 + xt)*xt*(yk^2 - yS^2)^2*
        \[Lambda]Rik - 2*xs^4*(2*xt*yk*\[Lambda]Lik + (1 + yk^2 - yS^2)*
          \[Lambda]Rik) + xs^3*((1 + yk^4 - 4*yS^2 + yS^4 - 
           2*yk^2*(-2 + yS^2))*\[Lambda]Rik - xt^2*(4*yk*\[Lambda]Lik + 
           \[Lambda]Rik) + xt*(4*yk^3*\[Lambda]Lik - 4*yk*(-2 + yS^2)*
            \[Lambda]Lik + \[Lambda]Rik - 2*yk^2*\[Lambda]Rik)) + 
       xs*((yk^2 - yS^2)^2*\[Lambda]Rik - xt*(yk^2 - yS^2)*
          (-4*yk*\[Lambda]Lik + 3*yk^2*\[Lambda]Rik - 3*yS^2*\[Lambda]Rik) + 
         xt^2*(-4*yk^3*\[Lambda]Lik + 4*yk*yS^2*\[Lambda]Lik + 
           yk^4*\[Lambda]Rik - 2*yk^2*yS^2*\[Lambda]Rik + 
           yS^2*(2 + yS^2)*\[Lambda]Rik)) + 
       xs^2*(-2*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4)*\[Lambda]Rik + 
         xt^2*(4*yk^3*\[Lambda]Lik - 4*yk*(-1 + yS^2)*\[Lambda]Lik + 
           \[Lambda]Rik - 6*yS^2*\[Lambda]Rik) + 
         xt*(-8*yk^3*\[Lambda]Lik + 4*yk*(-1 + 2*yS^2)*\[Lambda]Lik + 
           2*yk^4*\[Lambda]Rik + yk^2*(2 - 4*yS^2)*\[Lambda]Rik + 
           (-1 + 2*yS^4)*\[Lambda]Rik))))*\[Lambda]Rjk*CC[7])/
   (2*xs^3*xt^3*(-1 + xs + xt)^2) - 
  (mi^2*Qk*(Qk*(-1 + xs)*xs*(xs^4*\[Lambda]Rik - 2*(-1 + xt)*xt*(yk^2 - yS^2)*
        \[Lambda]Rik + xs^3*(2*(-1 + yk^2 - yS^2)*\[Lambda]Rik + 
         xt*(4*yk*\[Lambda]Lik + 3*\[Lambda]Rik)) + 
       xs^2*(\[Lambda]Rik - 4*yk^2*\[Lambda]Rik + 4*yS^2*\[Lambda]Rik + 
         4*xt^2*(yk*\[Lambda]Lik + \[Lambda]Rik) + 
         xt*(-8*yk*\[Lambda]Lik - 4*\[Lambda]Rik + 6*yk^2*\[Lambda]Rik - 
           6*yS^2*\[Lambda]Rik)) + xs*(2*(yk^2 - yS^2)*\[Lambda]Rik + 
         xt^2*(-4*yk*\[Lambda]Lik - 2*\[Lambda]Rik + 4*yk^2*\[Lambda]Rik - 
           4*yS^2*\[Lambda]Rik) + xt*(4*yk*\[Lambda]Lik + \[Lambda]Rik - 
           8*yk^2*\[Lambda]Rik + 8*yS^2*\[Lambda]Rik))) + 
     Qj*(-((-1 + xt)^2*xt*(yk^2 - yS^2)^2*\[Lambda]Rik) + 
       xs^4*(xt^2*(4*yk*\[Lambda]Lik - \[Lambda]Rik) + (yk^2 - yS^2)^2*
          \[Lambda]Rik + 2*xt*yk*(2*yk^2*\[Lambda]Lik - 2*yS^2*\[Lambda]Lik - 
           yk*\[Lambda]Rik)) + xs^3*(-3*(yk^2 - yS^2)^2*\[Lambda]Rik + 
         xt^3*(4*yk*\[Lambda]Lik + \[Lambda]Rik) + 
         xt^2*(8*yk^3*\[Lambda]Lik - 8*yk*(1 + yS^2)*\[Lambda]Lik + 
           \[Lambda]Rik - 6*yk^2*\[Lambda]Rik - 2*yS^2*\[Lambda]Rik) + 
         xt*(-12*yk^3*\[Lambda]Lik + 12*yk*yS^2*\[Lambda]Lik + 
           3*yk^4*\[Lambda]Rik + 3*yS^4*\[Lambda]Rik + yk^2*(4 - 6*yS^2)*
            \[Lambda]Rik)) + xs*(-1 + xt)*((yk^2 - yS^2)^2*\[Lambda]Rik - 
         4*xt*(yk^2 - yS^2)*(-(yk*\[Lambda]Lik) + yk^2*\[Lambda]Rik - 
           yS^2*\[Lambda]Rik) + xt^2*(-4*yk^3*\[Lambda]Lik + 
           4*yk*yS^2*\[Lambda]Lik + yk^4*\[Lambda]Rik - 
           2*yk^2*yS^2*\[Lambda]Rik + yS^2*(2 + yS^2)*\[Lambda]Rik)) + 
       xs^2*(3*(yk^2 - yS^2)^2*\[Lambda]Rik - xt*(-12*yk^3*\[Lambda]Lik + 
           12*yk*yS^2*\[Lambda]Lik + 7*yk^4*\[Lambda]Rik + 
           7*yS^4*\[Lambda]Rik + 2*yk^2*(1 - 7*yS^2)*\[Lambda]Rik) + 
         xt^3*(4*yk^3*\[Lambda]Lik - 4*yk*(1 + yS^2)*\[Lambda]Lik - 
           4*yk^2*\[Lambda]Rik - (1 + 2*yS^2)*\[Lambda]Rik) + 
         xt^2*(-16*yk^3*\[Lambda]Lik + 4*yk*(\[Lambda]Lik + 
             4*yS^2*\[Lambda]Lik) + 3*yk^4*\[Lambda]Rik - 6*yk^2*(-1 + yS^2)*
            \[Lambda]Rik + yS^2*(4 + 3*yS^2)*\[Lambda]Rik))))*\[Lambda]Rjk*
    CC[8])/(2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*(Qj - Qk)*(-(Qk*(-1 + xt)*xt*(-2*xs^3*(xt - yk^2 + yS^2)*
         \[Lambda]Rik + (-1 + xt)^2*xt*(xt - 2*yk^2 + 2*yS^2)*\[Lambda]Rik - 
        xs*(-1 + xt)*xt*(4*(-1 + xt)*yk*\[Lambda]Lik + 4*yk^2*\[Lambda]Rik - 
          (1 + xt + 4*yS^2)*\[Lambda]Rik) + xs^2*(-4*xt^2*yk*\[Lambda]Lik + 
          2*(-yk^2 + yS^2)*\[Lambda]Rik + 2*xt*(2*yk*\[Lambda]Lik + 
            \[Lambda]Rik)))) + Qj*(-1 + xs + xt)*
      (4*xs^3*xt*yS^2*\[Lambda]Rik + 
       xt*(xt^2 + yk^2 - yS^2 + xt*(-1 - yk^2 + yS^2))^2*\[Lambda]Rik - 
       xs*(-1 + xt)*(4*xt^3*yk*\[Lambda]Lik + (yk^2 - yS^2)^2*\[Lambda]Rik - 
         xt^2*(4*yk^3*\[Lambda]Lik - 4*yk*(-1 + yS^2)*\[Lambda]Lik + 
           \[Lambda]Rik - 2*yk^2*\[Lambda]Rik) - 2*xt*(yk^2 - yS^2)*
          (-2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - yS^2*\[Lambda]Rik)) - 
       xs^2*((yk^2 - yS^2)^2*\[Lambda]Rik + xt^3*(4*yk*\[Lambda]Lik + 
           \[Lambda]Rik) - xt^2*(4*yk^3*\[Lambda]Lik - 4*yk*(-1 + yS^2)*
            \[Lambda]Lik + \[Lambda]Rik - 2*yS^2*\[Lambda]Rik) - 
         xt*(-4*yk^3*\[Lambda]Lik + 4*yk*yS^2*\[Lambda]Lik + 
           yk^4*\[Lambda]Rik - 2*yk^2*yS^2*\[Lambda]Rik + yS^2*(-2 + yS^2)*
            \[Lambda]Rik))))*\[Lambda]Rjk*CC[9])/
   (2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*Qk*(-(Qk*(-1 + xt)*xt*((-1 + xt)^2*xt*(xt + 2*yk^2 - 2*yS^2)*
         \[Lambda]Rik - 2*xs^3*(2*xt + yk^2 - yS^2)*\[Lambda]Rik + 
        xs*(-1 + xt)*xt*(4*(-1 + xt)*yk*\[Lambda]Lik + 4*yk^2*\[Lambda]Rik + 
          (1 + xt - 4*yS^2)*\[Lambda]Rik) + 
        2*xs^2*(xt^2*(2*yk*\[Lambda]Lik - \[Lambda]Rik) + 
          (yk^2 - yS^2)*\[Lambda]Rik + xt*(-2*yk*\[Lambda]Lik + 
            2*\[Lambda]Rik)))) + 
     Qj*(-((-1 + xt)^3*xt*(yk^2 - yS^2)^2*\[Lambda]Rik) + 
       4*xs^4*xt*yk*((-1 + 2*xt)*\[Lambda]Lik - yk*\[Lambda]Rik) - 
       xs*(-1 + xt)^2*(-((yk^2 - yS^2)^2*\[Lambda]Rik) + 
         xt^2*(4*yk^3*\[Lambda]Lik - 4*yk*yS^2*\[Lambda]Lik - 
           2*yS^2*\[Lambda]Rik) + xt*(yk^2 - yS^2)*(-4*yk*\[Lambda]Lik + 
           3*yk^2*\[Lambda]Rik - (2 + 3*yS^2)*\[Lambda]Rik)) + 
       xs^2*(-1 + xt)*(xt^3*(4*yk*\[Lambda]Lik - \[Lambda]Rik) + 
         2*(yk^2 - yS^2)^2*\[Lambda]Rik + xt^2*(-8*yk^3*\[Lambda]Lik + 
           8*yk*(-1 + yS^2)*\[Lambda]Lik - 2*yk^2*\[Lambda]Rik + 
           2*(1 + 3*yS^2)*\[Lambda]Rik) + xt*(8*yk^3*\[Lambda]Lik + 
           yk*(4 - 8*yS^2)*\[Lambda]Lik - 3*yk^4*\[Lambda]Rik - 
           3*yS^2*(2 + yS^2)*\[Lambda]Rik + 2*yk^2*(4 + 3*yS^2)*
            \[Lambda]Rik)) + xs^3*(3*xt^3*(4*yk*\[Lambda]Lik - 
           \[Lambda]Rik) + (yk^2 - yS^2)^2*\[Lambda]Rik - 
         xt*(-4*yk^3*\[Lambda]Lik + 4*yk*(-2 + yS^2)*\[Lambda]Lik + 
           yk^4*\[Lambda]Rik + yS^2*(4 + yS^2)*\[Lambda]Rik - 
           2*yk^2*(5 + yS^2)*\[Lambda]Rik) + xt^2*(-4*yk^3*\[Lambda]Lik + 
           4*yk*(-5 + yS^2)*\[Lambda]Lik - 6*yk^2*\[Lambda]Rik + 
           (3 + 4*yS^2)*\[Lambda]Rik))))*\[Lambda]Rjk*CC[10])/
   (2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*Qk^2*(xs^5*\[Lambda]Rik + xt^4*(xt + 2*yk^2 - 2*yS^2)*\[Lambda]Rik + 
     xs^4*(4*xt*yk*\[Lambda]Lik + 3*xt*\[Lambda]Rik + 2*yk^2*\[Lambda]Rik - 
       2*yS^2*\[Lambda]Rik) + 4*xs^2*xt^2*(xt*yk*\[Lambda]Lik + 
       (-yk^2 + yS^2)*\[Lambda]Rik) + 2*xs^3*xt*
      (2*(yk^2 - yS^2)*\[Lambda]Rik + xt*(2*yk*\[Lambda]Lik + 
         \[Lambda]Rik)) + xs*xt^3*(4*(yk^2 - yS^2)*\[Lambda]Rik + 
       xt*(4*yk*\[Lambda]Lik + \[Lambda]Rik)))*\[Lambda]Rjk*CC[11])/
   (2*xs^3*xt^3*(xs + xt)) - 
  (mi^2*(Qj - Qk)^2*(xs^5*\[Lambda]Rik + xt^4*(xt - 2*yk^2 + 2*yS^2)*
      \[Lambda]Rik + xs^4*(2*(-yk^2 + yS^2)*\[Lambda]Rik + 
       xt*(-4*yk*\[Lambda]Lik + \[Lambda]Rik)) + 
     xs*xt^3*(4*(-yk^2 + yS^2)*\[Lambda]Rik + 
       xt*(-4*yk*\[Lambda]Lik + \[Lambda]Rik)) - 
     2*xs^3*xt*(2*(yk^2 - yS^2)*\[Lambda]Rik + 
       xt*(2*yk*\[Lambda]Lik + \[Lambda]Rik)) - 
     2*xs^2*xt^2*(2*(-yk^2 + yS^2)*\[Lambda]Rik + 
       xt*(2*yk*\[Lambda]Lik + \[Lambda]Rik)))*\[Lambda]Rjk*CC[12])/
   (2*xs^3*xt^3*(xs + xt)) + 
  (mi^4*(Qj - Qk)*Qk*((-1 + xt)^2*(yk^2 - yS^2)^2 + 
     xs^2*(xt^2 + (yk^2 - yS^2)^2 - 2*xt*(yk^2 + yS^2)) - 
     2*xs*((yk^2 - yS^2)^2 + xt^2*(yk^2 + yS^2) - 
       xt*(yk^2 + yk^4 + yS^2 - 2*yk^2*yS^2 + yS^4)))*
    ((-1 + xt)*xt*(yk^2 - yS^2)*\[Lambda]Rik + 
     xs^2*(4*xt*yk*\[Lambda]Lik + 3*xt*\[Lambda]Rik + yk^2*\[Lambda]Rik - 
       yS^2*\[Lambda]Rik) + xs*((-yk^2 + yS^2)*\[Lambda]Rik + 
       xt^2*(4*yk*\[Lambda]Lik + \[Lambda]Rik) + 
       2*xt*(-2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - 
         (1 + yS^2)*\[Lambda]Rik)))*\[Lambda]Rjk*DD[1])/
   (2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^4*(Qj - Qk)*Qk*(xs*(xt - (yk - yS)^2) - (-1 + xt)*(yk - yS)^2)*
    (-((-1 + xt)*(yk + yS)^2) + xs*(xt - (yk + yS)^2))*
    ((-1 + xt)*xt*(yk^2 - yS^2)*\[Lambda]Rik + 
     xs^2*(4*xt*yk*\[Lambda]Lik - xt*\[Lambda]Rik + yk^2*\[Lambda]Rik - 
       yS^2*\[Lambda]Rik) + xs*((-yk^2 + yS^2)*\[Lambda]Rik + 
       xt^2*(4*yk*\[Lambda]Lik + \[Lambda]Rik) + 
       2*xt*(-2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - yS^2*\[Lambda]Rik)))*
    \[Lambda]Rjk*DD[2])/(2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^4*(Qj - Qk)^2*(xs^3 + (-1 + xt)*(yk^2 - yS^2)^2 + 
     xs^2*(-1 + xt - 2*yk^2 + 2*yS^2) + xs*(yk^4 + yS^2*(-2 - 2*xt + yS^2) - 
       2*yk^2*(-1 + xt + yS^2)))*(xs^2*\[Lambda]Rik + 
     xt*(-yk^2 + yS^2)*\[Lambda]Rik - xs*((yk^2 - yS^2)*\[Lambda]Rik + 
       xt*(4*yk*\[Lambda]Lik + \[Lambda]Rik)))*\[Lambda]Rjk*DD[3])/
   (2*xs^3*xt^3) - (mi^4*Qk^2*(xs^3 + xs^2*(-1 + xt + 2*yk^2 - 2*yS^2) + 
     (-1 + xt)*(yk^2 - yS^2)^2 + xs*(yk^4 + yS^2*(2 - 2*xt + yS^2) - 
       2*yk^2*(1 + xt + yS^2)))*(xs^2*\[Lambda]Rik + 
     xs*(yk^2 - yS^2)*\[Lambda]Rik + xt*(yk^2 - yS^2)*\[Lambda]Rik + 
     xs*xt*(4*yk*\[Lambda]Lik + \[Lambda]Rik))*\[Lambda]Rjk*DD[4])/
   (2*xs^3*xt^3) + (mi^4*(Qj - Qk)^2*(xt^3 + (-1 + xs)*(yk^2 - yS^2)^2 + 
     xt^2*(-1 + xs - 2*yk^2 + 2*yS^2) + xt*(yk^4 + yS^2*(-2 - 2*xs + yS^2) - 
       2*yk^2*(-1 + xs + yS^2)))*(xt*(xt - yk^2 + yS^2)*\[Lambda]Rik - 
     xs*((yk^2 - yS^2)*\[Lambda]Rik + xt*(4*yk*\[Lambda]Lik + \[Lambda]Rik)))*
    \[Lambda]Rjk*DD[5])/(2*xs^3*xt^3) + 
  (mi^4*Qk^2*(-((-1 + xt)*xt*(xt + yk^2 - yS^2)^3*\[Lambda]Rik) + 
     xs^2*(-((yk^2 - yS^2)^3*\[Lambda]Rik) + 
       xt^3*(4*yk*\[Lambda]Lik + \[Lambda]Rik) - xt*(yk^2 - yS^2)*
        (4*yk^3*\[Lambda]Lik - 4*yk*yS^2*\[Lambda]Lik - yk^2*\[Lambda]Rik - 
         3*yS^2*\[Lambda]Rik) + xt^2*(8*yk^3*\[Lambda]Lik + 
         4*yk*(-1 + 2*yS^2)*\[Lambda]Lik + 3*yk^2*\[Lambda]Rik + 
         yS^2*\[Lambda]Rik)) - xs*(xt + yk^2 - yS^2)*
      (4*xt^3*yk*\[Lambda]Lik - (yk^2 - yS^2)^2*\[Lambda]Rik + 
       xt^2*(4*yk^3*\[Lambda]Lik - 4*yk*(1 + yS^2)*\[Lambda]Lik + 
         \[Lambda]Rik - 4*yS^2*\[Lambda]Rik) + 2*xt*(yk^2 - yS^2)*
        (-2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - (1 + yS^2)*\[Lambda]Rik)))*
    \[Lambda]Rjk*DD[6])/(2*xs^3*xt^3), 
 (2*Qj*\[Lambda]Ljk*(2*Qk*(-1 + 2*xt)*\[Lambda]Lik + 
     Qj*(1 - yk^2 + yS^2 - 2*xt*(1 + yk^2 - yS^2))*\[Lambda]Lik - 
     4*Qj*xt*yk*\[Lambda]Rik))/(xs*(-1 + xt)*xt^2) - 
  (2*Qj^2*yk^2*\[Lambda]Ljk*(-((1 + 2*xt)*yS^2*\[Lambda]Lik) + 
     yk^2*(\[Lambda]Lik + 2*xt*\[Lambda]Lik) + 4*xt*yk*\[Lambda]Rik)*BB[1])/
   (xs*(-1 + xt)*xt^2*(yk^2 - yS^2)) + 
  (2*Qj^2*yS^2*\[Lambda]Ljk*(-((1 + 2*xt)*yS^2*\[Lambda]Lik) + 
     yk^2*(\[Lambda]Lik + 2*xt*\[Lambda]Lik) + 4*xt*yk*\[Lambda]Rik)*BB[2])/
   (xs*(-1 + xt)*xt^2*(yk^2 - yS^2)) + 
  (2*Qj*\[Lambda]Ljk*(2*Qk*(-((-2 + xt)*(-1 + xt)^2) + xs^2*xt + 
       xs*(-2 + 3*xt - 2*xt^2))*\[Lambda]Lik + Qj*(-1 + xs + xt)*
      ((2 + 2*(-2 + xs)*yk^2 + 4*yS^2 - 2*xs*yS^2 + 
         2*xt^2*(-2*yk^2 + 2*yS^2 + xs*(1 + yk^2 - yS^2)) - 
         xt*(1 - 7*yk^2 + 7*yS^2 + 3*xs*(1 + yk^2 - yS^2)))*\[Lambda]Lik + 
       2*(-1 + xs)*(1 - 3*xt + 2*xt^2)*yk*\[Lambda]Rik))*BB[3])/
   ((-1 + xs)*xs*(-1 + xt)^2*xt^2*(-1 + xs + xt)) - 
  (4*Qj*(-(Qk*(-2*xs + 2*xs^2 + xt)) + Qj*(-1 + xs + xt)*(xs - yk^2 + yS^2))*
    \[Lambda]Lik*\[Lambda]Ljk*BB[4])/((-1 + xs)*xs*xt^2*(-1 + xs + xt)) - 
  (2*Qj*\[Lambda]Ljk*(2*Qk*(xs - 2*(-1 + xt)^2)*xt*\[Lambda]Lik + 
     Qj*(-1 + xs + xt)*(2*xt^2*\[Lambda]Lik + 3*yk^2*\[Lambda]Lik - 
       3*yS^2*\[Lambda]Lik + 2*yk*\[Lambda]Rik - 
       xt*((3 + 2*yk^2 - 2*yS^2)*\[Lambda]Lik + 2*yk*\[Lambda]Rik)))*BB[5])/
   (xs*(-1 + xt)^2*xt^2*(-1 + xs + xt)) - 
  (4*Qk^2*\[Lambda]Lik*\[Lambda]Ljk*BB[6])/(xs*xt^2) + 
  (4*(Qj - Qk)^2*\[Lambda]Lik*\[Lambda]Ljk*BB[7])/(xs*xt^2) - 
  (2*mi^2*Qk*\[Lambda]Ljk*(Qj*((-1 + xt)^2*(yk^2 - yS^2)^2 + 
       xs^2*(xt^2 - 2*xt*yS^2 + (yk^2 - yS^2)^2) + 
       2*xs*((-1 + xt)*yk^4 - 2*(-1 + xt)*yk^2*yS^2 + 
         yS^2*(xt - xt^2 - yS^2 + xt*yS^2)))*\[Lambda]Lik + 
     Qk*xs*(xs^3*\[Lambda]Lik + 2*xs^2*(-1 + xt + yk^2 - yS^2)*\[Lambda]Lik + 
       xs*(1 - 4*yk^2 + 4*yS^2 + 2*xt*(-1 + yk^2 - yS^2))*\[Lambda]Lik - 
       xs*xt*yk*\[Lambda]Rik - (-1 + xt)*(2*yk^2*\[Lambda]Lik - 
         2*yS^2*\[Lambda]Lik + xt*yk*\[Lambda]Rik)))*CC[1])/
   (xs*xt^3*(-1 + xs + xt)^2) + (2*mi^2*(Qj - Qk)*\[Lambda]Ljk*
    (Qj*(xs^4 + (-1 + xt)^2*(yk^2 - yS^2)^2 + 
       2*xs^3*(-1 + xt - yk^2 + yS^2) + xs^2*(1 + xt^2 + yk^4 - 4*yS^2 + 
         yS^4 - 2*yk^2*(-2 + yS^2) + xt*(-2 - 4*yk^2 + 2*yS^2)) + 
       2*xs*((-1 + xt)*yk^4 + (-1 + xt)*yS^2*(-1 + yS^2) - 
         yk^2*(1 + xt^2 - 2*yS^2 + 2*xt*(-1 + yS^2))))*\[Lambda]Lik - 
     Qk*xs*(xs^3*\[Lambda]Lik + 2*(-yk^2 + yS^2)*\[Lambda]Lik + 
       2*xs^2*(-1 + xt - yk^2 + yS^2)*\[Lambda]Lik + 
       xs*(1 + 2*xt^2 + 4*yk^2 - 4*yS^2 - 2*xt*(1 + yk^2 - yS^2))*
        \[Lambda]Lik + xs*xt*yk*\[Lambda]Rik + 
       xt*(2*yk^2*\[Lambda]Lik - 2*yS^2*\[Lambda]Lik - yk*\[Lambda]Rik) + 
       xt^2*(-\[Lambda]Lik + yk*\[Lambda]Rik)))*CC[2])/
   (xs*xt^3*(-1 + xs + xt)^2) + 
  (2*mi^2*Qk*\[Lambda]Ljk*
    (-(Qk*(-1 + xt)*xt*(xt^3 + (-1 + xs)*xt*(-1 + 4*yk^2 - 4*yS^2) + 
        2*(-1 + xs)^2*(yk^2 - yS^2) + 2*xt^2*(-1 + xs + yk^2 - yS^2))*
       \[Lambda]Lik) + Qj*(-((-1 + xt)^3*(yk^2 - yS^2)^2*\[Lambda]Lik) + 
       4*xs^3*xt*yk*(yk*\[Lambda]Lik + \[Lambda]Rik) - 
       xs^2*(-1 + xt)*(xt^2*\[Lambda]Lik + (yk^2 - yS^2)^2*\[Lambda]Lik - 
         xt*yk*(6*yk*\[Lambda]Lik + 7*\[Lambda]Rik)) + 
       xs*(-1 + xt)*(2*(yk^2 - yS^2)^2*\[Lambda]Lik + 
         xt^2*(\[Lambda]Lik + 2*yk^2*\[Lambda]Lik + 3*yk*\[Lambda]Rik) - 
         xt*(2*yk^4*\[Lambda]Lik + 2*yS^4*\[Lambda]Lik + yk^2*(2 - 4*yS^2)*
            \[Lambda]Lik + 3*yk*\[Lambda]Rik))))*CC[3])/
   (xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)^2) + 
  (2*mi^2*(Qj - Qk)*(-(Qk*(-1 + xt)*xt*(xt^3 - 2*(-1 + xs)^2*(yk^2 - yS^2) + 
        2*xt^2*(-1 + xs - yk^2 + yS^2) + xt*(1 + 2*xs^2 + 4*yk^2 - 4*yS^2 + 
          xs*(-2 - 4*yk^2 + 4*yS^2)))) + 
     Qj*(xt^5 - (-1 + xs)^2*(yk^2 - yS^2)^2 + 
       xt^4*(-3 + 2*xs - 2*yk^2 + 2*yS^2) + xt^3*(3 + xs^2 + yk^4 - 6*yS^2 + 
         yS^4 - 2*yk^2*(-3 + yS^2) + xs*(-4 - 4*yk^2 + 2*yS^2)) + 
       xt*((3 - 4*xs + xs^2)*yk^4 - (-1 + xs)*yS^2*(-2 + 4*xs^2 + 3*yS^2 - 
           xs*yS^2) - 2*(-1 + xs)*yk^2*(1 - 3*yS^2 + xs*(-1 + yS^2))) - 
       xt^2*(1 + 3*yk^4 - 6*yS^2 + 3*yS^4 - 6*yk^2*(-1 + yS^2) + 
         xs^2*(1 + 2*yk^2 + 4*yS^2) - 2*xs*(yk^4 - 2*yk^2*(-2 + yS^2) + 
           (-1 + yS^2)^2))))*\[Lambda]Lik*\[Lambda]Ljk*CC[4])/
   (xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)^2) + 
  (2*mi^2*Qk^2*\[Lambda]Ljk*(xs^3*\[Lambda]Lik + 
     xs^2*(-1 + xt + 2*yk^2 - 2*yS^2)*\[Lambda]Lik + 
     (xt^3 + xt^2*(-1 + 2*yk^2 - 2*yS^2) - 2*(yk^2 - yS^2)^2 + 
       2*xt*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2)))*\[Lambda]Lik + 
     xs*(xt^2*\[Lambda]Lik + 2*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2))*
        \[Lambda]Lik - 2*xt*(2*yS^2*\[Lambda]Lik + yk*\[Lambda]Rik)))*CC[5])/
   (xs^2*xt^3) - (2*mi^2*(Qj - Qk)^2*(xs^3 + xt^3 - 2*(yk^2 - yS^2)^2 + 
     xt^2*(-1 - 2*yk^2 + 2*yS^2) + xs^2*(-1 + xt - 2*yk^2 + 2*yS^2) + 
     2*xt*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4) + 
     xs*(xt^2 - 4*xt*yk^2 + 2*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4)))*
    \[Lambda]Lik*\[Lambda]Ljk*CC[6])/(xs^2*xt^3) + 
  (2*mi^2*(Qj - Qk)*(-(Qk*(-1 + xs)*xs*(xs^3 - 2*(-1 + xt)^2*(yk^2 - yS^2) + 
        2*xs^2*(-1 + xt - yk^2 + yS^2) + xs*(1 + 2*xt^2 + 4*yk^2 - 4*yS^2 + 
          xt*(-2 - 4*yk^2 + 4*yS^2)))) + 
     Qj*(xs^5 - (-1 + xt)^2*(yk^2 - yS^2)^2 + 
       xs^4*(-3 + 2*xt - 2*yk^2 + 2*yS^2) + xs^3*(3 + xt^2 + yk^4 - 6*yS^2 + 
         yS^4 - 2*yk^2*(-3 + yS^2) + xt*(-4 - 4*yk^2 + 2*yS^2)) + 
       xs*((3 - 4*xt + xt^2)*yk^4 - (-1 + xt)*yS^2*(-2 + 4*xt^2 + 3*yS^2 - 
           xt*yS^2) - 2*(-1 + xt)*yk^2*(1 - 3*yS^2 + xt*(-1 + yS^2))) - 
       xs^2*(1 + 3*yk^4 - 6*yS^2 + 3*yS^4 - 6*yk^2*(-1 + yS^2) + 
         xt^2*(1 + 2*yk^2 + 4*yS^2) - 2*xt*(yk^4 - 2*yk^2*(-2 + yS^2) + 
           (-1 + yS^2)^2))))*\[Lambda]Lik*\[Lambda]Ljk*CC[7])/
   (xs^2*xt^3*(-1 + xs + xt)^2) - 
  (2*mi^2*Qk*\[Lambda]Ljk*(Qk*(-1 + xs)*xs*
      (xs^3 + xt^2*(1 + 2*yk^2 - 2*yS^2) + 2*(yk^2 - yS^2) - 
       4*xt*(yk^2 - yS^2) + 2*xs^2*(-1 + xt + yk^2 - yS^2) + 
       xs*(1 - 4*yk^2 + 4*yS^2 + xt*(-2 + 4*yk^2 - 4*yS^2)))*\[Lambda]Lik + 
     Qj*(-((-1 + xt)^2*(yk^2 - yS^2)^2*\[Lambda]Lik) + 
       xs^3*(xt^2*\[Lambda]Lik + (yk^2 - yS^2)^2*\[Lambda]Lik - 
         xt*yk*(2*yk*\[Lambda]Lik + \[Lambda]Rik)) + 
       xs^2*(-3*(yk^2 - yS^2)^2*\[Lambda]Lik + 2*xt*(yk^4*\[Lambda]Lik + 
           yS^4*\[Lambda]Lik - 2*yk^2*(-1 + yS^2)*\[Lambda]Lik + 
           yk*\[Lambda]Rik) - xt^2*((2 + 6*yk^2)*\[Lambda]Lik + 
           3*yk*\[Lambda]Rik)) + xs*(3*(yk^2 - yS^2)^2*\[Lambda]Lik - 
         2*xt^3*yk*(2*yk*\[Lambda]Lik + \[Lambda]Rik) - 
         xt*(4*yk^4*\[Lambda]Lik + 4*yS^4*\[Lambda]Lik + yk^2*(2 - 8*yS^2)*
            \[Lambda]Lik + yk*\[Lambda]Rik) + 
         xt^2*((1 + yk^4 + yS^4 - 2*yk^2*(-3 + yS^2))*\[Lambda]Lik + 
           3*yk*\[Lambda]Rik))))*CC[8])/(xs^2*xt^3*(-1 + xs + xt)^2) + 
  (2*mi^2*(Qj - Qk)*\[Lambda]Ljk*
    (Qj*(xt^6 + (-1 + xs)^2*(yk^2 - yS^2)^2 + 
       2*xt^5*(-2 + xs - yk^2 + yS^2) + xt^4*(6 + xs^2 + yk^4 - 8*yS^2 + 
         yS^4 - 2*yk^2*(-4 + yS^2) + xs*(-6 - 4*yk^2 + 2*yS^2)) - 
       2*xt^3*(xs^2*(1 + yk^2) + 2*(1 + yk^4 - 3*yS^2 + yS^4 + 
           yk^2*(3 - 2*yS^2)) - xs*(3 + yk^4 - 2*yS^2 + yS^4 - 
           2*yk^2*(-3 + yS^2))) + xt^2*(1 + 6*yk^4 - 8*yS^2 + 6*yS^4 + 
         yk^2*(8 - 12*yS^2) + xs^2*(1 + yk^4 + 4*yS^2 + yS^4 - 
           2*yk^2*(-2 + yS^2)) - 2*xs*(1 + 3*yk^4 - yS^2 + 3*yS^4 - 
           6*yk^2*(-1 + yS^2))) + 2*xt*(-((2 - 3*xs + xs^2)*yk^4) + 
         (-1 + xs)*yS^2*(-1 + xs^2 + 2*yS^2 - xs*(1 + yS^2)) + 
         (-1 + xs)*yk^2*(1 - 4*yS^2 + xs*(-1 + 2*yS^2))))*\[Lambda]Lik - 
     Qk*(-1 + xt)^2*xt*(xt^3*\[Lambda]Lik + 2*xt^2*(-1 + xs - yk^2 + yS^2)*
        \[Lambda]Lik + xt*(1 + 2*xs^2 + 4*yk^2 - 4*yS^2 + 
         xs*(-3 - 2*yk^2 + 2*yS^2))*\[Lambda]Lik + xs*xt*yk*\[Lambda]Rik + 
       (-1 + xs)*(2*yk^2*\[Lambda]Lik - 2*yS^2*\[Lambda]Lik + 
         xs*yk*\[Lambda]Rik)))*CC[9])/(xs^2*(-1 + xt)*xt^3*
    (-1 + xs + xt)^2) - (2*mi^2*Qk*\[Lambda]Ljk*
    (Qk*(-1 + xt)^2*xt*(xt^3*\[Lambda]Lik + 2*xt^2*(-1 + xs + yk^2 - yS^2)*
        \[Lambda]Lik + xt*(1 - 4*yk^2 + 4*yS^2 + 2*xs*(-1 + yk^2 - yS^2))*
        \[Lambda]Lik - xs*xt*yk*\[Lambda]Rik - 
       (-1 + xs)*(2*yk^2*\[Lambda]Lik - 2*yS^2*\[Lambda]Lik + 
         xs*yk*\[Lambda]Rik)) + Qj*((-1 + xt)^4*(yk^2 - yS^2)^2*
        \[Lambda]Lik + 2*xs^3*xt*yk*(yk*\[Lambda]Lik + \[Lambda]Rik - 
         xt*\[Lambda]Rik) - 2*xs*(-1 + xt)^2*((yk^2 - yS^2)^2*\[Lambda]Lik + 
         xt^2*(yS^2*\[Lambda]Lik + yk*\[Lambda]Rik) - 
         xt*(yk^4*\[Lambda]Lik + yS^2*(1 + yS^2)*\[Lambda]Lik + 
           yk^2*(\[Lambda]Lik - 2*yS^2*\[Lambda]Lik) + yk*\[Lambda]Rik)) + 
       xs^2*(-1 + xt)*(xt^3*\[Lambda]Lik - (yk^2 - yS^2)^2*\[Lambda]Lik - 
         xt^2*(\[Lambda]Lik + 2*yS^2*\[Lambda]Lik + 4*yk*\[Lambda]Rik) + 
         xt*(yk^4*\[Lambda]Lik - 2*yk^2*(-2 + yS^2)*\[Lambda]Lik + 
           yS^2*(2 + yS^2)*\[Lambda]Lik + 4*yk*\[Lambda]Rik))))*CC[10])/
   (xs^2*(-1 + xt)*xt^3*(-1 + xs + xt)^2) + 
  (2*mi^2*Qk^2*(xs^2 + xt^2)*(xs + xt + 2*yk^2 - 2*yS^2)*\[Lambda]Lik*
    \[Lambda]Ljk*CC[11])/(xs^2*xt^3) - 
  (2*mi^2*(Qj - Qk)^2*(xs^2 + xt^2)*(xs + xt - 2*yk^2 + 2*yS^2)*\[Lambda]Lik*
    \[Lambda]Ljk*CC[12])/(xs^2*xt^3) + 
  (2*mi^4*(Qj - Qk)*Qk*(xs*(xt - (yk - yS)^2) - (-1 + xt)*(yk - yS)^2)*
    ((-1 + xt)*(yk^2 - yS^2) + xs*(xt + yk^2 - yS^2))*
    (-((-1 + xt)*(yk + yS)^2) + xs*(xt - (yk + yS)^2))*\[Lambda]Lik*
    \[Lambda]Ljk*DD[1])/(xs^2*xt^3*(-1 + xs + xt)^2) + 
  (2*mi^4*(Qj - Qk)*Qk*\[Lambda]Ljk*((-1 + xt)^3*(yk^2 - yS^2)^3*
      \[Lambda]Lik + xs^3*(xt^3*\[Lambda]Lik + (yk^2 - yS^2)^3*\[Lambda]Lik + 
       xt^2*(yk^2*\[Lambda]Lik - yS^2*\[Lambda]Lik + yk*\[Lambda]Rik) - 
       xt*(yk^2 - yS^2)*(3*yk^2*\[Lambda]Lik + yS^2*\[Lambda]Lik + 
         yk*\[Lambda]Rik)) - xs*(-1 + xt)*(3*(yk^2 - yS^2)^3*\[Lambda]Lik - 
       xt*(yk^2 - yS^2)*(3*yk^4*\[Lambda]Lik + yk^2*(3 - 6*yS^2)*
          \[Lambda]Lik + yS^2*(1 + 3*yS^2)*\[Lambda]Lik + yk*\[Lambda]Rik) + 
       xt^2*(3*yk^4*\[Lambda]Lik - yS^2*(1 + yS^2)*\[Lambda]Lik - 
         yk^2*(\[Lambda]Lik + 2*yS^2*\[Lambda]Lik) + yk^3*\[Lambda]Rik - 
         yk*yS^2*\[Lambda]Rik)) + xs^2*(-3*(yk^2 - yS^2)^3*\[Lambda]Lik + 
       xt^3*((-1 + yk^2 - yS^2)*\[Lambda]Lik + yk*\[Lambda]Rik) + 
       xt*(yk^2 - yS^2)*(3*yk^4*\[Lambda]Lik - 6*yk^2*(-1 + yS^2)*
          \[Lambda]Lik + yS^2*(2 + 3*yS^2)*\[Lambda]Lik + 
         2*yk*\[Lambda]Rik) + xt^2*(-6*yk^4*\[Lambda]Lik + 
         4*yk^2*yS^2*\[Lambda]Lik + 2*yS^2*(1 + yS^2)*\[Lambda]Lik - 
         2*yk^3*\[Lambda]Rik + yk*(-1 + 2*yS^2)*\[Lambda]Rik)))*DD[2])/
   (xs^2*xt^3*(-1 + xs + xt)^2) + (2*mi^4*(Qj - Qk)^2*(xs - yk^2 + yS^2)*
    (xs^3 + (-1 + xt)*(yk^2 - yS^2)^2 + xs^2*(-1 + xt - 2*yk^2 + 2*yS^2) + 
     xs*(yk^4 + yS^2*(-2 - 2*xt + yS^2) - 2*yk^2*(-1 + xt + yS^2)))*
    \[Lambda]Lik*\[Lambda]Ljk*DD[3])/(xs^2*xt^3) - 
  (2*mi^4*Qk^2*\[Lambda]Ljk*(xs^5*\[Lambda]Lik + 
     xs^4*(-2 + 2*xt + 3*yk^2 - 3*yS^2)*\[Lambda]Lik + 
     (-1 + xt)^2*(yk^2 - yS^2)^3*\[Lambda]Lik + 
     xs^3*((1 + xt^2 + 3*yk^4 + 6*yS^2 + 3*yS^4 + 2*xt*(-1 + yk^2 - 3*yS^2) - 
         6*yk^2*(1 + yS^2))*\[Lambda]Lik - xt*yk*\[Lambda]Rik) + 
     xs^2*(yk^6*\[Lambda]Lik - yS^2*(3 + 3*xt^2 + 6*yS^2 + yS^4 - 
         6*xt*(1 + yS^2))*\[Lambda]Lik + yk^4*(2*xt - 3*(2 + yS^2))*
        \[Lambda]Lik + yk^2*(-5*xt^2 - 2*xt*(1 + 4*yS^2) + 
         3*(1 + 4*yS^2 + yS^4))*\[Lambda]Lik - xt*yk^3*\[Lambda]Rik + 
       xt*yk*(1 - xt + yS^2)*\[Lambda]Rik) + 
     xs*(-((-3 + 2*yk^2 - 2*yS^2)*(yk^2 - yS^2)^2*\[Lambda]Lik) + 
       xt*(yk^2 - yS^2)*(2*yk^4*\[Lambda]Lik + 2*yS^2*(3 + yS^2)*
          \[Lambda]Lik - 2*yk^2*(\[Lambda]Lik + 2*yS^2*\[Lambda]Lik) + 
         yk*\[Lambda]Rik) + xt^2*(-(yk^4*\[Lambda]Lik) + 
         3*yS^4*\[Lambda]Lik - 2*yk^2*(-1 + yS^2)*\[Lambda]Lik - 
         yk^3*\[Lambda]Rik + yk*yS^2*\[Lambda]Rik)))*DD[4])/
   (xs^2*xt^3*(-1 + xs + xt)) + (2*mi^4*(Qj - Qk)^2*(xt - yk^2 + yS^2)*
    (xt^3 + (-1 + xs)*(yk^2 - yS^2)^2 + xt^2*(-1 + xs - 2*yk^2 + 2*yS^2) + 
     xt*(yk^4 + yS^2*(-2 - 2*xs + yS^2) - 2*yk^2*(-1 + xs + yS^2)))*
    \[Lambda]Lik*\[Lambda]Ljk*DD[5])/(xs^2*xt^3) - 
  (2*mi^4*Qk^2*\[Lambda]Ljk*(xt^5*\[Lambda]Lik + 
     xt^4*(-2 + 2*xs + 3*yk^2 - 3*yS^2)*\[Lambda]Lik + 
     (-1 + xs)^2*(yk^2 - yS^2)^3*\[Lambda]Lik + 
     xt^3*((1 + xs^2 + 3*yk^4 + 6*yS^2 + 3*yS^4 + 2*xs*(-1 + yk^2 - 3*yS^2) - 
         6*yk^2*(1 + yS^2))*\[Lambda]Lik - 3*xs*yk*\[Lambda]Rik) - 
     (-1 + xs)*xt*(yk^2 - yS^2)*(-2*yk^4*\[Lambda]Lik + 
       yS^2*(-3 + 3*xs - 2*yS^2)*\[Lambda]Lik + yk^2*(3 + xs + 4*yS^2)*
        \[Lambda]Lik + xs*yk*\[Lambda]Rik) + 
     xt^2*(yk^6*\[Lambda]Lik - yS^2*(3 + 3*xs^2 + 6*yS^2 + yS^4 - 
         6*xs*(1 + yS^2))*\[Lambda]Lik + yk^4*(2*xs - 3*(2 + yS^2))*
        \[Lambda]Lik + yk^2*(-5*xs^2 - 8*xs*yS^2 + 3*(1 + 4*yS^2 + yS^4))*
        \[Lambda]Lik - xs*yk^3*\[Lambda]Rik + xs*yk*(3 - 3*xs + yS^2)*
        \[Lambda]Rik))*DD[6])/(xs^2*xt^3*(-1 + xs + xt)), 
 (-2*Qj*(-2*Qk*xs*xt*(-2 + xs + xt) + Qj*(-((-1 + xt)*xt*(yk^2 - yS^2)) + 
       xs^2*(xt - yk^2 + 2*xt*yk^2 + yS^2 - 2*xt*yS^2) + 
       xs*(yk^2 - yS^2 + xt^2*(1 + 2*yk^2 - 2*yS^2) + 
         xt*(-2 - 4*yk^2 + 4*yS^2))))*\[Lambda]Rik*\[Lambda]Rjk)/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)) + 
  (2*Qj^2*(xs + xt - 2*xs*xt)*yk^2*\[Lambda]Rik*\[Lambda]Rjk*BB[1])/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2) + 
  (2*Qj^2*(-xt + xs*(-1 + 2*xt))*yS^2*\[Lambda]Rik*\[Lambda]Rjk*BB[2])/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2) + 
  (2*Qj*(-2*Qk*(2 - 2*xs + xs^2 - 2*xt + xt^2)*\[Lambda]Rik + 
     Qj*(4*(-1 + xs)*(-1 + xt)*(-2 + xs + xt)*yk*\[Lambda]Lik + 
       (-6 + 2*xs*(-2 + xt)^2 + 8*xt - 3*xt^2 + xs^2*(-3 + 2*xt))*yk^2*
        \[Lambda]Rik - (2 - 4*xt + xt^2 - 6*yS^2 + 8*xt*yS^2 - 3*xt^2*yS^2 + 
         xs^2*(1 - 3*yS^2 + 2*xt*(-1 + yS^2)) + 
         2*xs*(-2 + 4*yS^2 - 4*xt*(-1 + yS^2) + xt^2*(-1 + yS^2)))*
        \[Lambda]Rik))*\[Lambda]Rjk*BB[3])/((-1 + xs)^2*xs*(-1 + xt)^2*xt*
    (-1 + xs + xt)) + 
  (2*Qj*(2*Qk*xs*\[Lambda]Rik + Qj*((yk^2 - yS^2)*\[Lambda]Rik - 
       2*xs^2*(2*yk*\[Lambda]Lik + \[Lambda]Rik) + 
       xs*(4*yk*\[Lambda]Lik + \[Lambda]Rik)))*\[Lambda]Rjk*BB[4])/
   ((-1 + xs)^2*xs^2*xt*(-1 + xs + xt)) + 
  (2*Qj*(2*Qk*xt*\[Lambda]Rik + Qj*((yk^2 - yS^2)*\[Lambda]Rik - 
       2*xt^2*(2*yk*\[Lambda]Lik + \[Lambda]Rik) + 
       xt*(4*yk*\[Lambda]Lik + \[Lambda]Rik)))*\[Lambda]Rjk*BB[5])/
   (xs*(-1 + xt)^2*xt^2*(-1 + xs + xt)) + 
  (2*mi^2*Qk*(Qk*xs*(2*xt^2*yk*\[Lambda]Lik + (-1 + xs)*(yk^2 - yS^2)*
        \[Lambda]Rik + xt*(2*(-1 + xs)*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik + 
         (-1 + 2*xs - yS^2)*\[Lambda]Rik)) + 
     Qj*(-((-1 + xt)*xt*yk^2*\[Lambda]Rik) + 
       xs^2*(yS^2*\[Lambda]Rik - 2*xt*(yk*\[Lambda]Lik + \[Lambda]Rik)) + 
       xs*(-2*xt^2*yk*\[Lambda]Lik - yS^2*\[Lambda]Rik + 
         xt*(2*yk*\[Lambda]Lik + \[Lambda]Rik - yk^2*\[Lambda]Rik + 
           yS^2*\[Lambda]Rik))))*\[Lambda]Rjk*CC[1])/
   (xs*xt^2*(-1 + xs + xt)^3) - 
  (2*mi^2*(Qj - Qk)*(Qj*(-xs + xs^2 + xt - xt^2)*yS^2*\[Lambda]Rik + 
     Qk*xs*((-1 + xs)*(yk^2 - yS^2)*\[Lambda]Rik + 
       2*xt^2*(yk*\[Lambda]Lik + \[Lambda]Rik) + 
       xt*(2*(-1 + xs)*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - 
         (1 + yS^2)*\[Lambda]Rik)))*\[Lambda]Rjk*CC[2])/
   (xs*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^2*Qk*(Qk*xt*(2*xs^2*yk*\[Lambda]Lik + (-1 + xt)*(yk^2 - yS^2)*
        \[Lambda]Rik + xs*(2*(-1 + xt)*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik + 
         (-1 + 2*xt - yS^2)*\[Lambda]Rik)) + 
     Qj*((-1 + xt)*xt*yS^2*\[Lambda]Rik - xs^2*yk*(2*xt*\[Lambda]Lik + 
         yk*\[Lambda]Rik) + xs*(yk^2*\[Lambda]Rik - 
         2*xt^2*(yk*\[Lambda]Lik + \[Lambda]Rik) + 
         xt*(2*yk*\[Lambda]Lik + \[Lambda]Rik - yk^2*\[Lambda]Rik + 
           yS^2*\[Lambda]Rik))))*\[Lambda]Rjk*CC[3])/
   (xs^2*xt*(-1 + xs + xt)^3) + 
  (2*mi^2*(Qj - Qk)*(Qj*(-xs + xs^2 + xt - xt^2)*yS^2*\[Lambda]Rik + 
     Qk*xt*(-((-1 + xt)*(yk^2 - yS^2)*\[Lambda]Rik) - 
       2*xs^2*(yk*\[Lambda]Lik + \[Lambda]Rik) + 
       xs*(2*yk*\[Lambda]Lik - 2*xt*yk*\[Lambda]Lik + \[Lambda]Rik - 
         yk^2*\[Lambda]Rik + yS^2*\[Lambda]Rik)))*\[Lambda]Rjk*CC[4])/
   (xs^2*xt*(-1 + xs + xt)^3) - (4*mi^2*Qk^2*(xs + xt)*yk^2*\[Lambda]Rik*
    \[Lambda]Rjk*CC[5])/(xs^2*xt^2*(-1 + xs + xt)) + 
  (4*mi^2*(Qj - Qk)^2*(xs + xt)*yS^2*\[Lambda]Rik*\[Lambda]Rjk*CC[6])/
   (xs^2*xt^2*(-1 + xs + xt)) - 
  (2*mi^2*(Qj - Qk)*(Qj*(xs^4 + (-1 + xt)*xt + xs^3*(-3 + 6*xt) + 
       xs^2*(3 - 13*xt + 7*xt^2) + xs*(-1 + 8*xt - 8*xt^2 + 2*xt^3))*yS^2*
      \[Lambda]Rik + Qk*(-1 + xs)^2*xt*((-1 + xt)*(yk^2 - yS^2)*
        \[Lambda]Rik + 2*xs^2*(yk*\[Lambda]Lik + \[Lambda]Rik) + 
       xs*(2*(-1 + xt)*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - 
         (1 + yS^2)*\[Lambda]Rik)))*\[Lambda]Rjk*CC[7])/
   ((-1 + xs)*xs^2*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^2*Qk*(Qk*(-1 + xs)^2*xt*(2*xs^2*yk*\[Lambda]Lik + 
       (-1 + xt)*(yk^2 - yS^2)*\[Lambda]Rik + 
       xs*(2*(-1 + xt)*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik + 
         (-1 + 2*xt - yS^2)*\[Lambda]Rik)) + 
     Qj*((-1 + xt)*xt*yS^2*\[Lambda]Rik + xs^4*yk*(-2*xt*\[Lambda]Lik + 
         yk*\[Lambda]Rik) + xs^3*(-3*yk^2*\[Lambda]Rik - 
         2*xt^2*(yk*\[Lambda]Lik + \[Lambda]Rik) + 
         xt*(6*yk*\[Lambda]Lik + \[Lambda]Rik + 5*yk^2*\[Lambda]Rik + 
           yS^2*\[Lambda]Rik)) + xs*(-(yk^2*\[Lambda]Rik) + 
         2*xt^3*yk^2*\[Lambda]Rik - 2*xt^2*(yk*\[Lambda]Lik + \[Lambda]Rik + 
           3*yk^2*\[Lambda]Rik + yS^2*\[Lambda]Rik) + 
         xt*(2*yk*\[Lambda]Lik + \[Lambda]Rik + 5*yk^2*\[Lambda]Rik + 
           3*yS^2*\[Lambda]Rik)) + xs^2*(3*yk^2*\[Lambda]Rik + 
         xt^2*(4*yk*\[Lambda]Lik + 6*yk^2*\[Lambda]Rik + 
           (4 + yS^2)*\[Lambda]Rik) - xt*(6*yk*\[Lambda]Lik + 
           10*yk^2*\[Lambda]Rik + (2 + 3*yS^2)*\[Lambda]Rik))))*\[Lambda]Rjk*
    CC[8])/((-1 + xs)*xs^2*xt^2*(-1 + xs + xt)^3) - 
  (2*mi^2*(Qj - Qk)*(Qj*(2*xs^3*xt + (-1 + xt)^3*xt + 
       xs*(-1 + xt)^2*(-1 + 6*xt) + xs^2*(1 - 8*xt + 7*xt^2))*yS^2*
      \[Lambda]Rik + Qk*xs*(-1 + xt)^2*((-1 + xs)*(yk^2 - yS^2)*
        \[Lambda]Rik + 2*xt^2*(yk*\[Lambda]Lik + \[Lambda]Rik) + 
       xt*(2*(-1 + xs)*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - 
         (1 + yS^2)*\[Lambda]Rik)))*\[Lambda]Rjk*CC[9])/
   (xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^2*Qk*(Qk*xs*(-1 + xt)^2*(2*xt^2*yk*\[Lambda]Lik + 
       (-1 + xs)*(yk^2 - yS^2)*\[Lambda]Rik + 
       xt*(2*(-1 + xs)*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik + 
         (-1 + 2*xs - yS^2)*\[Lambda]Rik)) + 
     Qj*(2*xs^3*xt*yk^2*\[Lambda]Rik + (-1 + xt)^3*xt*yk^2*\[Lambda]Rik - 
       xs*(-1 + xt)^2*(2*xt^2*yk*\[Lambda]Lik + yS^2*\[Lambda]Rik - 
         xt*(2*yk*\[Lambda]Lik + \[Lambda]Rik + 5*yk^2*\[Lambda]Rik + 
           yS^2*\[Lambda]Rik)) - xs^2*(-1 + xt)*(yS^2*\[Lambda]Rik + 
         2*xt^2*(yk*\[Lambda]Lik + \[Lambda]Rik) - 
         xt*(2*yk*\[Lambda]Lik + 6*yk^2*\[Lambda]Rik + 
           (2 + yS^2)*\[Lambda]Rik))))*\[Lambda]Rjk*CC[10])/
   (xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^4*(Qj - Qk)*Qk*(-((-1 + xt)^2*xt*yk^2*(yk^2 - yS^2)*\[Lambda]Rik) + 
     xs^3*(yS^2*(-yk^2 + yS^2)*\[Lambda]Rik + 
       2*xt^2*(yk*\[Lambda]Lik + \[Lambda]Rik) - 
       xt*(2*yk^3*\[Lambda]Lik + 2*yk*yS^2*\[Lambda]Lik + 
         2*yk^2*\[Lambda]Rik + 3*yS^2*\[Lambda]Rik)) - 
     xs*(-1 + xt)*(yS^2*(-yk^2 + yS^2)*\[Lambda]Rik + 
       xt^2*yk*(2*yk^2*\[Lambda]Lik + 2*yS^2*\[Lambda]Lik - 
         yk*\[Lambda]Rik) + xt*(-2*yk^3*\[Lambda]Lik - 
         2*yk*yS^2*\[Lambda]Lik + 2*yk^4*\[Lambda]Rik - 
         yk^2*(1 + yS^2)*\[Lambda]Rik - yS^2*(1 + yS^2)*\[Lambda]Rik)) + 
     xs^2*(2*xt^3*yk*\[Lambda]Lik + 2*yS^2*(yk^2 - yS^2)*\[Lambda]Rik - 
       xt^2*(4*yk^3*\[Lambda]Lik + 2*yk*(\[Lambda]Lik + 
           2*yS^2*\[Lambda]Lik) + \[Lambda]Rik + yk^2*\[Lambda]Rik + 
         3*yS^2*\[Lambda]Rik) + xt*(4*yk^3*\[Lambda]Lik + 
         4*yk*yS^2*\[Lambda]Lik - yk^4*\[Lambda]Rik - yk^2*(-3 + yS^2)*
          \[Lambda]Rik + 2*yS^2*(2 + yS^2)*\[Lambda]Rik)))*\[Lambda]Rjk*
    DD[1])/(xs^2*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^4*(Qj - Qk)*Qk*((-1 + xt)^2*xt*yS^2*(-yk^2 + yS^2)*\[Lambda]Rik + 
     xs^3*yk*(2*xt^2*\[Lambda]Lik + yk*(-yk^2 + yS^2)*\[Lambda]Rik + 
       xt*(-2*yk^2*\[Lambda]Lik - 2*yS^2*\[Lambda]Lik + yk*\[Lambda]Rik)) + 
     xs^2*(2*yk^2*(yk^2 - yS^2)*\[Lambda]Rik + 
       2*xt^3*(yk*\[Lambda]Lik + \[Lambda]Rik) - 
       xt^2*(2*yk*\[Lambda]Lik + 4*yk^3*\[Lambda]Lik + 
         4*yk*yS^2*\[Lambda]Lik + \[Lambda]Rik + yk^2*\[Lambda]Rik + 
         3*yS^2*\[Lambda]Rik) + xt*(4*yk^3*\[Lambda]Lik + 
         4*yk*yS^2*\[Lambda]Lik - 2*yk^4*\[Lambda]Rik + yS^2*\[Lambda]Rik + 
         yk^2*yS^2*\[Lambda]Rik + yS^4*\[Lambda]Rik)) - 
     xs*(-1 + xt)*(yk^2*(-yk^2 + yS^2)*\[Lambda]Rik + 
       xt^2*(2*yk^3*\[Lambda]Lik + 2*yk*yS^2*\[Lambda]Lik + 
         2*yk^2*\[Lambda]Rik + 3*yS^2*\[Lambda]Rik) + 
       xt*(-2*yk^3*\[Lambda]Lik - 2*yk*yS^2*\[Lambda]Lik + 
         yk^4*\[Lambda]Rik + yk^2*(-1 + yS^2)*\[Lambda]Rik - 
         yS^2*(1 + 2*yS^2)*\[Lambda]Rik)))*\[Lambda]Rjk*DD[2])/
   (xs^2*xt^2*(-1 + xs + xt)^3) - 
  (2*mi^4*(Qj - Qk)^2*yS^2*(xs^2*\[Lambda]Rik + xt*(-yk^2 + yS^2)*
      \[Lambda]Rik - xs*((yk^2 - yS^2)*\[Lambda]Rik + 
       xt*(4*yk*\[Lambda]Lik + \[Lambda]Rik)))*\[Lambda]Rjk*DD[3])/
   (xs^2*xt^2*(-1 + xs + xt)) + 
  (2*mi^4*Qk^2*yk^2*(xs^2*\[Lambda]Rik + xt*(yk^2 - yS^2)*\[Lambda]Rik + 
     xs*(4*xt*yk*\[Lambda]Lik + 3*xt*\[Lambda]Rik + yk^2*\[Lambda]Rik - 
       yS^2*\[Lambda]Rik))*\[Lambda]Rjk*DD[4])/(xs^2*xt^2*(-1 + xs + xt)) + 
  (2*mi^4*(Qj - Qk)^2*yS^2*(xs*(yk^2 - yS^2)*\[Lambda]Rik - 
     xt*(xt - yk^2 + yS^2)*\[Lambda]Rik + 
     xs*xt*(4*yk*\[Lambda]Lik + \[Lambda]Rik))*\[Lambda]Rjk*DD[5])/
   (xs^2*xt^2*(-1 + xs + xt)) + 
  (2*mi^4*Qk^2*yk^2*(xt*(xt + yk^2 - yS^2)*\[Lambda]Rik + 
     xs*(4*xt*yk*\[Lambda]Lik + 3*xt*\[Lambda]Rik + yk^2*\[Lambda]Rik - 
       yS^2*\[Lambda]Rik))*\[Lambda]Rjk*DD[6])/(xs^2*xt^2*(-1 + xs + xt)), 
 (Qj*\[Lambda]Ljk*(2*Qk*(1 - 2*xt)*\[Lambda]Lik + 
     Qj*(-1 + yk^2 - yS^2 + xt*(2 + yk^2 - yS^2))*\[Lambda]Lik + 
     2*Qj*xt*yk*\[Lambda]Rik))/(xs*(-1 + xt)*xt^2) + 
  (Qj^2*yk^2*\[Lambda]Ljk*((1 + xt)*yk^2*\[Lambda]Lik - 
     (1 + xt)*yS^2*\[Lambda]Lik + 2*xt*yk*\[Lambda]Rik)*BB[1])/
   (xs*(-1 + xt)*xt^2*(yk^2 - yS^2)) + 
  (Qj^2*yS^2*\[Lambda]Ljk*(-((1 + xt)*yk^2*\[Lambda]Lik) + 
     (1 + xt)*yS^2*\[Lambda]Lik - 2*xt*yk*\[Lambda]Rik)*BB[2])/
   (xs*(-1 + xt)*xt^2*(yk^2 - yS^2)) + 
  (Qj*\[Lambda]Ljk*(2*Qk*(xs - xt)*\[Lambda]Lik + 
     Qj*(xt*(1 - yk^2 + yS^2)*\[Lambda]Lik + xs*(-1 + xt)*
        ((1 + yk^2 - yS^2)*\[Lambda]Lik + 2*yk*\[Lambda]Rik) - 
       xs^2*(xt*(1 + yk^2 - yS^2)*\[Lambda]Lik + 2*xt*yk*\[Lambda]Rik - 
         2*(yk^2*\[Lambda]Lik - yS^2*\[Lambda]Lik + yk*\[Lambda]Rik))))*
    BB[3])/((-1 + xs)*xs^2*(-1 + xt)*xt^2) + 
  (Qj*(2*Qk*xs - Qj*(xs - yk^2 + yS^2))*\[Lambda]Lik*\[Lambda]Ljk*BB[4])/
   ((-1 + xs)*xs*xt^2) + (Qj*\[Lambda]Ljk*(-2*Qk*xt^2*\[Lambda]Lik + 
     Qj*xt*(xt - yk^2 + yS^2)*\[Lambda]Lik + 
     Qj*xs*(xt*\[Lambda]Lik - 3*yk^2*\[Lambda]Lik + 3*yS^2*\[Lambda]Lik - 
       2*yk*\[Lambda]Rik))*BB[5])/(xs^2*(-1 + xt)*xt^2) + 
  (Qk^2*(-xs + xt)*\[Lambda]Lik*\[Lambda]Ljk*BB[6])/(xs^2*xt^2) + 
  ((Qj - Qk)^2*(xs - xt)*\[Lambda]Lik*\[Lambda]Ljk*BB[7])/(xs^2*xt^2) - 
  (mi^2*Qk*(Qk*(-1 + xs)*xs^2*(xs^2 + xt - 2*yk^2 + 2*xt*yk^2 + 2*yS^2 - 
       2*xt*yS^2 + xs*(-1 + xt + 2*yk^2 - 2*yS^2)) - 
     Qj*((-1 + xt)^2*xt*(yk^2 - yS^2)^2 + xs^3*(xt^2 + 2*xt*yk^2 - 
         (yk^2 - yS^2)^2) + xs^2*(xt^3 + 2*(yk^2 - yS^2)^2 - 
         xt*(yk^4 - 2*yk^2*(-2 + yS^2) + yS^2*(-2 + yS^2)) + 
         xt^2*(6*yk^2 - 2*(1 + yS^2))) + xs*(xt^3*(4*yk^2 - 2*yS^2) + 
         2*xt*(yk^2 - yS^2) - (yk^2 - yS^2)^2 + 
         xt^2*(yk^4 - 2*yk^2*(3 + yS^2) + yS^2*(4 + yS^2)))))*\[Lambda]Lik*
    \[Lambda]Ljk*CC[1])/(2*xs^2*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*(Qk^2*(-1 + xs)*xs^2*(xs^2 + 2*xt^2 + 2*(yk^2 - yS^2) + 
       xt*(-1 - 2*yk^2 + 2*yS^2) + xs*(-1 + 3*xt - 2*yk^2 + 2*yS^2)) - 
     Qj*Qk*(2*xs^5 - (-1 + xt)^2*xt*(yk^2 - yS^2)^2 + 
       xs^4*(-4 + 6*xt - 4*yk^2 + 4*yS^2) + xs^3*(2 + 5*xt^2 + yk^4 - 
         8*yS^2 + yS^4 - 2*yk^2*(-4 + yS^2) + xt*(-8 - 6*yk^2 + 4*yS^2)) + 
       xs*(-((-1 + xt^2)*yk^4) + 2*(-1 + xt^2)*yk^2*yS^2 - 
         (-1 + xt)*yS^2*(2*xt^2 + yS^2 + xt*yS^2)) + 
       xs^2*(xt^3 - 2*xt^2*(2 + yk^2 + yS^2) + xt*(2 + yk^4 - 4*yS^2 + yS^4 - 
           2*yk^2*(-3 + yS^2)) - 2*(yk^4 + yS^2*(-2 + yS^2) - 
           2*yk^2*(-1 + yS^2)))) + 
     Qj^2*(xs^5 - (-1 + xt)^2*xt*(yk^2 - yS^2)^2 + 
       xs^4*(-2 + 3*xt - 2*yk^2 + 2*yS^2) + xs^3*(1 + 3*xt^2 + yk^4 - 
         4*yS^2 + yS^4 - 2*yk^2*(-2 + yS^2) + xt*(-4 - 4*yk^2 + 2*yS^2)) + 
       xs*(-((-1 + xt^2)*yk^4) + 2*(-1 + xt^2)*yk^2*yS^2 - 
         (-1 + xt)*yS^2*(2*xt^2 + yS^2 + xt*yS^2)) + 
       xs^2*(xt^3 - 2*xt^2*(1 + yk^2 + yS^2) - 2*(yk^2 + yk^4 - yS^2 - 
           2*yk^2*yS^2 + yS^4) + xt*(yk^4 - 2*yk^2*(-2 + yS^2) + 
           (-1 + yS^2)^2))))*\[Lambda]Lik*\[Lambda]Ljk*CC[2])/
   (2*xs^2*xt^3*(-1 + xs + xt)^2) - 
  (mi^2*Qk*\[Lambda]Ljk*
    (-(Qk*(-1 + xt)*xt*((-1 + xt)^2*xt*(xt + 2*yk^2 - 2*yS^2) + 
        2*xs^3*(yk^2 - yS^2) + xs*(-1 + xt)*(3*xt^2 - 2*yk^2 + 2*yS^2 + 
          xt*(-1 + 4*yk^2 - 4*yS^2)) + 2*xs^2*(xt^2 - 2*yk^2 + 2*yS^2 + 
          xt*(-1 + 2*yk^2 - 2*yS^2)))*\[Lambda]Lik) + 
     Qj*(-((-1 + xt)^3*xt*(yk^2 - yS^2)^2*\[Lambda]Lik) - 
       xs*(-1 + xt)^2*((1 + xt)*yk^4 - 2*xt^2*yS^2 - 2*(1 + xt)*yk^2*yS^2 + 
         yS^4 + xt*yS^4)*\[Lambda]Lik + 4*xs^4*xt*yk*(yk*\[Lambda]Lik + 
         \[Lambda]Rik) - xs^3*(-1 + xt)*(xt^2*\[Lambda]Lik - 
         (yk^2 - yS^2)^2*\[Lambda]Lik - 2*xt*yk*(5*yk*\[Lambda]Lik + 
           4*\[Lambda]Rik)) - xs^2*(-1 + xt)*(xt^3*\[Lambda]Lik + 
         2*(yk^2 - yS^2)^2*\[Lambda]Lik - 2*xt^2*(3*yk^2*\[Lambda]Lik + 
           yS^2*\[Lambda]Lik + 2*yk*\[Lambda]Rik) + 
         xt*(-(yk^4*\[Lambda]Lik) - yS^4*\[Lambda]Lik + 2*yk^2*(3 + yS^2)*
            \[Lambda]Lik + 4*yk*\[Lambda]Rik))))*CC[3])/
   (2*xs^3*(-1 + xt)*xt^2*(-1 + xs + xt)^2) + 
  (mi^2*(-(Qk^2*(-1 + xt)*xt*(2*xs^3*(xt - yk^2 + yS^2) + 
        (-1 + xt)^2*xt*(xt - 2*yk^2 + 2*yS^2) + 
        4*xs^2*(xt^2 + yk^2 - yS^2 + xt*(-1 - yk^2 + yS^2)) + 
        xs*(-1 + xt)*(3*xt^2 + 2*(yk^2 - yS^2) + 
          xt*(-1 - 4*yk^2 + 4*yS^2)))) + 
     Qj^2*(4*xs^4*xt*yS^2 - (-1 + xt)^3*xt*(xt - yk^2 + yS^2)^2 - 
       xs*(-1 + xt)^2*(3*xt^3 + (yk^2 - yS^2)^2 + xt*(yk^2 - yS^2)^2 + 
         xt^2*(-1 - 4*yk^2 + 2*yS^2)) - xs^3*(xt^3 + (yk^2 - yS^2)^2 - 
         xt^2*(1 + 10*yS^2) - xt*(yk^4 - 10*yS^2 - 2*yk^2*yS^2 + yS^4)) + 
       xs^2*(-3*xt^4 + 2*(yk^2 - yS^2)^2 + xt^3*(5 + 2*yk^2 + 6*yS^2) - 
         3*xt*(yk^4 - 2*yS^2 - 2*yk^2*yS^2 + yS^4) + 
         xt^2*(-2 + yk^4 - 12*yS^2 + yS^4 - 2*yk^2*(1 + yS^2)))) - 
     Qj*Qk*(4*xs^4*xt*yS^2 - (-1 + xt)^3*xt*(2*xt^2 - 4*xt*(yk^2 - yS^2) + 
         (yk^2 - yS^2)^2) - xs*(-1 + xt)^2*(6*xt^3 + (yk^2 - yS^2)^2 + 
         xt^2*(-2 - 8*yk^2 + 6*yS^2) + xt*(yk^4 + yS^2*(-2 + yS^2) - 
           2*yk^2*(-1 + yS^2))) - xs^3*(3*xt^3 + (yk^2 - yS^2)^2 - 
         xt^2*(3 + 2*yk^2 + 8*yS^2) - xt*(yk^4 + yS^2*(-8 + yS^2) - 
           2*yk^2*(1 + yS^2))) + xs^2*(-7*xt^4 + 2*(yk^2 - yS^2)^2 + 
         xt^3*(13 + 6*yk^2 + 2*yS^2) + xt^2*(-6 + yk^4 - 4*yS^2 + yS^4 - 
           2*yk^2*(5 + yS^2)) + xt*(-3*yk^4 + 2*yS^2 - 3*yS^4 + 
           yk^2*(4 + 6*yS^2)))))*\[Lambda]Lik*\[Lambda]Ljk*CC[4])/
   (2*xs^3*(-1 + xt)*xt^2*(-1 + xs + xt)^2) + 
  (mi^2*Qk^2*(xs^4 + xs^3*(-1 + 2*yk^2 - 2*yS^2) + 
     xs*(-2*xt^3 + 4*xt*(yk^2 - yS^2) - 2*(yk^2 - yS^2)^2 + 
       xt^2*(1 - 2*yk^2 + 6*yS^2)) + 
     xs^2*(-2*xt^2 + xt*(1 - 6*yk^2 + 2*yS^2) + 
       2*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2))) + 
     xt*(-xt^3 + 2*(yk^2 - yS^2)^2 + xt^2*(1 - 2*yk^2 + 2*yS^2) - 
       2*xt*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2))))*\[Lambda]Lik*
    \[Lambda]Ljk*CC[5])/(2*xs^3*xt^3) - 
  (mi^2*(Qj - Qk)^2*(xs^4 + xs^3*(-1 + 2*xt - 2*yk^2 + 2*yS^2) + 
     xs*(-2*xt^3 - 2*(yk^2 - yS^2)^2 + xt^2*(1 + 2*yk^2 + 2*yS^2)) - 
     xs^2*(xt*(1 + 2*yk^2 + 2*yS^2) - 2*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + 
         yS^4)) + xt*(-xt^3 + xt^2*(1 + 2*yk^2 - 2*yS^2) + 
       2*(yk^2 - yS^2)^2 - 2*xt*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4)))*
    \[Lambda]Lik*\[Lambda]Ljk*CC[6])/(2*xs^3*xt^3) + 
  (mi^2*(Qk^2*(-1 + xs)*xs*(xs^4 - 2*(-1 + xt)*xt*(yk^2 - yS^2) + 
       xs*(xt - 2*yk^2 + 2*yS^2) + xs^3*(-2 + 3*xt - 2*yk^2 + 2*yS^2) + 
       xs^2*(1 + 2*xt^2 + 4*yk^2 - 4*yS^2 - 2*xt*(2 + yk^2 - yS^2))) + 
     Qj*Qk*(-2*xs^6 + xs^5*(6 - 6*xt + 4*yk^2 - 4*yS^2) - 
       (-1 + xt)^2*xt*(yk^2 - yS^2)^2 - xs^4*(6 + 5*xt^2 + yk^4 - 12*yS^2 + 
         yS^4 - 2*xt*(7 + 3*yk^2 - 2*yS^2) - 2*yk^2*(-6 + yS^2)) - 
       xs^3*(-2 + xt^3 - 3*yk^4 + 12*yS^2 - 3*yS^4 + 6*yk^2*(-2 + yS^2) - 
         xt^2*(7 + 2*yk^2 + 2*yS^2) + xt*(10 + yk^4 - 6*yS^2 + yS^4 - 
           2*yk^2*(-5 + yS^2))) + xs^2*(-3*yk^4 + 4*yS^2 - 3*yS^4 + 
         xt^3*(1 + 2*yS^2) + yk^2*(-4 + 6*yS^2) + 
         xt^2*(-2 + yk^4 - 2*yS^2 - 2*yk^2*yS^2 + yS^4) + 
         xt*(2 + yk^4 + yS^4 - 2*yk^2*(-1 + yS^2))) + 
       xs*((yk^2 - yS^2)^2 + xt^3*(yk^4 + 2*yS^2 - 2*yk^2*yS^2 + yS^4) - 
         xt^2*(3*yk^4 + 3*yS^4 + yk^2*(2 - 6*yS^2)) + 
         xt*(yk^4 + yS^2*(-2 + yS^2) - 2*yk^2*(-1 + yS^2)))) + 
     Qj^2*(xs^6 + (-1 + xt)^2*xt*(yk^2 - yS^2)^2 + 
       xs^5*(-3 + 3*xt - 2*yk^2 + 2*yS^2) + xs^4*(3 + 3*xt^2 + yk^4 - 
         6*yS^2 + yS^4 - 2*yk^2*(-3 + yS^2) + xt*(-7 - 4*yk^2 + 2*yS^2)) + 
       xs^3*(-1 + xt^3 - 3*yk^4 + 6*yS^2 - 3*yS^4 + 6*yk^2*(-1 + yS^2) - 
         xt^2*(5 + 2*yk^2 + 2*yS^2) + xt*(5 + yk^4 - 4*yS^2 + yS^4 - 
           2*yk^2*(-4 + yS^2))) - xs^2*(-3*yk^4 + 2*yS^2 - 3*yS^4 + 
         xt^3*(1 + 2*yS^2) + yk^2*(-2 + 6*yS^2) + 
         xt*(yk^4 - 2*yk^2*(-2 + yS^2) + (-1 + yS^2)^2) + 
         xt^2*(-2 + yk^4 + yS^4 - 2*yk^2*(1 + yS^2))) - 
       xs*((1 + xt - 3*xt^2 + xt^3)*yk^4 - 2*(1 + xt - 3*xt^2 + xt^3)*yk^2*
          yS^2 + (-1 + xt)*yS^2*(-yS^2 - 2*xt*yS^2 + xt^2*(2 + yS^2)))))*
    \[Lambda]Lik*\[Lambda]Ljk*CC[7])/(2*xs^3*xt^3*(-1 + xs + xt)^2) - 
  (mi^2*Qk*(Qk*(-1 + xs)*xs*(xs^4 + xs^3*(-2 + xt + 2*yk^2 - 2*yS^2) + 
       xs*(-xt + 2*xt^2 + 2*yk^2 - 2*yS^2) + 2*(-1 + xt)*xt*(yk^2 - yS^2) + 
       xs^2*(1 + 2*(-2 + xt)*yk^2 - 2*(-2 + xt)*yS^2)) - 
     Qj*(-((-1 + xt)^2*xt*(yk^2 - yS^2)^2) + 
       xs^4*(xt^2 + 2*xt*yk^2 - (yk^2 - yS^2)^2) + 
       xs^3*(xt^3 + xt^2*(-1 + 6*yk^2 - 2*yS^2) + 3*(yk^2 - yS^2)^2 - 
         xt*(yk^4 + yS^4 - 2*yk^2*(-2 + yS^2))) + 
       xs*((1 + xt - 3*xt^2 + xt^3)*yk^4 - 2*(1 + xt - 3*xt^2 + xt^3)*yk^2*
          yS^2 + (-1 + xt)*yS^2*(-yS^2 - 2*xt*yS^2 + xt^2*(2 + yS^2))) + 
       xs^2*(xt^3*(-1 + 4*yk^2 - 2*yS^2) - 3*(yk^2 - yS^2)^2 + 
         xt*(yk^4 + yS^4 - 2*yk^2*(-1 + yS^2)) + 
         xt^2*(yk^4 - 2*yk^2*(3 + yS^2) + yS^2*(4 + yS^2)))))*\[Lambda]Lik*
    \[Lambda]Ljk*CC[8])/(2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*(-(Qk^2*(-1 + xt)*xt*(2*xs^3*(xt - yk^2 + yS^2) + 
        (-1 + xt)^2*xt*(xt - 2*yk^2 + 2*yS^2) + 
        xs*xt*(1 + 3*xt^2 + 4*yk^2 - 4*yS^2 - 4*xt*(1 + yk^2 - yS^2)) + 
        2*xs^2*(2*xt^2 + yk^2 - yS^2 + xt*(-1 - 2*yk^2 + 2*yS^2)))) + 
     Qj^2*(4*xs^4*xt*yS^2 - (-1 + xt)^3*xt*(xt - yk^2 + yS^2)^2 - 
       xs*(-1 + xt)^2*(3*xt^3 + (yk^2 - yS^2)^2 + xt*(yk^2 - yS^2)^2 + 
         xt^2*(-1 - 4*yk^2 + 2*yS^2)) - xs^3*(xt^3 + (yk^2 - yS^2)^2 - 
         xt^2*(1 + 10*yS^2) - xt*(yk^4 - 6*yS^2 - 2*yk^2*yS^2 + yS^4)) + 
       xs^2*(-3*xt^4 + 2*(yk^2 - yS^2)^2 + xt^3*(5 + 2*yk^2 + 6*yS^2) + 
         xt*(-3*yk^4 + 2*yS^2 + 6*yk^2*yS^2 - 3*yS^4) + 
         xt^2*(-2 + yk^4 - 8*yS^2 + yS^4 - 2*yk^2*(1 + yS^2)))) - 
     Qj*Qk*(4*xs^4*xt*yS^2 - (-1 + xt)^3*xt*(2*xt^2 - 4*xt*(yk^2 - yS^2) + 
         (yk^2 - yS^2)^2) - xs*(-1 + xt)^2*(6*xt^3 + (yk^2 - yS^2)^2 + 
         xt*(yk^2 - yS^2)^2 + xt^2*(-2 - 8*yk^2 + 6*yS^2)) - 
       xs^3*(3*xt^3 + (yk^2 - yS^2)^2 - xt^2*(3 + 2*yk^2 + 8*yS^2) - 
         xt*(yk^4 + yS^2*(-4 + yS^2) - 2*yk^2*(1 + yS^2))) + 
       xs^2*(-7*xt^4 + 2*(yk^2 - yS^2)^2 + xt^3*(11 + 6*yk^2 + 2*yS^2) + 
         xt^2*(-4 + yk^4 - 2*yS^2 + yS^4 - 2*yk^2*(4 + yS^2)) + 
         xt*(-3*yk^4 - 3*yS^4 + yk^2*(2 + 6*yS^2)))))*\[Lambda]Lik*
    \[Lambda]Ljk*CC[9])/(2*xs^3*xt^3*(-1 + xs + xt)^2) - 
  (mi^2*Qk*\[Lambda]Ljk*
    (-(Qk*(-1 + xt)*xt*((-1 + xt)^2*xt*(xt + 2*yk^2 - 2*yS^2) + 
        2*xs^3*(yk^2 - yS^2) + xs*xt*(1 + 3*xt^2 - 4*yk^2 + 4*yS^2 + 
          4*xt*(-1 + yk^2 - yS^2)) + 2*xs^2*(xt^2 - yk^2 + yS^2 + 
          2*xt*(yk^2 - yS^2)))*\[Lambda]Lik) + 
     Qj*(-((-1 + xt)^3*xt*(yk^2 - yS^2)^2*\[Lambda]Lik) + 
       xs*(-1 + xt)^2*(2*xt^2*yS^2 - (yk^2 - yS^2)^2 - 
         xt*(yk^4 - 2*yk^2*(1 + yS^2) + yS^2*(2 + yS^2)))*\[Lambda]Lik + 
       4*xs^4*xt*yk*(yk*\[Lambda]Lik + \[Lambda]Rik) - 
       xs^2*(-1 + xt)*(xt^3*\[Lambda]Lik + 2*(yk^2 - yS^2)^2*\[Lambda]Lik - 
         2*xt^2*((1 + 3*yk^2 + yS^2)*\[Lambda]Lik + 2*yk*\[Lambda]Rik) + 
         xt*(-(yk^4*\[Lambda]Lik) + 2*yk^2*yS^2*\[Lambda]Lik - 
           yS^2*(-2 + yS^2)*\[Lambda]Lik + 4*yk*\[Lambda]Rik)) + 
       xs^3*(-(xt^3*\[Lambda]Lik) - (yk^2 - yS^2)^2*\[Lambda]Lik + 
         xt*(yk^4*\[Lambda]Lik + yS^4*\[Lambda]Lik - 2*yk^2*(3 + yS^2)*
            \[Lambda]Lik - 8*yk*\[Lambda]Rik) + 
         xt^2*(\[Lambda]Lik + 10*yk^2*\[Lambda]Lik + 8*yk*\[Lambda]Rik))))*
    CC[10])/(2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*Qk^2*(xs + xt)*(xs^3 - xs*xt^2 - xt^2*(xt + 2*yk^2 - 2*yS^2) - 
     xs^2*(xt - 2*yk^2 + 2*yS^2))*\[Lambda]Lik*\[Lambda]Ljk*CC[11])/
   (2*xs^3*xt^3) - (mi^2*(Qj - Qk)^2*(xs - xt)*(xs + xt)^2*
    (xs + xt - 2*yk^2 + 2*yS^2)*\[Lambda]Lik*\[Lambda]Ljk*CC[12])/
   (2*xs^3*xt^3) - (mi^4*(Qj - Qk)*Qk*(xs*(xt - (yk - yS)^2) - 
     (-1 + xt)*(yk - yS)^2)*((-1 + xt)*xt*(yk^2 - yS^2) + 
     xs*(-2*xt + xt^2 + yk^2 - yS^2) + xs^2*(xt - yk^2 + yS^2))*
    (-((-1 + xt)*(yk + yS)^2) + xs*(xt - (yk + yS)^2))*\[Lambda]Lik*
    \[Lambda]Ljk*DD[1])/(2*xs^3*xt^3*(-1 + xs + xt)^2) - 
  (mi^4*(Qj - Qk)*Qk*(xs*(xt - (yk - yS)^2) - (-1 + xt)*(yk - yS)^2)*
    ((-1 + xt)*xt*(yk^2 - yS^2) + xs*(xt^2 + yk^2 - yS^2) + 
     xs^2*(xt - yk^2 + yS^2))*(-((-1 + xt)*(yk + yS)^2) + 
     xs*(xt - (yk + yS)^2))*\[Lambda]Lik*\[Lambda]Ljk*DD[2])/
   (2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^4*(Qj - Qk)^2*(xs^2 + xt*(yk^2 - yS^2) + xs*(xt - yk^2 + yS^2))*
    (xs^3 + (-1 + xt)*(yk^2 - yS^2)^2 + xs^2*(-1 + xt - 2*yk^2 + 2*yS^2) + 
     xs*(yk^4 + yS^2*(-2 - 2*xt + yS^2) - 2*yk^2*(-1 + xt + yS^2)))*
    \[Lambda]Lik*\[Lambda]Ljk*DD[3])/(2*xs^3*xt^3) - 
  (mi^4*Qk^2*(xs - xt)*(xs + yk^2 - yS^2)*
    (xs^3 + xs^2*(-1 + xt + 2*yk^2 - 2*yS^2) + (-1 + xt)*(yk^2 - yS^2)^2 + 
     xs*(yk^4 + yS^2*(2 - 2*xt + yS^2) - 2*yk^2*(1 + xt + yS^2)))*
    \[Lambda]Lik*\[Lambda]Ljk*DD[4])/(2*xs^3*xt^3) - 
  (mi^4*(Qj - Qk)^2*(xs*(xt + yk^2 - yS^2) + xt*(xt - yk^2 + yS^2))*
    (xt^3 + (-1 + xs)*(yk^2 - yS^2)^2 + xt^2*(-1 + xs - 2*yk^2 + 2*yS^2) + 
     xt*(yk^4 + yS^2*(-2 - 2*xs + yS^2) - 2*yk^2*(-1 + xs + yS^2)))*
    \[Lambda]Lik*\[Lambda]Ljk*DD[5])/(2*xs^3*xt^3) + 
  (mi^4*Qk^2*\[Lambda]Ljk*((-1 + xt)*xt*(xt + yk^2 - yS^2)^3*\[Lambda]Lik + 
     xs*(2*xt^4 + xt^3*(-1 + 2*yk^2 - 6*yS^2) - xt*(yk^2 - yS^2)^2 + 
       (yk^2 - yS^2)^3 - xt^2*(yk^2 - yS^2)*(3 + 4*yS^2))*\[Lambda]Lik + 
     xs^2*(xt^3*\[Lambda]Lik - (yk^2 - yS^2)^3*\[Lambda]Lik + 
       xt*(5*yk^4 - 6*yk^2*yS^2 + yS^4)*\[Lambda]Lik - 
       xt^2*(9*yk^2*\[Lambda]Lik + 3*yS^2*\[Lambda]Lik + 4*yk*\[Lambda]Rik)))*
    DD[6])/(2*xs^3*xt^3), 
 (-2*Qj*(-2*Qk + Qj*(1 + yk^2 - yS^2))*\[Lambda]Rik*\[Lambda]Rjk)/
   (xs*(-1 + xt)*xt^2) + (2*Qj^2*yk^2*\[Lambda]Rik*\[Lambda]Rjk*BB[1])/
   (xs*xt^2 - xs*xt^3) + (2*Qj^2*yS^2*\[Lambda]Rik*\[Lambda]Rjk*BB[2])/
   (xs*(-1 + xt)*xt^2) + 
  (2*Qj*(2*Qk*(xs - xs^2 + (-1 + xt)^2)*xt*\[Lambda]Rik + 
     Qj*(-1 + xs)*(-1 + xs + xt)*(2*(-1 + xt)*yk*\[Lambda]Lik + 
       (-2 + xt)*yk^2*\[Lambda]Rik + (xt + 2*yS^2 - xt*yS^2)*\[Lambda]Rik))*
    \[Lambda]Rjk*BB[3])/((-1 + xs)*xs*(-1 + xt)^2*xt^2*(-1 + xs + xt)) - 
  (4*Qj*Qk*\[Lambda]Rik*\[Lambda]Rjk*BB[4])/((-1 + xs)*xs*xt*
    (-1 + xs + xt)) - (2*Qj*(-2*Qk*xs*xt*\[Lambda]Rik + 
     Qj*(-1 + xs + xt)*(2*(-1 + xt)*yk*\[Lambda]Lik - yk^2*\[Lambda]Rik + 
       (xt + yS^2)*\[Lambda]Rik))*\[Lambda]Rjk*BB[5])/
   (xs*(-1 + xt)^2*xt^2*(-1 + xs + xt)) + 
  (2*mi^2*Qk^2*yk*\[Lambda]Lik*\[Lambda]Rjk*CC[1])/(xt^2*(-1 + xs + xt)) - 
  (2*mi^2*(Qj - Qk)*Qk*((-1 + xs + xt)*yk*\[Lambda]Lik + xt*\[Lambda]Rik)*
    \[Lambda]Rjk*CC[2])/(xt^2*(-1 + xs + xt)^2) - 
  (2*mi^2*Qk*(Qj*(-1 + xs + xt)*yk*\[Lambda]Lik + Qj*xt*\[Lambda]Rik - 
     Qk*xt*\[Lambda]Rik)*\[Lambda]Rjk*CC[3])/(xs*xt*(-1 + xs + xt)^2) - 
  (4*mi^2*Qk^2*yk*\[Lambda]Lik*\[Lambda]Rjk*CC[5])/(xs*xt^2) + 
  (2*mi^2*Qk*(Qj*(1 + xs^2 - 3*xt + 2*xt^2 + xs*(-2 + 3*xt))*yk*
      \[Lambda]Lik - Qj*(-1 + xs)*xt*\[Lambda]Rik + 
     Qk*(-1 + xs)*xt*\[Lambda]Rik)*\[Lambda]Rjk*CC[8])/
   (xs*xt^2*(-1 + xs + xt)^2) - 
  (2*mi^2*(Qj - Qk)*(2*Qj*(-1 + xs + xt)^2*yS^2*\[Lambda]Rik + 
     Qk*(-1 + xt)^2*((-1 + xs + xt)*yk*\[Lambda]Lik + xt*\[Lambda]Rik))*
    \[Lambda]Rjk*CC[9])/(xs*(-1 + xt)*xt^2*(-1 + xs + xt)^2) + 
  (2*mi^2*Qk*yk*(Qk*(-1 + xt)^2*\[Lambda]Lik - 2*Qj*(-1 + xs + xt)*
      ((-1 + xt)*\[Lambda]Lik - yk*\[Lambda]Rik))*\[Lambda]Rjk*CC[10])/
   (xs*(-1 + xt)*xt^2*(-1 + xs + xt)) + 
  (2*mi^4*(Qj - Qk)*Qk*(xs^2*yk*(xt - yk^2 + yS^2)*\[Lambda]Lik - 
     (-1 + xt)*((-1 + xt)*yk^3*\[Lambda]Lik - (-1 + xt)*yk*yS^2*
        \[Lambda]Lik + xt*yk^2*\[Lambda]Rik + xt*yS^2*\[Lambda]Rik) + 
     xs*(2*yk*(yk^2 - yS^2)*\[Lambda]Lik + xt^2*(yk*\[Lambda]Lik + 
         \[Lambda]Rik) - xt*(yk*\[Lambda]Lik + 2*yk^3*\[Lambda]Lik - 
         2*yk*yS^2*\[Lambda]Lik + yk^2*\[Lambda]Rik + yS^2*\[Lambda]Rik)))*
    \[Lambda]Rjk*DD[2])/(xs*xt^2*(-1 + xs + xt)^2) + 
  (2*mi^4*Qk^2*yk*(xs^2*\[Lambda]Lik + (-1 + xt)*yk^2*\[Lambda]Lik - 
     (-1 + xt)*yS^2*\[Lambda]Lik + xs*(-1 + xt + yk^2 - yS^2)*\[Lambda]Lik + 
     2*xt*yk*\[Lambda]Rik)*\[Lambda]Rjk*DD[4])/(xs*xt^2*(-1 + xs + xt)) - 
  (2*mi^4*Qk^2*yk*(xt^2*\[Lambda]Lik - (-1 + xs)*(yk^2 - yS^2)*\[Lambda]Lik + 
     xt*(-1 + xs - yk^2 + yS^2)*\[Lambda]Lik - 2*xt*yk*\[Lambda]Rik)*
    \[Lambda]Rjk*DD[6])/(xs*xt^2*(-1 + xs + xt)), 
 (2*Qj*(xs - xt)*(-2*Qk*xs*xt + Qj*(-1 + xt)*(yk^2 - yS^2) + 
     Qj*xs*(xt + yk^2 - yS^2))*\[Lambda]Lik*\[Lambda]Ljk)/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)) + 
  (2*Qj^2*(xs - xt)*yk^2*\[Lambda]Lik*\[Lambda]Ljk*BB[1])/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2) + 
  (2*Qj^2*(-xs + xt)*yS^2*\[Lambda]Lik*\[Lambda]Ljk*BB[2])/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2) + 
  (2*Qj*(-2*xs + xs^2 - (-2 + xt)*xt)*(2*Qk + Qj*(-1 + yk^2 - yS^2))*
    \[Lambda]Lik*\[Lambda]Ljk*BB[3])/((-1 + xs)^2*xs*(-1 + xt)^2*xt*
    (-1 + xs + xt)) + 
  (2*Qj*(2*Qk*xs + Qj*(-yk^2 + yS^2 + xs*(-1 + 2*yk^2 - 2*yS^2)))*
    \[Lambda]Lik*\[Lambda]Ljk*BB[4])/((-1 + xs)^2*xs^2*xt*(-1 + xs + xt)) + 
  (2*Qj*(-2*Qk*xt + Qj*(xt + yk^2 - 2*xt*yk^2 - yS^2 + 2*xt*yS^2))*
    \[Lambda]Lik*\[Lambda]Ljk*BB[5])/(xs*(-1 + xt)^2*xt^2*(-1 + xs + xt)) + 
  (2*mi^2*Qk*(Qk*xs*((-1 + xs)*(yk^2 - yS^2) + xt*(1 + yk^2 - yS^2)) + 
     Qj*((-1 + xt)*xt*yk^2 - xs*yS^2 + xs^2*yS^2 + xs*xt*(-1 + yk^2 + yS^2)))*
    \[Lambda]Lik*\[Lambda]Ljk*CC[1])/(xs*xt^2*(-1 + xs + xt)^3) - 
  (2*mi^2*(Qj^2*(xs^2 + (-1 + xt)*xt + xs*(-1 + 2*xt))*yS^2 + 
     Qk^2*xs*(-((-1 + xs)*(yk^2 - yS^2)) + xt*(1 - yk^2 + yS^2)) + 
     Qj*Qk*(-((-1 + xt)*xt*yS^2) + xs^2*(yk^2 - 2*yS^2) + 
       xs*(-yk^2 + 2*yS^2 + xt*(-1 + yk^2 - 3*yS^2))))*\[Lambda]Lik*
    \[Lambda]Ljk*CC[2])/(xs*xt^2*(-1 + xs + xt)^3) - 
  (2*mi^2*Qk*(Qk*xt*((-1 + xt)*(yk^2 - yS^2) + xs*(1 + yk^2 - yS^2)) + 
     Qj*(-(xs*yk^2) + xs^2*yk^2 + (-1 + xt)*xt*yS^2 + 
       xs*xt*(-1 + yk^2 + yS^2)))*\[Lambda]Lik*\[Lambda]Ljk*CC[3])/
   (xs^2*xt*(-1 + xs + xt)^3) + 
  (2*mi^2*(Qj^2*(xs^2 + (-1 + xt)*xt + xs*(-1 + 2*xt))*yS^2 + 
     Qk^2*xt*(-((-1 + xt)*(yk^2 - yS^2)) + xs*(1 - yk^2 + yS^2)) + 
     Qj*Qk*(-(xs^2*yS^2) + (-1 + xt)*xt*(yk^2 - 2*yS^2) + 
       xs*(yS^2 + xt*(-1 + yk^2 - 3*yS^2))))*\[Lambda]Lik*\[Lambda]Ljk*CC[4])/
   (xs^2*xt*(-1 + xs + xt)^3) - (4*mi^2*Qk^2*(xs - xt)*yk^2*\[Lambda]Lik*
    \[Lambda]Ljk*CC[5])/(xs^2*xt^2*(-1 + xs + xt)) + 
  (4*mi^2*(Qj - Qk)^2*(xs - xt)*yS^2*\[Lambda]Lik*\[Lambda]Ljk*CC[6])/
   (xs^2*xt^2*(-1 + xs + xt)) - 
  (2*mi^2*(Qj - Qk)*(Qj*(xs^4 - (-1 + xt)*xt + xs^3*(-3 + 4*xt) + 
       xs^2*(3 - 7*xt + 5*xt^2) + xs*(-1 + 2*xt - 4*xt^2 + 2*xt^3))*yS^2 - 
     Qk*(-1 + xs)^2*xt*(xs*(-1 + yk^2 - yS^2) + (-1 + xt)*(yk^2 - yS^2)))*
    \[Lambda]Lik*\[Lambda]Ljk*CC[7])/((-1 + xs)*xs^2*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^2*Qk*(-(Qk*(-1 + xs)^2*xt*((-1 + xt)*(yk^2 - yS^2) + 
        xs*(1 + yk^2 - yS^2))) + Qj*(xs^4*yk^2 - (-1 + xt)*xt*yS^2 + 
       xs^3*(xt - 3*yk^2 + 5*xt*yk^2 - xt*yS^2) + 
       xs*(-yk^2 + 2*xt^3*yk^2 + xt*(1 + 5*yk^2 - 3*yS^2) + 
         xt^2*(-6*yk^2 + 2*yS^2)) + xs^2*(3*yk^2 + xt^2*(6*yk^2 - yS^2) + 
         xt*(-2 - 10*yk^2 + 3*yS^2))))*\[Lambda]Lik*\[Lambda]Ljk*CC[8])/
   ((-1 + xs)*xs^2*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^2*(Qj - Qk)*(Qj*(2*xs^3*xt + (-1 + xt)^3*xt + 
       xs*(-1 + xt)^2*(1 + 4*xt) + xs^2*(-1 - 4*xt + 5*xt^2))*yS^2 - 
     Qk*xs*(-1 + xt)^2*(xt*(-1 + yk^2 - yS^2) + (-1 + xs)*(yk^2 - yS^2)))*
    \[Lambda]Lik*\[Lambda]Ljk*CC[9])/(xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)^3) - 
  (2*mi^2*Qk*(-(Qk*xs*(-1 + xt)^2*((-1 + xs)*(yk^2 - yS^2) + 
        xt*(1 + yk^2 - yS^2))) + Qj*(2*xs^3*xt*yk^2 + (-1 + xt)^3*xt*yk^2 + 
       xs*(-1 + xt)^2*(xt + 5*xt*yk^2 + yS^2 - xt*yS^2) + 
       xs^2*(-1 + xt)*(6*xt*yk^2 + yS^2 - xt*yS^2)))*\[Lambda]Lik*
    \[Lambda]Ljk*CC[10])/(xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)^3) - 
  (2*mi^4*(Qj - Qk)*Qk*(-((-1 + xt)^2*xt*yk^2*(yk^2 - yS^2)) + 
     xs^3*yS^2*(xt + yk^2 - yS^2) + xs^2*(-2*yk^2*yS^2 + 2*yS^4 + 
       xt^2*(-1 + yk^2 + yS^2) + xt*(yk^2 - yk^4 + 3*yk^2*yS^2 - 2*yS^4)) + 
     xs*(xt^3*yk^2 + yS^2*(yk^2 - yS^2) + xt^2*(-2*yk^4 + yS^2 + 
         3*yk^2*yS^2 - yS^4) + xt*(2*yk^4 + yS^2*(-1 + 2*yS^2) - 
         yk^2*(1 + 4*yS^2))))*\[Lambda]Lik*\[Lambda]Ljk*DD[1])/
   (xs^2*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^4*(Qj - Qk)*Qk*((-1 + xt)^2*xt*yS^2*(yk^2 - yS^2) + 
     xs^3*yk^2*(xt - yk^2 + yS^2) + xs^2*(2*yk^2*(yk^2 - yS^2) + 
       xt^2*(-1 + yk^2 + yS^2) + xt*(-2*yk^4 + yS^2 + 3*yk^2*yS^2 - yS^4)) + 
     xs*(-yk^4 + xt^3*yS^2 + yk^2*yS^2 + xt^2*(-yk^4 - 2*yS^4 + 
         yk^2*(1 + 3*yS^2)) + xt*(2*yk^4 + yS^2*(-1 + 2*yS^2) - 
         yk^2*(1 + 4*yS^2))))*\[Lambda]Lik*\[Lambda]Ljk*DD[2])/
   (xs^2*xt^2*(-1 + xs + xt)^3) - 
  (2*mi^4*(Qj - Qk)^2*yS^2*(xs^2 + xt*(yk^2 - yS^2) + xs*(xt - yk^2 + yS^2))*
    \[Lambda]Lik*\[Lambda]Ljk*DD[3])/(xs^2*xt^2*(-1 + xs + xt)) + 
  (2*mi^4*Qk^2*yk^2*(xs^2 + xs*(xt + yk^2 - yS^2) + xt*(-yk^2 + yS^2))*
    \[Lambda]Lik*\[Lambda]Ljk*DD[4])/(xs^2*xt^2*(-1 + xs + xt)) + 
  (2*mi^4*(Qj - Qk)^2*yS^2*(xs*(xt + yk^2 - yS^2) + xt*(xt - yk^2 + yS^2))*
    \[Lambda]Lik*\[Lambda]Ljk*DD[5])/(xs^2*xt^2*(-1 + xs + xt)) - 
  (2*mi^4*Qk^2*yk^2*(xt*(xt + yk^2 - yS^2) + xs*(xt - yk^2 + yS^2))*
    \[Lambda]Lik*\[Lambda]Ljk*DD[6])/(xs^2*xt^2*(-1 + xs + xt))};



coef1b={(Qj*\[Lambda]Ljk*(2*Qk*\[Lambda]Lik + Qj*(-1 + (-3 + xs)*yk^2 - 
       (-3 + xs)*yS^2)*\[Lambda]Lik - 2*Qj*xs*yk*\[Lambda]Rik))/
   ((-1 + xs)*xs^2*xt) + (Qj^2*yk^2*\[Lambda]Ljk*
    ((-3 + xs)*yk^2*\[Lambda]Lik - (-3 + xs)*yS^2*\[Lambda]Lik - 
     2*xs*yk*\[Lambda]Rik)*BB[1])/((-1 + xs)*xs^2*xt*(yk^2 - yS^2)) + 
  (Qj^2*yS^2*\[Lambda]Ljk*(-((-3 + xs)*yk^2*\[Lambda]Lik) + 
     (-3 + xs)*yS^2*\[Lambda]Lik + 2*xs*yk*\[Lambda]Rik)*BB[2])/
   ((-1 + xs)*xs^2*xt*(yk^2 - yS^2)) - 
  (Qj*\[Lambda]Ljk*(-2*Qk*(xt^2 + xs^2*(1 - 2*xt^2) + 2*xs*xt*(-1 + xt^2))*
      \[Lambda]Lik + Qj*(xt^2*(\[Lambda]Lik + (1 - 2*xt)*yk^2*\[Lambda]Lik + 
         (-1 + 2*xt)*yS^2*\[Lambda]Lik - 2*(-1 + xt)*yk*\[Lambda]Rik) + 
       xs*(-1 + xt)*xt*((2 + xt*(1 + yk^2 - yS^2))*\[Lambda]Lik + 
         2*(-1 + xt)*yk*\[Lambda]Rik) + 
       xs^2*((1 - yk^2 + yS^2 + xt^2*(-3 + yk^2 - yS^2) + 
           xt*(1 + yk^2 - yS^2))*\[Lambda]Lik + 2*(-1 + xt)*xt*yk*
          \[Lambda]Rik)))*BB[3])/((-1 + xs)*xs^2*(-1 + xt)*xt^2*(xs + xt)) + 
  (Qj*\[Lambda]Ljk*(2*Qk*xs*(xs - 2*xt)*\[Lambda]Lik + 
     Qj*(-(xs^2*\[Lambda]Lik) + xs*(yk^2 - yS^2)*\[Lambda]Lik + 
       xt*(yk^2*\[Lambda]Lik - yS^2*\[Lambda]Lik - 2*yk*\[Lambda]Rik) + 
       xs*xt*(\[Lambda]Lik + 4*yk*\[Lambda]Rik)))*BB[4])/
   ((-1 + xs)*xs^2*xt^2) + (Qj*(2*Qk*xt - Qj*(xt - yk^2 + yS^2))*\[Lambda]Lik*
    \[Lambda]Ljk*BB[5])/(xs^2*(-1 + xt)*xt) - 
  (Qk^2*(xs - xt)^2*\[Lambda]Lik*\[Lambda]Ljk*BB[6])/(xs^2*xt^2*(xs + xt)) + 
  ((Qj - Qk)^2*(xs - xt)^2*\[Lambda]Lik*\[Lambda]Ljk*BB[7])/
   (xs^2*xt^2*(xs + xt)) - (mi^2*Qk*\[Lambda]Ljk*
    (Qk*(-1 + xs)*xs*(xs^4*\[Lambda]Lik + 2*(-1 + xt)^2*xt*(yk^2 - yS^2)*
        \[Lambda]Lik + xs^3*((-2 + xt + 2*yk^2 - 2*yS^2)*\[Lambda]Lik + 
         4*xt*yk*\[Lambda]Rik) + xs^2*((1 - 2*xt^2 - 4*yk^2 + 4*xt*yk^2 + 
           4*yS^2 - 4*xt*yS^2)*\[Lambda]Lik + 4*(-2 + xt)*xt*yk*
          \[Lambda]Rik) + xs*(2*(yk^2 - yS^2)*\[Lambda]Lik - 
         xt*(\[Lambda]Lik + 6*yk^2*\[Lambda]Lik - 6*yS^2*\[Lambda]Lik - 
           4*yk*\[Lambda]Rik) + 2*xt^2*(\[Lambda]Lik + 2*yk^2*\[Lambda]Lik - 
           2*yS^2*\[Lambda]Lik - 2*yk*\[Lambda]Rik))) + 
     Qj*(-((-1 + xt)^2*xt*(yk^2 - yS^2)^2*\[Lambda]Lik) + 
       xs^4*((yk^2 - yS^2)^2*\[Lambda]Lik + 4*xt*yk^3*\[Lambda]Rik + 
         xt^2*(\[Lambda]Lik - 4*yk*\[Lambda]Rik) - 
         2*xt*yS^2*(\[Lambda]Lik + 2*yk*\[Lambda]Rik)) - 
       xs^3*(3*(yk^2 - yS^2)^2*\[Lambda]Lik + xt^3*(\[Lambda]Lik + 
           12*yk*\[Lambda]Rik) - xt*(3*yk^4*\[Lambda]Lik - 
           6*yk^2*yS^2*\[Lambda]Lik + yS^2*(4 + 3*yS^2)*\[Lambda]Lik - 
           12*yk^3*\[Lambda]Rik + 12*yk*yS^2*\[Lambda]Rik) + 
         xt^2*((1 + 2*yk^2 + 2*yS^2)*\[Lambda]Lik - 
           4*yk*(3 + 2*yk^2 - 2*yS^2)*\[Lambda]Rik)) - 
       xs*(-1 + xt)*(-((yk^2 - yS^2)^2*\[Lambda]Lik) + 
         4*xt^3*yk*(yk*\[Lambda]Lik - \[Lambda]Rik) + 4*xt*(yk^2 - yS^2)*
          (yk^2*\[Lambda]Lik - yS^2*\[Lambda]Lik - yk*\[Lambda]Rik) - 
         xt^2*(yk^4*\[Lambda]Lik + yS^4*\[Lambda]Lik - 2*yk^2*(-1 + yS^2)*
            \[Lambda]Lik - 4*yk^3*\[Lambda]Rik + 4*yk*(-1 + yS^2)*
            \[Lambda]Rik)) + xs^2*(3*(yk^2 - yS^2)^2*\[Lambda]Lik - 
         8*xt^4*yk*\[Lambda]Rik - xt*(7*yk^4*\[Lambda]Lik - 
           14*yk^2*yS^2*\[Lambda]Lik + yS^2*(2 + 7*yS^2)*\[Lambda]Lik - 
           12*yk^3*\[Lambda]Rik + 12*yk*yS^2*\[Lambda]Rik) + 
         xt^3*(\[Lambda]Lik - 6*yk^2*\[Lambda]Lik + 4*yk*(5 + yk^2 - yS^2)*
            \[Lambda]Rik) + xt^2*(3*yk^4*\[Lambda]Lik + yk^2*(4 - 6*yS^2)*
            \[Lambda]Lik + yS^2*(2 + 3*yS^2)*\[Lambda]Lik - 
           16*yk^3*\[Lambda]Rik + 4*yk*(-3 + 4*yS^2)*\[Lambda]Rik))))*CC[1])/
   (2*(-1 + xs)*xs^2*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*(Qj - Qk)*\[Lambda]Ljk*
    (-(Qk*(-1 + xs)*xs*(xs^4*\[Lambda]Lik - 2*(-1 + xt)^2*xt*(yk^2 - yS^2)*
         \[Lambda]Lik + xs^3*((-2 + xt - 2*yk^2 + 2*yS^2)*\[Lambda]Lik - 
          4*xt*yk*\[Lambda]Rik) + xs^2*(\[Lambda]Lik - 4*(-1 + xt)*yk^2*
           \[Lambda]Lik + 4*(-1 + xt)*yS^2*\[Lambda]Lik - 
          4*(-2 + xt)*xt*yk*\[Lambda]Rik) + xs*(2*xt^3*\[Lambda]Lik + 
          xt*(-1 + 6*yk^2 - 6*yS^2)*\[Lambda]Lik + 2*(-yk^2 + yS^2)*
           \[Lambda]Lik - 4*xt*yk*\[Lambda]Rik + xt^2*(-4*yk^2*\[Lambda]Lik + 
            4*yS^2*\[Lambda]Lik + 4*yk*\[Lambda]Rik)))) + 
     Qj*(-1 + xs + xt)*(xs^5*\[Lambda]Lik - (-1 + xt)*xt*(yk^2 - yS^2)^2*
        \[Lambda]Lik - 2*xs^4*((1 + yk^2 - yS^2)*\[Lambda]Lik + 
         2*xt*yk*\[Lambda]Rik) + xs*((1 - 3*xt + xt^2)*yk^4*\[Lambda]Lik - 
         2*(1 - 3*xt + xt^2)*yk^2*yS^2*\[Lambda]Lik + 
         yS^2*(-4*xt^3 + yS^2 - 3*xt*yS^2 + xt^2*(2 + yS^2))*\[Lambda]Lik - 
         4*(-1 + xt)*xt*yk^3*\[Lambda]Rik + 4*(-1 + xt)*xt*yk*yS^2*
          \[Lambda]Rik) + xs^3*((1 + xt - xt^2 - 2*xt*yk^2 + yk^4 - 4*yS^2 + 
           yS^4 - 2*yk^2*(-2 + yS^2))*\[Lambda]Lik - 
         4*xt*yk*(-2 + xt - yk^2 + yS^2)*\[Lambda]Rik) + 
       xs^2*(-2*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4)*\[Lambda]Lik + 
         xt*((-1 + 2*yk^4 + 2*yS^4 + yk^2*(2 - 4*yS^2))*\[Lambda]Lik - 
           4*yk*(1 + 2*yk^2 - 2*yS^2)*\[Lambda]Rik) + 
         xt^2*(\[Lambda]Lik - 2*yS^2*\[Lambda]Lik + 4*yk*(1 + yk^2 - yS^2)*
            \[Lambda]Rik))))*CC[2])/(2*(-1 + xs)*xs^2*xt^3*
    (-1 + xs + xt)^2) - (mi^2*Qk*\[Lambda]Ljk*
    (Qk*(-1 + xt)*xt^2*((xt^2 - 2*yk^2 + 2*yS^2 + xt*(-1 + 2*yk^2 - 2*yS^2))*
        \[Lambda]Lik + 4*xs^2*yk*\[Lambda]Rik + 
       xs*((-1 + 3*xt + 2*yk^2 - 2*yS^2)*\[Lambda]Lik + 
         4*(-1 + xt)*yk*\[Lambda]Rik)) + 
     Qj*((-1 + xt)^2*xt*(yk^2 - yS^2)^2*\[Lambda]Lik + 
       xs^3*((yk^2 - yS^2)^2*\[Lambda]Lik + xt^2*(\[Lambda]Lik + 
           4*yk*\[Lambda]Rik) + xt*(4*yk^2*\[Lambda]Lik - 
           2*yS^2*\[Lambda]Lik + 4*yk^3*\[Lambda]Rik - 
           4*yk*yS^2*\[Lambda]Rik)) + xs*(-1 + xt)*
        (-((yk^2 - yS^2)^2*\[Lambda]Lik) + xt*(yk^2 - yS^2)*
          ((-2 + 3*yk^2 - 3*yS^2)*\[Lambda]Lik - 4*yk*\[Lambda]Rik) + 
         2*xt^2*(yk^2*\[Lambda]Lik - 2*yS^2*\[Lambda]Lik + 
           2*yk^3*\[Lambda]Rik - 2*yk*yS^2*\[Lambda]Rik)) + 
       xs^2*(-2*(yk^2 - yS^2)^2*\[Lambda]Lik + xt^3*(3*\[Lambda]Lik + 
           4*yk*\[Lambda]Rik) + xt*(3*yk^4*\[Lambda]Lik - 6*yk^2*(1 + yS^2)*
            \[Lambda]Lik + yS^2*(4 + 3*yS^2)*\[Lambda]Lik - 
           8*yk^3*\[Lambda]Rik + 8*yk*yS^2*\[Lambda]Rik) + 
         xt^2*((-2 + 6*yk^2 - 6*yS^2)*\[Lambda]Lik + 
           4*yk*(-1 + 2*yk^2 - 2*yS^2)*\[Lambda]Rik))))*CC[3])/
   (2*xs^3*xt^2*(-1 + xs + xt)^2) - 
  (mi^2*(Qj - Qk)*\[Lambda]Ljk*
    (Qk*(-1 + xt)*xt^2*((xt^2 + 2*(yk^2 - yS^2) + xt*(-1 - 2*yk^2 + 2*yS^2))*
        \[Lambda]Lik - 2*xs^2*(\[Lambda]Lik + 2*yk*\[Lambda]Rik) + 
       xs*((1 + xt - 2*yk^2 + 2*yS^2)*\[Lambda]Lik - 4*(-1 + xt)*yk*
          \[Lambda]Rik)) + Qj*(-1 + xs + xt)*
      (-((-1 + xt)*xt*(xt - yk^2 + yS^2)^2*\[Lambda]Lik) + 
       xs*((yk^2 - yS^2)^2*\[Lambda]Lik + 4*xt^3*yk*\[Lambda]Rik - 
         2*xt*(yk^2 - yS^2)*(yk^2*\[Lambda]Lik - yS^2*\[Lambda]Lik - 
           2*yk*\[Lambda]Rik) + xt^2*((-1 + 2*yk^2)*\[Lambda]Lik - 
           4*yk*(1 + yk^2 - yS^2)*\[Lambda]Rik)) + 
       xs^2*(-((yk^2 - yS^2)^2*\[Lambda]Lik) + 
         xt^2*(\[Lambda]Lik + 4*yk*\[Lambda]Rik) - 
         2*xt*(2*yk^3*\[Lambda]Rik + yS^2*(\[Lambda]Lik - 
             2*yk*\[Lambda]Rik)))))*CC[4])/(2*xs^3*xt^2*(-1 + xs + xt)^2) + 
  (mi^2*Qk^2*\[Lambda]Ljk*(xs^4*\[Lambda]Lik + 
     xt*(xt^3 + xt^2*(-1 + 2*yk^2 - 2*yS^2) - 2*(yk^2 - yS^2)^2 + 
       2*xt*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2)))*\[Lambda]Lik + 
     xs^3*((-1 + 2*yk^2 - 2*yS^2)*\[Lambda]Lik + 4*xt*yk*\[Lambda]Rik) + 
     xs*(-2*(yk^2 - yS^2)^2*\[Lambda]Lik + 4*xt*(yk^2 - yS^2)*
        ((-1 + yk^2 - yS^2)*\[Lambda]Lik - 2*yk*\[Lambda]Rik) + 
       2*xt^3*(\[Lambda]Lik + 2*yk*\[Lambda]Rik) + 
       xt^2*((-1 + 2*yk^2 - 6*yS^2)*\[Lambda]Lik + 
         4*yk*(-1 + 2*yk^2 - 2*yS^2)*\[Lambda]Rik)) + 
     xs^2*(2*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2))*\[Lambda]Lik + 
       8*xt^2*yk*\[Lambda]Rik + xt*((1 + 2*yk^2 - 6*yS^2)*\[Lambda]Lik + 
         4*yk*(-1 + 2*yk^2 - 2*yS^2)*\[Lambda]Rik)))*CC[5])/(2*xs^3*xt^3) - 
  (mi^2*(Qj - Qk)^2*\[Lambda]Ljk*(xs^4*\[Lambda]Lik + 
     xt*(xt^3 - 2*(yk^2 - yS^2)^2 + xt^2*(-1 - 2*yk^2 + 2*yS^2) + 
       2*xt*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4))*\[Lambda]Lik - 
     xs^3*(\[Lambda]Lik + 2*yk^2*\[Lambda]Lik - 2*yS^2*\[Lambda]Lik + 
       4*xt*yk*\[Lambda]Rik) + xs^2*(xt*(1 - 2*yk^2 - 2*yS^2)*\[Lambda]Lik + 
       2*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4)*\[Lambda]Lik + 
       4*xt*yk*(1 + 2*yk^2 - 2*yS^2)*\[Lambda]Rik - 
       2*xt^2*(\[Lambda]Lik + 4*yk*\[Lambda]Rik)) + 
     xs*(-2*(yk^2 - yS^2)^2*\[Lambda]Lik - 4*xt^3*yk*\[Lambda]Rik + 
       4*xt*(yk^2 - yS^2)*(yk^2*\[Lambda]Lik - yS^2*\[Lambda]Lik - 
         2*yk*\[Lambda]Rik) + xt^2*((1 - 2*yk^2 - 2*yS^2)*\[Lambda]Lik + 
         4*yk*(1 + 2*yk^2 - 2*yS^2)*\[Lambda]Rik)))*CC[6])/(2*xs^3*xt^3) + 
  (mi^2*(Qj - Qk)*\[Lambda]Ljk*
    (-(Qk*(-1 + xs)*xs*(xs^4*\[Lambda]Lik + 2*(-1 + xt)*xt^2*(yk^2 - yS^2)*
         \[Lambda]Lik + xs^3*((-2 + xt - 2*yk^2 + 2*yS^2)*\[Lambda]Lik - 
          4*xt*yk*\[Lambda]Rik) + xs^2*(\[Lambda]Lik - 4*(-1 + xt)*yk^2*
           \[Lambda]Lik + 4*(-1 + xt)*yS^2*\[Lambda]Lik - 
          4*(-2 + xt)*xt*yk*\[Lambda]Rik) - xs*(2*xt^3*\[Lambda]Lik + 
          2*(yk^2 - yS^2)*\[Lambda]Lik - 2*xt^2*(\[Lambda]Lik + 
            2*yk*\[Lambda]Rik) + xt*(\[Lambda]Lik - 4*yk^2*\[Lambda]Lik + 
            4*yS^2*\[Lambda]Lik + 4*yk*\[Lambda]Rik)))) + 
     Qj*(-1 + xs + xt)*(xs^5*\[Lambda]Lik - (-1 + xt)*xt*(yk^2 - yS^2)^2*
        \[Lambda]Lik - 2*xs^4*((1 + yk^2 - yS^2)*\[Lambda]Lik + 
         2*xt*yk*\[Lambda]Rik) + xs*((1 - 3*xt + xt^2)*yk^4*\[Lambda]Lik - 
         2*(1 - 3*xt + xt^2)*yk^2*yS^2*\[Lambda]Lik + 
         yS^2*(4*xt^3 + yS^2 - 3*xt*yS^2 + xt^2*(-2 + yS^2))*\[Lambda]Lik - 
         4*(-1 + xt)*xt*yk^3*\[Lambda]Rik + 4*(-1 + xt)*xt*yk*yS^2*
          \[Lambda]Rik) + xs^3*((1 + xt - xt^2 - 2*xt*yk^2 + yk^4 - 4*yS^2 + 
           yS^4 - 2*yk^2*(-2 + yS^2))*\[Lambda]Lik - 
         4*xt*yk*(-2 + xt - yk^2 + yS^2)*\[Lambda]Rik) + 
       xs^2*(-2*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4)*\[Lambda]Lik + 
         xt*((-1 + 2*yk^4 + 2*yS^4 + yk^2*(2 - 4*yS^2))*\[Lambda]Lik - 
           4*yk*(1 + 2*yk^2 - 2*yS^2)*\[Lambda]Rik) + 
         xt^2*(\[Lambda]Lik - 2*yS^2*\[Lambda]Lik + 4*yk*(1 + yk^2 - yS^2)*
            \[Lambda]Rik))))*CC[7])/(2*xs^3*xt^3*(-1 + xs + xt)^2) - 
  (mi^2*Qk*\[Lambda]Ljk*(Qk*(-1 + xs)*xs*(xs^4*\[Lambda]Lik - 
       2*(-1 + xt)*xt^2*(yk^2 - yS^2)*\[Lambda]Lik + 
       xs^3*((-2 + xt + 2*yk^2 - 2*yS^2)*\[Lambda]Lik + 
         4*xt*yk*\[Lambda]Rik) + xs^2*((1 - 2*xt^2 - 4*yk^2 + 4*xt*yk^2 + 
           4*yS^2 - 4*xt*yS^2)*\[Lambda]Lik + 4*(-2 + xt)*xt*yk*
          \[Lambda]Rik) - xs*(4*xt^3*\[Lambda]Lik + 2*(-yk^2 + yS^2)*
          \[Lambda]Lik + xt*(\[Lambda]Lik + 4*yk^2*\[Lambda]Lik - 
           4*yS^2*\[Lambda]Lik - 4*yk*\[Lambda]Rik) - 
         4*xt^2*(\[Lambda]Lik - yk*\[Lambda]Rik))) + 
     Qj*(-((-1 + xt)^2*xt*(yk^2 - yS^2)^2*\[Lambda]Lik) + 
       xs^4*((yk^2 - yS^2)^2*\[Lambda]Lik + 4*xt*yk^3*\[Lambda]Rik + 
         xt^2*(\[Lambda]Lik - 4*yk*\[Lambda]Rik) - 
         2*xt*yS^2*(\[Lambda]Lik + 2*yk*\[Lambda]Rik)) + 
       xs^3*(-3*(yk^2 - yS^2)^2*\[Lambda]Lik + 3*xt^3*(\[Lambda]Lik - 
           4*yk*\[Lambda]Rik) + xt*(3*yk^4*\[Lambda]Lik + 3*yS^2*(2 + yS^2)*
            \[Lambda]Lik - 2*yk^2*(\[Lambda]Lik + 3*yS^2*\[Lambda]Lik) - 
           12*yk^3*\[Lambda]Rik + 12*yk*yS^2*\[Lambda]Rik) + 
         xt^2*((-3 + 2*yk^2 - 6*yS^2)*\[Lambda]Lik + 
           4*yk*(3 + 2*yk^2 - 2*yS^2)*\[Lambda]Rik)) + 
       xs^2*(3*(yk^2 - yS^2)^2*\[Lambda]Lik - 8*xt^4*yk*\[Lambda]Rik + 
         xt*(-7*yk^4*\[Lambda]Lik + 2*yk^2*(2 + 7*yS^2)*\[Lambda]Lik - 
           yS^2*(6 + 7*yS^2)*\[Lambda]Lik + 12*yk^3*\[Lambda]Rik - 
           12*yk*yS^2*\[Lambda]Rik) + 
         xt^2*((2 + 3*yk^4 + 12*yS^2 + 3*yS^4 - 2*yk^2*(5 + 3*yS^2))*
            \[Lambda]Lik - 4*yk*(3 + 4*yk^2 - 4*yS^2)*\[Lambda]Rik) + 
         xt^3*((-3 + 6*yk^2 - 4*yS^2)*\[Lambda]Lik + 4*yk*(5 + yk^2 - yS^2)*
            \[Lambda]Rik)) + xs*(-1 + xt)*((yk^2 - yS^2)^2*\[Lambda]Lik + 
         4*xt^3*yk*(yk*\[Lambda]Lik + \[Lambda]Rik) - 2*xt*(yk^2 - yS^2)*
          ((-1 + 2*yk^2 - 2*yS^2)*\[Lambda]Lik - 2*yk*\[Lambda]Rik) + 
         xt^2*(yk^4*\[Lambda]Lik - 2*yk^2*(3 + yS^2)*\[Lambda]Lik + 
           yS^2*(4 + yS^2)*\[Lambda]Lik - 4*yk^3*\[Lambda]Rik + 
           4*yk*(-1 + yS^2)*\[Lambda]Rik))))*CC[8])/
   (2*xs^3*xt^3*(-1 + xs + xt)^2) - 
  (mi^2*(Qj - Qk)*\[Lambda]Ljk*
    (Qk*(-1 + xt)*xt*((-1 + xt)^2*xt*(xt - 2*yk^2 + 2*yS^2)*\[Lambda]Lik + 
       2*xs^2*((yk^2 - yS^2)*\[Lambda]Lik + xt^2*(\[Lambda]Lik - 
           2*yk*\[Lambda]Rik) + xt*(-2*yk^2*\[Lambda]Lik + 
           2*yS^2*\[Lambda]Lik + 2*yk*\[Lambda]Rik)) + 
       xs*(-1 + xt)*(2*(yk^2 - yS^2)*\[Lambda]Lik + 
         xt^2*(\[Lambda]Lik - 4*yk*\[Lambda]Rik) + 
         xt*(\[Lambda]Lik - 6*yk^2*\[Lambda]Lik + 6*yS^2*\[Lambda]Lik + 
           4*yk*\[Lambda]Rik))) + Qj*(-1 + xs + xt)*
      (-(xt*(xt^2 + yk^2 - yS^2 + xt*(-1 - yk^2 + yS^2))^2*\[Lambda]Lik) + 
       xs*(-1 + xt)*((yk^2 - yS^2)^2*\[Lambda]Lik + 4*xt^3*yk*\[Lambda]Rik - 
         2*xt*(yk^2 - yS^2)*(yk^2*\[Lambda]Lik - yS^2*\[Lambda]Lik - 
           2*yk*\[Lambda]Rik) + xt^2*((-1 + 2*yk^2)*\[Lambda]Lik - 
           4*yk*(1 + yk^2 - yS^2)*\[Lambda]Rik)) + 
       xs^2*((yk^2 - yS^2)^2*\[Lambda]Lik + xt^3*(\[Lambda]Lik + 
           4*yk*\[Lambda]Rik) - xt*(yk^4*\[Lambda]Lik - 
           2*yk^2*yS^2*\[Lambda]Lik + yS^2*(2 + yS^2)*\[Lambda]Lik - 
           4*yk^3*\[Lambda]Rik + 4*yk*yS^2*\[Lambda]Rik) + 
         xt^2*((-1 + 6*yS^2)*\[Lambda]Lik - 4*yk*(1 + yk^2 - yS^2)*
            \[Lambda]Rik))))*CC[9])/(2*xs^3*xt^3*(-1 + xs + xt)^2) - 
  (mi^2*Qk*\[Lambda]Ljk*(Qk*(-1 + xt)*xt*
      ((-1 + xt)^2*xt*(xt + 2*yk^2 - 2*yS^2)*\[Lambda]Lik + 
       2*xs^2*(xt*(-1 + 2*yk^2 - 2*yS^2)*\[Lambda]Lik + 
         (-yk^2 + yS^2)*\[Lambda]Lik - 2*xt*yk*\[Lambda]Rik + 
         2*xt^2*(\[Lambda]Lik + yk*\[Lambda]Rik)) + 
       xs*(-1 + xt)*(xt*(-1 + 6*yk^2 - 6*yS^2)*\[Lambda]Lik + 
         2*(-yk^2 + yS^2)*\[Lambda]Lik - 4*xt*yk*\[Lambda]Rik + 
         xt^2*(3*\[Lambda]Lik + 4*yk*\[Lambda]Rik))) + 
     Qj*((-1 + xt)^3*xt*(yk^2 - yS^2)^2*\[Lambda]Lik + 
       xs*(-1 + xt)^2*(-((yk^2 - yS^2)^2*\[Lambda]Lik) + 
         xt*(yk^2 - yS^2)*(3*yk^2*\[Lambda]Lik - 3*yS^2*\[Lambda]Lik - 
           4*yk*\[Lambda]Rik) + 2*xt^2*yk*(-(yk*\[Lambda]Lik) + 
           2*yk^2*\[Lambda]Rik - 2*yS^2*\[Lambda]Rik)) + 
       xs^3*(-((yk^2 - yS^2)^2*\[Lambda]Lik) + 
         xt^3*(\[Lambda]Lik + 4*yk*\[Lambda]Rik) + 
         xt*(yk^4*\[Lambda]Lik - 2*yk^2*yS^2*\[Lambda]Lik + 
           yS^2*(2 + yS^2)*\[Lambda]Lik - 4*yk^3*\[Lambda]Rik + 
           4*yk*yS^2*\[Lambda]Rik) - xt^2*((1 + 4*yk^2 + 2*yS^2)*
            \[Lambda]Lik + 4*yk*(1 - yk^2 + yS^2)*\[Lambda]Rik)) - 
       xs^2*(-1 + xt)*(2*(yk^2 - yS^2)^2*\[Lambda]Lik + 
         xt^3*(\[Lambda]Lik - 4*yk*\[Lambda]Rik) - 
         xt*(3*yk^4*\[Lambda]Lik - 6*yk^2*yS^2*\[Lambda]Lik + 
           yS^2*(2 + 3*yS^2)*\[Lambda]Lik - 8*yk^3*\[Lambda]Rik + 
           8*yk*yS^2*\[Lambda]Rik) + xt^2*(6*yk^2*\[Lambda]Lik + 
           2*yS^2*\[Lambda]Lik - 8*yk^3*\[Lambda]Rik + 
           4*yk*(\[Lambda]Rik + 2*yS^2*\[Lambda]Rik)))))*CC[10])/
   (2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*Qk^2*\[Lambda]Ljk*(xs^5*\[Lambda]Lik + xt^4*(xt + 2*yk^2 - 2*yS^2)*
      \[Lambda]Lik + 4*xs^3*xt*(yk^2*\[Lambda]Lik - yS^2*\[Lambda]Lik + 
       xt*yk*\[Lambda]Rik) + xs*xt^3*(3*xt*\[Lambda]Lik + 
       4*yk^2*\[Lambda]Lik - 4*yS^2*\[Lambda]Lik + 4*xt*yk*\[Lambda]Rik) + 
     2*xs^2*xt^2*(2*(-yk^2 + yS^2)*\[Lambda]Lik + 
       xt*(\[Lambda]Lik + 2*yk*\[Lambda]Rik)) + 
     xs^4*(2*(yk^2 - yS^2)*\[Lambda]Lik + 
       xt*(\[Lambda]Lik + 4*yk*\[Lambda]Rik)))*CC[11])/
   (2*xs^3*xt^3*(xs + xt)) - (mi^2*(Qj - Qk)^2*\[Lambda]Ljk*
    (xs^5*\[Lambda]Lik + xt^4*(xt - 2*yk^2 + 2*yS^2)*\[Lambda]Lik + 
     xs^4*(2*(-yk^2 + yS^2)*\[Lambda]Lik + 
       xt*(\[Lambda]Lik - 4*yk*\[Lambda]Rik)) + 
     xs*xt^3*(4*(-yk^2 + yS^2)*\[Lambda]Lik + 
       xt*(\[Lambda]Lik - 4*yk*\[Lambda]Rik)) - 
     2*xs^3*xt*(2*(yk^2 - yS^2)*\[Lambda]Lik + 
       xt*(\[Lambda]Lik + 2*yk*\[Lambda]Rik)) - 
     2*xs^2*xt^2*(2*(-yk^2 + yS^2)*\[Lambda]Lik + 
       xt*(\[Lambda]Lik + 2*yk*\[Lambda]Rik)))*CC[12])/
   (2*xs^3*xt^3*(xs + xt)) + 
  (mi^4*(Qj - Qk)*Qk*(xs*(xt - (yk - yS)^2) - (-1 + xt)*(yk - yS)^2)*
    (-((-1 + xt)*(yk + yS)^2) + xs*(xt - (yk + yS)^2))*\[Lambda]Ljk*
    ((-1 + xt)*xt*(yk^2 - yS^2)*\[Lambda]Lik + 
     xs^2*((yk^2 - yS^2)*\[Lambda]Lik + 
       xt*(\[Lambda]Lik + 4*yk*\[Lambda]Rik)) - 
     xs*((yk^2 - yS^2)*\[Lambda]Lik + xt^2*(\[Lambda]Lik - 
         4*yk*\[Lambda]Rik) + xt*(-2*yk^2*\[Lambda]Lik + 
         2*yS^2*\[Lambda]Lik + 4*yk*\[Lambda]Rik)))*DD[1])/
   (2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^4*(Qj - Qk)*Qk*((-1 + xt)^2*(yk^2 - yS^2)^2 + 
     xs^2*(xt^2 + (yk^2 - yS^2)^2 - 2*xt*(yk^2 + yS^2)) - 
     2*xs*((yk^2 - yS^2)^2 + xt^2*(yk^2 + yS^2) - 
       xt*(yk^2 + yk^4 + yS^2 - 2*yk^2*yS^2 + yS^4)))*\[Lambda]Ljk*
    ((-1 + xt)*xt*(yk^2 - yS^2)*\[Lambda]Lik + 
     xs^2*((yk^2 - yS^2)*\[Lambda]Lik + 
       xt*(\[Lambda]Lik + 4*yk*\[Lambda]Rik)) + 
     xs*(2*xt*(-1 + yk^2 - yS^2)*\[Lambda]Lik + (-yk^2 + yS^2)*\[Lambda]Lik - 
       4*xt*yk*\[Lambda]Rik + xt^2*(3*\[Lambda]Lik + 4*yk*\[Lambda]Rik)))*
    DD[2])/(2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^4*(Qj - Qk)^2*(xs^3 + (-1 + xt)*(yk^2 - yS^2)^2 + 
     xs^2*(-1 + xt - 2*yk^2 + 2*yS^2) + xs*(yk^4 + yS^2*(-2 - 2*xt + yS^2) - 
       2*yk^2*(-1 + xt + yS^2)))*\[Lambda]Ljk*(xs^2*\[Lambda]Lik + 
     xt*(-yk^2 + yS^2)*\[Lambda]Lik - xs*((yk^2 - yS^2)*\[Lambda]Lik + 
       xt*(\[Lambda]Lik + 4*yk*\[Lambda]Rik)))*DD[3])/(2*xs^3*xt^3) - 
  (mi^4*Qk^2*\[Lambda]Ljk*(xs^5*\[Lambda]Lik + (-1 + xt)*xt*(yk^2 - yS^2)^3*
      \[Lambda]Lik + xs^4*((-1 + 3*yk^2 - 3*yS^2)*\[Lambda]Lik + 
       4*xt*yk*\[Lambda]Rik) + xs*(yk^2 - yS^2)*
      (-((yk^2 - yS^2)^2*\[Lambda]Lik) + xt*(yk^2 - yS^2)*
        ((-3 + 2*yk^2 - 2*yS^2)*\[Lambda]Lik - 4*yk*\[Lambda]Rik) + 
       xt^2*(-(yk^2*\[Lambda]Lik) - 3*yS^2*\[Lambda]Lik + 
         4*yk^3*\[Lambda]Rik - 4*yk*yS^2*\[Lambda]Rik)) - 
     xs^2*(-((-3 + yk^2 - yS^2)*(yk^2 - yS^2)^2*\[Lambda]Lik) - 
       xt*(yk^2 - yS^2)*((-1 + 2*yk^2 - 6*yS^2)*\[Lambda]Lik + 
         4*yk*(-2 + yk^2 - yS^2)*\[Lambda]Rik) + 
       xt^2*(3*yk^2*\[Lambda]Lik + yS^2*\[Lambda]Lik + 8*yk^3*\[Lambda]Rik + 
         4*yk*(-1 + 2*yS^2)*\[Lambda]Rik)) - 
     xs^3*(-3*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2))*\[Lambda]Lik + 
       xt^2*(\[Lambda]Lik + 4*yk*\[Lambda]Rik) + 
       xt*((-1 + 4*yS^2)*\[Lambda]Lik + 4*yk*(1 - 2*yk^2 + 2*yS^2)*
          \[Lambda]Rik)))*DD[4])/(2*xs^3*xt^3) + 
  (mi^4*(Qj - Qk)^2*(xt^3 + (-1 + xs)*(yk^2 - yS^2)^2 + 
     xt^2*(-1 + xs - 2*yk^2 + 2*yS^2) + xt*(yk^4 + yS^2*(-2 - 2*xs + yS^2) - 
       2*yk^2*(-1 + xs + yS^2)))*\[Lambda]Ljk*
    (xt*(xt - yk^2 + yS^2)*\[Lambda]Lik - 
     xs*((yk^2 - yS^2)*\[Lambda]Lik + xt*(\[Lambda]Lik + 4*yk*\[Lambda]Rik)))*
    DD[5])/(2*xs^3*xt^3) - 
  (mi^4*Qk^2*(xt^3 + xt^2*(-1 + xs + 2*yk^2 - 2*yS^2) + 
     (-1 + xs)*(yk^2 - yS^2)^2 + xt*(yk^4 + yS^2*(2 - 2*xs + yS^2) - 
       2*yk^2*(1 + xs + yS^2)))*\[Lambda]Ljk*(xs*(yk^2 - yS^2)*\[Lambda]Lik + 
     xt*(xt + yk^2 - yS^2)*\[Lambda]Lik + 
     xs*xt*(\[Lambda]Lik + 4*yk*\[Lambda]Rik))*DD[6])/(2*xs^3*xt^3), 
 (2*Qj*(2*Qk*(-1 + 2*xs)*\[Lambda]Rik + Qj*(1 - yk^2 + yS^2)*\[Lambda]Rik - 
     2*Qj*xs*(2*yk*\[Lambda]Lik + \[Lambda]Rik + yk^2*\[Lambda]Rik - 
       yS^2*\[Lambda]Rik))*\[Lambda]Rjk)/((-1 + xs)*xs^2*xt) - 
  (2*Qj^2*yk^2*((yk^2 - yS^2)*\[Lambda]Rik + 
     2*xs*(2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - yS^2*\[Lambda]Rik))*
    \[Lambda]Rjk*BB[1])/((-1 + xs)*xs^2*xt*(yk^2 - yS^2)) + 
  (2*Qj^2*yS^2*((yk^2 - yS^2)*\[Lambda]Rik + 
     2*xs*(2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - yS^2*\[Lambda]Rik))*
    \[Lambda]Rjk*BB[2])/((-1 + xs)*xs^2*xt*(yk^2 - yS^2)) + 
  (2*Qj*(-2*Qk*(xs^3 + 2*xs^2*(-2 + xt) + 2*(-1 + xt) - 
       xs*(-5 + 3*xt + xt^2))*\[Lambda]Rik + Qj*(-1 + xs + xt)*
      (2*(1 - 3*xs + 2*xs^2)*(-1 + xt)*yk*\[Lambda]Lik + 
       (xs*(7 - 3*xt) + 2*(-2 + xt) + 2*xs^2*(-2 + xt))*yk^2*\[Lambda]Rik - 
       (-2 + 2*(-2 + xt)*yS^2 + xs*(1 + 7*yS^2 - 3*xt*(-1 + yS^2)) + 
         2*xs^2*(-2*yS^2 + xt*(-1 + yS^2)))*\[Lambda]Rik))*\[Lambda]Rjk*
    BB[3])/((-1 + xs)^2*xs^2*(-1 + xt)*xt*(-1 + xs + xt)) - 
  (2*Qj*(2*Qk*xs*(-2 + 4*xs - 2*xs^2 + xt)*\[Lambda]Rik + 
     Qj*(-1 + xs + xt)*(-2*(-1 + xs)*yk*\[Lambda]Lik + 
       (3 - 2*xs)*yk^2*\[Lambda]Rik + (2*xs^2 - 3*yS^2 + xs*(-3 + 2*yS^2))*
        \[Lambda]Rik))*\[Lambda]Rjk*BB[4])/((-1 + xs)^2*xs^2*xt*
    (-1 + xs + xt)) - (4*Qj*(-(Qk*(xs + 2*(-1 + xt)*xt)) + 
     Qj*(-1 + xs + xt)*(xt - yk^2 + yS^2))*\[Lambda]Rik*\[Lambda]Rjk*BB[5])/
   (xs^2*(-1 + xt)*xt*(-1 + xs + xt)) - 
  (4*Qk^2*\[Lambda]Rik*\[Lambda]Rjk*BB[6])/(xs^2*xt) + 
  (4*(Qj - Qk)^2*\[Lambda]Rik*\[Lambda]Rjk*BB[7])/(xs^2*xt) - 
  (2*mi^2*Qk*(Qk*(-1 + xs)*xs*(xs^3 + xs*(-1 + xt)*(-1 + 4*yk^2 - 4*yS^2) + 
       2*(-1 + xt)^2*(yk^2 - yS^2) + 2*xs^2*(-1 + xt + yk^2 - yS^2))*
      \[Lambda]Rik + Qj*(-((-1 + xt)^2*(yk^2 - yS^2)^2*\[Lambda]Rik) + 
       xs^3*(xt^2*\[Lambda]Rik + (yk^2 - yS^2)^2*\[Lambda]Rik - 
         xt*(3*yk*\[Lambda]Lik + \[Lambda]Rik + 2*yk^2*\[Lambda]Rik)) + 
       xs^2*(-3*(yk^2 - yS^2)^2*\[Lambda]Rik - xt^2*(7*yk*\[Lambda]Lik + 
           \[Lambda]Rik + 6*yk^2*\[Lambda]Rik) + 
         xt*(6*yk*\[Lambda]Lik + \[Lambda]Rik + 2*yk^4*\[Lambda]Rik + 
           2*yS^4*\[Lambda]Rik - 4*yk^2*(-1 + yS^2)*\[Lambda]Rik)) - 
       xs*(-1 + xt)*(3*(yk^2 - yS^2)^2*\[Lambda]Rik + 
         4*xt^2*yk*(\[Lambda]Lik + yk*\[Lambda]Rik) - 
         xt*(3*yk*\[Lambda]Lik + yk^4*\[Lambda]Rik + yS^4*\[Lambda]Rik - 
           2*yk^2*(-1 + yS^2)*\[Lambda]Rik))))*\[Lambda]Rjk*CC[1])/
   ((-1 + xs)*xs^2*xt^2*(-1 + xs + xt)^2) + 
  (2*mi^2*(Qj - Qk)*(-(Qk*(-1 + xs)*xs*(xs^3 - 2*(-1 + xt)^2*(yk^2 - yS^2) + 
        2*xs^2*(-1 + xt - yk^2 + yS^2) + xs*(1 + 2*xt^2 + 4*yk^2 - 4*yS^2 + 
          xt*(-2 - 4*yk^2 + 4*yS^2)))) + 
     Qj*(xs^5 - (-1 + xt)^2*(yk^2 - yS^2)^2 + 
       xs^4*(-3 + 2*xt - 2*yk^2 + 2*yS^2) + xs^3*(3 + xt^2 + yk^4 - 6*yS^2 + 
         yS^4 - 2*yk^2*(-3 + yS^2) + xt*(-4 - 4*yk^2 + 2*yS^2)) + 
       xs*((3 - 4*xt + xt^2)*yk^4 - (-1 + xt)*yS^2*(-2 + 4*xt^2 + 3*yS^2 - 
           xt*yS^2) - 2*(-1 + xt)*yk^2*(1 - 3*yS^2 + xt*(-1 + yS^2))) - 
       xs^2*(1 + 3*yk^4 - 6*yS^2 + 3*yS^4 - 6*yk^2*(-1 + yS^2) + 
         xt^2*(1 + 2*yk^2 + 4*yS^2) - 2*xt*(yk^4 - 2*yk^2*(-2 + yS^2) + 
           (-1 + yS^2)^2))))*\[Lambda]Rik*\[Lambda]Rjk*CC[2])/
   ((-1 + xs)*xs^2*xt^2*(-1 + xs + xt)^2) - 
  (2*mi^2*Qk*(Qj*((-1 + xt)^2*(yk^2 - yS^2)^2 + 
       xs^2*(xt^2 - 2*xt*yS^2 + (yk^2 - yS^2)^2) + 
       2*xs*((-1 + xt)*yk^4 - 2*(-1 + xt)*yk^2*yS^2 + 
         yS^2*(xt - xt^2 - yS^2 + xt*yS^2)))*\[Lambda]Rik + 
     Qk*xt*(-(xs^2*yk*\[Lambda]Lik) + (-1 + xt)^2*(xt + 2*yk^2 - 2*yS^2)*
        \[Lambda]Rik + xs*(-1 + xt)*(-(yk*\[Lambda]Lik) + 
         2*yk^2*\[Lambda]Rik + 2*(xt - yS^2)*\[Lambda]Rik)))*\[Lambda]Rjk*
    CC[3])/(xs^3*xt*(-1 + xs + xt)^2) + 
  (2*mi^2*(Qj - Qk)*(Qj*(xt^4 + (-1 + xs)^2*(yk^2 - yS^2)^2 + 
       2*xt^3*(-1 + xs - yk^2 + yS^2) + xt^2*(1 + xs^2 + yk^4 - 4*yS^2 + 
         yS^4 - 2*yk^2*(-2 + yS^2) + xs*(-2 - 4*yk^2 + 2*yS^2)) + 
       2*xt*((-1 + xs)*yk^4 + (-1 + xs)*yS^2*(-1 + yS^2) - 
         yk^2*(1 + xs^2 - 2*yS^2 + 2*xs*(-1 + yS^2))))*\[Lambda]Rik - 
     Qk*xt*((-1 + xt)^2*(xt - 2*yk^2 + 2*yS^2)*\[Lambda]Rik + 
       xs^2*(yk*\[Lambda]Lik + (-1 + 2*xt)*\[Lambda]Rik) + 
       xs*(-1 + xt)*(yk*\[Lambda]Lik - 2*yk^2*\[Lambda]Rik + 
         2*(xt + yS^2)*\[Lambda]Rik)))*\[Lambda]Rjk*CC[4])/
   (xs^3*xt*(-1 + xs + xt)^2) + 
  (2*mi^2*Qk^2*(xs^3*\[Lambda]Rik + xs^2*(-1 + xt + 2*yk^2 - 2*yS^2)*
      \[Lambda]Rik + (xt^3 + xt^2*(-1 + 2*yk^2 - 2*yS^2) - 
       2*(yk^2 - yS^2)^2 + 2*xt*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2)))*
      \[Lambda]Rik + xs*(xt^2*\[Lambda]Rik + 
       2*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2))*\[Lambda]Rik - 
       2*xt*(yk*\[Lambda]Lik + 2*yS^2*\[Lambda]Rik)))*\[Lambda]Rjk*CC[5])/
   (xs^3*xt^2) - (2*mi^2*(Qj - Qk)^2*(xs^3 + xt^3 - 2*(yk^2 - yS^2)^2 + 
     xt^2*(-1 - 2*yk^2 + 2*yS^2) + xs^2*(-1 + xt - 2*yk^2 + 2*yS^2) + 
     2*xt*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4) + 
     xs*(xt^2 - 4*xt*yk^2 + 2*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4)))*
    \[Lambda]Rik*\[Lambda]Rjk*CC[6])/(xs^3*xt^2) + 
  (2*mi^2*(Qj - Qk)*(Qj*(xs^6 + (-1 + xt)^2*(yk^2 - yS^2)^2 + 
       2*xs^5*(-2 + xt - yk^2 + yS^2) + xs^4*(6 + xt^2 + yk^4 - 8*yS^2 + 
         yS^4 - 2*yk^2*(-4 + yS^2) + xt*(-6 - 4*yk^2 + 2*yS^2)) - 
       2*xs^3*(xt^2*(1 + yk^2) + 2*(1 + yk^4 - 3*yS^2 + yS^4 + 
           yk^2*(3 - 2*yS^2)) - xt*(3 + yk^4 - 2*yS^2 + yS^4 - 
           2*yk^2*(-3 + yS^2))) + xs^2*(1 + 6*yk^4 - 8*yS^2 + 6*yS^4 + 
         yk^2*(8 - 12*yS^2) + xt^2*(1 + yk^4 + 4*yS^2 + yS^4 - 
           2*yk^2*(-2 + yS^2)) - 2*xt*(1 + 3*yk^4 - yS^2 + 3*yS^4 - 
           6*yk^2*(-1 + yS^2))) + 2*xs*(-((2 - 3*xt + xt^2)*yk^4) + 
         (-1 + xt)*yS^2*(-1 + xt^2 + 2*yS^2 - xt*(1 + yS^2)) + 
         (-1 + xt)*yk^2*(1 - 4*yS^2 + xt*(-1 + 2*yS^2))))*\[Lambda]Rik - 
     Qk*(-1 + xs)^2*xs*((-1 + xs)^2*(xs - 2*yk^2 + 2*yS^2)*\[Lambda]Rik + 
       xt^2*(yk*\[Lambda]Lik + 2*xs*\[Lambda]Rik) + 
       xt*((-1 + xs)*yk*\[Lambda]Lik - 2*(-1 + xs)*yk^2*\[Lambda]Rik + 
         (-3*xs + 2*xs^2 - 2*yS^2 + 2*xs*yS^2)*\[Lambda]Rik)))*\[Lambda]Rjk*
    CC[7])/((-1 + xs)*xs^3*xt^2*(-1 + xs + xt)^2) - 
  (2*mi^2*Qk*(Qk*(-1 + xs)^2*xs*(-(xt^2*yk*\[Lambda]Lik) + 
       (-1 + xs)^2*(xs + 2*yk^2 - 2*yS^2)*\[Lambda]Rik + 
       (-1 + xs)*xt*(-(yk*\[Lambda]Lik) + 2*yk^2*\[Lambda]Rik + 
         2*(xs - yS^2)*\[Lambda]Rik)) + 
     Qj*((-1 + xt)^2*(yk^2 - yS^2)^2*\[Lambda]Rik + 
       xs^4*(xt^2*\[Lambda]Rik + (yk^2 - yS^2)^2*\[Lambda]Rik - 
         2*xt*(yk*\[Lambda]Lik + yS^2*\[Lambda]Rik)) + 
       xs^2*(-2*xt^3*yk*\[Lambda]Lik + 6*(yk^2 - yS^2)^2*\[Lambda]Rik - 
         2*xt*(3*yk*\[Lambda]Lik + 3*yk^4*\[Lambda]Rik + yk^2*(2 - 6*yS^2)*
            \[Lambda]Rik + 3*yS^2*(1 + yS^2)*\[Lambda]Rik) + 
         xt^2*(8*yk*\[Lambda]Lik + yk^4*\[Lambda]Rik - 2*yk^2*(-2 + yS^2)*
            \[Lambda]Rik + (1 + 4*yS^2 + yS^4)*\[Lambda]Rik)) + 
       2*xs*(-1 + xt)*(2*(yk^2 - yS^2)^2*\[Lambda]Rik + 
         xt^2*yk*(\[Lambda]Lik + yk*\[Lambda]Rik) - 
         xt*(yk*\[Lambda]Lik + yk^4*\[Lambda]Rik + yS^2*(1 + yS^2)*
            \[Lambda]Rik + yk^2*(\[Lambda]Rik - 2*yS^2*\[Lambda]Rik))) - 
       2*xs^3*(2*(yk^2 - yS^2)^2*\[Lambda]Rik + 
         xt^2*(2*yk*\[Lambda]Lik + \[Lambda]Rik + yS^2*\[Lambda]Rik) - 
         xt*(3*yk*\[Lambda]Lik + yk^4*\[Lambda]Rik + yS^2*(3 + yS^2)*
            \[Lambda]Rik + yk^2*(\[Lambda]Rik - 2*yS^2*\[Lambda]Rik)))))*
    \[Lambda]Rjk*CC[8])/((-1 + xs)*xs^3*xt^2*(-1 + xs + xt)^2) + 
  (2*mi^2*(Qj - Qk)*(-(Qk*(-1 + xt)*xt*(xt^3 - 2*(-1 + xs)^2*(yk^2 - yS^2) + 
        2*xt^2*(-1 + xs - yk^2 + yS^2) + xt*(1 + 2*xs^2 + 4*yk^2 - 4*yS^2 + 
          xs*(-2 - 4*yk^2 + 4*yS^2)))) + 
     Qj*(xt^5 - (-1 + xs)^2*(yk^2 - yS^2)^2 + 
       xt^4*(-3 + 2*xs - 2*yk^2 + 2*yS^2) + xt^3*(3 + xs^2 + yk^4 - 6*yS^2 + 
         yS^4 - 2*yk^2*(-3 + yS^2) + xs*(-4 - 4*yk^2 + 2*yS^2)) + 
       xt*((3 - 4*xs + xs^2)*yk^4 - (-1 + xs)*yS^2*(-2 + 4*xs^2 + 3*yS^2 - 
           xs*yS^2) - 2*(-1 + xs)*yk^2*(1 - 3*yS^2 + xs*(-1 + yS^2))) - 
       xt^2*(1 + 3*yk^4 - 6*yS^2 + 3*yS^4 - 6*yk^2*(-1 + yS^2) + 
         xs^2*(1 + 2*yk^2 + 4*yS^2) - 2*xs*(yk^4 - 2*yk^2*(-2 + yS^2) + 
           (-1 + yS^2)^2))))*\[Lambda]Rik*\[Lambda]Rjk*CC[9])/
   (xs^3*xt^2*(-1 + xs + xt)^2) + 
  (2*mi^2*Qk*(-(Qk*(-1 + xt)*xt*(xs^2*(1 + 2*yk^2 - 2*yS^2) + 
        (-1 + xt)^2*(xt + 2*yk^2 - 2*yS^2) + 2*xs*(xt^2 - 2*yk^2 + 2*yS^2 + 
          xt*(-1 + 2*yk^2 - 2*yS^2)))*\[Lambda]Rik) + 
     Qj*(-((-1 + xt)^3*(yk^2 - yS^2)^2*\[Lambda]Rik) + 
       2*xs^3*xt*yk*(\[Lambda]Lik + 2*yk*\[Lambda]Rik) + 
       xs*(-1 + xt)^2*(-2*(yk^2 - yS^2)^2*\[Lambda]Rik + 
         xt*yk*(\[Lambda]Lik + 2*yk*\[Lambda]Rik)) - 
       xs^2*(-1 + xt)*(xt^2*\[Lambda]Rik + (yk^2 - yS^2)^2*\[Lambda]Rik - 
         xt*(3*yk*\[Lambda]Lik + \[Lambda]Rik + 6*yk^2*\[Lambda]Rik))))*
    \[Lambda]Rjk*CC[10])/(xs^3*xt^2*(-1 + xs + xt)^2) + 
  (2*mi^2*Qk^2*(xs^2 + xt^2)*(xs + xt + 2*yk^2 - 2*yS^2)*\[Lambda]Rik*
    \[Lambda]Rjk*CC[11])/(xs^3*xt^2) - 
  (2*mi^2*(Qj - Qk)^2*(xs^2 + xt^2)*(xs + xt - 2*yk^2 + 2*yS^2)*\[Lambda]Rik*
    \[Lambda]Rjk*CC[12])/(xs^3*xt^2) + 
  (2*mi^4*(Qj - Qk)*Qk*((-1 + xt)^3*(yk^2 - yS^2)^3*\[Lambda]Rik - 
     xs*(-1 + xt)^2*(yk^2 - yS^2)*(-3*(yk^2 - yS^2)^2*\[Lambda]Rik + 
       xt*(yk*\[Lambda]Lik + 3*yk^2*\[Lambda]Rik + yS^2*\[Lambda]Rik)) + 
     xs^3*(xt^3*\[Lambda]Rik + (yk^2 - yS^2)^3*\[Lambda]Rik + 
       xt^2*(yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - (1 + yS^2)*\[Lambda]Rik) + 
       xt*(-(yk^3*\[Lambda]Lik) + yk*yS^2*\[Lambda]Lik - 
         3*yk^4*\[Lambda]Rik + yS^2*(1 + yS^2)*\[Lambda]Rik + 
         yk^2*(\[Lambda]Rik + 2*yS^2*\[Lambda]Rik))) + 
     xs^2*(-1 + xt)*(3*(yk^2 - yS^2)^3*\[Lambda]Rik + 
       xt^2*(yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - yS^2*\[Lambda]Rik) + 
       xt*(-2*yk^3*\[Lambda]Lik + 2*yk*yS^2*\[Lambda]Lik - 
         6*yk^4*\[Lambda]Rik + yS^2*(1 + 2*yS^2)*\[Lambda]Rik + 
         yk^2*(\[Lambda]Rik + 4*yS^2*\[Lambda]Rik))))*\[Lambda]Rjk*DD[1])/
   (xs^3*xt^2*(-1 + xs + xt)^2) + 
  (2*mi^4*(Qj - Qk)*Qk*(xs*(xt - (yk - yS)^2) - (-1 + xt)*(yk - yS)^2)*
    ((-1 + xt)*(yk^2 - yS^2) + xs*(xt + yk^2 - yS^2))*
    (-((-1 + xt)*(yk + yS)^2) + xs*(xt - (yk + yS)^2))*\[Lambda]Rik*
    \[Lambda]Rjk*DD[2])/(xs^3*xt^2*(-1 + xs + xt)^2) + 
  (2*mi^4*(Qj - Qk)^2*(xs - yk^2 + yS^2)*(xs^3 + (-1 + xt)*(yk^2 - yS^2)^2 + 
     xs^2*(-1 + xt - 2*yk^2 + 2*yS^2) + xs*(yk^4 + yS^2*(-2 - 2*xt + yS^2) - 
       2*yk^2*(-1 + xt + yS^2)))*\[Lambda]Rik*\[Lambda]Rjk*DD[3])/
   (xs^3*xt^2) - (2*mi^4*Qk^2*(xs^5*\[Lambda]Rik + 
     xs^4*(-2 + 2*xt + 3*yk^2 - 3*yS^2)*\[Lambda]Rik + 
     (-1 + xt)^2*(yk^2 - yS^2)^3*\[Lambda]Rik - xs*(-1 + xt)*(yk^2 - yS^2)*
      ((-2*yk^4 - yS^2*(3 + 2*yS^2) + yk^2*(3 + 4*yS^2))*\[Lambda]Rik + 
       xt*(yk*\[Lambda]Lik + yk^2*\[Lambda]Rik + 3*yS^2*\[Lambda]Rik)) + 
     xs^2*((yk^6 - 3*yk^4*(2 + yS^2) + 3*yk^2*(1 + 4*yS^2 + yS^4) - 
         yS^2*(3 + 6*yS^2 + yS^4))*\[Lambda]Rik - 
       xt^2*(3*yk*\[Lambda]Lik + 5*yk^2*\[Lambda]Rik + 3*yS^2*\[Lambda]Rik) + 
       xt*(-(yk^3*\[Lambda]Lik) + yk*(3 + yS^2)*\[Lambda]Lik + 
         2*yk^4*\[Lambda]Rik - 8*yk^2*yS^2*\[Lambda]Rik + 
         6*yS^2*(1 + yS^2)*\[Lambda]Rik)) + 
     xs^3*(xt^2*\[Lambda]Rik + (1 + 3*yk^4 + 6*yS^2 + 3*yS^4 - 
         6*yk^2*(1 + yS^2))*\[Lambda]Rik + xt*(-3*yk*\[Lambda]Lik + 
         2*yk^2*\[Lambda]Rik - 2*(\[Lambda]Rik + 3*yS^2*\[Lambda]Rik))))*
    \[Lambda]Rjk*DD[4])/(xs^3*xt^2*(-1 + xs + xt)) + 
  (2*mi^4*(Qj - Qk)^2*(xt - yk^2 + yS^2)*(xt^3 + (-1 + xs)*(yk^2 - yS^2)^2 + 
     xt^2*(-1 + xs - 2*yk^2 + 2*yS^2) + xt*(yk^4 + yS^2*(-2 - 2*xs + yS^2) - 
       2*yk^2*(-1 + xs + yS^2)))*\[Lambda]Rik*\[Lambda]Rjk*DD[5])/
   (xs^3*xt^2) - (2*mi^4*Qk^2*((-1 + xt)^2*(xt + yk^2 - yS^2)^3*
      \[Lambda]Rik + xs*(xt^2 - yk^2 + yS^2 + xt*(-1 + yk^2 - yS^2))*
      (2*xt^2*\[Lambda]Rik + 2*(yk^2 - yS^2)^2*\[Lambda]Rik - 
       xt*(yk*\[Lambda]Lik + 4*yS^2*\[Lambda]Rik)) + 
     xs^2*(xt^3*\[Lambda]Rik + (yk^2 - yS^2)^3*\[Lambda]Rik - 
       xt^2*(yk*\[Lambda]Lik + 5*yk^2*\[Lambda]Rik + 3*yS^2*\[Lambda]Rik) + 
       xt*(-(yk^3*\[Lambda]Lik) + yk*yS^2*\[Lambda]Lik - yk^4*\[Lambda]Rik + 
         3*yS^4*\[Lambda]Rik - 2*yk^2*(-1 + yS^2)*\[Lambda]Rik)))*
    \[Lambda]Rjk*DD[6])/(xs^3*xt^2*(-1 + xs + xt)), 
 (-2*Qj*(-2*Qk*xs*xt*(-2 + xs + xt) + Qj*(-((-1 + xt)*xt*(yk^2 - yS^2)) + 
       xs^2*(xt - yk^2 + 2*xt*yk^2 + yS^2 - 2*xt*yS^2) + 
       xs*(yk^2 - yS^2 + xt^2*(1 + 2*yk^2 - 2*yS^2) + 
         xt*(-2 - 4*yk^2 + 4*yS^2))))*\[Lambda]Lik*\[Lambda]Ljk)/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)) + 
  (2*Qj^2*(xs + xt - 2*xs*xt)*yk^2*\[Lambda]Lik*\[Lambda]Ljk*BB[1])/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2) + 
  (2*Qj^2*(-xt + xs*(-1 + 2*xt))*yS^2*\[Lambda]Lik*\[Lambda]Ljk*BB[2])/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2) + 
  (2*Qj*\[Lambda]Ljk*(-2*Qk*(2 - 2*xs + xs^2 - 2*xt + xt^2)*\[Lambda]Lik + 
     Qj*(-2 - 6*yk^2 + 6*yS^2 + xt*(4 + 8*yk^2 - 8*yS^2) + 
       xt^2*(-1 - 3*yk^2 + 3*yS^2) + xs^2*(-1 - 3*yk^2 + 3*yS^2 + 
         2*xt*(1 + yk^2 - yS^2)) + 2*xs*(2 + 4*yk^2 - 4*yS^2 - 
         4*xt*(1 + yk^2 - yS^2) + xt^2*(1 + yk^2 - yS^2)))*\[Lambda]Lik + 
     4*Qj*(-1 + xs)*(-1 + xt)*(-2 + xs + xt)*yk*\[Lambda]Rik)*BB[3])/
   ((-1 + xs)^2*xs*(-1 + xt)^2*xt*(-1 + xs + xt)) + 
  (2*Qj*\[Lambda]Ljk*(2*Qk*xs*\[Lambda]Lik + 
     Qj*((yk^2 - yS^2)*\[Lambda]Lik - 2*xs^2*(\[Lambda]Lik + 
         2*yk*\[Lambda]Rik) + xs*(\[Lambda]Lik + 4*yk*\[Lambda]Rik)))*BB[4])/
   ((-1 + xs)^2*xs^2*xt*(-1 + xs + xt)) + 
  (2*Qj*\[Lambda]Ljk*(2*Qk*xt*\[Lambda]Lik + 
     Qj*((yk^2 - yS^2)*\[Lambda]Lik - 2*xt^2*(\[Lambda]Lik + 
         2*yk*\[Lambda]Rik) + xt*(\[Lambda]Lik + 4*yk*\[Lambda]Rik)))*BB[5])/
   (xs*(-1 + xt)^2*xt^2*(-1 + xs + xt)) + 
  (2*mi^2*Qk*\[Lambda]Ljk*(Qk*xs*((-1 + xs)*(yk^2 - yS^2)*\[Lambda]Lik + 
       xt*(-1 + 2*xs + yk^2 - yS^2)*\[Lambda]Lik + 2*(-1 + xs)*xt*yk*
        \[Lambda]Rik + 2*xt^2*yk*\[Lambda]Rik) + 
     Qj*(-((-1 + xt)*xt*yk^2*\[Lambda]Lik) + 
       xs*(-(yS^2*\[Lambda]Lik) + xt*(1 - yk^2 + yS^2)*\[Lambda]Lik + 
         2*xt*yk*\[Lambda]Rik - 2*xt^2*yk*\[Lambda]Rik) + 
       xs^2*(yS^2*\[Lambda]Lik - 2*xt*(\[Lambda]Lik + yk*\[Lambda]Rik))))*
    CC[1])/(xs*xt^2*(-1 + xs + xt)^3) - 
  (2*mi^2*(Qj - Qk)*\[Lambda]Ljk*(Qj*(-xs + xs^2 + xt - xt^2)*yS^2*
      \[Lambda]Lik + Qk*xs*(xt*(-1 + yk^2 - yS^2)*\[Lambda]Lik + 
       (-1 + xs)*(yk^2 - yS^2)*\[Lambda]Lik + 2*(-1 + xs)*xt*yk*
        \[Lambda]Rik + 2*xt^2*(\[Lambda]Lik + yk*\[Lambda]Rik)))*CC[2])/
   (xs*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^2*Qk*\[Lambda]Ljk*(Qk*xt*((-1 + xt)*(yk^2 - yS^2)*\[Lambda]Lik + 
       xs*(-1 + 2*xt + yk^2 - yS^2)*\[Lambda]Lik + 2*xs^2*yk*\[Lambda]Rik + 
       2*xs*(-1 + xt)*yk*\[Lambda]Rik) + Qj*((-1 + xt)*xt*yS^2*\[Lambda]Lik - 
       xs^2*yk*(yk*\[Lambda]Lik + 2*xt*\[Lambda]Rik) + 
       xs*(yk^2*\[Lambda]Lik + xt*(1 - yk^2 + yS^2)*\[Lambda]Lik + 
         2*xt*yk*\[Lambda]Rik - 2*xt^2*(\[Lambda]Lik + yk*\[Lambda]Rik))))*
    CC[3])/(xs^2*xt*(-1 + xs + xt)^3) + 
  (2*mi^2*(Qj - Qk)*\[Lambda]Ljk*(Qj*(-xs + xs^2 + xt - xt^2)*yS^2*
      \[Lambda]Lik + Qk*xt*(-((-1 + xt)*(yk^2 - yS^2)*\[Lambda]Lik) + 
       xs*(1 - yk^2 + yS^2)*\[Lambda]Lik - 2*xs*(-1 + xt)*yk*\[Lambda]Rik - 
       2*xs^2*(\[Lambda]Lik + yk*\[Lambda]Rik)))*CC[4])/
   (xs^2*xt*(-1 + xs + xt)^3) - (4*mi^2*Qk^2*(xs + xt)*yk^2*\[Lambda]Lik*
    \[Lambda]Ljk*CC[5])/(xs^2*xt^2*(-1 + xs + xt)) + 
  (4*mi^2*(Qj - Qk)^2*(xs + xt)*yS^2*\[Lambda]Lik*\[Lambda]Ljk*CC[6])/
   (xs^2*xt^2*(-1 + xs + xt)) - (2*mi^2*(Qj - Qk)*\[Lambda]Ljk*
    (Qj*(xs^4 + (-1 + xt)*xt + xs^3*(-3 + 6*xt) + xs^2*(3 - 13*xt + 7*xt^2) + 
       xs*(-1 + 8*xt - 8*xt^2 + 2*xt^3))*yS^2*\[Lambda]Lik + 
     Qk*(-1 + xs)^2*xt*(xs*(-1 + yk^2 - yS^2)*\[Lambda]Lik + 
       (-1 + xt)*(yk^2 - yS^2)*\[Lambda]Lik + 2*xs*(-1 + xt)*yk*
        \[Lambda]Rik + 2*xs^2*(\[Lambda]Lik + yk*\[Lambda]Rik)))*CC[7])/
   ((-1 + xs)*xs^2*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^2*Qk*\[Lambda]Ljk*(Qk*(-1 + xs)^2*xt*
      ((-1 + xt)*(yk^2 - yS^2)*\[Lambda]Lik + xs*(-1 + 2*xt + yk^2 - yS^2)*
        \[Lambda]Lik + 2*xs^2*yk*\[Lambda]Rik + 2*xs*(-1 + xt)*yk*
        \[Lambda]Rik) + Qj*((-1 + xt)*xt*yS^2*\[Lambda]Lik + 
       xs^4*yk*(yk*\[Lambda]Lik - 2*xt*\[Lambda]Rik) + 
       xs^3*(-3*yk^2*\[Lambda]Lik + xt*(1 + 5*yk^2 + yS^2)*\[Lambda]Lik + 
         6*xt*yk*\[Lambda]Rik - 2*xt^2*(\[Lambda]Lik + yk*\[Lambda]Rik)) + 
       xs*(-(yk^2*\[Lambda]Lik) + 2*xt^3*yk^2*\[Lambda]Lik - 
         2*xt^2*((1 + 3*yk^2 + yS^2)*\[Lambda]Lik + yk*\[Lambda]Rik) + 
         xt*(\[Lambda]Lik + 5*yk^2*\[Lambda]Lik + 3*yS^2*\[Lambda]Lik + 
           2*yk*\[Lambda]Rik)) + xs^2*(3*yk^2*\[Lambda]Lik + 
         xt^2*((4 + 6*yk^2 + yS^2)*\[Lambda]Lik + 4*yk*\[Lambda]Rik) - 
         xt*((2 + 10*yk^2 + 3*yS^2)*\[Lambda]Lik + 6*yk*\[Lambda]Rik))))*
    CC[8])/((-1 + xs)*xs^2*xt^2*(-1 + xs + xt)^3) - 
  (2*mi^2*(Qj - Qk)*\[Lambda]Ljk*
    (Qj*(2*xs^3*xt + (-1 + xt)^3*xt + xs*(-1 + xt)^2*(-1 + 6*xt) + 
       xs^2*(1 - 8*xt + 7*xt^2))*yS^2*\[Lambda]Lik + 
     Qk*xs*(-1 + xt)^2*(xt*(-1 + yk^2 - yS^2)*\[Lambda]Lik + 
       (-1 + xs)*(yk^2 - yS^2)*\[Lambda]Lik + 2*(-1 + xs)*xt*yk*
        \[Lambda]Rik + 2*xt^2*(\[Lambda]Lik + yk*\[Lambda]Rik)))*CC[9])/
   (xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^2*Qk*\[Lambda]Ljk*(Qk*xs*(-1 + xt)^2*
      ((-1 + xs)*(yk^2 - yS^2)*\[Lambda]Lik + xt*(-1 + 2*xs + yk^2 - yS^2)*
        \[Lambda]Lik + 2*(-1 + xs)*xt*yk*\[Lambda]Rik + 
       2*xt^2*yk*\[Lambda]Rik) + Qj*(2*xs^3*xt*yk^2*\[Lambda]Lik + 
       (-1 + xt)^3*xt*yk^2*\[Lambda]Lik - xs*(-1 + xt)^2*
        (yS^2*\[Lambda]Lik + 2*xt^2*yk*\[Lambda]Rik - 
         xt*((1 + 5*yk^2 + yS^2)*\[Lambda]Lik + 2*yk*\[Lambda]Rik)) - 
       xs^2*(-1 + xt)*(yS^2*\[Lambda]Lik + 2*xt^2*(\[Lambda]Lik + 
           yk*\[Lambda]Rik) - xt*((2 + 6*yk^2 + yS^2)*\[Lambda]Lik + 
           2*yk*\[Lambda]Rik))))*CC[10])/(xs^2*(-1 + xt)*xt^2*
    (-1 + xs + xt)^3) + (2*mi^4*(Qj - Qk)*Qk*\[Lambda]Ljk*
    (-((-1 + xt)^2*xt*yk^2*(yk^2 - yS^2)*\[Lambda]Lik) - 
     xs*(-1 + xt)*(yS^2*(-yk^2 + yS^2)*\[Lambda]Lik + 
       xt^2*yk*(-(yk*\[Lambda]Lik) + 2*yk^2*\[Lambda]Rik + 
         2*yS^2*\[Lambda]Rik) + xt*(2*yk^4*\[Lambda]Lik - 
         yk^2*(1 + yS^2)*\[Lambda]Lik - yS^2*(1 + yS^2)*\[Lambda]Lik - 
         2*yk^3*\[Lambda]Rik - 2*yk*yS^2*\[Lambda]Rik)) + 
     xs^3*(yS^2*(-yk^2 + yS^2)*\[Lambda]Lik + 
       2*xt^2*(\[Lambda]Lik + yk*\[Lambda]Rik) - 
       xt*(2*yk^2*\[Lambda]Lik + 3*yS^2*\[Lambda]Lik + 2*yk^3*\[Lambda]Rik + 
         2*yk*yS^2*\[Lambda]Rik)) + xs^2*(2*yS^2*(yk^2 - yS^2)*\[Lambda]Lik + 
       2*xt^3*yk*\[Lambda]Rik + xt*(-(yk^4*\[Lambda]Lik) - 
         yk^2*(-3 + yS^2)*\[Lambda]Lik + 2*yS^2*(2 + yS^2)*\[Lambda]Lik + 
         4*yk^3*\[Lambda]Rik + 4*yk*yS^2*\[Lambda]Rik) - 
       xt^2*((1 + yk^2 + 3*yS^2)*\[Lambda]Lik + 2*yk*(1 + 2*yk^2 + 2*yS^2)*
          \[Lambda]Rik)))*DD[1])/(xs^2*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^4*(Qj - Qk)*Qk*\[Lambda]Ljk*((-1 + xt)^2*xt*yS^2*(-yk^2 + yS^2)*
      \[Lambda]Lik + xs^3*yk*(yk*(-yk^2 + yS^2)*\[Lambda]Lik + 
       2*xt^2*\[Lambda]Rik + xt*(yk*\[Lambda]Lik - 2*yk^2*\[Lambda]Rik - 
         2*yS^2*\[Lambda]Rik)) - xs*(-1 + xt)*
      (yk^2*(-yk^2 + yS^2)*\[Lambda]Lik + xt*(yk^4*\[Lambda]Lik + 
         yk^2*(-1 + yS^2)*\[Lambda]Lik - yS^2*(1 + 2*yS^2)*\[Lambda]Lik - 
         2*yk^3*\[Lambda]Rik - 2*yk*yS^2*\[Lambda]Rik) + 
       xt^2*(2*yk^2*\[Lambda]Lik + 3*yS^2*\[Lambda]Lik + 
         2*yk^3*\[Lambda]Rik + 2*yk*yS^2*\[Lambda]Rik)) + 
     xs^2*(2*yk^2*(yk^2 - yS^2)*\[Lambda]Lik + 
       2*xt^3*(\[Lambda]Lik + yk*\[Lambda]Rik) + 
       xt*(-2*yk^4*\[Lambda]Lik + yk^2*yS^2*\[Lambda]Lik + 
         yS^2*(1 + yS^2)*\[Lambda]Lik + 4*yk^3*\[Lambda]Rik + 
         4*yk*yS^2*\[Lambda]Rik) - xt^2*((1 + yk^2 + 3*yS^2)*\[Lambda]Lik + 
         2*yk*(1 + 2*yk^2 + 2*yS^2)*\[Lambda]Rik)))*DD[2])/
   (xs^2*xt^2*(-1 + xs + xt)^3) - (2*mi^4*(Qj - Qk)^2*yS^2*\[Lambda]Ljk*
    (xs^2*\[Lambda]Lik + xt*(-yk^2 + yS^2)*\[Lambda]Lik - 
     xs*((yk^2 - yS^2)*\[Lambda]Lik + xt*(\[Lambda]Lik + 4*yk*\[Lambda]Rik)))*
    DD[3])/(xs^2*xt^2*(-1 + xs + xt)) + 
  (2*mi^4*Qk^2*yk^2*\[Lambda]Ljk*(xs^2*\[Lambda]Lik + 
     xt*(yk^2 - yS^2)*\[Lambda]Lik + xs*(3*xt*\[Lambda]Lik + 
       yk^2*\[Lambda]Lik - yS^2*\[Lambda]Lik + 4*xt*yk*\[Lambda]Rik))*DD[4])/
   (xs^2*xt^2*(-1 + xs + xt)) + (2*mi^4*(Qj - Qk)^2*yS^2*\[Lambda]Ljk*
    (xs*(yk^2 - yS^2)*\[Lambda]Lik - xt*(xt - yk^2 + yS^2)*\[Lambda]Lik + 
     xs*xt*(\[Lambda]Lik + 4*yk*\[Lambda]Rik))*DD[5])/
   (xs^2*xt^2*(-1 + xs + xt)) + (2*mi^4*Qk^2*yk^2*\[Lambda]Ljk*
    (xt*(xt + yk^2 - yS^2)*\[Lambda]Lik + 
     xs*(3*xt*\[Lambda]Lik + yk^2*\[Lambda]Lik - yS^2*\[Lambda]Lik + 
       4*xt*yk*\[Lambda]Rik))*DD[6])/(xs^2*xt^2*(-1 + xs + xt)), 
 (Qj*(2*Qk*(1 - 2*xs)*\[Lambda]Rik + Qj*(-1 + yk^2 - yS^2)*\[Lambda]Rik + 
     Qj*xs*(2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - 
       (-2 + yS^2)*\[Lambda]Rik))*\[Lambda]Rjk)/((-1 + xs)*xs^2*xt) + 
  (Qj^2*yk^2*((yk^2 - yS^2)*\[Lambda]Rik + 
     xs*(2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - yS^2*\[Lambda]Rik))*
    \[Lambda]Rjk*BB[1])/((-1 + xs)*xs^2*xt*(yk^2 - yS^2)) + 
  (Qj^2*yS^2*((-yk^2 + yS^2)*\[Lambda]Rik + 
     xs*(-2*yk*\[Lambda]Lik - yk^2*\[Lambda]Rik + yS^2*\[Lambda]Rik))*
    \[Lambda]Rjk*BB[2])/((-1 + xs)*xs^2*xt*(yk^2 - yS^2)) + 
  (Qj*(2*Qk*(-xs + xt)*\[Lambda]Rik + Qj*(xs*(1 - yk^2 + yS^2)*\[Lambda]Rik + 
       (-1 + xs)*xt*(2*yk*\[Lambda]Lik + \[Lambda]Rik + yk^2*\[Lambda]Rik - 
         yS^2*\[Lambda]Rik) - xt^2*(2*(-1 + xs)*yk*\[Lambda]Lik + 
         (-2 + xs)*yk^2*\[Lambda]Rik + (xs + 2*yS^2 - xs*yS^2)*
          \[Lambda]Rik)))*\[Lambda]Rjk*BB[3])/((-1 + xs)*xs^2*(-1 + xt)*
    xt^2) + (Qj*(-2*Qk*xs^2*\[Lambda]Rik + Qj*xs*(xs - yk^2 + yS^2)*
      \[Lambda]Rik + Qj*xt*(-2*yk*\[Lambda]Lik - 3*yk^2*\[Lambda]Rik + 
       (xs + 3*yS^2)*\[Lambda]Rik))*\[Lambda]Rjk*BB[4])/
   ((-1 + xs)*xs^2*xt^2) + (Qj*(2*Qk*xt - Qj*(xt - yk^2 + yS^2))*\[Lambda]Rik*
    \[Lambda]Rjk*BB[5])/(xs^2*(-1 + xt)*xt) + 
  (Qk^2*(xs - xt)*\[Lambda]Rik*\[Lambda]Rjk*BB[6])/(xs^2*xt^2) - 
  ((Qj - Qk)^2*(xs - xt)*\[Lambda]Rik*\[Lambda]Rjk*BB[7])/(xs^2*xt^2) + 
  (mi^2*Qk*(Qk*(-1 + xs)*xs*(xs^4 + xs^3*(-2 + 3*xt + 2*yk^2 - 2*yS^2) + 
       2*(-1 + xt)^2*xt*(yk^2 - yS^2) + xs^2*(1 + 2*xt^2 - 4*yk^2 + 4*yS^2 + 
         4*xt*(-1 + yk^2 - yS^2)) + xs*(xt^2*(-2 + 4*yk^2 - 4*yS^2) + 
         2*(yk^2 - yS^2) + xt*(1 - 6*yk^2 + 6*yS^2)))*\[Lambda]Rik + 
     Qj*((-1 + xt)^2*xt*(yk^2 - yS^2)^2*\[Lambda]Rik + 
       xs^4*(xt^2 - 2*xt*yS^2 + (yk^2 - yS^2)^2)*\[Lambda]Rik + 
       xs^3*(xt^3*\[Lambda]Rik - 3*(yk^2 - yS^2)^2*\[Lambda]Rik + 
         xt*(yk^4 + 4*yS^2 - 2*yk^2*yS^2 + yS^4)*\[Lambda]Rik - 
         xt^2*(4*yk*\[Lambda]Lik + \[Lambda]Rik + 6*yk^2*\[Lambda]Rik + 
           2*yS^2*\[Lambda]Rik)) - xs*(-1 + xt)*
        (-((yk^2 - yS^2)^2*\[Lambda]Rik) - 2*xt*(yk^2 - yS^2)^2*
          \[Lambda]Rik + 4*xt^3*yk*(\[Lambda]Lik + yk*\[Lambda]Rik) + 
         xt^2*(-4*yk*\[Lambda]Lik + yk^4*\[Lambda]Rik + yS^4*\[Lambda]Rik - 
           2*yk^2*(3 + yS^2)*\[Lambda]Rik)) - 
       xs^2*(-3*(yk^2 - yS^2)^2*\[Lambda]Rik + 
         xt*(yk^4 + 2*yS^2 - 2*yk^2*yS^2 + yS^4)*\[Lambda]Rik + 
         xt^3*(8*yk*\[Lambda]Lik + \[Lambda]Rik + 10*yk^2*\[Lambda]Rik) + 
         xt^2*(-8*yk*\[Lambda]Lik + yk^4*\[Lambda]Rik + yS^2*(-2 + yS^2)*
            \[Lambda]Rik - 2*yk^2*(6 + yS^2)*\[Lambda]Rik))))*\[Lambda]Rjk*
    CC[1])/(2*(-1 + xs)*xs^2*xt^3*(-1 + xs + xt)^2) - 
  (mi^2*(Qk^2*(-1 + xs)*xs*(xs^4 - 2*(-1 + xt)^2*xt*(yk^2 - yS^2) + 
       xs^3*(-2 + 3*xt - 2*yk^2 + 2*yS^2) + xs^2*(1 + 4*xt^2 + 4*yk^2 - 
         4*yS^2 - 4*xt*(1 + yk^2 - yS^2)) + xs*(2*xt^3 - 2*yk^2 + 2*yS^2 + 
         xt*(1 + 6*yk^2 - 6*yS^2) - 4*xt^2*(1 + yk^2 - yS^2))) + 
     Qj^2*(xs^6 + (-1 + xt)^2*xt*(yk^2 - yS^2)^2 + 
       xs^5*(-3 + 3*xt - 2*yk^2 + 2*yS^2) + xs^4*(3 + 3*xt^2 + yk^4 - 
         6*yS^2 + yS^4 - 2*yk^2*(-3 + yS^2) + xt*(-7 - 4*yk^2 + 2*yS^2)) - 
       xs*(-1 + xt)*((-1 - 2*xt + xt^2)*yk^4 - 2*(-1 - 2*xt + xt^2)*yk^2*
          yS^2 + yS^2*(4*xt^3 - yS^2 - 2*xt*yS^2 + xt^2*(-6 + yS^2))) + 
       xs^3*(-1 + xt^3 - 3*yk^4 + 6*yS^2 - 3*yS^4 + 6*yk^2*(-1 + yS^2) - 
         xt^2*(5 + 2*yk^2 + 6*yS^2) + xt*(5 + yk^4 - 4*yS^2 + yS^4 - 
           2*yk^2*(-4 + yS^2))) - xs^2*(-3*yk^4 + 2*yS^2 - 3*yS^4 + 
         yk^2*(-2 + 6*yS^2) + xt^3*(1 + 10*yS^2) + 
         xt*(yk^4 - 2*yk^2*(-2 + yS^2) + (-1 + yS^2)^2) + 
         xt^2*(-2 + yk^4 - 12*yS^2 + yS^4 - 2*yk^2*(1 + yS^2)))) + 
     Qj*Qk*(-2*xs^6 + xs^5*(6 - 6*xt + 4*yk^2 - 4*yS^2) - 
       (-1 + xt)^2*xt*(yk^2 - yS^2)^2 - xs^4*(6 + 7*xt^2 + yk^4 - 12*yS^2 + 
         yS^4 - 2*xt*(7 + 4*yk^2 - 3*yS^2) - 2*yk^2*(-6 + yS^2)) + 
       xs^3*(2 - 3*xt^3 + 3*yk^4 - 12*yS^2 + 3*yS^4 - 6*yk^2*(-2 + yS^2) + 
         xt^2*(13 + 6*yk^2 + 2*yS^2) - xt*(10 + yk^4 - 14*yS^2 + yS^4 - 
           2*yk^2*(-9 + yS^2))) + xs^2*(-3*yk^4 + 4*yS^2 - 3*yS^4 + 
         yk^2*(-4 + 6*yS^2) + xt^3*(3 + 2*yk^2 + 8*yS^2) + 
         xt*(2 + yk^4 - 10*yS^2 + yS^4 - 2*yk^2*(-6 + yS^2)) + 
         xt^2*(-6 + yk^4 - 4*yS^2 + yS^4 - 2*yk^2*(5 + yS^2))) + 
       xs*(4*xt^4*yS^2 + (yk^2 - yS^2)^2 + xt^3*(yk^4 + yS^2*(-8 + yS^2) - 
           2*yk^2*(1 + yS^2)) + xt*(yk^4 - 2*yk^2*(1 + yS^2) + 
           yS^2*(2 + yS^2)) + xt^2*(-3*yk^4 + 2*yS^2 - 3*yS^4 + 
           yk^2*(4 + 6*yS^2)))))*\[Lambda]Rik*\[Lambda]Rjk*CC[2])/
   (2*(-1 + xs)*xs^2*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*Qk*(-(Qk*(-1 + xt)*xt^2*(xt^2 - 2*yk^2 + 2*yS^2 + 
        xt*(-1 + 2*yk^2 - 2*yS^2) + xs*(1 + xt + 2*yk^2 - 2*yS^2))) + 
     Qj*(-((-1 + xt)^2*xt*(yk^2 - yS^2)^2) + 
       xs^3*(xt^2 + xt*(4*yk^2 - 2*yS^2) + (yk^2 - yS^2)^2) + 
       xs*(2*xt^3*yk^2 + 2*xt*(yk^2 - yS^2) + (yk^2 - yS^2)^2 - 
         xt^2*(yk^4 - 2*yk^2*(-2 + yS^2) + yS^2*(-2 + yS^2))) + 
       xs^2*(xt^3 - 2*(yk^2 - yS^2)^2 + xt^2*(6*yk^2 - 2*(1 + yS^2)) + 
         xt*(yk^4 - 2*yk^2*(3 + yS^2) + yS^2*(4 + yS^2)))))*\[Lambda]Rik*
    \[Lambda]Rjk*CC[3])/(2*xs^3*xt^2*(-1 + xs + xt)^2) + 
  (mi^2*(Qk^2*(-1 + xt)*xt^2*(2*xs^2 + xt^2 + 2*(yk^2 - yS^2) + 
       xt*(-1 - 2*yk^2 + 2*yS^2) + xs*(-1 + 3*xt - 2*yk^2 + 2*yS^2)) - 
     Qj*Qk*(xs^3*(xt^2 - 2*xt*yS^2 - (yk^2 - yS^2)^2) + 
       (-1 + xt)^2*xt*(2*xt^2 - 4*xt*(yk^2 - yS^2) + (yk^2 - yS^2)^2) + 
       xs^2*(5*xt^3 + 2*(yk^2 - yS^2)^2 - 2*xt^2*(2 + yk^2 + yS^2) - 
         xt*(yk^4 - 2*yS^2 - 2*yk^2*yS^2 + yS^4)) + 
       xs*(6*xt^4 - (yk^2 - yS^2)^2 + xt^3*(-8 - 6*yk^2 + 4*yS^2) + 
         xt^2*(2 + yk^4 - 4*yS^2 + yS^4 - 2*yk^2*(-3 + yS^2)))) + 
     Qj^2*(xs^3*(xt^2 - 2*xt*yS^2 - (yk^2 - yS^2)^2) + 
       xt*(xt^2 + yk^2 - yS^2 + xt*(-1 - yk^2 + yS^2))^2 + 
       xs^2*(3*xt^3 + 2*(yk^2 - yS^2)^2 - 2*xt^2*(1 + yk^2 + yS^2) - 
         xt*(yk^4 - 2*yS^2 - 2*yk^2*yS^2 + yS^4)) + 
       xs*(3*xt^4 - (yk^2 - yS^2)^2 + xt^3*(-4 - 4*yk^2 + 2*yS^2) + 
         xt^2*(yk^4 - 2*yk^2*(-2 + yS^2) + (-1 + yS^2)^2))))*\[Lambda]Rik*
    \[Lambda]Rjk*CC[4])/(2*xs^3*xt^2*(-1 + xs + xt)^2) - 
  (mi^2*Qk^2*(xs^4 + xs^3*(-1 + 2*xt + 2*yk^2 - 2*yS^2) + 
     xs*(xt^2*(-1 + 6*yk^2 - 2*yS^2) - 4*xt*(yk^2 - yS^2) - 
       2*(yk^2 - yS^2)^2) + xs^2*(2*xt^2 + xt*(-1 + 2*yk^2 - 6*yS^2) + 
       2*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2))) + 
     xt*(-xt^3 + 2*(yk^2 - yS^2)^2 + xt^2*(1 - 2*yk^2 + 2*yS^2) - 
       2*xt*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2))))*\[Lambda]Rik*
    \[Lambda]Rjk*CC[5])/(2*xs^3*xt^3) + 
  (mi^2*(Qj - Qk)^2*(xs^4 + xs^3*(-1 + 2*xt - 2*yk^2 + 2*yS^2) + 
     xs*(-2*xt^3 - 2*(yk^2 - yS^2)^2 + xt^2*(1 + 2*yk^2 + 2*yS^2)) - 
     xs^2*(xt*(1 + 2*yk^2 + 2*yS^2) - 2*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + 
         yS^4)) + xt*(-xt^3 + xt^2*(1 + 2*yk^2 - 2*yS^2) + 
       2*(yk^2 - yS^2)^2 - 2*xt*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4)))*
    \[Lambda]Rik*\[Lambda]Rjk*CC[6])/(2*xs^3*xt^3) - 
  (mi^2*(Qk^2*(-1 + xs)*xs*(xs^4 - 2*(-1 + xt)*xt^2*(yk^2 - yS^2) + 
       xs*(1 - 2*xt + 2*xt^2)*(xt - 2*yk^2 + 2*yS^2) + 
       xs^3*(-2 + 3*xt - 2*yk^2 + 2*yS^2) + xs^2*(1 + 4*xt^2 + 4*yk^2 - 
         4*yS^2 - 4*xt*(1 + yk^2 - yS^2))) + 
     Qj^2*(xs^6 + (-1 + xt)^2*xt*(yk^2 - yS^2)^2 + 
       xs^5*(-3 + 3*xt - 2*yk^2 + 2*yS^2) + xs^4*(3 + 3*xt^2 + yk^4 - 
         6*yS^2 + yS^4 - 2*yk^2*(-3 + yS^2) + xt*(-7 - 4*yk^2 + 2*yS^2)) + 
       xs^3*(-1 + xt^3 - 3*yk^4 + 6*yS^2 - 3*yS^4 + 6*yk^2*(-1 + yS^2) - 
         xt^2*(5 + 2*yk^2 + 6*yS^2) + xt*(5 + yk^4 - 4*yS^2 + yS^4 - 
           2*yk^2*(-4 + yS^2))) - xs*(-1 + xt)*((-1 - 2*xt + xt^2)*yk^4 - 
         2*(-1 - 2*xt + xt^2)*yk^2*yS^2 + yS^2*(4*xt^3 - yS^2 - 2*xt*yS^2 + 
           xt^2*(-2 + yS^2))) - xs^2*(-3*yk^4 + 2*yS^2 - 3*yS^4 + 
         yk^2*(-2 + 6*yS^2) + xt^3*(1 + 10*yS^2) + 
         xt*(yk^4 - 2*yk^2*(-2 + yS^2) + (-1 + yS^2)^2) + 
         xt^2*(-2 + yk^4 - 8*yS^2 + yS^4 - 2*yk^2*(1 + yS^2)))) + 
     Qj*Qk*(-2*xs^6 + xs^5*(6 - 6*xt + 4*yk^2 - 4*yS^2) - 
       (-1 + xt)^2*xt*(yk^2 - yS^2)^2 - xs^4*(6 + 7*xt^2 + yk^4 - 12*yS^2 + 
         yS^4 - 2*xt*(7 + 4*yk^2 - 3*yS^2) - 2*yk^2*(-6 + yS^2)) + 
       xs^3*(2 - 3*xt^3 + 3*yk^4 - 12*yS^2 + 3*yS^4 - 6*yk^2*(-2 + yS^2) + 
         xt^2*(11 + 6*yk^2 + 2*yS^2) - xt*(10 + yk^4 - 12*yS^2 + yS^4 - 
           2*yk^2*(-8 + yS^2))) + xs^2*(-3*yk^4 + 4*yS^2 - 3*yS^4 + 
         yk^2*(-4 + 6*yS^2) + xt^3*(3 + 2*yk^2 + 8*yS^2) + 
         xt*(2 + yk^4 - 6*yS^2 + yS^4 - 2*yk^2*(-4 + yS^2)) + 
         xt^2*(-4 + yk^4 - 2*yS^2 + yS^4 - 2*yk^2*(4 + yS^2))) + 
       xs*(4*xt^4*yS^2 + (yk^2 - yS^2)^2 + xt*(yk^2 - yS^2)^2 + 
         xt^3*(yk^4 + yS^2*(-4 + yS^2) - 2*yk^2*(1 + yS^2)) + 
         xt^2*(-3*yk^4 - 3*yS^4 + yk^2*(2 + 6*yS^2)))))*\[Lambda]Rik*
    \[Lambda]Rjk*CC[7])/(2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*Qk*(Qk*(-1 + xs)*xs*(xs^4 + xs^3*(-2 + 3*xt + 2*yk^2 - 2*yS^2) + 
       2*(-1 + xt)*xt^2*(yk^2 - yS^2) + xs^2*(1 + 2*xt^2 - 4*yk^2 + 4*yS^2 + 
         4*xt*(-1 + yk^2 - yS^2)) + xs*(2*(yk^2 - yS^2) + 
         4*xt^2*(yk^2 - yS^2) + xt*(1 - 4*yk^2 + 4*yS^2)))*\[Lambda]Rik + 
     Qj*((-1 + xt)^2*xt*(yk^2 - yS^2)^2*\[Lambda]Rik + 
       xs^4*(xt^2 - 2*xt*yS^2 + (yk^2 - yS^2)^2)*\[Lambda]Rik - 
       xs*(-1 + xt)*(-((yk^2 - yS^2)^2*\[Lambda]Rik) - 
         2*xt*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4)*\[Lambda]Rik + 
         4*xt^3*yk*(\[Lambda]Lik + yk*\[Lambda]Rik) + 
         xt^2*(-4*yk*\[Lambda]Lik + yk^4*\[Lambda]Rik + yS^4*\[Lambda]Rik - 
           2*yk^2*(1 + yS^2)*\[Lambda]Rik)) + 
       xs^3*(xt^3*\[Lambda]Rik - 3*(yk^2 - yS^2)^2*\[Lambda]Rik + 
         xt*(yk^4 - 2*yk^2*(1 + yS^2) + yS^2*(6 + yS^2))*\[Lambda]Rik - 
         xt^2*(4*yk*\[Lambda]Lik + 6*yk^2*\[Lambda]Rik + 
           (3 + 2*yS^2)*\[Lambda]Rik)) - 
       xs^2*(-3*(yk^2 - yS^2)^2*\[Lambda]Rik + xt*(yk^4 - 2*yk^2*(2 + yS^2) + 
           yS^2*(6 + yS^2))*\[Lambda]Rik + xt^3*(8*yk*\[Lambda]Lik + 
           \[Lambda]Rik + 10*yk^2*\[Lambda]Rik) + 
         xt^2*(-8*yk*\[Lambda]Lik + yk^4*\[Lambda]Rik - 2*yk^2*(3 + yS^2)*
            \[Lambda]Rik + (-2 - 4*yS^2 + yS^4)*\[Lambda]Rik))))*\[Lambda]Rjk*
    CC[8])/(2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*(Qk^2*(-1 + xt)*xt*(2*xs^2*(xt^2 - yk^2 + yS^2) + 
       (-1 + xt)^2*xt*(xt - 2*yk^2 + 2*yS^2) + xs*(-1 + xt)*
        (3*xt^2 - 2*yk^2 + 2*yS^2 + xt*(-1 - 2*yk^2 + 2*yS^2))) + 
     Qj^2*((-1 + xt)^3*xt*(xt - yk^2 + yS^2)^2 + xs*(-1 + xt)^2*
        (3*xt^3 + (yk^2 - yS^2)^2 + xt*(yk^2 - yS^2)^2 + 
         xt^2*(-1 - 4*yk^2 + 2*yS^2)) + xs^3*(xt^3 + (yk^2 - yS^2)^2 - 
         xt^2*(1 + 2*yS^2) - xt*(yk^4 + 2*yS^2 - 2*yk^2*yS^2 + yS^4)) + 
       xs^2*(3*xt^4 - 2*(yk^2 - yS^2)^2 - xt^3*(5 + 2*yk^2 + 2*yS^2) + 
         xt*(3*yk^4 + 2*yS^2 - 6*yk^2*yS^2 + 3*yS^4) - 
         xt^2*(-2 + yk^4 + yS^4 - 2*yk^2*(1 + yS^2)))) + 
     Qj*Qk*(-((-1 + xt)^3*xt*(2*xt^2 - 4*xt*(yk^2 - yS^2) + 
          (yk^2 - yS^2)^2)) + xs^3*(-xt^3 - (yk^2 - yS^2)^2 + 
         xt^2*(1 + 2*yS^2) + xt*(yk^4 + 2*yS^2 - 2*yk^2*yS^2 + yS^4)) + 
       xs^2*(-5*xt^4 + 2*(yk^2 - yS^2)^2 + xt^3*(7 + 2*yk^2 + 2*yS^2) + 
         xt^2*(-2 + yk^4 - 2*yS^2 - 2*yk^2*yS^2 + yS^4) - 
         xt*(3*yk^4 + 3*yS^4 + yk^2*(2 - 6*yS^2))) - 
       xs*(-1 + xt)^2*(6*xt^3 + (yk^2 - yS^2)^2 + 
         xt^2*(-2 - 6*yk^2 + 4*yS^2) + xt*(yk^4 - 2*yk^2*(1 + yS^2) + 
           yS^2*(2 + yS^2)))))*\[Lambda]Rik*\[Lambda]Rjk*CC[9])/
   (2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*Qk*(-(Qk*(-1 + xt)*xt*((-1 + xt)^2*xt*(xt + 2*yk^2 - 2*yS^2) + 
        xs*(-1 + xt^2)*(xt + 2*yk^2 - 2*yS^2) + 2*xs^2*(xt + yk^2 - yS^2))) + 
     Qj*(-((-1 + xt)^3*xt*(yk^2 - yS^2)^2) + xs*(-1 + xt)^2*
        (2*xt^2*yk^2 - (yk^2 - yS^2)^2 - xt*(yk^2 - yS^2)^2) + 
       xs^3*(xt^3 + xt^2*(-1 + 4*yk^2 - 2*yS^2) - (yk^2 - yS^2)^2 + 
         xt*(yk^4 + 2*yS^2 - 2*yk^2*yS^2 + yS^4)) + 
       xs^2*(xt^4 + xt^3*(-1 + 6*yk^2 - 2*yS^2) + 2*(yk^2 - yS^2)^2 - 
         xt*(3*yk^4 + 2*yS^2 - 6*yk^2*yS^2 + 3*yS^4) + 
         xt^2*(yk^4 - 2*yk^2*(3 + yS^2) + yS^2*(4 + yS^2)))))*\[Lambda]Rik*
    \[Lambda]Rjk*CC[10])/(2*xs^3*xt^3*(-1 + xs + xt)^2) - 
  (mi^2*Qk^2*(xs + xt)*(xs^3 + xs*xt^2 + xs^2*(xt + 2*yk^2 - 2*yS^2) - 
     xt^2*(xt + 2*yk^2 - 2*yS^2))*\[Lambda]Rik*\[Lambda]Rjk*CC[11])/
   (2*xs^3*xt^3) + (mi^2*(Qj - Qk)^2*(xs - xt)*(xs + xt)^2*
    (xs + xt - 2*yk^2 + 2*yS^2)*\[Lambda]Rik*\[Lambda]Rjk*CC[12])/
   (2*xs^3*xt^3) - (mi^4*(Qj - Qk)*Qk*(xs*(xt - (yk - yS)^2) - 
     (-1 + xt)*(yk - yS)^2)*(-((-1 + xt)*xt*(yk^2 - yS^2)) + 
     xs^2*(xt + yk^2 - yS^2) + xs*(xt^2 - yk^2 + yS^2))*
    (-((-1 + xt)*(yk + yS)^2) + xs*(xt - (yk + yS)^2))*\[Lambda]Rik*
    \[Lambda]Rjk*DD[1])/(2*xs^3*xt^3*(-1 + xs + xt)^2) - 
  (mi^4*(Qj - Qk)*Qk*(xs*(xt - (yk - yS)^2) - (-1 + xt)*(yk - yS)^2)*
    (-((-1 + xt)*xt*(yk^2 - yS^2)) + xs^2*(xt + yk^2 - yS^2) + 
     xs*(-2*xt + xt^2 - yk^2 + yS^2))*(-((-1 + xt)*(yk + yS)^2) + 
     xs*(xt - (yk + yS)^2))*\[Lambda]Rik*\[Lambda]Rjk*DD[2])/
   (2*xs^3*xt^3*(-1 + xs + xt)^2) - 
  (mi^4*(Qj - Qk)^2*(xs^2 + xt*(yk^2 - yS^2) + xs*(xt - yk^2 + yS^2))*
    (xs^3 + (-1 + xt)*(yk^2 - yS^2)^2 + xs^2*(-1 + xt - 2*yk^2 + 2*yS^2) + 
     xs*(yk^4 + yS^2*(-2 - 2*xt + yS^2) - 2*yk^2*(-1 + xt + yS^2)))*
    \[Lambda]Rik*\[Lambda]Rjk*DD[3])/(2*xs^3*xt^3) + 
  (mi^4*Qk^2*(xs^5*\[Lambda]Rik + xs^4*(-1 + 2*xt + 3*yk^2 - 3*yS^2)*
      \[Lambda]Rik - (-1 + xt)*xt*(yk^2 - yS^2)^3*\[Lambda]Rik - 
     xs*(yk^2 - yS^2)*(xt*(yk^2 - yS^2) + (yk^2 - yS^2)^2 + 
       xt^2*(-5*yk^2 + yS^2))*\[Lambda]Rik + 
     xs^3*(xt^2 + xt*(-1 + 2*yk^2 - 6*yS^2) + 
       3*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2)))*\[Lambda]Rik - 
     xs^2*(-((-3 + yk^2 - yS^2)*(yk^2 - yS^2)^2*\[Lambda]Rik) + 
       xt*(yk^2 - yS^2)*(3 + 4*yS^2)*\[Lambda]Rik + 
       xt^2*(4*yk*\[Lambda]Lik + 9*yk^2*\[Lambda]Rik + 3*yS^2*\[Lambda]Rik)))*
    \[Lambda]Rjk*DD[4])/(2*xs^3*xt^3) + 
  (mi^4*(Qj - Qk)^2*(xs*(xt + yk^2 - yS^2) + xt*(xt - yk^2 + yS^2))*
    (xt^3 + (-1 + xs)*(yk^2 - yS^2)^2 + xt^2*(-1 + xs - 2*yk^2 + 2*yS^2) + 
     xt*(yk^4 + yS^2*(-2 - 2*xs + yS^2) - 2*yk^2*(-1 + xs + yS^2)))*
    \[Lambda]Rik*\[Lambda]Rjk*DD[5])/(2*xs^3*xt^3) + 
  (mi^4*Qk^2*(xs - xt)*(xt + yk^2 - yS^2)*
    (xt^3 + xt^2*(-1 + xs + 2*yk^2 - 2*yS^2) + (-1 + xs)*(yk^2 - yS^2)^2 + 
     xt*(yk^4 + yS^2*(2 - 2*xs + yS^2) - 2*yk^2*(1 + xs + yS^2)))*
    \[Lambda]Rik*\[Lambda]Rjk*DD[6])/(2*xs^3*xt^3), 
 (-2*Qj*(-2*Qk + Qj*(1 + yk^2 - yS^2))*\[Lambda]Lik*\[Lambda]Ljk)/
   ((-1 + xs)*xs^2*xt) + (2*Qj^2*yk^2*\[Lambda]Lik*\[Lambda]Ljk*BB[1])/
   (xs^2*xt - xs^3*xt) + (2*Qj^2*yS^2*\[Lambda]Lik*\[Lambda]Ljk*BB[2])/
   ((-1 + xs)*xs^2*xt) + 
  (2*Qj*\[Lambda]Ljk*(2*Qk*xs*(1 - 2*xs + xs^2 + xt - xt^2)*\[Lambda]Lik + 
     Qj*(-1 + xt)*(-1 + xs + xt)*(xs*(1 + yk^2 - yS^2)*\[Lambda]Lik + 
       2*xs*yk*\[Lambda]Rik - 2*(yk^2*\[Lambda]Lik - yS^2*\[Lambda]Lik + 
         yk*\[Lambda]Rik)))*BB[3])/((-1 + xs)^2*xs^2*(-1 + xt)*xt*
    (-1 + xs + xt)) - (2*Qj*\[Lambda]Ljk*(-2*Qk*xs*xt*\[Lambda]Lik + 
     Qj*(-1 + xs + xt)*(-(yk^2*\[Lambda]Lik) + yS^2*\[Lambda]Lik - 
       2*yk*\[Lambda]Rik + xs*(\[Lambda]Lik + 2*yk*\[Lambda]Rik)))*BB[4])/
   ((-1 + xs)^2*xs^2*xt*(-1 + xs + xt)) - 
  (4*Qj*Qk*\[Lambda]Lik*\[Lambda]Ljk*BB[5])/(xs*(-1 + xt)*xt*
    (-1 + xs + xt)) - (2*mi^2*Qk*\[Lambda]Ljk*(-(Qk*xs*\[Lambda]Lik) + 
     Qj*(-1 + xt)*yk*\[Lambda]Rik + Qj*xs*(\[Lambda]Lik + yk*\[Lambda]Rik))*
    CC[1])/(xs*xt*(-1 + xs + xt)^2) + 
  (2*mi^2*Qk^2*yk*\[Lambda]Ljk*\[Lambda]Rik*CC[3])/(xs^2*(-1 + xs + xt)) - 
  (2*mi^2*(Qj - Qk)*Qk*\[Lambda]Ljk*((-1 + xt)*yk*\[Lambda]Rik + 
     xs*(\[Lambda]Lik + yk*\[Lambda]Rik))*CC[4])/(xs^2*(-1 + xs + xt)^2) - 
  (4*mi^2*Qk^2*yk*\[Lambda]Ljk*\[Lambda]Rik*CC[5])/(xs^2*xt) - 
  (2*mi^2*(Qj - Qk)*\[Lambda]Ljk*(2*Qj*(-1 + xs + xt)^2*yS^2*\[Lambda]Lik + 
     Qk*(-1 + xs)^2*((-1 + xt)*yk*\[Lambda]Rik + 
       xs*(\[Lambda]Lik + yk*\[Lambda]Rik)))*CC[7])/
   ((-1 + xs)*xs^2*xt*(-1 + xs + xt)^2) + 
  (2*mi^2*Qk*yk*\[Lambda]Ljk*(Qk*(-1 + xs)^2*\[Lambda]Rik + 
     2*Qj*(-1 + xs + xt)*(yk*\[Lambda]Lik + \[Lambda]Rik - xs*\[Lambda]Rik))*
    CC[8])/((-1 + xs)*xs^2*xt*(-1 + xs + xt)) + 
  (2*mi^2*Qk*\[Lambda]Ljk*(Qk*xs*(-1 + xt)*\[Lambda]Lik + 
     Qj*(2*xs^2*yk*\[Lambda]Rik + (-1 + xt)^2*yk*\[Lambda]Rik - 
       xs*(-1 + xt)*(\[Lambda]Lik - 3*yk*\[Lambda]Rik)))*CC[10])/
   (xs^2*xt*(-1 + xs + xt)^2) + (2*mi^4*(Qj - Qk)*Qk*\[Lambda]Ljk*
    (-((-1 + xt)^2*yk*(yk^2 - yS^2)*\[Lambda]Rik) + 
     xs*(-1 + xt)*(-(yk^2*\[Lambda]Lik) - yS^2*\[Lambda]Lik - 
       2*yk^3*\[Lambda]Rik + yk*(xt + 2*yS^2)*\[Lambda]Rik) + 
     xs^2*(-(yk^2*\[Lambda]Lik) - yS^2*\[Lambda]Lik - yk^3*\[Lambda]Rik + 
       yk*yS^2*\[Lambda]Rik + xt*(\[Lambda]Lik + yk*\[Lambda]Rik)))*DD[1])/
   (xs^2*xt*(-1 + xs + xt)^2) + (2*mi^4*Qk^2*yk*\[Lambda]Ljk*
    (-(xs^2*\[Lambda]Rik) + (-1 + xt)*(yk^2 - yS^2)*\[Lambda]Rik + 
     xs*(2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - (-1 + xt + yS^2)*
        \[Lambda]Rik))*DD[4])/(xs^2*xt*(-1 + xs + xt)) + 
  (2*mi^4*Qk^2*yk*\[Lambda]Ljk*((xt^2 - yk^2 + yS^2 + xt*(-1 + yk^2 - yS^2))*
      \[Lambda]Rik + xs*(2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik + 
       (xt - yS^2)*\[Lambda]Rik))*DD[6])/(xs^2*xt*(-1 + xs + xt)), 
 (-2*Qj*(xs - xt)*(-2*Qk*xs*xt + Qj*(-1 + xt)*(yk^2 - yS^2) + 
     Qj*xs*(xt + yk^2 - yS^2))*\[Lambda]Rik*\[Lambda]Rjk)/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)) + 
  (2*Qj^2*(-xs + xt)*yk^2*\[Lambda]Rik*\[Lambda]Rjk*BB[1])/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2) + 
  (2*Qj^2*(xs - xt)*yS^2*\[Lambda]Rik*\[Lambda]Rjk*BB[2])/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2) - 
  (2*Qj*(-2*xs + xs^2 - (-2 + xt)*xt)*(2*Qk + Qj*(-1 + yk^2 - yS^2))*
    \[Lambda]Rik*\[Lambda]Rjk*BB[3])/((-1 + xs)^2*xs*(-1 + xt)^2*xt*
    (-1 + xs + xt)) + 
  (2*Qj*(-2*Qk*xs + Qj*(xs + yk^2 - 2*xs*yk^2 - yS^2 + 2*xs*yS^2))*
    \[Lambda]Rik*\[Lambda]Rjk*BB[4])/((-1 + xs)^2*xs^2*xt*(-1 + xs + xt)) + 
  (2*Qj*(2*Qk*xt + Qj*(-yk^2 + yS^2 + xt*(-1 + 2*yk^2 - 2*yS^2)))*
    \[Lambda]Rik*\[Lambda]Rjk*BB[5])/(xs*(-1 + xt)^2*xt^2*(-1 + xs + xt)) - 
  (2*mi^2*Qk*(Qk*xs*((-1 + xs)*(yk^2 - yS^2) + xt*(1 + yk^2 - yS^2)) + 
     Qj*((-1 + xt)*xt*yk^2 - xs*yS^2 + xs^2*yS^2 + xs*xt*(-1 + yk^2 + yS^2)))*
    \[Lambda]Rik*\[Lambda]Rjk*CC[1])/(xs*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^2*(Qj^2*(xs^2 + (-1 + xt)*xt + xs*(-1 + 2*xt))*yS^2 + 
     Qk^2*xs*(-((-1 + xs)*(yk^2 - yS^2)) + xt*(1 - yk^2 + yS^2)) + 
     Qj*Qk*(-((-1 + xt)*xt*yS^2) + xs^2*(yk^2 - 2*yS^2) + 
       xs*(-yk^2 + 2*yS^2 + xt*(-1 + yk^2 - 3*yS^2))))*\[Lambda]Rik*
    \[Lambda]Rjk*CC[2])/(xs*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^2*Qk*(Qk*xt*((-1 + xt)*(yk^2 - yS^2) + xs*(1 + yk^2 - yS^2)) + 
     Qj*(-(xs*yk^2) + xs^2*yk^2 + (-1 + xt)*xt*yS^2 + 
       xs*xt*(-1 + yk^2 + yS^2)))*\[Lambda]Rik*\[Lambda]Rjk*CC[3])/
   (xs^2*xt*(-1 + xs + xt)^3) - 
  (2*mi^2*(Qj^2*(xs^2 + (-1 + xt)*xt + xs*(-1 + 2*xt))*yS^2 + 
     Qk^2*xt*(-((-1 + xt)*(yk^2 - yS^2)) + xs*(1 - yk^2 + yS^2)) + 
     Qj*Qk*(-(xs^2*yS^2) + (-1 + xt)*xt*(yk^2 - 2*yS^2) + 
       xs*(yS^2 + xt*(-1 + yk^2 - 3*yS^2))))*\[Lambda]Rik*\[Lambda]Rjk*CC[4])/
   (xs^2*xt*(-1 + xs + xt)^3) + (4*mi^2*Qk^2*(xs - xt)*yk^2*\[Lambda]Rik*
    \[Lambda]Rjk*CC[5])/(xs^2*xt^2*(-1 + xs + xt)) - 
  (4*mi^2*(Qj - Qk)^2*(xs - xt)*yS^2*\[Lambda]Rik*\[Lambda]Rjk*CC[6])/
   (xs^2*xt^2*(-1 + xs + xt)) + 
  (2*mi^2*(Qj - Qk)*(Qj*(xs^4 - (-1 + xt)*xt + xs^3*(-3 + 4*xt) + 
       xs^2*(3 - 7*xt + 5*xt^2) + xs*(-1 + 2*xt - 4*xt^2 + 2*xt^3))*yS^2 - 
     Qk*(-1 + xs)^2*xt*(xs*(-1 + yk^2 - yS^2) + (-1 + xt)*(yk^2 - yS^2)))*
    \[Lambda]Rik*\[Lambda]Rjk*CC[7])/((-1 + xs)*xs^2*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^2*Qk*(Qk*(-1 + xs)^2*xt*((-1 + xt)*(yk^2 - yS^2) + 
       xs*(1 + yk^2 - yS^2)) - Qj*(xs^4*yk^2 - (-1 + xt)*xt*yS^2 + 
       xs^3*(xt - 3*yk^2 + 5*xt*yk^2 - xt*yS^2) + 
       xs*(-yk^2 + 2*xt^3*yk^2 + xt*(1 + 5*yk^2 - 3*yS^2) + 
         xt^2*(-6*yk^2 + 2*yS^2)) + xs^2*(3*yk^2 + xt^2*(6*yk^2 - yS^2) + 
         xt*(-2 - 10*yk^2 + 3*yS^2))))*\[Lambda]Rik*\[Lambda]Rjk*CC[8])/
   ((-1 + xs)*xs^2*xt^2*(-1 + xs + xt)^3) - 
  (2*mi^2*(Qj - Qk)*(Qj*(2*xs^3*xt + (-1 + xt)^3*xt + 
       xs*(-1 + xt)^2*(1 + 4*xt) + xs^2*(-1 - 4*xt + 5*xt^2))*yS^2 - 
     Qk*xs*(-1 + xt)^2*(xt*(-1 + yk^2 - yS^2) + (-1 + xs)*(yk^2 - yS^2)))*
    \[Lambda]Rik*\[Lambda]Rjk*CC[9])/(xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^2*Qk*(-(Qk*xs*(-1 + xt)^2*((-1 + xs)*(yk^2 - yS^2) + 
        xt*(1 + yk^2 - yS^2))) + Qj*(2*xs^3*xt*yk^2 + (-1 + xt)^3*xt*yk^2 + 
       xs*(-1 + xt)^2*(xt + 5*xt*yk^2 + yS^2 - xt*yS^2) + 
       xs^2*(-1 + xt)*(6*xt*yk^2 + yS^2 - xt*yS^2)))*\[Lambda]Rik*
    \[Lambda]Rjk*CC[10])/(xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^4*(Qj - Qk)*Qk*(-((-1 + xt)^2*xt*yk^2*(yk^2 - yS^2)) + 
     xs^3*yS^2*(xt + yk^2 - yS^2) + xs^2*(-2*yk^2*yS^2 + 2*yS^4 + 
       xt^2*(-1 + yk^2 + yS^2) + xt*(yk^2 - yk^4 + 3*yk^2*yS^2 - 2*yS^4)) + 
     xs*(xt^3*yk^2 + yS^2*(yk^2 - yS^2) + xt^2*(-2*yk^4 + yS^2 + 
         3*yk^2*yS^2 - yS^4) + xt*(2*yk^4 + yS^2*(-1 + 2*yS^2) - 
         yk^2*(1 + 4*yS^2))))*\[Lambda]Rik*\[Lambda]Rjk*DD[1])/
   (xs^2*xt^2*(-1 + xs + xt)^3) - 
  (2*mi^4*(Qj - Qk)*Qk*((-1 + xt)^2*xt*yS^2*(yk^2 - yS^2) + 
     xs^3*yk^2*(xt - yk^2 + yS^2) + xs^2*(2*yk^2*(yk^2 - yS^2) + 
       xt^2*(-1 + yk^2 + yS^2) + xt*(-2*yk^4 + yS^2 + 3*yk^2*yS^2 - yS^4)) + 
     xs*(-yk^4 + xt^3*yS^2 + yk^2*yS^2 + xt^2*(-yk^4 - 2*yS^4 + 
         yk^2*(1 + 3*yS^2)) + xt*(2*yk^4 + yS^2*(-1 + 2*yS^2) - 
         yk^2*(1 + 4*yS^2))))*\[Lambda]Rik*\[Lambda]Rjk*DD[2])/
   (xs^2*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^4*(Qj - Qk)^2*yS^2*(xs^2 + xt*(yk^2 - yS^2) + xs*(xt - yk^2 + yS^2))*
    \[Lambda]Rik*\[Lambda]Rjk*DD[3])/(xs^2*xt^2*(-1 + xs + xt)) - 
  (2*mi^4*Qk^2*yk^2*(xs^2 + xs*(xt + yk^2 - yS^2) + xt*(-yk^2 + yS^2))*
    \[Lambda]Rik*\[Lambda]Rjk*DD[4])/(xs^2*xt^2*(-1 + xs + xt)) - 
  (2*mi^4*(Qj - Qk)^2*yS^2*(xs*(xt + yk^2 - yS^2) + xt*(xt - yk^2 + yS^2))*
    \[Lambda]Rik*\[Lambda]Rjk*DD[5])/(xs^2*xt^2*(-1 + xs + xt)) + 
  (2*mi^4*Qk^2*yk^2*(xt*(xt + yk^2 - yS^2) + xs*(xt - yk^2 + yS^2))*
    \[Lambda]Rik*\[Lambda]Rjk*DD[6])/(xs^2*xt^2*(-1 + xs + xt)), 
 (Qj*(2*Qk*\[Lambda]Rik + Qj*(-1 - 3*yk^2 + 3*yS^2)*\[Lambda]Rik + 
     Qj*xs*(-2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - yS^2*\[Lambda]Rik))*
    \[Lambda]Rjk)/((-1 + xs)*xs^2*xt) + 
  (Qj^2*yk^2*(3*(-yk^2 + yS^2)*\[Lambda]Rik + 
     xs*(-2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - yS^2*\[Lambda]Rik))*
    \[Lambda]Rjk*BB[1])/((-1 + xs)*xs^2*xt*(yk^2 - yS^2)) + 
  (Qj^2*yS^2*(3*(yk^2 - yS^2)*\[Lambda]Rik + 
     xs*(2*yk*\[Lambda]Lik - yk^2*\[Lambda]Rik + yS^2*\[Lambda]Rik))*
    \[Lambda]Rjk*BB[2])/((-1 + xs)*xs^2*xt*(yk^2 - yS^2)) - 
  (Qj*(-2*Qk*(xt^2 + xs^2*(1 - 2*xt^2) + 2*xs*xt*(-1 + xt^2))*\[Lambda]Rik + 
     Qj*(xs*(-1 + xt)*xt*(2*(-1 + xt)*yk*\[Lambda]Lik + 
         xt*yk^2*\[Lambda]Rik + (2 + xt - xt*yS^2)*\[Lambda]Rik) + 
       xt^2*(-2*(-1 + xt)*yk*\[Lambda]Lik + \[Lambda]Rik + 
         (-1 + 2*xt)*yS^2*\[Lambda]Rik + yk^2*(\[Lambda]Rik - 
           2*xt*\[Lambda]Rik)) + xs^2*((1 - yk^2 + yS^2)*\[Lambda]Rik + 
         xt*(-2*yk*\[Lambda]Lik + \[Lambda]Rik + yk^2*\[Lambda]Rik - 
           yS^2*\[Lambda]Rik) + xt^2*(2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - 
           (3 + yS^2)*\[Lambda]Rik))))*\[Lambda]Rjk*BB[3])/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2*(xs + xt)) + 
  (Qj*(2*Qk*xs*(xs - 2*xt)*\[Lambda]Rik - Qj*xs*(xs - yk^2 + yS^2)*
      \[Lambda]Rik + Qj*xt*(2*(-1 + 2*xs)*yk*\[Lambda]Lik + 
       yk^2*\[Lambda]Rik + (xs - yS^2)*\[Lambda]Rik))*\[Lambda]Rjk*BB[4])/
   ((-1 + xs)*xs^2*xt^2) + (Qj*(2*Qk*xt - Qj*(xt - yk^2 + yS^2))*\[Lambda]Rik*
    \[Lambda]Rjk*BB[5])/(xs^2*(-1 + xt)*xt) - 
  (Qk^2*(xs - xt)^2*\[Lambda]Rik*\[Lambda]Rjk*BB[6])/(xs^2*xt^2*(xs + xt)) + 
  ((Qj - Qk)^2*(xs - xt)^2*\[Lambda]Rik*\[Lambda]Rjk*BB[7])/
   (xs^2*xt^2*(xs + xt)) - 
  (mi^2*Qk*(Qk*(-1 + xs)*xs*(xs^4*\[Lambda]Rik + 2*(-1 + xt)^2*xt*
        (yk^2 - yS^2)*\[Lambda]Rik + 
       xs^3*(2*(-1 + yk^2 - yS^2)*\[Lambda]Rik + 
         xt*(4*yk*\[Lambda]Lik + \[Lambda]Rik)) + 
       xs*(2*(yk^2 - yS^2)*\[Lambda]Rik - xt*(-4*yk*\[Lambda]Lik + 
           \[Lambda]Rik + 6*yk^2*\[Lambda]Rik - 6*yS^2*\[Lambda]Rik) + 
         2*xt^2*(-2*yk*\[Lambda]Lik + \[Lambda]Rik + 2*yk^2*\[Lambda]Rik - 
           2*yS^2*\[Lambda]Rik)) + 
       xs^2*(xt^2*(4*yk*\[Lambda]Lik - 2*\[Lambda]Rik) + \[Lambda]Rik - 
         4*yk^2*\[Lambda]Rik + 4*yS^2*\[Lambda]Rik + 
         4*xt*(-2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - 
           yS^2*\[Lambda]Rik))) - 
     Qj*((-1 + xt)^2*xt*(yk^2 - yS^2)^2*\[Lambda]Rik + 
       xs^4*(xt^2*(4*yk*\[Lambda]Lik - \[Lambda]Rik) - (yk^2 - yS^2)^2*
          \[Lambda]Rik + xt*(-4*yk^3*\[Lambda]Lik + 4*yk*yS^2*\[Lambda]Lik + 
           2*yS^2*\[Lambda]Rik)) + xs*(-1 + xt)*
        (-((yk^2 - yS^2)^2*\[Lambda]Rik) + 4*xt^3*yk*(-\[Lambda]Lik + 
           yk*\[Lambda]Rik) + 4*xt*(yk^2 - yS^2)*(-(yk*\[Lambda]Lik) + 
           yk^2*\[Lambda]Rik - yS^2*\[Lambda]Rik) - 
         xt^2*(-4*yk^3*\[Lambda]Lik + 4*yk*(-1 + yS^2)*\[Lambda]Lik + 
           yk^4*\[Lambda]Rik + yS^4*\[Lambda]Rik - 2*yk^2*(-1 + yS^2)*
            \[Lambda]Rik)) + xs^3*(3*(yk^2 - yS^2)^2*\[Lambda]Rik + 
         xt^3*(12*yk*\[Lambda]Lik + \[Lambda]Rik) + 
         xt^2*(-8*yk^3*\[Lambda]Lik + 4*yk*(-3 + 2*yS^2)*\[Lambda]Lik + 
           \[Lambda]Rik + 2*yk^2*\[Lambda]Rik + 2*yS^2*\[Lambda]Rik) - 
         xt*(-12*yk^3*\[Lambda]Lik + 12*yk*yS^2*\[Lambda]Lik + 
           3*yk^4*\[Lambda]Rik - 6*yk^2*yS^2*\[Lambda]Rik + 
           yS^2*(4 + 3*yS^2)*\[Lambda]Rik)) + 
       xs^2*(8*xt^4*yk*\[Lambda]Lik - 3*(yk^2 - yS^2)^2*\[Lambda]Rik - 
         xt^3*(4*yk^3*\[Lambda]Lik - 4*yk*(-5 + yS^2)*\[Lambda]Lik + 
           \[Lambda]Rik - 6*yk^2*\[Lambda]Rik) - 
         xt^2*(-16*yk^3*\[Lambda]Lik + 4*yk*(-3 + 4*yS^2)*\[Lambda]Lik + 
           3*yk^4*\[Lambda]Rik + yk^2*(4 - 6*yS^2)*\[Lambda]Rik + 
           yS^2*(2 + 3*yS^2)*\[Lambda]Rik) + xt*(-12*yk^3*\[Lambda]Lik + 
           12*yk*yS^2*\[Lambda]Lik + 7*yk^4*\[Lambda]Rik - 
           14*yk^2*yS^2*\[Lambda]Rik + yS^2*(2 + 7*yS^2)*\[Lambda]Rik))))*
    \[Lambda]Rjk*CC[1])/(2*(-1 + xs)*xs^2*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*(Qj - Qk)*(-(Qk*(-1 + xs)*xs*(xs^4*\[Lambda]Rik - 
        2*(-1 + xt)^2*xt*(yk^2 - yS^2)*\[Lambda]Rik + 
        xs^2*(8*xt*yk*\[Lambda]Lik - 4*xt^2*yk*\[Lambda]Lik + \[Lambda]Rik + 
          4*yk^2*\[Lambda]Rik - 4*xt*yk^2*\[Lambda]Rik - 
          4*yS^2*\[Lambda]Rik + 4*xt*yS^2*\[Lambda]Rik) + 
        xs^3*(-2*(1 + yk^2 - yS^2)*\[Lambda]Rik + 
          xt*(-4*yk*\[Lambda]Lik + \[Lambda]Rik)) + 
        xs*(2*xt^3*\[Lambda]Rik + 2*(-yk^2 + yS^2)*\[Lambda]Rik + 
          xt*(-4*yk*\[Lambda]Lik - \[Lambda]Rik + 6*yk^2*\[Lambda]Rik - 
            6*yS^2*\[Lambda]Rik) + xt^2*(4*yk*\[Lambda]Lik - 
            4*yk^2*\[Lambda]Rik + 4*yS^2*\[Lambda]Rik)))) + 
     Qj*(-1 + xs + xt)*(xs^5*\[Lambda]Rik - (-1 + xt)*xt*(yk^2 - yS^2)^2*
        \[Lambda]Rik - 2*xs^4*(2*xt*yk*\[Lambda]Lik + (1 + yk^2 - yS^2)*
          \[Lambda]Rik) + xs^3*((1 + yk^4 - 4*yS^2 + yS^4 - 
           2*yk^2*(-2 + yS^2))*\[Lambda]Rik - xt^2*(4*yk*\[Lambda]Lik + 
           \[Lambda]Rik) + xt*(4*yk^3*\[Lambda]Lik - 4*yk*(-2 + yS^2)*
            \[Lambda]Lik + \[Lambda]Rik - 2*yk^2*\[Lambda]Rik)) + 
       xs*(-4*xt^3*yS^2*\[Lambda]Rik + (yk^2 - yS^2)^2*\[Lambda]Rik - 
         xt*(yk^2 - yS^2)*(-4*yk*\[Lambda]Lik + 3*yk^2*\[Lambda]Rik - 
           3*yS^2*\[Lambda]Rik) + xt^2*(-4*yk^3*\[Lambda]Lik + 
           4*yk*yS^2*\[Lambda]Lik + yk^4*\[Lambda]Rik - 
           2*yk^2*yS^2*\[Lambda]Rik + yS^2*(2 + yS^2)*\[Lambda]Rik)) + 
       xs^2*(-2*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4)*\[Lambda]Rik + 
         xt^2*(4*yk^3*\[Lambda]Lik - 4*yk*(-1 + yS^2)*\[Lambda]Lik + 
           \[Lambda]Rik - 2*yS^2*\[Lambda]Rik) + 
         xt*(-8*yk^3*\[Lambda]Lik + 4*yk*(-1 + 2*yS^2)*\[Lambda]Lik + 
           2*yk^4*\[Lambda]Rik + yk^2*(2 - 4*yS^2)*\[Lambda]Rik + 
           (-1 + 2*yS^4)*\[Lambda]Rik))))*\[Lambda]Rjk*CC[2])/
   (2*(-1 + xs)*xs^2*xt^3*(-1 + xs + xt)^2) - 
  (mi^2*Qk*(Qk*(-1 + xt)*xt^2*(4*xs^2*yk*\[Lambda]Lik + 
       (xt^2 - 2*yk^2 + 2*yS^2 + xt*(-1 + 2*yk^2 - 2*yS^2))*\[Lambda]Rik + 
       xs*(4*(-1 + xt)*yk*\[Lambda]Lik + 2*yk^2*\[Lambda]Rik + 
         (-1 + 3*xt - 2*yS^2)*\[Lambda]Rik)) + 
     Qj*((-1 + xt)^2*xt*(yk^2 - yS^2)^2*\[Lambda]Rik + 
       xs*(-1 + xt)*(-((yk^2 - yS^2)^2*\[Lambda]Rik) + xt*(yk^2 - yS^2)*
          (-4*yk*\[Lambda]Lik - 2*\[Lambda]Rik + 3*yk^2*\[Lambda]Rik - 
           3*yS^2*\[Lambda]Rik) + 2*xt^2*(2*yk^3*\[Lambda]Lik - 
           2*yk*yS^2*\[Lambda]Lik + yk^2*\[Lambda]Rik - 
           2*yS^2*\[Lambda]Rik)) + xs^3*((yk^2 - yS^2)^2*\[Lambda]Rik + 
         xt^2*(4*yk*\[Lambda]Lik + \[Lambda]Rik) + 
         xt*(4*yk^3*\[Lambda]Lik - 4*yk*yS^2*\[Lambda]Lik + 
           4*yk^2*\[Lambda]Rik - 2*yS^2*\[Lambda]Rik)) + 
       xs^2*(-2*(yk^2 - yS^2)^2*\[Lambda]Rik + xt^3*(4*yk*\[Lambda]Lik + 
           3*\[Lambda]Rik) + xt*(-8*yk^3*\[Lambda]Lik + 
           8*yk*yS^2*\[Lambda]Lik + 3*yk^4*\[Lambda]Rik - 6*yk^2*(1 + yS^2)*
            \[Lambda]Rik + yS^2*(4 + 3*yS^2)*\[Lambda]Rik) + 
         xt^2*(8*yk^3*\[Lambda]Lik - 4*yk*(\[Lambda]Lik + 
             2*yS^2*\[Lambda]Lik) + 6*yk^2*\[Lambda]Rik - 
           2*(\[Lambda]Rik + 3*yS^2*\[Lambda]Rik)))))*\[Lambda]Rjk*CC[3])/
   (2*xs^3*xt^2*(-1 + xs + xt)^2) - 
  (mi^2*(Qj - Qk)*(Qk*(-1 + xt)*xt^2*
      ((xt^2 + 2*(yk^2 - yS^2) + xt*(-1 - 2*yk^2 + 2*yS^2))*\[Lambda]Rik - 
       2*xs^2*(2*yk*\[Lambda]Lik + \[Lambda]Rik) + 
       xs*(-4*(-1 + xt)*yk*\[Lambda]Lik - 2*yk^2*\[Lambda]Rik + 
         (1 + xt + 2*yS^2)*\[Lambda]Rik)) + Qj*(-1 + xs + xt)*
      (-((-1 + xt)*xt*(xt - yk^2 + yS^2)^2*\[Lambda]Rik) + 
       xs*(4*xt^3*yk*\[Lambda]Lik + (yk^2 - yS^2)^2*\[Lambda]Rik - 
         xt^2*(4*yk^3*\[Lambda]Lik - 4*yk*(-1 + yS^2)*\[Lambda]Lik + 
           \[Lambda]Rik - 2*yk^2*\[Lambda]Rik) - 2*xt*(yk^2 - yS^2)*
          (-2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - yS^2*\[Lambda]Rik)) + 
       xs^2*(-((yk^2 - yS^2)^2*\[Lambda]Rik) + xt^2*(4*yk*\[Lambda]Lik + 
           \[Lambda]Rik) - 2*xt*(2*yk^3*\[Lambda]Lik - 
           2*yk*yS^2*\[Lambda]Lik + yS^2*\[Lambda]Rik))))*\[Lambda]Rjk*CC[4])/
   (2*xs^3*xt^2*(-1 + xs + xt)^2) + 
  (mi^2*Qk^2*(xs^4*\[Lambda]Rik + xt*(xt^3 + xt^2*(-1 + 2*yk^2 - 2*yS^2) - 
       2*(yk^2 - yS^2)^2 + 2*xt*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2)))*
      \[Lambda]Rik + xs^3*(4*xt*yk*\[Lambda]Lik + (-1 + 2*yk^2 - 2*yS^2)*
        \[Lambda]Rik) + xs^2*(8*xt^2*yk*\[Lambda]Lik + 
       2*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2))*\[Lambda]Rik + 
       xt*(8*yk^3*\[Lambda]Lik - 4*yk*(\[Lambda]Lik + 2*yS^2*\[Lambda]Lik) + 
         \[Lambda]Rik + 2*yk^2*\[Lambda]Rik - 6*yS^2*\[Lambda]Rik)) + 
     xs*(-2*(yk^2 - yS^2)^2*\[Lambda]Rik + 2*xt^3*(2*yk*\[Lambda]Lik + 
         \[Lambda]Rik) + 4*xt*(yk^2 - yS^2)*(-2*yk*\[Lambda]Lik + 
         yk^2*\[Lambda]Rik - (1 + yS^2)*\[Lambda]Rik) + 
       xt^2*(8*yk^3*\[Lambda]Lik - 4*yk*(\[Lambda]Lik + 
           2*yS^2*\[Lambda]Lik) + 2*yk^2*\[Lambda]Rik - 
         (1 + 6*yS^2)*\[Lambda]Rik)))*\[Lambda]Rjk*CC[5])/(2*xs^3*xt^3) - 
  (mi^2*(Qj - Qk)^2*(xs^4*\[Lambda]Rik + xt*(xt^3 - 2*(yk^2 - yS^2)^2 + 
       xt^2*(-1 - 2*yk^2 + 2*yS^2) + 2*xt*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + 
         yS^4))*\[Lambda]Rik - xs^3*(4*xt*yk*\[Lambda]Lik + \[Lambda]Rik + 
       2*yk^2*\[Lambda]Rik - 2*yS^2*\[Lambda]Rik) + 
     xs^2*(2*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4)*\[Lambda]Rik - 
       2*xt^2*(4*yk*\[Lambda]Lik + \[Lambda]Rik) + 
       xt*(8*yk^3*\[Lambda]Lik + yk*(4 - 8*yS^2)*\[Lambda]Lik + 
         \[Lambda]Rik - 2*yk^2*\[Lambda]Rik - 2*yS^2*\[Lambda]Rik)) + 
     xs*(-4*xt^3*yk*\[Lambda]Lik - 2*(yk^2 - yS^2)^2*\[Lambda]Rik + 
       xt^2*(8*yk^3*\[Lambda]Lik + yk*(4 - 8*yS^2)*\[Lambda]Lik + 
         \[Lambda]Rik - 2*yk^2*\[Lambda]Rik - 2*yS^2*\[Lambda]Rik) + 
       4*xt*(yk^2 - yS^2)*(-2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - 
         yS^2*\[Lambda]Rik)))*\[Lambda]Rjk*CC[6])/(2*xs^3*xt^3) + 
  (mi^2*(Qj - Qk)*(-(Qk*(-1 + xs)*xs*(xs^4*\[Lambda]Rik + 
        2*(-1 + xt)*xt^2*(yk^2 - yS^2)*\[Lambda]Rik + 
        xs^2*(8*xt*yk*\[Lambda]Lik - 4*xt^2*yk*\[Lambda]Lik + \[Lambda]Rik + 
          4*yk^2*\[Lambda]Rik - 4*xt*yk^2*\[Lambda]Rik - 
          4*yS^2*\[Lambda]Rik + 4*xt*yS^2*\[Lambda]Rik) + 
        xs^3*(-2*(1 + yk^2 - yS^2)*\[Lambda]Rik + 
          xt*(-4*yk*\[Lambda]Lik + \[Lambda]Rik)) - 
        xs*(2*xt^3*\[Lambda]Rik + 2*(yk^2 - yS^2)*\[Lambda]Rik - 
          2*xt^2*(2*yk*\[Lambda]Lik + \[Lambda]Rik) + 
          xt*(4*yk*\[Lambda]Lik + \[Lambda]Rik - 4*yk^2*\[Lambda]Rik + 
            4*yS^2*\[Lambda]Rik)))) + Qj*(-1 + xs + xt)*
      (xs^5*\[Lambda]Rik - (-1 + xt)*xt*(yk^2 - yS^2)^2*\[Lambda]Rik - 
       2*xs^4*(2*xt*yk*\[Lambda]Lik + (1 + yk^2 - yS^2)*\[Lambda]Rik) + 
       xs^3*((1 + yk^4 - 4*yS^2 + yS^4 - 2*yk^2*(-2 + yS^2))*\[Lambda]Rik - 
         xt^2*(4*yk*\[Lambda]Lik + \[Lambda]Rik) + 
         xt*(4*yk^3*\[Lambda]Lik - 4*yk*(-2 + yS^2)*\[Lambda]Lik + 
           \[Lambda]Rik - 2*yk^2*\[Lambda]Rik)) + 
       xs*(4*xt^3*yS^2*\[Lambda]Rik + (yk^2 - yS^2)^2*\[Lambda]Rik - 
         xt*(yk^2 - yS^2)*(-4*yk*\[Lambda]Lik + 3*yk^2*\[Lambda]Rik - 
           3*yS^2*\[Lambda]Rik) + xt^2*(-4*yk^3*\[Lambda]Lik + 
           4*yk*yS^2*\[Lambda]Lik + yk^4*\[Lambda]Rik - 
           2*yk^2*yS^2*\[Lambda]Rik + yS^2*(-2 + yS^2)*\[Lambda]Rik)) + 
       xs^2*(-2*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4)*\[Lambda]Rik + 
         xt^2*(4*yk^3*\[Lambda]Lik - 4*yk*(-1 + yS^2)*\[Lambda]Lik + 
           \[Lambda]Rik - 2*yS^2*\[Lambda]Rik) + 
         xt*(-8*yk^3*\[Lambda]Lik + 4*yk*(-1 + 2*yS^2)*\[Lambda]Lik + 
           2*yk^4*\[Lambda]Rik + yk^2*(2 - 4*yS^2)*\[Lambda]Rik + 
           (-1 + 2*yS^4)*\[Lambda]Rik))))*\[Lambda]Rjk*CC[7])/
   (2*xs^3*xt^3*(-1 + xs + xt)^2) - 
  (mi^2*Qk*(Qk*(-1 + xs)*xs*(xs^4*\[Lambda]Rik - 2*(-1 + xt)*xt^2*
        (yk^2 - yS^2)*\[Lambda]Rik + 
       xs^3*(2*(-1 + yk^2 - yS^2)*\[Lambda]Rik + 
         xt*(4*yk*\[Lambda]Lik + \[Lambda]Rik)) - 
       xs*(4*xt^2*(yk*\[Lambda]Lik - \[Lambda]Rik) + 4*xt^3*\[Lambda]Rik + 
         2*(-yk^2 + yS^2)*\[Lambda]Rik + xt*(-4*yk*\[Lambda]Lik + 
           \[Lambda]Rik + 4*yk^2*\[Lambda]Rik - 4*yS^2*\[Lambda]Rik)) + 
       xs^2*(xt^2*(4*yk*\[Lambda]Lik - 2*\[Lambda]Rik) + \[Lambda]Rik - 
         4*yk^2*\[Lambda]Rik + 4*yS^2*\[Lambda]Rik + 
         4*xt*(-2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - 
           yS^2*\[Lambda]Rik))) + 
     Qj*(-((-1 + xt)^2*xt*(yk^2 - yS^2)^2*\[Lambda]Rik) + 
       xs^4*((yk^2 - yS^2)^2*\[Lambda]Rik + xt^2*(-4*yk*\[Lambda]Lik + 
           \[Lambda]Rik) + xt*(4*yk^3*\[Lambda]Lik - 4*yk*yS^2*\[Lambda]Lik - 
           2*yS^2*\[Lambda]Rik)) + xs*(-1 + xt)*
        ((yk^2 - yS^2)^2*\[Lambda]Rik + 4*xt^3*yk*(\[Lambda]Lik + 
           yk*\[Lambda]Rik) + xt^2*(-4*yk^3*\[Lambda]Lik + 
           4*yk*(-1 + yS^2)*\[Lambda]Lik + yk^4*\[Lambda]Rik - 
           2*yk^2*(3 + yS^2)*\[Lambda]Rik + yS^2*(4 + yS^2)*\[Lambda]Rik) - 
         2*xt*(yk^2 - yS^2)*(-2*yk*\[Lambda]Lik + 2*yk^2*\[Lambda]Rik - 
           (1 + 2*yS^2)*\[Lambda]Rik)) + xs^2*(-8*xt^4*yk*\[Lambda]Lik + 
         3*(yk^2 - yS^2)^2*\[Lambda]Rik + xt^3*(4*yk^3*\[Lambda]Lik - 
           4*yk*(-5 + yS^2)*\[Lambda]Lik + 6*yk^2*\[Lambda]Rik - 
           (3 + 4*yS^2)*\[Lambda]Rik) + xt*(12*yk^3*\[Lambda]Lik - 
           12*yk*yS^2*\[Lambda]Lik - 7*yk^4*\[Lambda]Rik + 
           2*yk^2*(2 + 7*yS^2)*\[Lambda]Rik - yS^2*(6 + 7*yS^2)*
            \[Lambda]Rik) + xt^2*(-16*yk^3*\[Lambda]Lik + 4*yk*(-3 + 4*yS^2)*
            \[Lambda]Lik + 3*yk^4*\[Lambda]Rik - 2*yk^2*(5 + 3*yS^2)*
            \[Lambda]Rik + (2 + 12*yS^2 + 3*yS^4)*\[Lambda]Rik)) + 
       xs^3*(-3*(yk^2 - yS^2)^2*\[Lambda]Rik + 3*xt^3*(-4*yk*\[Lambda]Lik + 
           \[Lambda]Rik) + xt^2*(8*yk^3*\[Lambda]Lik + 4*yk*(3 - 2*yS^2)*
            \[Lambda]Lik + 2*yk^2*\[Lambda]Rik - 
           3*(\[Lambda]Rik + 2*yS^2*\[Lambda]Rik)) + 
         xt*(-12*yk^3*\[Lambda]Lik + 12*yk*yS^2*\[Lambda]Lik + 
           3*yk^4*\[Lambda]Rik + 3*yS^2*(2 + yS^2)*\[Lambda]Rik - 
           2*yk^2*(\[Lambda]Rik + 3*yS^2*\[Lambda]Rik)))))*\[Lambda]Rjk*
    CC[8])/(2*xs^3*xt^3*(-1 + xs + xt)^2) - 
  (mi^2*(Qj - Qk)*(Qk*(-1 + xt)*xt*((-1 + xt)^2*xt*(xt - 2*yk^2 + 2*yS^2)*
        \[Lambda]Rik - xs*(-1 + xt)*(xt^2*(4*yk*\[Lambda]Lik - 
           \[Lambda]Rik) + 2*(-yk^2 + yS^2)*\[Lambda]Rik + 
         xt*(-4*yk*\[Lambda]Lik - \[Lambda]Rik + 6*yk^2*\[Lambda]Rik - 
           6*yS^2*\[Lambda]Rik)) + xs^2*(2*(yk^2 - yS^2)*\[Lambda]Rik + 
         xt^2*(-4*yk*\[Lambda]Lik + 2*\[Lambda]Rik) + 
         xt*(4*yk*\[Lambda]Lik - 4*yk^2*\[Lambda]Rik + 
           4*yS^2*\[Lambda]Rik))) + Qj*(-1 + xs + xt)*
      (-(xt*(xt^2 + yk^2 - yS^2 + xt*(-1 - yk^2 + yS^2))^2*\[Lambda]Rik) + 
       xs*(-1 + xt)*(4*xt^3*yk*\[Lambda]Lik + (yk^2 - yS^2)^2*\[Lambda]Rik - 
         xt^2*(4*yk^3*\[Lambda]Lik - 4*yk*(-1 + yS^2)*\[Lambda]Lik + 
           \[Lambda]Rik - 2*yk^2*\[Lambda]Rik) - 2*xt*(yk^2 - yS^2)*
          (-2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - yS^2*\[Lambda]Rik)) + 
       xs^2*((yk^2 - yS^2)^2*\[Lambda]Rik + xt^3*(4*yk*\[Lambda]Lik + 
           \[Lambda]Rik) - xt^2*(4*yk^3*\[Lambda]Lik - 4*yk*(-1 + yS^2)*
            \[Lambda]Lik + \[Lambda]Rik - 6*yS^2*\[Lambda]Rik) - 
         xt*(-4*yk^3*\[Lambda]Lik + 4*yk*yS^2*\[Lambda]Lik + 
           yk^4*\[Lambda]Rik - 2*yk^2*yS^2*\[Lambda]Rik + 
           yS^2*(2 + yS^2)*\[Lambda]Rik))))*\[Lambda]Rjk*CC[9])/
   (2*xs^3*xt^3*(-1 + xs + xt)^2) - 
  (mi^2*Qk*(Qk*(-1 + xt)*xt*((-1 + xt)^2*xt*(xt + 2*yk^2 - 2*yS^2)*
        \[Lambda]Rik + xs*(-1 + xt)*(2*(-yk^2 + yS^2)*\[Lambda]Rik + 
         xt^2*(4*yk*\[Lambda]Lik + 3*\[Lambda]Rik) + 
         xt*(-4*yk*\[Lambda]Lik - \[Lambda]Rik + 6*yk^2*\[Lambda]Rik - 
           6*yS^2*\[Lambda]Rik)) + 2*xs^2*((-yk^2 + yS^2)*\[Lambda]Rik + 
         2*xt^2*(yk*\[Lambda]Lik + \[Lambda]Rik) + 
         xt*(-2*yk*\[Lambda]Lik - \[Lambda]Rik + 2*yk^2*\[Lambda]Rik - 
           2*yS^2*\[Lambda]Rik))) + 
     Qj*((-1 + xt)^3*xt*(yk^2 - yS^2)^2*\[Lambda]Rik + 
       xs*(-1 + xt)^2*(-((yk^2 - yS^2)^2*\[Lambda]Rik) + 
         2*xt^2*yk*(2*yk^2*\[Lambda]Lik - 2*yS^2*\[Lambda]Lik - 
           yk*\[Lambda]Rik) + xt*(yk^2 - yS^2)*(-4*yk*\[Lambda]Lik + 
           3*yk^2*\[Lambda]Rik - 3*yS^2*\[Lambda]Rik)) + 
       xs^3*(-((yk^2 - yS^2)^2*\[Lambda]Rik) + xt^3*(4*yk*\[Lambda]Lik + 
           \[Lambda]Rik) + xt*(-4*yk^3*\[Lambda]Lik + 
           4*yk*yS^2*\[Lambda]Lik + yk^4*\[Lambda]Rik - 
           2*yk^2*yS^2*\[Lambda]Rik + yS^2*(2 + yS^2)*\[Lambda]Rik) + 
         xt^2*(4*yk^3*\[Lambda]Lik - 4*yk*(1 + yS^2)*\[Lambda]Lik - 
           4*yk^2*\[Lambda]Rik - (1 + 2*yS^2)*\[Lambda]Rik)) + 
       xs^2*(-1 + xt)*(xt^3*(4*yk*\[Lambda]Lik - \[Lambda]Rik) - 
         2*(yk^2 - yS^2)^2*\[Lambda]Rik + xt^2*(8*yk^3*\[Lambda]Lik - 
           4*yk*(\[Lambda]Lik + 2*yS^2*\[Lambda]Lik) - 6*yk^2*\[Lambda]Rik - 
           2*yS^2*\[Lambda]Rik) + xt*(-8*yk^3*\[Lambda]Lik + 
           8*yk*yS^2*\[Lambda]Lik + 3*yk^4*\[Lambda]Rik - 
           6*yk^2*yS^2*\[Lambda]Rik + yS^2*(2 + 3*yS^2)*\[Lambda]Rik))))*
    \[Lambda]Rjk*CC[10])/(2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*Qk^2*(xs^5*\[Lambda]Rik + xt^4*(xt + 2*yk^2 - 2*yS^2)*\[Lambda]Rik + 
     xs*xt^3*(4*xt*yk*\[Lambda]Lik + 3*xt*\[Lambda]Rik + 
       4*yk^2*\[Lambda]Rik - 4*yS^2*\[Lambda]Rik) + 
     4*xs^3*xt*(xt*yk*\[Lambda]Lik + (yk^2 - yS^2)*\[Lambda]Rik) + 
     2*xs^2*xt^2*(2*(-yk^2 + yS^2)*\[Lambda]Rik + 
       xt*(2*yk*\[Lambda]Lik + \[Lambda]Rik)) + 
     xs^4*(2*(yk^2 - yS^2)*\[Lambda]Rik + xt*(4*yk*\[Lambda]Lik + 
         \[Lambda]Rik)))*\[Lambda]Rjk*CC[11])/(2*xs^3*xt^3*(xs + xt)) - 
  (mi^2*(Qj - Qk)^2*(xs^5*\[Lambda]Rik + xt^4*(xt - 2*yk^2 + 2*yS^2)*
      \[Lambda]Rik + xs^4*(2*(-yk^2 + yS^2)*\[Lambda]Rik + 
       xt*(-4*yk*\[Lambda]Lik + \[Lambda]Rik)) + 
     xs*xt^3*(4*(-yk^2 + yS^2)*\[Lambda]Rik + 
       xt*(-4*yk*\[Lambda]Lik + \[Lambda]Rik)) - 
     2*xs^3*xt*(2*(yk^2 - yS^2)*\[Lambda]Rik + 
       xt*(2*yk*\[Lambda]Lik + \[Lambda]Rik)) - 
     2*xs^2*xt^2*(2*(-yk^2 + yS^2)*\[Lambda]Rik + 
       xt*(2*yk*\[Lambda]Lik + \[Lambda]Rik)))*\[Lambda]Rjk*CC[12])/
   (2*xs^3*xt^3*(xs + xt)) + 
  (mi^4*(Qj - Qk)*Qk*(xs*(xt - (yk - yS)^2) - (-1 + xt)*(yk - yS)^2)*
    (-((-1 + xt)*(yk + yS)^2) + xs*(xt - (yk + yS)^2))*
    ((-1 + xt)*xt*(yk^2 - yS^2)*\[Lambda]Rik + 
     xs^2*((yk^2 - yS^2)*\[Lambda]Rik + xt*(4*yk*\[Lambda]Lik + 
         \[Lambda]Rik)) + xs*(xt^2*(4*yk*\[Lambda]Lik - \[Lambda]Rik) + 
       (-yk^2 + yS^2)*\[Lambda]Rik + 2*xt*(-2*yk*\[Lambda]Lik + 
         yk^2*\[Lambda]Rik - yS^2*\[Lambda]Rik)))*\[Lambda]Rjk*DD[1])/
   (2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^4*(Qj - Qk)*Qk*((-1 + xt)^2*(yk^2 - yS^2)^2 + 
     xs^2*(xt^2 + (yk^2 - yS^2)^2 - 2*xt*(yk^2 + yS^2)) - 
     2*xs*((yk^2 - yS^2)^2 + xt^2*(yk^2 + yS^2) - 
       xt*(yk^2 + yk^4 + yS^2 - 2*yk^2*yS^2 + yS^4)))*
    ((-1 + xt)*xt*(yk^2 - yS^2)*\[Lambda]Rik + 
     xs^2*((yk^2 - yS^2)*\[Lambda]Rik + xt*(4*yk*\[Lambda]Lik + 
         \[Lambda]Rik)) + xs*((-yk^2 + yS^2)*\[Lambda]Rik + 
       xt^2*(4*yk*\[Lambda]Lik + 3*\[Lambda]Rik) + 
       2*xt*(-2*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - 
         (1 + yS^2)*\[Lambda]Rik)))*\[Lambda]Rjk*DD[2])/
   (2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^4*(Qj - Qk)^2*(xs^3 + (-1 + xt)*(yk^2 - yS^2)^2 + 
     xs^2*(-1 + xt - 2*yk^2 + 2*yS^2) + xs*(yk^4 + yS^2*(-2 - 2*xt + yS^2) - 
       2*yk^2*(-1 + xt + yS^2)))*(xs^2*\[Lambda]Rik + 
     xt*(-yk^2 + yS^2)*\[Lambda]Rik - xs*((yk^2 - yS^2)*\[Lambda]Rik + 
       xt*(4*yk*\[Lambda]Lik + \[Lambda]Rik)))*\[Lambda]Rjk*DD[3])/
   (2*xs^3*xt^3) - (mi^4*Qk^2*(xs^5*\[Lambda]Rik + 
     (-1 + xt)*xt*(yk^2 - yS^2)^3*\[Lambda]Rik + 
     xs^4*(4*xt*yk*\[Lambda]Lik + (-1 + 3*yk^2 - 3*yS^2)*\[Lambda]Rik) + 
     xs*(yk^2 - yS^2)*(-((yk^2 - yS^2)^2*\[Lambda]Rik) + 
       xt^2*(4*yk^3*\[Lambda]Lik - 4*yk*yS^2*\[Lambda]Lik - 
         yk^2*\[Lambda]Rik - 3*yS^2*\[Lambda]Rik) + xt*(yk^2 - yS^2)*
        (-4*yk*\[Lambda]Lik + 2*yk^2*\[Lambda]Rik - (3 + 2*yS^2)*
          \[Lambda]Rik)) - xs^3*(-3*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2))*
        \[Lambda]Rik + xt^2*(4*yk*\[Lambda]Lik + \[Lambda]Rik) + 
       xt*(-8*yk^3*\[Lambda]Lik + 4*yk*(\[Lambda]Lik + 2*yS^2*\[Lambda]Lik) + 
         (-1 + 4*yS^2)*\[Lambda]Rik)) - 
     xs^2*(-((-3 + yk^2 - yS^2)*(yk^2 - yS^2)^2*\[Lambda]Rik) + 
       xt^2*(8*yk^3*\[Lambda]Lik + 4*yk*(-1 + 2*yS^2)*\[Lambda]Lik + 
         3*yk^2*\[Lambda]Rik + yS^2*\[Lambda]Rik) - xt*(yk^2 - yS^2)*
        (4*yk^3*\[Lambda]Lik - 4*yk*(2 + yS^2)*\[Lambda]Lik + 
         2*yk^2*\[Lambda]Rik - (1 + 6*yS^2)*\[Lambda]Rik)))*\[Lambda]Rjk*
    DD[4])/(2*xs^3*xt^3) + 
  (mi^4*(Qj - Qk)^2*(xt^3 + (-1 + xs)*(yk^2 - yS^2)^2 + 
     xt^2*(-1 + xs - 2*yk^2 + 2*yS^2) + xt*(yk^4 + yS^2*(-2 - 2*xs + yS^2) - 
       2*yk^2*(-1 + xs + yS^2)))*(xt*(xt - yk^2 + yS^2)*\[Lambda]Rik - 
     xs*((yk^2 - yS^2)*\[Lambda]Rik + xt*(4*yk*\[Lambda]Lik + \[Lambda]Rik)))*
    \[Lambda]Rjk*DD[5])/(2*xs^3*xt^3) - 
  (mi^4*Qk^2*(xt^3 + xt^2*(-1 + xs + 2*yk^2 - 2*yS^2) + 
     (-1 + xs)*(yk^2 - yS^2)^2 + xt*(yk^4 + yS^2*(2 - 2*xs + yS^2) - 
       2*yk^2*(1 + xs + yS^2)))*(xs*(yk^2 - yS^2)*\[Lambda]Rik + 
     xt*(xt + yk^2 - yS^2)*\[Lambda]Rik + 
     xs*xt*(4*yk*\[Lambda]Lik + \[Lambda]Rik))*\[Lambda]Rjk*DD[6])/
   (2*xs^3*xt^3), (2*Qj*\[Lambda]Ljk*(2*Qk*(-1 + 2*xs)*\[Lambda]Lik + 
     Qj*(1 - yk^2 + yS^2 - 2*xs*(1 + yk^2 - yS^2))*\[Lambda]Lik - 
     4*Qj*xs*yk*\[Lambda]Rik))/((-1 + xs)*xs^2*xt) - 
  (2*Qj^2*yk^2*\[Lambda]Ljk*(-((1 + 2*xs)*yS^2*\[Lambda]Lik) + 
     yk^2*(\[Lambda]Lik + 2*xs*\[Lambda]Lik) + 4*xs*yk*\[Lambda]Rik)*BB[1])/
   ((-1 + xs)*xs^2*xt*(yk^2 - yS^2)) + 
  (2*Qj^2*yS^2*\[Lambda]Ljk*(-((1 + 2*xs)*yS^2*\[Lambda]Lik) + 
     yk^2*(\[Lambda]Lik + 2*xs*\[Lambda]Lik) + 4*xs*yk*\[Lambda]Rik)*BB[2])/
   ((-1 + xs)*xs^2*xt*(yk^2 - yS^2)) + 
  (2*Qj*\[Lambda]Ljk*(-2*Qk*(xs^3 + 2*xs^2*(-2 + xt) + 2*(-1 + xt) - 
       xs*(-5 + 3*xt + xt^2))*\[Lambda]Lik + Qj*(-1 + xs + xt)*
      ((2 + 2*(-2 + xt)*yk^2 + 4*yS^2 - 2*xt*yS^2 + 
         2*xs^2*(-2*yk^2 + 2*yS^2 + xt*(1 + yk^2 - yS^2)) - 
         xs*(1 - 7*yk^2 + 7*yS^2 + 3*xt*(1 + yk^2 - yS^2)))*\[Lambda]Lik + 
       2*(1 - 3*xs + 2*xs^2)*(-1 + xt)*yk*\[Lambda]Rik))*BB[3])/
   ((-1 + xs)^2*xs^2*(-1 + xt)*xt*(-1 + xs + xt)) - 
  (2*Qj*\[Lambda]Ljk*(2*Qk*xs*(-2 + 4*xs - 2*xs^2 + xt)*\[Lambda]Lik + 
     Qj*(-1 + xs + xt)*(2*xs^2*\[Lambda]Lik + 3*yk^2*\[Lambda]Lik - 
       3*yS^2*\[Lambda]Lik + 2*yk*\[Lambda]Rik - 
       xs*((3 + 2*yk^2 - 2*yS^2)*\[Lambda]Lik + 2*yk*\[Lambda]Rik)))*BB[4])/
   ((-1 + xs)^2*xs^2*xt*(-1 + xs + xt)) - 
  (4*Qj*(-(Qk*(xs + 2*(-1 + xt)*xt)) + Qj*(-1 + xs + xt)*(xt - yk^2 + yS^2))*
    \[Lambda]Lik*\[Lambda]Ljk*BB[5])/(xs^2*(-1 + xt)*xt*(-1 + xs + xt)) - 
  (4*Qk^2*\[Lambda]Lik*\[Lambda]Ljk*BB[6])/(xs^2*xt) + 
  (4*(Qj - Qk)^2*\[Lambda]Lik*\[Lambda]Ljk*BB[7])/(xs^2*xt) - 
  (2*mi^2*Qk*\[Lambda]Ljk*(Qk*(-1 + xs)*xs*
      (xs^3 + xs*(-1 + xt)*(-1 + 4*yk^2 - 4*yS^2) + 
       2*(-1 + xt)^2*(yk^2 - yS^2) + 2*xs^2*(-1 + xt + yk^2 - yS^2))*
      \[Lambda]Lik + Qj*(-((-1 + xt)^2*(yk^2 - yS^2)^2*\[Lambda]Lik) + 
       xs^3*(xt^2*\[Lambda]Lik + (yk^2 - yS^2)^2*\[Lambda]Lik - 
         xt*(\[Lambda]Lik + 2*yk^2*\[Lambda]Lik + 3*yk*\[Lambda]Rik)) - 
       xs*(-1 + xt)*(3*(yk^2 - yS^2)^2*\[Lambda]Lik + 
         4*xt^2*yk*(yk*\[Lambda]Lik + \[Lambda]Rik) - 
         xt*(yk^4*\[Lambda]Lik + yS^4*\[Lambda]Lik - 2*yk^2*(-1 + yS^2)*
            \[Lambda]Lik + 3*yk*\[Lambda]Rik)) + 
       xs^2*(-3*(yk^2 - yS^2)^2*\[Lambda]Lik + xt*(1 + 2*yk^4 + 2*yS^4 - 
           4*yk^2*(-1 + yS^2))*\[Lambda]Lik + 6*xt*yk*\[Lambda]Rik - 
         xt^2*(\[Lambda]Lik + 6*yk^2*\[Lambda]Lik + 7*yk*\[Lambda]Rik))))*
    CC[1])/((-1 + xs)*xs^2*xt^2*(-1 + xs + xt)^2) + 
  (2*mi^2*(Qj - Qk)*(-(Qk*(-1 + xs)*xs*(xs^3 - 2*(-1 + xt)^2*(yk^2 - yS^2) + 
        2*xs^2*(-1 + xt - yk^2 + yS^2) + xs*(1 + 2*xt^2 + 4*yk^2 - 4*yS^2 + 
          xt*(-2 - 4*yk^2 + 4*yS^2)))) + 
     Qj*(xs^5 - (-1 + xt)^2*(yk^2 - yS^2)^2 + 
       xs^4*(-3 + 2*xt - 2*yk^2 + 2*yS^2) + xs^3*(3 + xt^2 + yk^4 - 6*yS^2 + 
         yS^4 - 2*yk^2*(-3 + yS^2) + xt*(-4 - 4*yk^2 + 2*yS^2)) + 
       xs*((3 - 4*xt + xt^2)*yk^4 - (-1 + xt)*yS^2*(-2 + 4*xt^2 + 3*yS^2 - 
           xt*yS^2) - 2*(-1 + xt)*yk^2*(1 - 3*yS^2 + xt*(-1 + yS^2))) - 
       xs^2*(1 + 3*yk^4 - 6*yS^2 + 3*yS^4 - 6*yk^2*(-1 + yS^2) + 
         xt^2*(1 + 2*yk^2 + 4*yS^2) - 2*xt*(yk^4 - 2*yk^2*(-2 + yS^2) + 
           (-1 + yS^2)^2))))*\[Lambda]Lik*\[Lambda]Ljk*CC[2])/
   ((-1 + xs)*xs^2*xt^2*(-1 + xs + xt)^2) - 
  (2*mi^2*Qk*\[Lambda]Ljk*(Qj*((-1 + xt)^2*(yk^2 - yS^2)^2 + 
       xs^2*(xt^2 - 2*xt*yS^2 + (yk^2 - yS^2)^2) + 
       2*xs*((-1 + xt)*yk^4 - 2*(-1 + xt)*yk^2*yS^2 + 
         yS^2*(xt - xt^2 - yS^2 + xt*yS^2)))*\[Lambda]Lik + 
     Qk*xt*(xt^3*\[Lambda]Lik + 2*xt^2*(-1 + xs + yk^2 - yS^2)*\[Lambda]Lik + 
       xt*(1 - 4*yk^2 + 4*yS^2 + 2*xs*(-1 + yk^2 - yS^2))*\[Lambda]Lik - 
       xs*xt*yk*\[Lambda]Rik - (-1 + xs)*(2*yk^2*\[Lambda]Lik - 
         2*yS^2*\[Lambda]Lik + xs*yk*\[Lambda]Rik)))*CC[3])/
   (xs^3*xt*(-1 + xs + xt)^2) + (2*mi^2*(Qj - Qk)*\[Lambda]Ljk*
    (Qj*(xt^4 + (-1 + xs)^2*(yk^2 - yS^2)^2 + 
       2*xt^3*(-1 + xs - yk^2 + yS^2) + xt^2*(1 + xs^2 + yk^4 - 4*yS^2 + 
         yS^4 - 2*yk^2*(-2 + yS^2) + xs*(-2 - 4*yk^2 + 2*yS^2)) + 
       2*xt*((-1 + xs)*yk^4 + (-1 + xs)*yS^2*(-1 + yS^2) - 
         yk^2*(1 + xs^2 - 2*yS^2 + 2*xs*(-1 + yS^2))))*\[Lambda]Lik - 
     Qk*xt*((-1 + xt)^2*(xt - 2*yk^2 + 2*yS^2)*\[Lambda]Lik + 
       xs^2*((-1 + 2*xt)*\[Lambda]Lik + yk*\[Lambda]Rik) + 
       xs*(-1 + xt)*(2*xt*\[Lambda]Lik - 2*yk^2*\[Lambda]Lik + 
         2*yS^2*\[Lambda]Lik + yk*\[Lambda]Rik)))*CC[4])/
   (xs^3*xt*(-1 + xs + xt)^2) + 
  (2*mi^2*Qk^2*\[Lambda]Ljk*(xs^3*\[Lambda]Lik + 
     xs^2*(-1 + xt + 2*yk^2 - 2*yS^2)*\[Lambda]Lik + 
     (xt^3 + xt^2*(-1 + 2*yk^2 - 2*yS^2) - 2*(yk^2 - yS^2)^2 + 
       2*xt*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2)))*\[Lambda]Lik + 
     xs*(xt^2*\[Lambda]Lik + 2*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2))*
        \[Lambda]Lik - 2*xt*(2*yS^2*\[Lambda]Lik + yk*\[Lambda]Rik)))*CC[5])/
   (xs^3*xt^2) - (2*mi^2*(Qj - Qk)^2*(xs^3 + xt^3 - 2*(yk^2 - yS^2)^2 + 
     xt^2*(-1 - 2*yk^2 + 2*yS^2) + xs^2*(-1 + xt - 2*yk^2 + 2*yS^2) + 
     2*xt*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4) + 
     xs*(xt^2 - 4*xt*yk^2 + 2*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4)))*
    \[Lambda]Lik*\[Lambda]Ljk*CC[6])/(xs^3*xt^2) + 
  (2*mi^2*(Qj - Qk)*\[Lambda]Ljk*
    (Qj*(xs^6 + (-1 + xt)^2*(yk^2 - yS^2)^2 + 
       2*xs^5*(-2 + xt - yk^2 + yS^2) + xs^4*(6 + xt^2 + yk^4 - 8*yS^2 + 
         yS^4 - 2*yk^2*(-4 + yS^2) + xt*(-6 - 4*yk^2 + 2*yS^2)) - 
       2*xs^3*(xt^2*(1 + yk^2) + 2*(1 + yk^4 - 3*yS^2 + yS^4 + 
           yk^2*(3 - 2*yS^2)) - xt*(3 + yk^4 - 2*yS^2 + yS^4 - 
           2*yk^2*(-3 + yS^2))) + xs^2*(1 + 6*yk^4 - 8*yS^2 + 6*yS^4 + 
         yk^2*(8 - 12*yS^2) + xt^2*(1 + yk^4 + 4*yS^2 + yS^4 - 
           2*yk^2*(-2 + yS^2)) - 2*xt*(1 + 3*yk^4 - yS^2 + 3*yS^4 - 
           6*yk^2*(-1 + yS^2))) + 2*xs*(-((2 - 3*xt + xt^2)*yk^4) + 
         (-1 + xt)*yS^2*(-1 + xt^2 + 2*yS^2 - xt*(1 + yS^2)) + 
         (-1 + xt)*yk^2*(1 - 4*yS^2 + xt*(-1 + 2*yS^2))))*\[Lambda]Lik - 
     Qk*(-1 + xs)^2*xs*(xs^3*\[Lambda]Lik + 2*xs^2*(-1 + xt - yk^2 + yS^2)*
        \[Lambda]Lik + xs*(1 + 2*xt^2 + 4*yk^2 - 4*yS^2 + 
         xt*(-3 - 2*yk^2 + 2*yS^2))*\[Lambda]Lik + xs*xt*yk*\[Lambda]Rik + 
       (-1 + xt)*(2*yk^2*\[Lambda]Lik - 2*yS^2*\[Lambda]Lik + 
         xt*yk*\[Lambda]Rik)))*CC[7])/((-1 + xs)*xs^3*xt^2*
    (-1 + xs + xt)^2) - (2*mi^2*Qk*\[Lambda]Ljk*
    (Qk*(-1 + xs)^2*xs*(xs^3*\[Lambda]Lik + 2*xs^2*(-1 + xt + yk^2 - yS^2)*
        \[Lambda]Lik + xs*(1 - 4*yk^2 + 4*yS^2 + 2*xt*(-1 + yk^2 - yS^2))*
        \[Lambda]Lik - xs*xt*yk*\[Lambda]Rik - 
       (-1 + xt)*(2*yk^2*\[Lambda]Lik - 2*yS^2*\[Lambda]Lik + 
         xt*yk*\[Lambda]Rik)) + Qj*((-1 + xt)^2*(yk^2 - yS^2)^2*
        \[Lambda]Lik + xs^4*(xt^2*\[Lambda]Lik + (yk^2 - yS^2)^2*
          \[Lambda]Lik - 2*xt*(yS^2*\[Lambda]Lik + yk*\[Lambda]Rik)) + 
       2*xs*(-1 + xt)*(2*(yk^2 - yS^2)^2*\[Lambda]Lik + 
         xt^2*yk*(yk*\[Lambda]Lik + \[Lambda]Rik) - 
         xt*(yk^4*\[Lambda]Lik + yS^2*(1 + yS^2)*\[Lambda]Lik + 
           yk^2*(\[Lambda]Lik - 2*yS^2*\[Lambda]Lik) + yk*\[Lambda]Rik)) - 
       2*xs^3*(2*(yk^2 - yS^2)^2*\[Lambda]Lik + 
         xt^2*(\[Lambda]Lik + yS^2*\[Lambda]Lik + 2*yk*\[Lambda]Rik) - 
         xt*(yk^4*\[Lambda]Lik + yS^2*(3 + yS^2)*\[Lambda]Lik + 
           yk^2*(\[Lambda]Lik - 2*yS^2*\[Lambda]Lik) + 3*yk*\[Lambda]Rik)) + 
       xs^2*(6*(yk^2 - yS^2)^2*\[Lambda]Lik - 2*xt^3*yk*\[Lambda]Rik - 
         2*xt*(3*yk^4*\[Lambda]Lik + yk^2*(2 - 6*yS^2)*\[Lambda]Lik + 
           3*yS^2*(1 + yS^2)*\[Lambda]Lik + 3*yk*\[Lambda]Rik) + 
         xt^2*((1 + yk^4 + 4*yS^2 + yS^4 - 2*yk^2*(-2 + yS^2))*\[Lambda]Lik + 
           8*yk*\[Lambda]Rik))))*CC[8])/((-1 + xs)*xs^3*xt^2*
    (-1 + xs + xt)^2) + 
  (2*mi^2*(Qj - Qk)*(-(Qk*(-1 + xt)*xt*(xt^3 - 2*(-1 + xs)^2*(yk^2 - yS^2) + 
        2*xt^2*(-1 + xs - yk^2 + yS^2) + xt*(1 + 2*xs^2 + 4*yk^2 - 4*yS^2 + 
          xs*(-2 - 4*yk^2 + 4*yS^2)))) + 
     Qj*(xt^5 - (-1 + xs)^2*(yk^2 - yS^2)^2 + 
       xt^4*(-3 + 2*xs - 2*yk^2 + 2*yS^2) + xt^3*(3 + xs^2 + yk^4 - 6*yS^2 + 
         yS^4 - 2*yk^2*(-3 + yS^2) + xs*(-4 - 4*yk^2 + 2*yS^2)) + 
       xt*((3 - 4*xs + xs^2)*yk^4 - (-1 + xs)*yS^2*(-2 + 4*xs^2 + 3*yS^2 - 
           xs*yS^2) - 2*(-1 + xs)*yk^2*(1 - 3*yS^2 + xs*(-1 + yS^2))) - 
       xt^2*(1 + 3*yk^4 - 6*yS^2 + 3*yS^4 - 6*yk^2*(-1 + yS^2) + 
         xs^2*(1 + 2*yk^2 + 4*yS^2) - 2*xs*(yk^4 - 2*yk^2*(-2 + yS^2) + 
           (-1 + yS^2)^2))))*\[Lambda]Lik*\[Lambda]Ljk*CC[9])/
   (xs^3*xt^2*(-1 + xs + xt)^2) + 
  (2*mi^2*Qk*\[Lambda]Ljk*(-(Qk*(-1 + xt)*xt*(xs^2*(1 + 2*yk^2 - 2*yS^2) + 
        (-1 + xt)^2*(xt + 2*yk^2 - 2*yS^2) + 2*xs*(xt^2 - 2*yk^2 + 2*yS^2 + 
          xt*(-1 + 2*yk^2 - 2*yS^2)))*\[Lambda]Lik) + 
     Qj*(-((-1 + xt)^3*(yk^2 - yS^2)^2*\[Lambda]Lik) + 
       2*xs^3*xt*yk*(2*yk*\[Lambda]Lik + \[Lambda]Rik) + 
       xs*(-1 + xt)^2*(-2*(yk^2 - yS^2)^2*\[Lambda]Lik + 
         xt*yk*(2*yk*\[Lambda]Lik + \[Lambda]Rik)) - 
       xs^2*(-1 + xt)*(xt^2*\[Lambda]Lik + (yk^2 - yS^2)^2*\[Lambda]Lik - 
         xt*(\[Lambda]Lik + 6*yk^2*\[Lambda]Lik + 3*yk*\[Lambda]Rik))))*
    CC[10])/(xs^3*xt^2*(-1 + xs + xt)^2) + 
  (2*mi^2*Qk^2*(xs^2 + xt^2)*(xs + xt + 2*yk^2 - 2*yS^2)*\[Lambda]Lik*
    \[Lambda]Ljk*CC[11])/(xs^3*xt^2) - 
  (2*mi^2*(Qj - Qk)^2*(xs^2 + xt^2)*(xs + xt - 2*yk^2 + 2*yS^2)*\[Lambda]Lik*
    \[Lambda]Ljk*CC[12])/(xs^3*xt^2) + 
  (2*mi^4*(Qj - Qk)*Qk*\[Lambda]Ljk*((-1 + xt)^3*(yk^2 - yS^2)^3*
      \[Lambda]Lik - xs*(-1 + xt)^2*(yk^2 - yS^2)*
      (-3*(yk^2 - yS^2)^2*\[Lambda]Lik + xt*(3*yk^2*\[Lambda]Lik + 
         yS^2*\[Lambda]Lik + yk*\[Lambda]Rik)) + 
     xs^3*(xt^3*\[Lambda]Lik + (yk^2 - yS^2)^3*\[Lambda]Lik + 
       xt^2*((-1 + yk^2 - yS^2)*\[Lambda]Lik + yk*\[Lambda]Rik) + 
       xt*(-3*yk^4*\[Lambda]Lik + yS^2*(1 + yS^2)*\[Lambda]Lik + 
         yk^2*(\[Lambda]Lik + 2*yS^2*\[Lambda]Lik) - yk^3*\[Lambda]Rik + 
         yk*yS^2*\[Lambda]Rik)) + xs^2*(-1 + xt)*
      (3*(yk^2 - yS^2)^3*\[Lambda]Lik + xt^2*(yk^2*\[Lambda]Lik - 
         yS^2*\[Lambda]Lik + yk*\[Lambda]Rik) + 
       xt*(-6*yk^4*\[Lambda]Lik + yS^2*(1 + 2*yS^2)*\[Lambda]Lik + 
         yk^2*(\[Lambda]Lik + 4*yS^2*\[Lambda]Lik) - 2*yk^3*\[Lambda]Rik + 
         2*yk*yS^2*\[Lambda]Rik)))*DD[1])/(xs^3*xt^2*(-1 + xs + xt)^2) + 
  (2*mi^4*(Qj - Qk)*Qk*(xs*(xt - (yk - yS)^2) - (-1 + xt)*(yk - yS)^2)*
    ((-1 + xt)*(yk^2 - yS^2) + xs*(xt + yk^2 - yS^2))*
    (-((-1 + xt)*(yk + yS)^2) + xs*(xt - (yk + yS)^2))*\[Lambda]Lik*
    \[Lambda]Ljk*DD[2])/(xs^3*xt^2*(-1 + xs + xt)^2) + 
  (2*mi^4*(Qj - Qk)^2*(xs - yk^2 + yS^2)*(xs^3 + (-1 + xt)*(yk^2 - yS^2)^2 + 
     xs^2*(-1 + xt - 2*yk^2 + 2*yS^2) + xs*(yk^4 + yS^2*(-2 - 2*xt + yS^2) - 
       2*yk^2*(-1 + xt + yS^2)))*\[Lambda]Lik*\[Lambda]Ljk*DD[3])/
   (xs^3*xt^2) - (2*mi^4*Qk^2*\[Lambda]Ljk*(xs^5*\[Lambda]Lik + 
     xs^4*(-2 + 2*xt + 3*yk^2 - 3*yS^2)*\[Lambda]Lik + 
     (-1 + xt)^2*(yk^2 - yS^2)^3*\[Lambda]Lik + 
     xs^3*((1 + xt^2 + 3*yk^4 + 6*yS^2 + 3*yS^4 + 2*xt*(-1 + yk^2 - 3*yS^2) - 
         6*yk^2*(1 + yS^2))*\[Lambda]Lik - 3*xt*yk*\[Lambda]Rik) - 
     xs*(-1 + xt)*(yk^2 - yS^2)*(-2*yk^4*\[Lambda]Lik + 
       yS^2*(-3 + 3*xt - 2*yS^2)*\[Lambda]Lik + yk^2*(3 + xt + 4*yS^2)*
        \[Lambda]Lik + xt*yk*\[Lambda]Rik) + 
     xs^2*(yk^6*\[Lambda]Lik - yS^2*(3 + 3*xt^2 + 6*yS^2 + yS^4 - 
         6*xt*(1 + yS^2))*\[Lambda]Lik + yk^4*(2*xt - 3*(2 + yS^2))*
        \[Lambda]Lik + yk^2*(-5*xt^2 - 8*xt*yS^2 + 3*(1 + 4*yS^2 + yS^4))*
        \[Lambda]Lik - xt*yk^3*\[Lambda]Rik + xt*yk*(3 - 3*xt + yS^2)*
        \[Lambda]Rik))*DD[4])/(xs^3*xt^2*(-1 + xs + xt)) + 
  (2*mi^4*(Qj - Qk)^2*(xt - yk^2 + yS^2)*(xt^3 + (-1 + xs)*(yk^2 - yS^2)^2 + 
     xt^2*(-1 + xs - 2*yk^2 + 2*yS^2) + xt*(yk^4 + yS^2*(-2 - 2*xs + yS^2) - 
       2*yk^2*(-1 + xs + yS^2)))*\[Lambda]Lik*\[Lambda]Ljk*DD[5])/
   (xs^3*xt^2) - (2*mi^4*Qk^2*\[Lambda]Ljk*(xt^5*\[Lambda]Lik + 
     xt^4*(-2 + 2*xs + 3*yk^2 - 3*yS^2)*\[Lambda]Lik + 
     (-1 + xs)^2*(yk^2 - yS^2)^3*\[Lambda]Lik + 
     xt^3*((1 + xs^2 + 3*yk^4 + 6*yS^2 + 3*yS^4 + 2*xs*(-1 + yk^2 - 3*yS^2) - 
         6*yk^2*(1 + yS^2))*\[Lambda]Lik - xs*yk*\[Lambda]Rik) + 
     xt^2*(yk^6*\[Lambda]Lik - yS^2*(3 + 3*xs^2 + 6*yS^2 + yS^4 - 
         6*xs*(1 + yS^2))*\[Lambda]Lik + yk^4*(2*xs - 3*(2 + yS^2))*
        \[Lambda]Lik + yk^2*(-5*xs^2 - 2*xs*(1 + 4*yS^2) + 
         3*(1 + 4*yS^2 + yS^4))*\[Lambda]Lik - xs*yk^3*\[Lambda]Rik + 
       xs*yk*(1 - xs + yS^2)*\[Lambda]Rik) + 
     xt*(-((-3 + 2*yk^2 - 2*yS^2)*(yk^2 - yS^2)^2*\[Lambda]Lik) + 
       xs*(yk^2 - yS^2)*(2*yk^4*\[Lambda]Lik + 2*yS^2*(3 + yS^2)*
          \[Lambda]Lik - 2*yk^2*(\[Lambda]Lik + 2*yS^2*\[Lambda]Lik) + 
         yk*\[Lambda]Rik) + xs^2*(-(yk^4*\[Lambda]Lik) + 
         3*yS^4*\[Lambda]Lik - 2*yk^2*(-1 + yS^2)*\[Lambda]Lik - 
         yk^3*\[Lambda]Rik + yk*yS^2*\[Lambda]Rik)))*DD[6])/
   (xs^3*xt^2*(-1 + xs + xt)), 
 (-2*Qj*(-2*Qk*xs*xt*(-2 + xs + xt) + Qj*(-((-1 + xt)*xt*(yk^2 - yS^2)) + 
       xs^2*(xt - yk^2 + 2*xt*yk^2 + yS^2 - 2*xt*yS^2) + 
       xs*(yk^2 - yS^2 + xt^2*(1 + 2*yk^2 - 2*yS^2) + 
         xt*(-2 - 4*yk^2 + 4*yS^2))))*\[Lambda]Rik*\[Lambda]Rjk)/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)) + 
  (2*Qj^2*(xs + xt - 2*xs*xt)*yk^2*\[Lambda]Rik*\[Lambda]Rjk*BB[1])/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2) + 
  (2*Qj^2*(-xt + xs*(-1 + 2*xt))*yS^2*\[Lambda]Rik*\[Lambda]Rjk*BB[2])/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2) + 
  (2*Qj*(-2*Qk*(2 - 2*xs + xs^2 - 2*xt + xt^2)*\[Lambda]Rik + 
     Qj*(4*(-1 + xs)*(-1 + xt)*(-2 + xs + xt)*yk*\[Lambda]Lik + 
       (-6 + 2*xs*(-2 + xt)^2 + 8*xt - 3*xt^2 + xs^2*(-3 + 2*xt))*yk^2*
        \[Lambda]Rik - (2 - 4*xt + xt^2 - 6*yS^2 + 8*xt*yS^2 - 3*xt^2*yS^2 + 
         xs^2*(1 - 3*yS^2 + 2*xt*(-1 + yS^2)) + 
         2*xs*(-2 + 4*yS^2 - 4*xt*(-1 + yS^2) + xt^2*(-1 + yS^2)))*
        \[Lambda]Rik))*\[Lambda]Rjk*BB[3])/((-1 + xs)^2*xs*(-1 + xt)^2*xt*
    (-1 + xs + xt)) + 
  (2*Qj*(2*Qk*xs*\[Lambda]Rik + Qj*((yk^2 - yS^2)*\[Lambda]Rik - 
       2*xs^2*(2*yk*\[Lambda]Lik + \[Lambda]Rik) + 
       xs*(4*yk*\[Lambda]Lik + \[Lambda]Rik)))*\[Lambda]Rjk*BB[4])/
   ((-1 + xs)^2*xs^2*xt*(-1 + xs + xt)) + 
  (2*Qj*(2*Qk*xt*\[Lambda]Rik + Qj*((yk^2 - yS^2)*\[Lambda]Rik - 
       2*xt^2*(2*yk*\[Lambda]Lik + \[Lambda]Rik) + 
       xt*(4*yk*\[Lambda]Lik + \[Lambda]Rik)))*\[Lambda]Rjk*BB[5])/
   (xs*(-1 + xt)^2*xt^2*(-1 + xs + xt)) + 
  (2*mi^2*Qk*(Qk*xs*(2*xt^2*yk*\[Lambda]Lik + (-1 + xs)*(yk^2 - yS^2)*
        \[Lambda]Rik + xt*(2*(-1 + xs)*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik + 
         (-1 + 2*xs - yS^2)*\[Lambda]Rik)) + 
     Qj*(-((-1 + xt)*xt*yk^2*\[Lambda]Rik) + 
       xs^2*(yS^2*\[Lambda]Rik - 2*xt*(yk*\[Lambda]Lik + \[Lambda]Rik)) + 
       xs*(-2*xt^2*yk*\[Lambda]Lik - yS^2*\[Lambda]Rik + 
         xt*(2*yk*\[Lambda]Lik + \[Lambda]Rik - yk^2*\[Lambda]Rik + 
           yS^2*\[Lambda]Rik))))*\[Lambda]Rjk*CC[1])/
   (xs*xt^2*(-1 + xs + xt)^3) - 
  (2*mi^2*(Qj - Qk)*(Qj*(-xs + xs^2 + xt - xt^2)*yS^2*\[Lambda]Rik + 
     Qk*xs*((-1 + xs)*(yk^2 - yS^2)*\[Lambda]Rik + 
       2*xt^2*(yk*\[Lambda]Lik + \[Lambda]Rik) + 
       xt*(2*(-1 + xs)*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - 
         (1 + yS^2)*\[Lambda]Rik)))*\[Lambda]Rjk*CC[2])/
   (xs*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^2*Qk*(Qk*xt*(2*xs^2*yk*\[Lambda]Lik + (-1 + xt)*(yk^2 - yS^2)*
        \[Lambda]Rik + xs*(2*(-1 + xt)*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik + 
         (-1 + 2*xt - yS^2)*\[Lambda]Rik)) + 
     Qj*((-1 + xt)*xt*yS^2*\[Lambda]Rik - xs^2*yk*(2*xt*\[Lambda]Lik + 
         yk*\[Lambda]Rik) + xs*(yk^2*\[Lambda]Rik - 
         2*xt^2*(yk*\[Lambda]Lik + \[Lambda]Rik) + 
         xt*(2*yk*\[Lambda]Lik + \[Lambda]Rik - yk^2*\[Lambda]Rik + 
           yS^2*\[Lambda]Rik))))*\[Lambda]Rjk*CC[3])/
   (xs^2*xt*(-1 + xs + xt)^3) + 
  (2*mi^2*(Qj - Qk)*(Qj*(-xs + xs^2 + xt - xt^2)*yS^2*\[Lambda]Rik + 
     Qk*xt*(-((-1 + xt)*(yk^2 - yS^2)*\[Lambda]Rik) - 
       2*xs^2*(yk*\[Lambda]Lik + \[Lambda]Rik) + 
       xs*(2*yk*\[Lambda]Lik - 2*xt*yk*\[Lambda]Lik + \[Lambda]Rik - 
         yk^2*\[Lambda]Rik + yS^2*\[Lambda]Rik)))*\[Lambda]Rjk*CC[4])/
   (xs^2*xt*(-1 + xs + xt)^3) - (4*mi^2*Qk^2*(xs + xt)*yk^2*\[Lambda]Rik*
    \[Lambda]Rjk*CC[5])/(xs^2*xt^2*(-1 + xs + xt)) + 
  (4*mi^2*(Qj - Qk)^2*(xs + xt)*yS^2*\[Lambda]Rik*\[Lambda]Rjk*CC[6])/
   (xs^2*xt^2*(-1 + xs + xt)) - 
  (2*mi^2*(Qj - Qk)*(Qj*(xs^4 + (-1 + xt)*xt + xs^3*(-3 + 6*xt) + 
       xs^2*(3 - 13*xt + 7*xt^2) + xs*(-1 + 8*xt - 8*xt^2 + 2*xt^3))*yS^2*
      \[Lambda]Rik + Qk*(-1 + xs)^2*xt*((-1 + xt)*(yk^2 - yS^2)*
        \[Lambda]Rik + 2*xs^2*(yk*\[Lambda]Lik + \[Lambda]Rik) + 
       xs*(2*(-1 + xt)*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - 
         (1 + yS^2)*\[Lambda]Rik)))*\[Lambda]Rjk*CC[7])/
   ((-1 + xs)*xs^2*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^2*Qk*(Qk*(-1 + xs)^2*xt*(2*xs^2*yk*\[Lambda]Lik + 
       (-1 + xt)*(yk^2 - yS^2)*\[Lambda]Rik + 
       xs*(2*(-1 + xt)*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik + 
         (-1 + 2*xt - yS^2)*\[Lambda]Rik)) + 
     Qj*((-1 + xt)*xt*yS^2*\[Lambda]Rik + xs^4*yk*(-2*xt*\[Lambda]Lik + 
         yk*\[Lambda]Rik) + xs^3*(-3*yk^2*\[Lambda]Rik - 
         2*xt^2*(yk*\[Lambda]Lik + \[Lambda]Rik) + 
         xt*(6*yk*\[Lambda]Lik + \[Lambda]Rik + 5*yk^2*\[Lambda]Rik + 
           yS^2*\[Lambda]Rik)) + xs*(-(yk^2*\[Lambda]Rik) + 
         2*xt^3*yk^2*\[Lambda]Rik - 2*xt^2*(yk*\[Lambda]Lik + \[Lambda]Rik + 
           3*yk^2*\[Lambda]Rik + yS^2*\[Lambda]Rik) + 
         xt*(2*yk*\[Lambda]Lik + \[Lambda]Rik + 5*yk^2*\[Lambda]Rik + 
           3*yS^2*\[Lambda]Rik)) + xs^2*(3*yk^2*\[Lambda]Rik + 
         xt^2*(4*yk*\[Lambda]Lik + 6*yk^2*\[Lambda]Rik + 
           (4 + yS^2)*\[Lambda]Rik) - xt*(6*yk*\[Lambda]Lik + 
           10*yk^2*\[Lambda]Rik + (2 + 3*yS^2)*\[Lambda]Rik))))*\[Lambda]Rjk*
    CC[8])/((-1 + xs)*xs^2*xt^2*(-1 + xs + xt)^3) - 
  (2*mi^2*(Qj - Qk)*(Qj*(2*xs^3*xt + (-1 + xt)^3*xt + 
       xs*(-1 + xt)^2*(-1 + 6*xt) + xs^2*(1 - 8*xt + 7*xt^2))*yS^2*
      \[Lambda]Rik + Qk*xs*(-1 + xt)^2*((-1 + xs)*(yk^2 - yS^2)*
        \[Lambda]Rik + 2*xt^2*(yk*\[Lambda]Lik + \[Lambda]Rik) + 
       xt*(2*(-1 + xs)*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik - 
         (1 + yS^2)*\[Lambda]Rik)))*\[Lambda]Rjk*CC[9])/
   (xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^2*Qk*(Qk*xs*(-1 + xt)^2*(2*xt^2*yk*\[Lambda]Lik + 
       (-1 + xs)*(yk^2 - yS^2)*\[Lambda]Rik + 
       xt*(2*(-1 + xs)*yk*\[Lambda]Lik + yk^2*\[Lambda]Rik + 
         (-1 + 2*xs - yS^2)*\[Lambda]Rik)) + 
     Qj*(2*xs^3*xt*yk^2*\[Lambda]Rik + (-1 + xt)^3*xt*yk^2*\[Lambda]Rik - 
       xs*(-1 + xt)^2*(2*xt^2*yk*\[Lambda]Lik + yS^2*\[Lambda]Rik - 
         xt*(2*yk*\[Lambda]Lik + \[Lambda]Rik + 5*yk^2*\[Lambda]Rik + 
           yS^2*\[Lambda]Rik)) - xs^2*(-1 + xt)*(yS^2*\[Lambda]Rik + 
         2*xt^2*(yk*\[Lambda]Lik + \[Lambda]Rik) - 
         xt*(2*yk*\[Lambda]Lik + 6*yk^2*\[Lambda]Rik + 
           (2 + yS^2)*\[Lambda]Rik))))*\[Lambda]Rjk*CC[10])/
   (xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^4*(Qj - Qk)*Qk*(-((-1 + xt)^2*xt*yk^2*(yk^2 - yS^2)*\[Lambda]Rik) + 
     xs^3*(yS^2*(-yk^2 + yS^2)*\[Lambda]Rik + 
       2*xt^2*(yk*\[Lambda]Lik + \[Lambda]Rik) - 
       xt*(2*yk^3*\[Lambda]Lik + 2*yk*yS^2*\[Lambda]Lik + 
         2*yk^2*\[Lambda]Rik + 3*yS^2*\[Lambda]Rik)) - 
     xs*(-1 + xt)*(yS^2*(-yk^2 + yS^2)*\[Lambda]Rik + 
       xt^2*yk*(2*yk^2*\[Lambda]Lik + 2*yS^2*\[Lambda]Lik - 
         yk*\[Lambda]Rik) + xt*(-2*yk^3*\[Lambda]Lik - 
         2*yk*yS^2*\[Lambda]Lik + 2*yk^4*\[Lambda]Rik - 
         yk^2*(1 + yS^2)*\[Lambda]Rik - yS^2*(1 + yS^2)*\[Lambda]Rik)) + 
     xs^2*(2*xt^3*yk*\[Lambda]Lik + 2*yS^2*(yk^2 - yS^2)*\[Lambda]Rik - 
       xt^2*(4*yk^3*\[Lambda]Lik + 2*yk*(\[Lambda]Lik + 
           2*yS^2*\[Lambda]Lik) + \[Lambda]Rik + yk^2*\[Lambda]Rik + 
         3*yS^2*\[Lambda]Rik) + xt*(4*yk^3*\[Lambda]Lik + 
         4*yk*yS^2*\[Lambda]Lik - yk^4*\[Lambda]Rik - yk^2*(-3 + yS^2)*
          \[Lambda]Rik + 2*yS^2*(2 + yS^2)*\[Lambda]Rik)))*\[Lambda]Rjk*
    DD[1])/(xs^2*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^4*(Qj - Qk)*Qk*((-1 + xt)^2*xt*yS^2*(-yk^2 + yS^2)*\[Lambda]Rik + 
     xs^3*yk*(2*xt^2*\[Lambda]Lik + yk*(-yk^2 + yS^2)*\[Lambda]Rik + 
       xt*(-2*yk^2*\[Lambda]Lik - 2*yS^2*\[Lambda]Lik + yk*\[Lambda]Rik)) + 
     xs^2*(2*yk^2*(yk^2 - yS^2)*\[Lambda]Rik + 
       2*xt^3*(yk*\[Lambda]Lik + \[Lambda]Rik) - 
       xt^2*(2*yk*\[Lambda]Lik + 4*yk^3*\[Lambda]Lik + 
         4*yk*yS^2*\[Lambda]Lik + \[Lambda]Rik + yk^2*\[Lambda]Rik + 
         3*yS^2*\[Lambda]Rik) + xt*(4*yk^3*\[Lambda]Lik + 
         4*yk*yS^2*\[Lambda]Lik - 2*yk^4*\[Lambda]Rik + yS^2*\[Lambda]Rik + 
         yk^2*yS^2*\[Lambda]Rik + yS^4*\[Lambda]Rik)) - 
     xs*(-1 + xt)*(yk^2*(-yk^2 + yS^2)*\[Lambda]Rik + 
       xt^2*(2*yk^3*\[Lambda]Lik + 2*yk*yS^2*\[Lambda]Lik + 
         2*yk^2*\[Lambda]Rik + 3*yS^2*\[Lambda]Rik) + 
       xt*(-2*yk^3*\[Lambda]Lik - 2*yk*yS^2*\[Lambda]Lik + 
         yk^4*\[Lambda]Rik + yk^2*(-1 + yS^2)*\[Lambda]Rik - 
         yS^2*(1 + 2*yS^2)*\[Lambda]Rik)))*\[Lambda]Rjk*DD[2])/
   (xs^2*xt^2*(-1 + xs + xt)^3) - 
  (2*mi^4*(Qj - Qk)^2*yS^2*(xs^2*\[Lambda]Rik + xt*(-yk^2 + yS^2)*
      \[Lambda]Rik - xs*((yk^2 - yS^2)*\[Lambda]Rik + 
       xt*(4*yk*\[Lambda]Lik + \[Lambda]Rik)))*\[Lambda]Rjk*DD[3])/
   (xs^2*xt^2*(-1 + xs + xt)) + 
  (2*mi^4*Qk^2*yk^2*(xs^2*\[Lambda]Rik + xt*(yk^2 - yS^2)*\[Lambda]Rik + 
     xs*(4*xt*yk*\[Lambda]Lik + 3*xt*\[Lambda]Rik + yk^2*\[Lambda]Rik - 
       yS^2*\[Lambda]Rik))*\[Lambda]Rjk*DD[4])/(xs^2*xt^2*(-1 + xs + xt)) + 
  (2*mi^4*(Qj - Qk)^2*yS^2*(xs*(yk^2 - yS^2)*\[Lambda]Rik - 
     xt*(xt - yk^2 + yS^2)*\[Lambda]Rik + 
     xs*xt*(4*yk*\[Lambda]Lik + \[Lambda]Rik))*\[Lambda]Rjk*DD[5])/
   (xs^2*xt^2*(-1 + xs + xt)) + 
  (2*mi^4*Qk^2*yk^2*(xt*(xt + yk^2 - yS^2)*\[Lambda]Rik + 
     xs*(4*xt*yk*\[Lambda]Lik + 3*xt*\[Lambda]Rik + yk^2*\[Lambda]Rik - 
       yS^2*\[Lambda]Rik))*\[Lambda]Rjk*DD[6])/(xs^2*xt^2*(-1 + xs + xt)), 
 (Qj*\[Lambda]Ljk*(2*Qk*(1 - 2*xs)*\[Lambda]Lik + 
     Qj*(-1 + yk^2 - yS^2 + xs*(2 + yk^2 - yS^2))*\[Lambda]Lik + 
     2*Qj*xs*yk*\[Lambda]Rik))/((-1 + xs)*xs^2*xt) + 
  (Qj^2*yk^2*\[Lambda]Ljk*((1 + xs)*yk^2*\[Lambda]Lik - 
     (1 + xs)*yS^2*\[Lambda]Lik + 2*xs*yk*\[Lambda]Rik)*BB[1])/
   ((-1 + xs)*xs^2*xt*(yk^2 - yS^2)) + 
  (Qj^2*yS^2*\[Lambda]Ljk*(-((1 + xs)*yk^2*\[Lambda]Lik) + 
     (1 + xs)*yS^2*\[Lambda]Lik - 2*xs*yk*\[Lambda]Rik)*BB[2])/
   ((-1 + xs)*xs^2*xt*(yk^2 - yS^2)) + 
  (Qj*\[Lambda]Ljk*(2*Qk*(-xs + xt)*\[Lambda]Lik + 
     Qj*(xt*(-1 + (-1 + 2*xt)*yk^2 + yS^2 - 2*xt*yS^2)*\[Lambda]Lik + 
       xs*(1 - yk^2 + yS^2 + xt*(1 + yk^2 - yS^2) + xt^2*(-1 - yk^2 + yS^2))*
        \[Lambda]Lik + 2*(-1 + xt)*xt*yk*\[Lambda]Rik - 
       2*xs*(-1 + xt)*xt*yk*\[Lambda]Rik))*BB[3])/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2) + 
  (Qj*\[Lambda]Ljk*(-2*Qk*xs^2*\[Lambda]Lik + 
     Qj*(xs^2*\[Lambda]Lik + xs*(xt - yk^2 + yS^2)*\[Lambda]Lik + 
       xt*(-3*yk^2*\[Lambda]Lik + 3*yS^2*\[Lambda]Lik - 2*yk*\[Lambda]Rik)))*
    BB[4])/((-1 + xs)*xs^2*xt^2) + (Qj*(2*Qk*xt - Qj*(xt - yk^2 + yS^2))*
    \[Lambda]Lik*\[Lambda]Ljk*BB[5])/(xs^2*(-1 + xt)*xt) + 
  (Qk^2*(xs - xt)*\[Lambda]Lik*\[Lambda]Ljk*BB[6])/(xs^2*xt^2) - 
  ((Qj - Qk)^2*(xs - xt)*\[Lambda]Lik*\[Lambda]Ljk*BB[7])/(xs^2*xt^2) + 
  (mi^2*Qk*\[Lambda]Ljk*(Qk*(-1 + xs)*xs*
      (xs^4 + xs^3*(-2 + 3*xt + 2*yk^2 - 2*yS^2) + 2*(-1 + xt)^2*xt*
        (yk^2 - yS^2) + xs^2*(1 + 2*xt^2 - 4*yk^2 + 4*yS^2 + 
         4*xt*(-1 + yk^2 - yS^2)) + xs*(xt^2*(-2 + 4*yk^2 - 4*yS^2) + 
         2*(yk^2 - yS^2) + xt*(1 - 6*yk^2 + 6*yS^2)))*\[Lambda]Lik + 
     Qj*((-1 + xt)^2*xt*(yk^2 - yS^2)^2*\[Lambda]Lik + 
       xs^4*(xt^2 - 2*xt*yS^2 + (yk^2 - yS^2)^2)*\[Lambda]Lik - 
       xs*(-1 + xt)*(-((yk^2 - yS^2)^2*\[Lambda]Lik) - 2*xt*(yk^2 - yS^2)^2*
          \[Lambda]Lik + 4*xt^3*yk*(yk*\[Lambda]Lik + \[Lambda]Rik) + 
         xt^2*(yk^4*\[Lambda]Lik + yS^4*\[Lambda]Lik - 2*yk^2*(3 + yS^2)*
            \[Lambda]Lik - 4*yk*\[Lambda]Rik)) + 
       xs^3*(xt^3*\[Lambda]Lik - 3*(yk^2 - yS^2)^2*\[Lambda]Lik + 
         xt*(yk^4 + 4*yS^2 - 2*yk^2*yS^2 + yS^4)*\[Lambda]Lik - 
         xt^2*(\[Lambda]Lik + 6*yk^2*\[Lambda]Lik + 2*yS^2*\[Lambda]Lik + 
           4*yk*\[Lambda]Rik)) - xs^2*(-3*(yk^2 - yS^2)^2*\[Lambda]Lik + 
         xt*(yk^4 + 2*yS^2 - 2*yk^2*yS^2 + yS^4)*\[Lambda]Lik + 
         xt^2*(yk^4*\[Lambda]Lik + yS^2*(-2 + yS^2)*\[Lambda]Lik - 
           2*yk^2*(6 + yS^2)*\[Lambda]Lik - 8*yk*\[Lambda]Rik) + 
         xt^3*(\[Lambda]Lik + 10*yk^2*\[Lambda]Lik + 8*yk*\[Lambda]Rik))))*
    CC[1])/(2*(-1 + xs)*xs^2*xt^3*(-1 + xs + xt)^2) - 
  (mi^2*(Qk^2*(-1 + xs)*xs*(xs^4 - 2*(-1 + xt)^2*xt*(yk^2 - yS^2) + 
       xs^3*(-2 + 3*xt - 2*yk^2 + 2*yS^2) + xs^2*(1 + 4*xt^2 + 4*yk^2 - 
         4*yS^2 - 4*xt*(1 + yk^2 - yS^2)) + xs*(2*xt^3 - 2*yk^2 + 2*yS^2 + 
         xt*(1 + 6*yk^2 - 6*yS^2) - 4*xt^2*(1 + yk^2 - yS^2))) + 
     Qj^2*(xs^6 + (-1 + xt)^2*xt*(yk^2 - yS^2)^2 + 
       xs^5*(-3 + 3*xt - 2*yk^2 + 2*yS^2) + xs^4*(3 + 3*xt^2 + yk^4 - 
         6*yS^2 + yS^4 - 2*yk^2*(-3 + yS^2) + xt*(-7 - 4*yk^2 + 2*yS^2)) - 
       xs*(-1 + xt)*((-1 - 2*xt + xt^2)*yk^4 - 2*(-1 - 2*xt + xt^2)*yk^2*
          yS^2 + yS^2*(4*xt^3 - yS^2 - 2*xt*yS^2 + xt^2*(-6 + yS^2))) + 
       xs^3*(-1 + xt^3 - 3*yk^4 + 6*yS^2 - 3*yS^4 + 6*yk^2*(-1 + yS^2) - 
         xt^2*(5 + 2*yk^2 + 6*yS^2) + xt*(5 + yk^4 - 4*yS^2 + yS^4 - 
           2*yk^2*(-4 + yS^2))) - xs^2*(-3*yk^4 + 2*yS^2 - 3*yS^4 + 
         yk^2*(-2 + 6*yS^2) + xt^3*(1 + 10*yS^2) + 
         xt*(yk^4 - 2*yk^2*(-2 + yS^2) + (-1 + yS^2)^2) + 
         xt^2*(-2 + yk^4 - 12*yS^2 + yS^4 - 2*yk^2*(1 + yS^2)))) + 
     Qj*Qk*(-2*xs^6 + xs^5*(6 - 6*xt + 4*yk^2 - 4*yS^2) - 
       (-1 + xt)^2*xt*(yk^2 - yS^2)^2 - xs^4*(6 + 7*xt^2 + yk^4 - 12*yS^2 + 
         yS^4 - 2*xt*(7 + 4*yk^2 - 3*yS^2) - 2*yk^2*(-6 + yS^2)) + 
       xs^3*(2 - 3*xt^3 + 3*yk^4 - 12*yS^2 + 3*yS^4 - 6*yk^2*(-2 + yS^2) + 
         xt^2*(13 + 6*yk^2 + 2*yS^2) - xt*(10 + yk^4 - 14*yS^2 + yS^4 - 
           2*yk^2*(-9 + yS^2))) + xs^2*(-3*yk^4 + 4*yS^2 - 3*yS^4 + 
         yk^2*(-4 + 6*yS^2) + xt^3*(3 + 2*yk^2 + 8*yS^2) + 
         xt*(2 + yk^4 - 10*yS^2 + yS^4 - 2*yk^2*(-6 + yS^2)) + 
         xt^2*(-6 + yk^4 - 4*yS^2 + yS^4 - 2*yk^2*(5 + yS^2))) + 
       xs*(4*xt^4*yS^2 + (yk^2 - yS^2)^2 + xt^3*(yk^4 + yS^2*(-8 + yS^2) - 
           2*yk^2*(1 + yS^2)) + xt*(yk^4 - 2*yk^2*(1 + yS^2) + 
           yS^2*(2 + yS^2)) + xt^2*(-3*yk^4 + 2*yS^2 - 3*yS^4 + 
           yk^2*(4 + 6*yS^2)))))*\[Lambda]Lik*\[Lambda]Ljk*CC[2])/
   (2*(-1 + xs)*xs^2*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*Qk*(-(Qk*(-1 + xt)*xt^2*(xt^2 - 2*yk^2 + 2*yS^2 + 
        xt*(-1 + 2*yk^2 - 2*yS^2) + xs*(1 + xt + 2*yk^2 - 2*yS^2))) + 
     Qj*(-((-1 + xt)^2*xt*(yk^2 - yS^2)^2) + 
       xs^3*(xt^2 + xt*(4*yk^2 - 2*yS^2) + (yk^2 - yS^2)^2) + 
       xs*(2*xt^3*yk^2 + 2*xt*(yk^2 - yS^2) + (yk^2 - yS^2)^2 - 
         xt^2*(yk^4 - 2*yk^2*(-2 + yS^2) + yS^2*(-2 + yS^2))) + 
       xs^2*(xt^3 - 2*(yk^2 - yS^2)^2 + xt^2*(6*yk^2 - 2*(1 + yS^2)) + 
         xt*(yk^4 - 2*yk^2*(3 + yS^2) + yS^2*(4 + yS^2)))))*\[Lambda]Lik*
    \[Lambda]Ljk*CC[3])/(2*xs^3*xt^2*(-1 + xs + xt)^2) + 
  (mi^2*(Qk^2*(-1 + xt)*xt^2*(2*xs^2 + xt^2 + 2*(yk^2 - yS^2) + 
       xt*(-1 - 2*yk^2 + 2*yS^2) + xs*(-1 + 3*xt - 2*yk^2 + 2*yS^2)) - 
     Qj*Qk*(xs^3*(xt^2 - 2*xt*yS^2 - (yk^2 - yS^2)^2) + 
       (-1 + xt)^2*xt*(2*xt^2 - 4*xt*(yk^2 - yS^2) + (yk^2 - yS^2)^2) + 
       xs^2*(5*xt^3 + 2*(yk^2 - yS^2)^2 - 2*xt^2*(2 + yk^2 + yS^2) - 
         xt*(yk^4 - 2*yS^2 - 2*yk^2*yS^2 + yS^4)) + 
       xs*(6*xt^4 - (yk^2 - yS^2)^2 + xt^3*(-8 - 6*yk^2 + 4*yS^2) + 
         xt^2*(2 + yk^4 - 4*yS^2 + yS^4 - 2*yk^2*(-3 + yS^2)))) + 
     Qj^2*(xs^3*(xt^2 - 2*xt*yS^2 - (yk^2 - yS^2)^2) + 
       xt*(xt^2 + yk^2 - yS^2 + xt*(-1 - yk^2 + yS^2))^2 + 
       xs^2*(3*xt^3 + 2*(yk^2 - yS^2)^2 - 2*xt^2*(1 + yk^2 + yS^2) - 
         xt*(yk^4 - 2*yS^2 - 2*yk^2*yS^2 + yS^4)) + 
       xs*(3*xt^4 - (yk^2 - yS^2)^2 + xt^3*(-4 - 4*yk^2 + 2*yS^2) + 
         xt^2*(yk^4 - 2*yk^2*(-2 + yS^2) + (-1 + yS^2)^2))))*\[Lambda]Lik*
    \[Lambda]Ljk*CC[4])/(2*xs^3*xt^2*(-1 + xs + xt)^2) - 
  (mi^2*Qk^2*(xs^4 + xs^3*(-1 + 2*xt + 2*yk^2 - 2*yS^2) + 
     xs*(xt^2*(-1 + 6*yk^2 - 2*yS^2) - 4*xt*(yk^2 - yS^2) - 
       2*(yk^2 - yS^2)^2) + xs^2*(2*xt^2 + xt*(-1 + 2*yk^2 - 6*yS^2) + 
       2*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2))) + 
     xt*(-xt^3 + 2*(yk^2 - yS^2)^2 + xt^2*(1 - 2*yk^2 + 2*yS^2) - 
       2*xt*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2))))*\[Lambda]Lik*
    \[Lambda]Ljk*CC[5])/(2*xs^3*xt^3) + 
  (mi^2*(Qj - Qk)^2*(xs^4 + xs^3*(-1 + 2*xt - 2*yk^2 + 2*yS^2) + 
     xs*(-2*xt^3 - 2*(yk^2 - yS^2)^2 + xt^2*(1 + 2*yk^2 + 2*yS^2)) - 
     xs^2*(xt*(1 + 2*yk^2 + 2*yS^2) - 2*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + 
         yS^4)) + xt*(-xt^3 + xt^2*(1 + 2*yk^2 - 2*yS^2) + 
       2*(yk^2 - yS^2)^2 - 2*xt*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4)))*
    \[Lambda]Lik*\[Lambda]Ljk*CC[6])/(2*xs^3*xt^3) - 
  (mi^2*(Qk^2*(-1 + xs)*xs*(xs^4 - 2*(-1 + xt)*xt^2*(yk^2 - yS^2) + 
       xs*(1 - 2*xt + 2*xt^2)*(xt - 2*yk^2 + 2*yS^2) + 
       xs^3*(-2 + 3*xt - 2*yk^2 + 2*yS^2) + xs^2*(1 + 4*xt^2 + 4*yk^2 - 
         4*yS^2 - 4*xt*(1 + yk^2 - yS^2))) + 
     Qj^2*(xs^6 + (-1 + xt)^2*xt*(yk^2 - yS^2)^2 + 
       xs^5*(-3 + 3*xt - 2*yk^2 + 2*yS^2) + xs^4*(3 + 3*xt^2 + yk^4 - 
         6*yS^2 + yS^4 - 2*yk^2*(-3 + yS^2) + xt*(-7 - 4*yk^2 + 2*yS^2)) + 
       xs^3*(-1 + xt^3 - 3*yk^4 + 6*yS^2 - 3*yS^4 + 6*yk^2*(-1 + yS^2) - 
         xt^2*(5 + 2*yk^2 + 6*yS^2) + xt*(5 + yk^4 - 4*yS^2 + yS^4 - 
           2*yk^2*(-4 + yS^2))) - xs*(-1 + xt)*((-1 - 2*xt + xt^2)*yk^4 - 
         2*(-1 - 2*xt + xt^2)*yk^2*yS^2 + yS^2*(4*xt^3 - yS^2 - 2*xt*yS^2 + 
           xt^2*(-2 + yS^2))) - xs^2*(-3*yk^4 + 2*yS^2 - 3*yS^4 + 
         yk^2*(-2 + 6*yS^2) + xt^3*(1 + 10*yS^2) + 
         xt*(yk^4 - 2*yk^2*(-2 + yS^2) + (-1 + yS^2)^2) + 
         xt^2*(-2 + yk^4 - 8*yS^2 + yS^4 - 2*yk^2*(1 + yS^2)))) + 
     Qj*Qk*(-2*xs^6 + xs^5*(6 - 6*xt + 4*yk^2 - 4*yS^2) - 
       (-1 + xt)^2*xt*(yk^2 - yS^2)^2 - xs^4*(6 + 7*xt^2 + yk^4 - 12*yS^2 + 
         yS^4 - 2*xt*(7 + 4*yk^2 - 3*yS^2) - 2*yk^2*(-6 + yS^2)) + 
       xs^3*(2 - 3*xt^3 + 3*yk^4 - 12*yS^2 + 3*yS^4 - 6*yk^2*(-2 + yS^2) + 
         xt^2*(11 + 6*yk^2 + 2*yS^2) - xt*(10 + yk^4 - 12*yS^2 + yS^4 - 
           2*yk^2*(-8 + yS^2))) + xs^2*(-3*yk^4 + 4*yS^2 - 3*yS^4 + 
         yk^2*(-4 + 6*yS^2) + xt^3*(3 + 2*yk^2 + 8*yS^2) + 
         xt*(2 + yk^4 - 6*yS^2 + yS^4 - 2*yk^2*(-4 + yS^2)) + 
         xt^2*(-4 + yk^4 - 2*yS^2 + yS^4 - 2*yk^2*(4 + yS^2))) + 
       xs*(4*xt^4*yS^2 + (yk^2 - yS^2)^2 + xt*(yk^2 - yS^2)^2 + 
         xt^3*(yk^4 + yS^2*(-4 + yS^2) - 2*yk^2*(1 + yS^2)) + 
         xt^2*(-3*yk^4 - 3*yS^4 + yk^2*(2 + 6*yS^2)))))*\[Lambda]Lik*
    \[Lambda]Ljk*CC[7])/(2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*Qk*\[Lambda]Ljk*(Qk*(-1 + xs)*xs*
      (xs^4 + xs^3*(-2 + 3*xt + 2*yk^2 - 2*yS^2) + 2*(-1 + xt)*xt^2*
        (yk^2 - yS^2) + xs^2*(1 + 2*xt^2 - 4*yk^2 + 4*yS^2 + 
         4*xt*(-1 + yk^2 - yS^2)) + xs*(2*(yk^2 - yS^2) + 
         4*xt^2*(yk^2 - yS^2) + xt*(1 - 4*yk^2 + 4*yS^2)))*\[Lambda]Lik + 
     Qj*((-1 + xt)^2*xt*(yk^2 - yS^2)^2*\[Lambda]Lik + 
       xs^4*(xt^2 - 2*xt*yS^2 + (yk^2 - yS^2)^2)*\[Lambda]Lik - 
       xs*(-1 + xt)*(-((yk^2 - yS^2)^2*\[Lambda]Lik) - 
         2*xt*(yk^2 + yk^4 - yS^2 - 2*yk^2*yS^2 + yS^4)*\[Lambda]Lik + 
         4*xt^3*yk*(yk*\[Lambda]Lik + \[Lambda]Rik) + 
         xt^2*(yk^4*\[Lambda]Lik + yS^4*\[Lambda]Lik - 2*yk^2*(1 + yS^2)*
            \[Lambda]Lik - 4*yk*\[Lambda]Rik)) + 
       xs^3*(xt^3*\[Lambda]Lik - 3*(yk^2 - yS^2)^2*\[Lambda]Lik + 
         xt*(yk^4 - 2*yk^2*(1 + yS^2) + yS^2*(6 + yS^2))*\[Lambda]Lik - 
         xt^2*((3 + 6*yk^2 + 2*yS^2)*\[Lambda]Lik + 4*yk*\[Lambda]Rik)) - 
       xs^2*(-3*(yk^2 - yS^2)^2*\[Lambda]Lik + xt*(yk^4 - 2*yk^2*(2 + yS^2) + 
           yS^2*(6 + yS^2))*\[Lambda]Lik + 
         xt^2*((-2 + yk^4 - 4*yS^2 + yS^4 - 2*yk^2*(3 + yS^2))*\[Lambda]Lik - 
           8*yk*\[Lambda]Rik) + xt^3*(\[Lambda]Lik + 10*yk^2*\[Lambda]Lik + 
           8*yk*\[Lambda]Rik))))*CC[8])/(2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*(Qk^2*(-1 + xt)*xt*(2*xs^2*(xt^2 - yk^2 + yS^2) + 
       (-1 + xt)^2*xt*(xt - 2*yk^2 + 2*yS^2) + xs*(-1 + xt)*
        (3*xt^2 - 2*yk^2 + 2*yS^2 + xt*(-1 - 2*yk^2 + 2*yS^2))) + 
     Qj^2*((-1 + xt)^3*xt*(xt - yk^2 + yS^2)^2 + xs*(-1 + xt)^2*
        (3*xt^3 + (yk^2 - yS^2)^2 + xt*(yk^2 - yS^2)^2 + 
         xt^2*(-1 - 4*yk^2 + 2*yS^2)) + xs^3*(xt^3 + (yk^2 - yS^2)^2 - 
         xt^2*(1 + 2*yS^2) - xt*(yk^4 + 2*yS^2 - 2*yk^2*yS^2 + yS^4)) + 
       xs^2*(3*xt^4 - 2*(yk^2 - yS^2)^2 - xt^3*(5 + 2*yk^2 + 2*yS^2) + 
         xt*(3*yk^4 + 2*yS^2 - 6*yk^2*yS^2 + 3*yS^4) - 
         xt^2*(-2 + yk^4 + yS^4 - 2*yk^2*(1 + yS^2)))) + 
     Qj*Qk*(-((-1 + xt)^3*xt*(2*xt^2 - 4*xt*(yk^2 - yS^2) + 
          (yk^2 - yS^2)^2)) + xs^3*(-xt^3 - (yk^2 - yS^2)^2 + 
         xt^2*(1 + 2*yS^2) + xt*(yk^4 + 2*yS^2 - 2*yk^2*yS^2 + yS^4)) + 
       xs^2*(-5*xt^4 + 2*(yk^2 - yS^2)^2 + xt^3*(7 + 2*yk^2 + 2*yS^2) + 
         xt^2*(-2 + yk^4 - 2*yS^2 - 2*yk^2*yS^2 + yS^4) - 
         xt*(3*yk^4 + 3*yS^4 + yk^2*(2 - 6*yS^2))) - 
       xs*(-1 + xt)^2*(6*xt^3 + (yk^2 - yS^2)^2 + 
         xt^2*(-2 - 6*yk^2 + 4*yS^2) + xt*(yk^4 - 2*yk^2*(1 + yS^2) + 
           yS^2*(2 + yS^2)))))*\[Lambda]Lik*\[Lambda]Ljk*CC[9])/
   (2*xs^3*xt^3*(-1 + xs + xt)^2) + 
  (mi^2*Qk*(-(Qk*(-1 + xt)*xt*((-1 + xt)^2*xt*(xt + 2*yk^2 - 2*yS^2) + 
        xs*(-1 + xt^2)*(xt + 2*yk^2 - 2*yS^2) + 2*xs^2*(xt + yk^2 - yS^2))) + 
     Qj*(-((-1 + xt)^3*xt*(yk^2 - yS^2)^2) + xs*(-1 + xt)^2*
        (2*xt^2*yk^2 - (yk^2 - yS^2)^2 - xt*(yk^2 - yS^2)^2) + 
       xs^3*(xt^3 + xt^2*(-1 + 4*yk^2 - 2*yS^2) - (yk^2 - yS^2)^2 + 
         xt*(yk^4 + 2*yS^2 - 2*yk^2*yS^2 + yS^4)) + 
       xs^2*(xt^4 + xt^3*(-1 + 6*yk^2 - 2*yS^2) + 2*(yk^2 - yS^2)^2 - 
         xt*(3*yk^4 + 2*yS^2 - 6*yk^2*yS^2 + 3*yS^4) + 
         xt^2*(yk^4 - 2*yk^2*(3 + yS^2) + yS^2*(4 + yS^2)))))*\[Lambda]Lik*
    \[Lambda]Ljk*CC[10])/(2*xs^3*xt^3*(-1 + xs + xt)^2) - 
  (mi^2*Qk^2*(xs + xt)*(xs^3 + xs*xt^2 + xs^2*(xt + 2*yk^2 - 2*yS^2) - 
     xt^2*(xt + 2*yk^2 - 2*yS^2))*\[Lambda]Lik*\[Lambda]Ljk*CC[11])/
   (2*xs^3*xt^3) + (mi^2*(Qj - Qk)^2*(xs - xt)*(xs + xt)^2*
    (xs + xt - 2*yk^2 + 2*yS^2)*\[Lambda]Lik*\[Lambda]Ljk*CC[12])/
   (2*xs^3*xt^3) - (mi^4*(Qj - Qk)*Qk*(xs*(xt - (yk - yS)^2) - 
     (-1 + xt)*(yk - yS)^2)*(-((-1 + xt)*xt*(yk^2 - yS^2)) + 
     xs^2*(xt + yk^2 - yS^2) + xs*(xt^2 - yk^2 + yS^2))*
    (-((-1 + xt)*(yk + yS)^2) + xs*(xt - (yk + yS)^2))*\[Lambda]Lik*
    \[Lambda]Ljk*DD[1])/(2*xs^3*xt^3*(-1 + xs + xt)^2) - 
  (mi^4*(Qj - Qk)*Qk*(xs*(xt - (yk - yS)^2) - (-1 + xt)*(yk - yS)^2)*
    (-((-1 + xt)*xt*(yk^2 - yS^2)) + xs^2*(xt + yk^2 - yS^2) + 
     xs*(-2*xt + xt^2 - yk^2 + yS^2))*(-((-1 + xt)*(yk + yS)^2) + 
     xs*(xt - (yk + yS)^2))*\[Lambda]Lik*\[Lambda]Ljk*DD[2])/
   (2*xs^3*xt^3*(-1 + xs + xt)^2) - 
  (mi^4*(Qj - Qk)^2*(xs^2 + xt*(yk^2 - yS^2) + xs*(xt - yk^2 + yS^2))*
    (xs^3 + (-1 + xt)*(yk^2 - yS^2)^2 + xs^2*(-1 + xt - 2*yk^2 + 2*yS^2) + 
     xs*(yk^4 + yS^2*(-2 - 2*xt + yS^2) - 2*yk^2*(-1 + xt + yS^2)))*
    \[Lambda]Lik*\[Lambda]Ljk*DD[3])/(2*xs^3*xt^3) + 
  (mi^4*Qk^2*\[Lambda]Ljk*(xs^5*\[Lambda]Lik + 
     xs^4*(-1 + 2*xt + 3*yk^2 - 3*yS^2)*\[Lambda]Lik - 
     (-1 + xt)*xt*(yk^2 - yS^2)^3*\[Lambda]Lik - xs*(yk^2 - yS^2)*
      (xt*(yk^2 - yS^2) + (yk^2 - yS^2)^2 + xt^2*(-5*yk^2 + yS^2))*
      \[Lambda]Lik + xs^3*(xt^2 + xt*(-1 + 2*yk^2 - 6*yS^2) + 
       3*(yk^4 + yS^2 + yS^4 - yk^2*(1 + 2*yS^2)))*\[Lambda]Lik - 
     xs^2*(-((-3 + yk^2 - yS^2)*(yk^2 - yS^2)^2*\[Lambda]Lik) + 
       xt*(yk^2 - yS^2)*(3 + 4*yS^2)*\[Lambda]Lik + 
       xt^2*(9*yk^2*\[Lambda]Lik + 3*yS^2*\[Lambda]Lik + 4*yk*\[Lambda]Rik)))*
    DD[4])/(2*xs^3*xt^3) + (mi^4*(Qj - Qk)^2*(xs*(xt + yk^2 - yS^2) + 
     xt*(xt - yk^2 + yS^2))*(xt^3 + (-1 + xs)*(yk^2 - yS^2)^2 + 
     xt^2*(-1 + xs - 2*yk^2 + 2*yS^2) + xt*(yk^4 + yS^2*(-2 - 2*xs + yS^2) - 
       2*yk^2*(-1 + xs + yS^2)))*\[Lambda]Lik*\[Lambda]Ljk*DD[5])/
   (2*xs^3*xt^3) + (mi^4*Qk^2*(xs - xt)*(xt + yk^2 - yS^2)*
    (xt^3 + xt^2*(-1 + xs + 2*yk^2 - 2*yS^2) + (-1 + xs)*(yk^2 - yS^2)^2 + 
     xt*(yk^4 + yS^2*(2 - 2*xs + yS^2) - 2*yk^2*(1 + xs + yS^2)))*
    \[Lambda]Lik*\[Lambda]Ljk*DD[6])/(2*xs^3*xt^3), 
 (-2*Qj*(-2*Qk + Qj*(1 + yk^2 - yS^2))*\[Lambda]Rik*\[Lambda]Rjk)/
   ((-1 + xs)*xs^2*xt) + (2*Qj^2*yk^2*\[Lambda]Rik*\[Lambda]Rjk*BB[1])/
   (xs^2*xt - xs^3*xt) + (2*Qj^2*yS^2*\[Lambda]Rik*\[Lambda]Rjk*BB[2])/
   ((-1 + xs)*xs^2*xt) + 
  (2*Qj*(2*Qk*xs*(1 - 2*xs + xs^2 + xt - xt^2)*\[Lambda]Rik + 
     Qj*(-1 + xt)*(-1 + xs + xt)*(2*(-1 + xs)*yk*\[Lambda]Lik + 
       (-2 + xs)*yk^2*\[Lambda]Rik + (xs + 2*yS^2 - xs*yS^2)*\[Lambda]Rik))*
    \[Lambda]Rjk*BB[3])/((-1 + xs)^2*xs^2*(-1 + xt)*xt*(-1 + xs + xt)) - 
  (2*Qj*(-2*Qk*xs*xt*\[Lambda]Rik + Qj*(-1 + xs + xt)*
      (2*(-1 + xs)*yk*\[Lambda]Lik - yk^2*\[Lambda]Rik + 
       (xs + yS^2)*\[Lambda]Rik))*\[Lambda]Rjk*BB[4])/
   ((-1 + xs)^2*xs^2*xt*(-1 + xs + xt)) - 
  (4*Qj*Qk*\[Lambda]Rik*\[Lambda]Rjk*BB[5])/(xs*(-1 + xt)*xt*
    (-1 + xs + xt)) - (2*mi^2*Qk*(Qj*(-1 + xs + xt)*yk*\[Lambda]Lik + 
     Qj*xs*\[Lambda]Rik - Qk*xs*\[Lambda]Rik)*\[Lambda]Rjk*CC[1])/
   (xs*xt*(-1 + xs + xt)^2) + (2*mi^2*Qk^2*yk*\[Lambda]Lik*\[Lambda]Rjk*
    CC[3])/(xs^2*(-1 + xs + xt)) - 
  (2*mi^2*(Qj - Qk)*Qk*((-1 + xs + xt)*yk*\[Lambda]Lik + xs*\[Lambda]Rik)*
    \[Lambda]Rjk*CC[4])/(xs^2*(-1 + xs + xt)^2) - 
  (4*mi^2*Qk^2*yk*\[Lambda]Lik*\[Lambda]Rjk*CC[5])/(xs^2*xt) - 
  (2*mi^2*(Qj - Qk)*(2*Qj*(-1 + xs + xt)^2*yS^2*\[Lambda]Rik + 
     Qk*(-1 + xs)^2*((-1 + xs + xt)*yk*\[Lambda]Lik + xs*\[Lambda]Rik))*
    \[Lambda]Rjk*CC[7])/((-1 + xs)*xs^2*xt*(-1 + xs + xt)^2) + 
  (2*mi^2*Qk*yk*(Qk*(-1 + xs)^2*\[Lambda]Lik - 2*Qj*(-1 + xs + xt)*
      ((-1 + xs)*\[Lambda]Lik - yk*\[Lambda]Rik))*\[Lambda]Rjk*CC[8])/
   ((-1 + xs)*xs^2*xt*(-1 + xs + xt)) + 
  (2*mi^2*Qk*(Qj*(2*xs^2 + 3*xs*(-1 + xt) + (-1 + xt)^2)*yk*\[Lambda]Lik - 
     Qj*xs*(-1 + xt)*\[Lambda]Rik + Qk*xs*(-1 + xt)*\[Lambda]Rik)*
    \[Lambda]Rjk*CC[10])/(xs^2*xt*(-1 + xs + xt)^2) + 
  (2*mi^4*(Qj - Qk)*Qk*(-((-1 + xt)^2*yk*(yk^2 - yS^2)*\[Lambda]Lik) + 
     xs*(-1 + xt)*(xt*yk*\[Lambda]Lik - 2*yk^3*\[Lambda]Lik + 
       2*yk*yS^2*\[Lambda]Lik - yk^2*\[Lambda]Rik - yS^2*\[Lambda]Rik) + 
     xs^2*(-(yk^3*\[Lambda]Lik) + yk*yS^2*\[Lambda]Lik - yk^2*\[Lambda]Rik - 
       yS^2*\[Lambda]Rik + xt*(yk*\[Lambda]Lik + \[Lambda]Rik)))*\[Lambda]Rjk*
    DD[1])/(xs^2*xt*(-1 + xs + xt)^2) - 
  (2*mi^4*Qk^2*yk*(xs^2*\[Lambda]Lik - (-1 + xt)*(yk^2 - yS^2)*\[Lambda]Lik + 
     xs*(-1 + xt - yk^2 + yS^2)*\[Lambda]Lik - 2*xs*yk*\[Lambda]Rik)*
    \[Lambda]Rjk*DD[4])/(xs^2*xt*(-1 + xs + xt)) + 
  (2*mi^4*Qk^2*yk*(xt^2*\[Lambda]Lik + (-1 + xs)*yk^2*\[Lambda]Lik - 
     (-1 + xs)*yS^2*\[Lambda]Lik + xt*(-1 + xs + yk^2 - yS^2)*\[Lambda]Lik + 
     2*xs*yk*\[Lambda]Rik)*\[Lambda]Rjk*DD[6])/(xs^2*xt*(-1 + xs + xt)), 
 (-2*Qj*(xs - xt)*(-2*Qk*xs*xt + Qj*(-1 + xt)*(yk^2 - yS^2) + 
     Qj*xs*(xt + yk^2 - yS^2))*\[Lambda]Lik*\[Lambda]Ljk)/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)) + 
  (2*Qj^2*(-xs + xt)*yk^2*\[Lambda]Lik*\[Lambda]Ljk*BB[1])/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2) + 
  (2*Qj^2*(xs - xt)*yS^2*\[Lambda]Lik*\[Lambda]Ljk*BB[2])/
   ((-1 + xs)*xs^2*(-1 + xt)*xt^2) - 
  (2*Qj*(-2*xs + xs^2 - (-2 + xt)*xt)*(2*Qk + Qj*(-1 + yk^2 - yS^2))*
    \[Lambda]Lik*\[Lambda]Ljk*BB[3])/((-1 + xs)^2*xs*(-1 + xt)^2*xt*
    (-1 + xs + xt)) + 
  (2*Qj*(-2*Qk*xs + Qj*(xs + yk^2 - 2*xs*yk^2 - yS^2 + 2*xs*yS^2))*
    \[Lambda]Lik*\[Lambda]Ljk*BB[4])/((-1 + xs)^2*xs^2*xt*(-1 + xs + xt)) + 
  (2*Qj*(2*Qk*xt + Qj*(-yk^2 + yS^2 + xt*(-1 + 2*yk^2 - 2*yS^2)))*
    \[Lambda]Lik*\[Lambda]Ljk*BB[5])/(xs*(-1 + xt)^2*xt^2*(-1 + xs + xt)) - 
  (2*mi^2*Qk*(Qk*xs*((-1 + xs)*(yk^2 - yS^2) + xt*(1 + yk^2 - yS^2)) + 
     Qj*((-1 + xt)*xt*yk^2 - xs*yS^2 + xs^2*yS^2 + xs*xt*(-1 + yk^2 + yS^2)))*
    \[Lambda]Lik*\[Lambda]Ljk*CC[1])/(xs*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^2*(Qj^2*(xs^2 + (-1 + xt)*xt + xs*(-1 + 2*xt))*yS^2 + 
     Qk^2*xs*(-((-1 + xs)*(yk^2 - yS^2)) + xt*(1 - yk^2 + yS^2)) + 
     Qj*Qk*(-((-1 + xt)*xt*yS^2) + xs^2*(yk^2 - 2*yS^2) + 
       xs*(-yk^2 + 2*yS^2 + xt*(-1 + yk^2 - 3*yS^2))))*\[Lambda]Lik*
    \[Lambda]Ljk*CC[2])/(xs*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^2*Qk*(Qk*xt*((-1 + xt)*(yk^2 - yS^2) + xs*(1 + yk^2 - yS^2)) + 
     Qj*(-(xs*yk^2) + xs^2*yk^2 + (-1 + xt)*xt*yS^2 + 
       xs*xt*(-1 + yk^2 + yS^2)))*\[Lambda]Lik*\[Lambda]Ljk*CC[3])/
   (xs^2*xt*(-1 + xs + xt)^3) - 
  (2*mi^2*(Qj^2*(xs^2 + (-1 + xt)*xt + xs*(-1 + 2*xt))*yS^2 + 
     Qk^2*xt*(-((-1 + xt)*(yk^2 - yS^2)) + xs*(1 - yk^2 + yS^2)) + 
     Qj*Qk*(-(xs^2*yS^2) + (-1 + xt)*xt*(yk^2 - 2*yS^2) + 
       xs*(yS^2 + xt*(-1 + yk^2 - 3*yS^2))))*\[Lambda]Lik*\[Lambda]Ljk*CC[4])/
   (xs^2*xt*(-1 + xs + xt)^3) + (4*mi^2*Qk^2*(xs - xt)*yk^2*\[Lambda]Lik*
    \[Lambda]Ljk*CC[5])/(xs^2*xt^2*(-1 + xs + xt)) - 
  (4*mi^2*(Qj - Qk)^2*(xs - xt)*yS^2*\[Lambda]Lik*\[Lambda]Ljk*CC[6])/
   (xs^2*xt^2*(-1 + xs + xt)) + 
  (2*mi^2*(Qj - Qk)*(Qj*(xs^4 - (-1 + xt)*xt + xs^3*(-3 + 4*xt) + 
       xs^2*(3 - 7*xt + 5*xt^2) + xs*(-1 + 2*xt - 4*xt^2 + 2*xt^3))*yS^2 - 
     Qk*(-1 + xs)^2*xt*(xs*(-1 + yk^2 - yS^2) + (-1 + xt)*(yk^2 - yS^2)))*
    \[Lambda]Lik*\[Lambda]Ljk*CC[7])/((-1 + xs)*xs^2*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^2*Qk*(Qk*(-1 + xs)^2*xt*((-1 + xt)*(yk^2 - yS^2) + 
       xs*(1 + yk^2 - yS^2)) - Qj*(xs^4*yk^2 - (-1 + xt)*xt*yS^2 + 
       xs^3*(xt - 3*yk^2 + 5*xt*yk^2 - xt*yS^2) + 
       xs*(-yk^2 + 2*xt^3*yk^2 + xt*(1 + 5*yk^2 - 3*yS^2) + 
         xt^2*(-6*yk^2 + 2*yS^2)) + xs^2*(3*yk^2 + xt^2*(6*yk^2 - yS^2) + 
         xt*(-2 - 10*yk^2 + 3*yS^2))))*\[Lambda]Lik*\[Lambda]Ljk*CC[8])/
   ((-1 + xs)*xs^2*xt^2*(-1 + xs + xt)^3) - 
  (2*mi^2*(Qj - Qk)*(Qj*(2*xs^3*xt + (-1 + xt)^3*xt + 
       xs*(-1 + xt)^2*(1 + 4*xt) + xs^2*(-1 - 4*xt + 5*xt^2))*yS^2 - 
     Qk*xs*(-1 + xt)^2*(xt*(-1 + yk^2 - yS^2) + (-1 + xs)*(yk^2 - yS^2)))*
    \[Lambda]Lik*\[Lambda]Ljk*CC[9])/(xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^2*Qk*(-(Qk*xs*(-1 + xt)^2*((-1 + xs)*(yk^2 - yS^2) + 
        xt*(1 + yk^2 - yS^2))) + Qj*(2*xs^3*xt*yk^2 + (-1 + xt)^3*xt*yk^2 + 
       xs*(-1 + xt)^2*(xt + 5*xt*yk^2 + yS^2 - xt*yS^2) + 
       xs^2*(-1 + xt)*(6*xt*yk^2 + yS^2 - xt*yS^2)))*\[Lambda]Lik*
    \[Lambda]Ljk*CC[10])/(xs^2*(-1 + xt)*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^4*(Qj - Qk)*Qk*(-((-1 + xt)^2*xt*yk^2*(yk^2 - yS^2)) + 
     xs^3*yS^2*(xt + yk^2 - yS^2) + xs^2*(-2*yk^2*yS^2 + 2*yS^4 + 
       xt^2*(-1 + yk^2 + yS^2) + xt*(yk^2 - yk^4 + 3*yk^2*yS^2 - 2*yS^4)) + 
     xs*(xt^3*yk^2 + yS^2*(yk^2 - yS^2) + xt^2*(-2*yk^4 + yS^2 + 
         3*yk^2*yS^2 - yS^4) + xt*(2*yk^4 + yS^2*(-1 + 2*yS^2) - 
         yk^2*(1 + 4*yS^2))))*\[Lambda]Lik*\[Lambda]Ljk*DD[1])/
   (xs^2*xt^2*(-1 + xs + xt)^3) - 
  (2*mi^4*(Qj - Qk)*Qk*((-1 + xt)^2*xt*yS^2*(yk^2 - yS^2) + 
     xs^3*yk^2*(xt - yk^2 + yS^2) + xs^2*(2*yk^2*(yk^2 - yS^2) + 
       xt^2*(-1 + yk^2 + yS^2) + xt*(-2*yk^4 + yS^2 + 3*yk^2*yS^2 - yS^4)) + 
     xs*(-yk^4 + xt^3*yS^2 + yk^2*yS^2 + xt^2*(-yk^4 - 2*yS^4 + 
         yk^2*(1 + 3*yS^2)) + xt*(2*yk^4 + yS^2*(-1 + 2*yS^2) - 
         yk^2*(1 + 4*yS^2))))*\[Lambda]Lik*\[Lambda]Ljk*DD[2])/
   (xs^2*xt^2*(-1 + xs + xt)^3) + 
  (2*mi^4*(Qj - Qk)^2*yS^2*(xs^2 + xt*(yk^2 - yS^2) + xs*(xt - yk^2 + yS^2))*
    \[Lambda]Lik*\[Lambda]Ljk*DD[3])/(xs^2*xt^2*(-1 + xs + xt)) - 
  (2*mi^4*Qk^2*yk^2*(xs^2 + xs*(xt + yk^2 - yS^2) + xt*(-yk^2 + yS^2))*
    \[Lambda]Lik*\[Lambda]Ljk*DD[4])/(xs^2*xt^2*(-1 + xs + xt)) - 
  (2*mi^4*(Qj - Qk)^2*yS^2*(xs*(xt + yk^2 - yS^2) + xt*(xt - yk^2 + yS^2))*
    \[Lambda]Lik*\[Lambda]Ljk*DD[5])/(xs^2*xt^2*(-1 + xs + xt)) + 
  (2*mi^4*Qk^2*yk^2*(xt*(xt + yk^2 - yS^2) + xs*(xt - yk^2 + yS^2))*
    \[Lambda]Lik*\[Lambda]Ljk*DD[6])/(xs^2*xt^2*(-1 + xs + xt))};


(* ::Text:: *)
(* (*tau contribution to the squared amplitude for the decay t->c\[Gamma]\[Gamma]. The muon contribution is negligible given the constraints on the corresponding couplings but it can be added straightforwardly if required, just uncomment the commented parts although the evaluation can be too slow. The notation used for the LQ couplings is that used in Reference arXiv:2203.10111*)*)


EvalM2=Compile[{{xxs,_Real},{xxt,_Real},{mS,_Real},{YRL2\[Tau],_Real},
{YLR2\[Tau],_Real},{YRL3\[Tau],_Real},{YLR3\[Tau],_Real}(*,{YRL2\[Mu],_Real},
{YLR2\[Mu],_Real},{YRL3\[Mu],_Real},{YLR3\[Mu],_Real}*)},
yyS=mS/mi;
(*tau contribution*)
sust1={xs->xxs,xt->xxt,xu->1-xxs-xxt,yS->yyS};
sust\[Tau]={yk->y\[Tau],\[Lambda]Lik->YRL3\[Tau],\[Lambda]Ljk->YRL2\[Tau],\[Lambda]Rik->YLR3\[Tau],\[Lambda]Rjk->YLR2\[Tau]};
sust=Union[sust1,sust\[Tau]];
{BB1[1],BB1[2],BB[1],BB[2],BB[3],BB[4],BB[5],BB[6],BB[7],CC[1],CC[2],CC[3],CC[4],CC[5],CC[6],CC[7],CC[8],CC[9],CC[10],CC[11],CC[12],DD[1],DD[2],DD[3],DD[4],DD[5],DD[6]}=funesc/.sust;
{F1,F2,F3,F4,F5,F6,FF1,FF2,FF3,FF4,FF5,FF6}=coef1/.sust;
{H1,H2,H3,H4,H5,H6,HH1,HH2,HH3,HH4,HH5,HH6}=coef1b/.sust;
(*muon contribution
sust\[Mu]={yk\[Rule]y\[Mu],\[Lambda]Lik\[Rule]YRL3\[Mu],\[Lambda]Ljk\[Rule]YRL2\[Mu],\[Lambda]Rik\[Rule]YLR3\[Mu],\[Lambda]Rjk\[Rule]YLR2\[Mu]};
sust=Union[sust1,sust\[Mu]];
{BB1[1],BB1[2],BB[1],BB[2],BB[3],BB[4],BB[5],BB[6],BB[7],CC[1],CC[2],CC[3],CC[4],CC[5],CC[6],CC[7],CC[8],CC[9],CC[10],CC[11],CC[12],DD[1],DD[2],DD[3],DD[4],DD[5],DD[6]}=funesc/.sust;
{F1,F2,F3,F4,F5,F6,FF1,FF2,FF3,FF4,FF5,FF6}={F1,F2,F3,F4,F5,F6,FF1,FF2,FF3,FF4,FF5,FF6}+(coef1/.sust);
{H1,H2,H3,H4,H5,H6,HH1,HH2,HH3,HH4,HH5,HH6}={H1,H2,H3,H4,H5,H6,HH1,HH2,HH3,HH4,HH5,HH6}+(coef1b/.sust);*)
Nc=3;
fac=\[Alpha] Nc/(4Pi);
res=fac^2 M2/.sust1;
res]


(* ::Text:: *)
(* (*Auxiliary function to numerical evaluation of the squared amplitude for the decay t->c\[Gamma]\[Gamma] over the t hat variable. Uncomment to add the muon contribution*)*)


 Fun[mS_,YRL2\[Tau]_,YLR2\[Tau]_,YRL3\[Tau]_,YLR3\[Tau]_,(*,{YRL2\[Mu],_Real},
{YLR2\[Mu],_Real},{YRL3\[Mu],_Real},{YLR3\[Mu],_Real}*)xs_,limsup_?NumberQ]:=(
liminf=0.1;
NIntegrate[EvalM2[xs,xt,mS,YRL2\[Tau],YLR2\[Tau],YRL3\[Tau],YLR3\[Tau](*,YRL2\[Mu],YLR2\[Mu],YRL3\[Mu],YLR3\[Mu]*)],{xt,0.1,limsup},Method->{"GlobalAdaptive","SymbolicProcessing"->0},AccuracyGoal->5]
)


(* ::Text:: *)
(* (*Decay width for the process t->c\[Gamma]\[Gamma]. Uncomment to add the muon contribution*)*)


 DecayWidth=Compile[{{mS,_Real},{YRL2\[Tau],_Real},{YLR2\[Tau],_Real},{YRL3\[Tau],_Real},{YLR3\[Tau],_Real}(*,{YRL2\[Mu],_Real},
{YLR2\[Mu],_Real},{YRL3\[Mu],_Real},{YLR3\[Mu],_Real}*)},
mi/(256Pi^3)NIntegrate[Fun[mS,YRL2\[Tau],YLR2\[Tau],YRL3\[Tau],YLR3\[Tau](*,YRL2\[Mu],YLR2\[Mu],YRL3\[Mu],YLR3\[Mu]*),xs,0.9-xs],{xs,0.1,0.8},Method->{"GlobalAdaptive","SymbolicProcessing"->0},AccuracyGoal->5]]
