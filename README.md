# TopQuarkDecay
Mathematica Code to compute the rare decay of the top quark into a charm quark and two photons in scalar leptoquark models as presented in manuscript https://arxiv.org/abs/2208.05064
Mathematica code used for the calculations presented in the
manuscript Rare decay t->c gamma gamama via scalar leptoquark
doublets by A. Bolanos, R. Sanchez-Velez and G. Tavares-Velasco
(https://arxiv.org/abs/2208.05064), which are the authors of this
code. All the mathematical expressions used in the calculation are
contained in the above manuscript. Please give the proper citations
if you use this code. The code is useful to obtain bounds on the
couplings of a charge 5/3 scalar leptoquark with the muon (tau) and
the c or t quarks from the current constraints on the muon anomaly
and the tau decay tau -> mu gamma. The results can be feed into the
code that evaluates the leptoquark contribution to the decay width
of the rare top quark decay t->c gamma gamma. The following files
are available:

BoundsDefinitions.wl: code with the definitions used for the
calculation of the leptoquark contributions to the muon anomaly and
the branching ratio for the tau decay tau -> mu gamma, along with
the definition of the numerical values of several constants used
through the calculation. This code requires que LoopTools
(https://feynarts.de/looptools/) package for the numerical
evaluation of Passarino-Veltman scalar functions.

DecayWidthDefinitions.wl: code with the definitions used for the
calculation of the leptoquark contributions to the branching ratio
for the top quark decay t -> c gamma gamma, along with the
definition of the numerical values of several constants used
through the calculation. This code requires que LoopTools
(https://feynarts.de/looptools/) package for the numerical
evaluation of Passarino-Veltman scalar functions.

BoundsSampleProgram.nb: sample code to compute the allowed values
of the leptoquark couplings from the current constraints on the
muon anomaly and the tau decay tau -> mu gamma.

BRTopRareDecay.nb: code to obtain the plot of Fig. 7 of manuscript
https://arxiv.org/abs/2208.05064.

BoundPlots.nb: code to obtain the plots of Figs. 4, 5 and 6 of
manuscript https://arxiv.org/abs/2208.05064.

All files run under Mathematica 12.
