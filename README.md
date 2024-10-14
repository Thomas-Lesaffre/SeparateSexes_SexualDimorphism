# The joint evolution of separate sexes and sexual dimorphism

This repository contains the code files pertaining to our article _The joint evolution of separate sexes and sexual dimorphism_ which has been accepted in _Journal of Evolutionary Biology_.

Authors: Thomas Lesaffre, John R. Pannell and Charles Mullon

Affiliation: Department of Ecology and Evolution, University of Lausanne (Switzerland)

Contact: thomas.lesaffre[at]unil.ch / thomaslesaffre.evolbiol[at]gmail.com

The _MathematicaNotebooks_ folder contains the three _Mathematica_ notebooks that are mentioned in the Appendix.

The _SimulationPrograms_ folder contains four simulation programs coded in C++11 (see ReadMe file in the folder for details):
* _Simulations_HaploidRecombinationSuppression_: program described in Appendix B.1.3, where recombination evolves jointly with allelic effects at both loci. The program described in Appendix A.3 can be obtained from this program by setting the mutation rate parameter <code>par.ur</code> and the initial recombination rate <code>par.r0</code> to zero.
* _Simulations_HaploidConditionalExpression_: program described in Appendix B.2.3 for the evolution of sex allocation-dependent expression of the conflict trait.
* _Simulations_DiploidBaseline_: program described in Appendix E.1.2.1 for the joint evolution of dominance and allelic effects at the loci encoding sex allocation and the conflict trait in diploids.
