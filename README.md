Description
===========

Oligonucleotides properties calculation packege.

Installation
============

from Repository on GitHub:

`pip install git+https://github.com/a1exsan/oligoMass.git#egg=oligoMass`

from PyPi:

`pip install oligoMass -U`

Getting started
==============

####1) Calculate average and monoisotopic molecular weight of any deoxyoligonucleotide consists of 
####five type bases: dA, dT, dG, dC, dU

`from oligoMass.molmassOligo import oligoNASequence`

`oligos = oligoNASequence('ATGCCGGTuugtU')`
`print(f"molecular formula: {oligos.getMolecularFormula()}")`
`print(f"MW: {oligoNASequence('ATGCC').getAvgMass()}")`
`print(f"monoMW: {oligoNASequence('ATGCC').getMonoMass()}")`

####2) Calculate average and monoisotopic molecular weight of some deoxyoligonucleotide modifications 
####such as: LNA(+), ribo(r), 2'-O-Me (m) or Phosphorothioated DNA bases as entered (*)

`print(f" MW: {oligoNASequence('+ATrG+CCmG*GTuu+gtU').getMonoMass()}")`
