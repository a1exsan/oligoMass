Description
===========

Calculation of oligonucleotides weights and some properties.

Installation
============

from Repository on GitHub:

`pip install git+https://github.com/a1exsan/oligoMass.git#egg=oligoMass`

from PyPi:

`pip install oligoMass -U`
https://pypi.org/project/oligoMass/

Getting started
===============

#### 1) Calculate average and monoisotopic molecular weight of any deoxyoligonucleotide sequense consists of 
#### five type bases: dA, dT, dG, dC, dU

`from oligoMass.molmassOligo import oligoNASequence`

`oligos = oligoNASequence('ATGCCGGTuugtU')`

`print(f"molecular formula: {oligos.getMolecularFormula()}")`

`print(f"MW: {oligoNASequence('ATGCC').getAvgMass()}")`

`print(f"monoMW: {oligoNASequence('ATGCC').getMonoMass()}")`

#### 2) Calculate average and monoisotopic molecular weight of some deoxyoligonucleotide modifications 
#### such as: LNA(+), ribo(r), 2'-O-Me (m) or Phosphorothioated DNA bases as entered (*)

`print(f" MW: {oligoNASequence('+ATrG+CCmG*GTuu+gtU').getMonoMass()}")`

#### 3) getPrefix(index) and getSuffix(index) methods can return prefix and suffix of the sequence 

`print(oligoNASequence('+ATrG+CCmG*GTuu+gtU').getPrefix(index=2)())`

`# +AT`

`print(oligoNASequence('+ATrG+CCmG*GTuu+gtU').getSuffix(index=4)())`

`# CmG*GTuu+gtU`

#### 4)Calc extinction coef:

`from oligoMass import dna`

`extinction = dna.get_simple_ssdna_extinction("ATGCTT", dna.get_extinction_dict())`

`print(extinction)`

`# 55100`