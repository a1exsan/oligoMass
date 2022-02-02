import molmass as mm
import dna

class oligoModifications():
    def __init__(self):
        self.modRead = False
        self.mod_formula = None
        self.mod_list = None

    def _add_mod(self, letter):
        pass

class oligoNAModifications(oligoModifications):
    def __init__(self):
        super().__init__()

        self.__set_modifications()

    def __set_modifications(self):
        self.mod_alphabet = '+ * m r'.split(' ')

        self.mod_formula = {}
        self.mod_formula['+'] = ['CO', '']
        self.mod_formula['*'] = ['S', 'O']
        self.mod_formula['m'] = ['OCH2', '']
        self.mod_formula['r'] = ['O', '']

        self.mod_list = {}
        for k in self.mod_formula.keys():
            self.mod_list[k] = 0

    def _add_mod(self, letter):
        if letter in self.mod_alphabet:
            self.modRead = False
            self.mod_list[letter] += 1

    def _get_mod_formula(self, formula):
        f_mod, f_ = '', ''
        for k in self.mod_list.keys():
            if self.mod_list[k] > 0:
                f_mod += f'({self.mod_formula[k][0]}){self.mod_list[k]}'
                if self.mod_formula[k][1] != '':
                    f_ += f'({self.mod_formula[k][1]}){self.mod_list[k]}'
        if f_ != '':
            return (mm.Formula(formula) + mm.Formula(f_mod) - mm.Formula(f_)).empirical
        else:
            return (mm.Formula(formula) + mm.Formula(f_mod)).empirical

class oligoSequence():
    def __init__(self, sequence):
        self.init_seq = sequence
        self.seq = None
        self.alphabet = None

    def sequence_parser(self):
        return self.init_seq

class oligoNASequence(oligoSequence):
    def __init__(self, sequence):
        super().__init__(sequence)
        self.alphabet = 'A G C T U a g c t u'.split(' ')
        self.modifications = oligoNAModifications()
        self.sequence_parser()

        self.dnaDB = dna.deoxynusleosideDB()

    def sequence_parser(self):
        seq_list = list(self.init_seq)
        self.seq = ''
        for letter in seq_list:
            if letter in self.alphabet and not self.modifications.modRead:
                self.seq += letter.upper()
            else:
                self.modifications._add_mod(letter)
        return self.seq

    def getMolecularFormula(self):
        f = ''
        for n in self.seq:
            f += self.dnaDB(n).seqformula
        for n in range(len(self.seq) - 1):
            f += 'HPO2'
        f += 'H2'
        return self.modifications._get_mod_formula(f)

    def getMonoMass(self):
        return mm.Formula(self.getMolecularFormula()).isotope

    def getAvgMass(self):
        return mm.Formula(self.getMolecularFormula()).mass


def test():
    olig1 = oligoNASequence('+a+c*grT*mTTuUrurU')
    print(olig1.seq)
    print(olig1.modifications.mod_list)
    print(olig1.getMolecularFormula())
    print(olig1.getAvgMass())


    sequence = '+A*AAUut+tcrgmG*'

    olig2 = oligoNASequence(sequence)
    print(olig2.seq)
    print(olig2.modifications.mod_list)
    print(olig2.getMolecularFormula())
    #print(olig2.getMonoMass())
    print(olig2.getAvgMass())

    seq = dna.oligoSeq(sequence)
    print(seq.getBruttoFormula())
    print(seq.getMolMass())

    print(seq.getMolMass() - olig2.getAvgMass())



if __name__ == '__main__':
    test()