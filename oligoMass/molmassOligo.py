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
        self.br_alphabet = '/ [ ] { }'.split(' ')
        self.ex_mod = {}
        self.ex_mod_key = ''

        self.mod_formula = {}
        self.mod_formula['+'] = ['CO', '']
        self.mod_formula['*'] = ['S', 'O']
        self.mod_formula['m'] = ['OCH2', '']
        self.mod_formula['r'] = ['O', '']

        self.mod_list = {}
        for k in self.mod_formula.keys():
            self.mod_list[k] = 0

    def _add_mod(self, letter):
        if not self.modRead:
            if letter in self.mod_alphabet:
                self.modRead = False
                self.mod_list[letter] += 1

            elif letter in self.br_alphabet:
                if not self.modRead:
                    self.modRead = True
                    self.ex_mod_key = ''
        else:
            if letter in self.br_alphabet:
                self.modRead = False
                if self.ex_mod_key in list(self.ex_mod.keys()):
                    self.ex_mod[self.ex_mod_key] += 1
                else:
                    self.ex_mod[self.ex_mod_key] = 1
            else:
                self.ex_mod_key += letter


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
        self.is_mixed = False
        self.alphabet = 'A G C T U a g c t u R Y M K ' \
                        'S W H B V D N'.split(' ')
        self.mixed_alphabet = 'R Y M K ' \
                        'S W H B V D N'.split(' ')

        self.modifications = oligoNAModifications()
        self.sequence_parser()

        self.dnaDB = dna.deoxynusleosideDB()

    def set_mixed_na(self):
        self.mixed_na = {}
        self.mixed_na['R'], self.mixed_na['Y'] = ['A', 'G'], ['C', 'T']
        self.mixed_na['M'], self.mixed_na['K'] = ['A', 'C'], ['G', 'T']
        self.mixed_na['S'], self.mixed_na['W'] = ['C', 'G'], ['A', 'T']
        self.mixed_na['H'], self.mixed_na['B'] = ['A', 'C', 'T'], ['C', 'G', 'T']
        self.mixed_na['V'], self.mixed_na['D'] = ['A', 'C', 'G'], ['A', 'G', 'T']
        self.mixed_na['N'] = ['A', 'C', 'T', 'G']

    def sequence_parser(self):
        seq_list = list(self.init_seq)
        self.seq = ''
        for letter in seq_list:
            if letter in self.alphabet and not self.modifications.modRead:
                if letter in self.mixed_alphabet:
                    self.is_mixed = True
                self.seq += letter.upper()
            else:
                self.modifications._add_mod(letter)
        return self.seq

    def getMolecularFormula(self):
        f = ''
        if not self.is_mixed:
            for n in self.seq:
                f += self.dnaDB(n).seqformula
            f += f'(HPO2){len(self.seq) - 1}H2'
            return self.modifications._get_mod_formula(f)
        else:
            return f

    def getMonoMass(self):
        if not self.is_mixed:
            return mm.Formula(self.getMolecularFormula()).isotope
        else:
            return 0

    def getAvgMass(self):
        if not self.is_mixed:
            return mm.Formula(self.getMolecularFormula()).mass
        else:
            return 0


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

    o1 = oligoNASequence('GT}test_mod{AR')
    print(o1.getMolecularFormula())
    print(o1.getAvgMass())
    print(o1.modifications.ex_mod)



if __name__ == '__main__':
    test()