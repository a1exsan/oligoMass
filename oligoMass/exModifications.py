import pandas as pd


class exModifDB():
    def __init__(self):
        self.data = None

    def get_mod_properties(self, mod_code):
        pass

class exModifDBtest(exModifDB):
    def __init__(self):
        super().__init__()

        self.data = {'code': [], 'mass': [], 'ext_cf': []}
        self.data['code'].append('5Phos')
        self.data['mass'].append(80)
        self.data['ext_cf'].append(0)

        self.data = pd.DataFrame(self.data)
        self.data = self.data.set_index('code')

    def get_mod_properties(self, mod_code):
        if mod_code in self.data['code']:
            for d in self.data:
                if mod_code == '5Phos':
                    pass
