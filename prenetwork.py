import numpy as np
import re


# given human readable reaction network output as np arrays zetas, rate_table
# species is list, reaction is list of tuples as (reactants, products, rate) as string, string, float.
# rate_indices is a N x K matrix indicating the rate and the indices of the corresponding concentrations
class preNetwork:
    def __init__(self, species, reactions):
        n, k = len(species), len(reactions)
        self.zeta_table = np.zeros((k, n), dtype = np.int32)
        self.rate_table = np.zeros((k, n), dtype = np.uint16)
        self.rates = np.zeros(k)
        self.species_dict = {item: index for index, item in enumerate(species)}

        for reaction_num in range(k):
            reactants, products, rate = reactions[reaction_num]
            self.rates[reaction_num] = rate
            
            # find reactants, add to zeta
            self.addSideToZeta(reactants, reaction_num, sign = -1)
            self.addSideToZeta(products, reaction_num, sign = 1)

    def addSideToZeta(self, side, reaction_num, sign = -1):
        zeta_indices = self.parse(side)
        for species_index, amount in zeta_indices:
            self.zeta_table[reaction_num, species_index] = amount * sign
            self.rate_table[reaction_num, species_index] = amount * (sign == -1)

    def parse(self, string):
        zeta_indices = []
        complexes_in_side = string.replace(" ","").split("+")
        for complex in complexes_in_side:
            zeta_indices.append(self.complex_to_index(complex))
        return zeta_indices

    def complex_to_index(self, complex):
        numbers_as_strings = re.findall(r'\d+', complex)
        amount = 0
        if numbers_as_strings:
            amount = int(numbers_as_strings[0])
        else:
            amount = 1
        species_index = self.species_dict[re.sub(r'[^a-zA-Z]', '', complex)]
        return species_index, amount
    
    def get_tables(self):
        return self.zeta_table, self.rate_table, self.rates