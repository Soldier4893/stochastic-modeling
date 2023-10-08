from numpy import random

class CompartmentManager:
    def __init__(self):
        self.d = {}
        self.l = []
    
    def add(self, element):
        # enforce unique items
        # if self.d[element]:
        #     return
        self.l.append(element)
        self.d[element] = len(self.l)-1

    def remove(self, element):
        i = self.d[element]
        del self.d[element]
        if i != len(self.l):
            self.l[i] = self.l.pop()
            self.d[self.l[i]] = i
        return element

    def sample(self):
        return random.choice(self.l, p = self.generate_distribution)

    # normalize before returning
    def generate_distribution(self):
        pass