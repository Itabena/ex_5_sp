import numpy as np
import matplotlib.pyplot as plt

class particle:
    """
    :param s: spin of the particle ,val- 1 or -1
    :param epsilon: energy of the particle
    """
    def __init__(self,s,epsilon):
      self.s=s
      self.epsilon=epsilon
      return
    def flip(self):
        self.s = -self.s
        return self.s
class lattice:
    def __init__(self,n):
        self.n=n
        self.lattice = np.zeros((n,n))
        return
    def Scatter(self):
        for i in range(self.n):
            for j in range(self.n):
                self.lattice[i,j] = particle(np.random.randint(2)*2-1)
        return self.lattice
    def flip(self,i,j):
        self.lattice[i,j].flip()
        return self.lattice




