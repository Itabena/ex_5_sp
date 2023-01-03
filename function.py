import numpy as np
import matplotlib.pyplot as plt

#region - Assiting Functions
def Calc_energy(lattice,i,j,eta,h):
       p=lattice[i][j].s
       p_right=lattice[(i+1)%len(lattice)][j].s
       p_left=lattice[(i-1)%len(lattice)][j].s
       p_down=lattice[i][(j-1)%len(lattice)].s
       p_up=lattice[i][(j+1)%len(lattice)].s
       epsilon = -p*eta - h *p* (p_up + p_down + p_left + p_right)
       return epsilon
def Choose_to_flip(lattice,i,j,eta,h):
    confirm=0
    epsilon=Calc_energy(lattice,i,j,mu,B,J)
    P=1/(1+np.exp(-2*epsilon)) #assuming Tkb=1
    tmp=np.random.binomial(1,P)
    if tmp==1:
        lattice.flip(i,j)
        confirm=1
    return lattice , confirm
def Calc_net_magenetization(lattice):
    M=0
    for i in range(len(lattice)):
        for j in range(len(lattice)):
            M+=lattice[i][j].s
    return M
def Calc_net_energy(lattice):
    E=0
    for i in range(len(lattice)):
        for j in range(len(lattice)):
            E+=Calc_energy(lattice,i,j,eta,h)
    return E
def Calc_specific_heat(U_mean_tot,U_mean_sq_tot,n):
    cv=(1/n**2)*(U_mean_sq_tot-(U_mean_tot**2))
    return cv
def Sweep(lattice,eta,h):
    flips=0
    for i,j in np.ndindex(lattice.n,lattice.n):
        lattice,confirm = Choose_to_flip(lattice,i,j,eta,h)
        flips+=confirm
    return lattice , flips

#endregion

#region -Classes
class particle:
    """
    :param s: spin of the particle ,val- 1 or -1
    :param epsilon: energy of the particle
    """
    def __init__(self,s):
      self.s=s
      return
    def flip(self):
        """
        :return: flip the spin of the particle
        """
        self.s = -self.s
        return self.s
class lattice:
    def __init__(self,n,E=0,M=0):
        """
        :param n: The size of the lattice
        :param E: The total energy of the lattice
        :param M: The total magnetization of the lattice
        """
        self.n=n
        self.lattice = [[particle(1) for i in range(n)] for j in range(n)]
        self.E=E
        self.M=M
        return
    def Update(self):
        """
        Updates the net energy and magnetization of the lattice
        :return: M total magnetization of the lattice , E total energy of the lattice
        """
        self.M=Calc_net_magenetization(self.lattice)
        self.E=Calc_net_energy(self.lattice)
        return self.M,self.E
    def Scatter(self):
        """
        Creating the initial random scatter of the particles
        :return: lattice form of the scattered particles
        """
        for i in range(self.n):
            for j in range(self.n):
                self.lattice[i][j] = particle(np.random.randint(2)*2-1)
        self.Update()
        return self.lattice
    def flip(self,i,j):
        """
        Flips the spin of the particle in the lattice
        :param i: The column of the particle to flip
        :param j: The row of the particle to flip
        :return: The lattice with the flipped particle
        """
        self.lattice[i][j].flip()
        return self.lattice

#endregion

#region - Main
def Simulate_spin_resrviour(n,eta,h,nsweep,Kinit,limit):

    delta=limit+1
    latt=lattice(n)
    latt.Scatter()
    latt.Update()
    U_sum_tot=0
    U_sum_sq_tot=0
    M_sum_tot=0
    M_sum_sq_tot=0
    flips=0
    sweep_cnt=0
    K=Kinit
    mean_M_half=0
    mean_M=0
    mean_M_s=0
    mean_U=0
    mean_U_s=0

    clr=0



    while delta>limit:

        if K == 10**8:
            print("This is too long for me, I quit!!!")
            break
        latt,flipstmp= Sweep(latt,eta,h)
        flips+=flipstmp
        sweep_cnt+=1

        if sweep_cnt%nsweep==0:
            latt.update()
            U_sum_tot+=latt.E
            U_sum_sq_tot+=latt.E**2
            M_sum_tot+=latt.M
            M_sum_sq_tot+=latt.M**2

        if Kinit<=flips and clr==0:
            U_sum_tot=0
            U_sum_sq_tot=0
            M_sum_tot=0
            M_sum_sq_tot=0
            sweep_cnt=0
            K=round(Kinit/2)
            flips=0
            clr=1

        if   K<=flips:
            nsample=sweep_cnt/nsweep
            mean_M=M_sum_tot/nsample
            mean_M_s=M_sum_sq_tot/nsample
            mean_U=U_sum_tot/nsample
            mean_U_s=U_sum_sq_tot/nsample
            cv = Calc_specific_heat(mean_U, mea_U_s,lattice.n)
            delta=abs(mean_M-mean_M_half)/abs(mean_M_half)
            U_sum_tot = 0
            U_sum_sq_tot = 0
            M_sum_tot = 0
            M_sum_sq_tot = 0
            sweep_cnt = 0
            K = 2*K
            flips = 0
            mean_M_half=mean_M


    return mean_M,mean_M_s,mean_U,mean_U_s,cv
#endregion
#hi yoav



#example of how to use the class
latt=lattice(n)
latt.Scatter()
latt2=lattice(n)
latt2.Scatter()
w=latt.lattice[1][1].s
y=latt2.lattice[1][1].s
print(w,y)
a=latt.lattice[1][1].flip()
b=latt2.flip(1,1)
print(a,b[1][1].s)



