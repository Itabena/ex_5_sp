import numpy as np
import matplotlib.pyplot as plt
import time
import json


#paths
result_path='C:\\Users\\ita-no\\My Drive\\BCs Physics\\3rd year\\Introduction to statistical physics\\Numerical ex\\ex4\\results\\'
data_path='C:\\Users\\ita-no\\My Drive\\BCs Physics\\3rd year\\Introduction to statistical physics\\Numerical ex\\ex4\\data\\'
#region - Assiting Functions
def Timer(t_start,t_end):
  if 1<((t_end-t_start)/60)<=60:
   print('total time ran:',"{:.2f}".format((t_end-t_start)/60),'min')
  elif ((t_end-t_start)/60)>60:
   print('total time ran:',"{:.2f}".format((t_end-t_start)/60/60),'hours')
  else:
    print('total time ran:',"{:.2f}".format(t_end-t_start),'sec')
  return
def Calc_energy(lattice,i,j,eta,h):
       latt=lattice.lattice
       p=latt[i][j].s
       p_right=latt[(i+1)%len(latt)][j].s
       p_left=latt[(i-1)%len(latt)][j].s
       p_down=latt[i][(j-1)%len(latt)].s
       p_up=latt[i][(j+1)%len(latt)].s
       epsilon = -p*h - eta *p* (p_up + p_down + p_left + p_right)
       return epsilon
def Choose_to_flip(lattice,i,j,eta,h):
    confirm=0
    epsilon=Calc_energy(lattice,i,j,eta,h)
    P=1/(1+np.exp(-2*epsilon)) #assuming Tkb=1
    tmp=np.random.binomial(1,P)
    if tmp==1:
        lattice.flip(i,j)
        confirm=1
    return lattice , confirm
def Calc_net_magenetization(lattice):
    M=0
    latt=lattice.lattice
    for i in range(len(latt)):
        for j in range(len(latt)):
            M += latt[i][j].s
    return M
def Calc_net_energy(lattice,eta,h):
    E=0
    latt=lattice.lattice
    for i in range(len(latt)):
        for j in range(len(latt)):
            E+=Calc_energy(lattice,i,j,eta,h)
    return E
def Calc_specific_heat(U_mean_tot,U_mean_sq_tot,n):
    kB=1
    cv=(1/n**2)*kB*(U_mean_sq_tot-(U_mean_tot**2))
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
    def Update(self,eta,h):
        """
        Updates the net energy and magnetization of the lattice
        :return: M total magnetization of the lattice , E total energy of the lattice
        """
        self.M=Calc_net_magenetization(self)
        self.E=Calc_net_energy(self,eta,h)
        return self.M,self.E
    def Scatter(self):
        """
        Creating the initial random scatter of the particles
        :return: lattice form of the scattered particles
        """
        for i in range(self.n):
            for j in range(self.n):
                self.lattice[i][j] = particle(np.random.randint(2)*2-1)
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

#region - Simulate
def Simulate_spin_resrviour(n,eta,h,nsweep,Kinit,delta):


    delta_M= delta + 1

    latt=lattice(n)
    latt.Scatter()
    latt.Update(eta,h)


    U_sum_tot=0
    U_sum_sq_tot=0
    M_sum_tot=0
    M_sum_sq_tot=0
    flips=0
    sweep_cnt=0
    K=Kinit
    mean_M_half=1
    mean_M=0
    mean_M_s=0
    mean_U=0
    mean_U_s=0

    clr=0

    tic=time.time()
    toc=time.time()
    ttt=0

    while delta_M>delta:
        # print('K is:',K)

        if K > 10**8:
            print('we got to',delta_M, '||', delta)
            print("This is too long for me, I quit!!!")
            break
        latt,flipstmp= Sweep(latt,eta,h)
        flips+=flipstmp
        sweep_cnt+=1
        if sweep_cnt%nsweep==0 or sweep_cnt==1:
            latt.Update(eta,h)
            U_sum_tot+=latt.E
            U_sum_sq_tot+=latt.E**2
            M_sum_tot+=latt.M
            M_sum_sq_tot+=latt.M**2

        if (Kinit/2)<=flips and clr==0:
            print('Clearing the system')
            U_sum_tot=0
            U_sum_sq_tot=0
            M_sum_tot=0
            M_sum_sq_tot=0
            sweep_cnt=0
            K=round(Kinit/2)
            flips=0
            clr=1

        if   K <= flips and clr==1:

##clock

            # Timer(tic,toc)

##checks for the code:
            # print(delta_M,'||',delta)
            # print(K , '||', flips)

            nsample=(sweep_cnt/nsweep)+1
            mean_M=M_sum_tot/nsample
            mean_M_s=M_sum_sq_tot/nsample
            mean_U=U_sum_tot/nsample
            mean_U_s=U_sum_sq_tot/nsample
            cv = Calc_specific_heat(mean_U, mean_U_s,latt.n)

            delta_M=np.array(abs(mean_M-mean_M_half))/np.array(abs(mean_M))


            U_sum_tot = 0
            U_sum_sq_tot = 0
            M_sum_tot = 0
            M_sum_sq_tot = 0
            sweep_cnt = 0
            K = 2*K
            toc=time.time()
            flips = 0
            mean_M_half=mean_M


    return mean_M,mean_M_s,mean_U,mean_U_s,cv
#endregion


#region- Main
#put this before the main
t_start=time.time()
#put the main here
eta2=np.arange(0.1,0.85,0.05)
eta1=np.arange(0.42,0.465,0.005)
eta3=np.arange(0.46,0.8,0.05)

M_list=[]
M_s_list=[]
U_list=[]
U_s_list=[]
cv_list=[]
# Open a file to write the results to
save_dict=dict()

# with open('results.txt', 'w') as f:
# Iterate over the values of eta1
print('Starting')
for h in [0]:
 cnt=0
 for i in range(0,len(eta1)):
        cnt+=1
        # Run the simulation and save the results
        tik=time.time()
        M,M_s,U,U_s,cv=Simulate_spin_resrviour(32, eta1[i], h, 5, 10000, 10**(-3))
        tok=time.time()
        Timer(tik,tok)
        print('Eta value:',"{:.2f}".format(eta1[i]),'M is: ',"{:.1f}".format(M),'simulation time was',)
        M_list.append(M)
        M_s_list.append(M_s)
        U_list.append(U)
        U_s_list.append(U_s)
        cv_list.append(cv)
        if cnt%5 == 0 or i==len(eta1)-1:
            save_dict['M']=M_list
            save_dict['Ms'] = M_s_list
            save_dict['U'] = U_list
            save_dict['Us'] = U_s_list
            save_dict['Cv'] = cv_list
            save_dict['eta'] = eta1[cnt-5:cnt].tolist()
            with open(data_path +'h= '+str(h)+'_field_eta_vs_data_until_HR_'+ "{:.2f}".format(eta1[i])+'.txt','wt') as file:
               json.dump(save_dict,file,indent=4)
            save_dict.clear()
            M_list=[]
            M_s_list=[]
            U_list=[]
            U_s_list=[]
            cv_list=[]
        # Write the results to the file
        # f.write(f"{M}, {M_s}, {U}, {U_s}, {cv}\n")

#
#     # Iterate over the values of eta2
#     for i in range(0,len(eta2)-1):
#         # Run the simulation and save the results
#         M,M_s,U,U_s,cv=Simulate_spin_resrviour(32, eta2[i], 0, 5, 10000, 10^(-3))
#         M_list.append(M)
#         M_s_list.append(M_s)
#         U_list.append(U)
#         U_s_list.append(U_s)
#         cv_list.append(cv)
#         # Write the results to the file
#         f.write(f"{M}, {M_s}, {U}, {U_s}, {cv}\n")
#
#     # Iterate over the values of eta3
#     for i in range(0,len(eta3)-1):
#         # Run the simulation and save the results
#         M,M_s,U,U_s,cv=Simulate_spin_resrviour(32, eta3[i], 0, 5, 10000, 10^(-3))
#         M_list.append(M)
#         M_s_list.append(M_s)
#         U_list.append(U)
#         U_s_list.append(U_s)
#         cv_list.append(cv)
#         # Write the results to the file
#         f.write(f"{M}, {M_s}, {U}, {U_s}, {cv}\n")
#
# print("M: " , M_list)
# print("M_s: " , M_s_list)
# print("U: " , U_list)
# print("U_s: " , U_s_list)
# print("cv: " , cv_list)
# #put this after the main
t_end=time.time()
Timer(t_start,t_end)
# print(M,M_s,U,U_s,cv)
#endregion


