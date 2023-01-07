import numpy as np
import matplotlib.pyplot as plt
import time
import json
from Projects.Bcs_physics.Numeric_statistical_mechanics.ex_5.Simulation import Calc_specific_heat

#region - Paths
result_path='C:\\Users\\ita-no\\My Drive\\BCs Physics\\3rd year\\Introduction to statistical physics\\Numerical ex\\ex4\\results\\'
data_path='C:\\Users\\ita-no\\My Drive\\BCs Physics\\3rd year\\Introduction to statistical physics\\Numerical ex\\ex4\\data\\'
#endregion
# region -Importiong simulated data
h=[0,0.1,0.5,1]
#region - h=0.1
M_h_01=[]
M_s_h_01=[]
U_h_01=[]
U_s_h_01=[]
cv_h_01=[]
eta_h_01=[]
for i in [0.30,0.55,0.75]:
  with open(data_path+'h= '+str(h[1])+'_field_eta_vs_data_until_'+ "{:.2f}".format(i)+'.txt','rt') as file:
     tmp_dict=json.load(file)
  M_h_01.extend(tmp_dict["M"])
  M_s_h_01.extend(tmp_dict["Ms"])
  U_h_01.extend(tmp_dict["U"])
  U_s_h_01.extend(tmp_dict["Us"])
  cv_h_01.extend(tmp_dict["Cv"])
  eta_h_01.extend(tmp_dict["eta"])
#endregion
tmp_dict={}
#region - h=0.5
M_h_05=[]
M_s_h_05=[]
U_h_05=[]
U_s_h_05=[]
cv_h_05=[]
eta_h_05=[]
for i in [0.30,0.55,0.75]:
  with open(data_path+'h= '+str(h[2])+'_field_eta_vs_data_until_'+ "{:.2f}".format(i)+'.txt','rt') as file:
     tmp_dict=json.load(file)
  M_h_05.extend(tmp_dict["M"])
  M_s_h_05.extend(tmp_dict["Ms"])
  U_h_05.extend(tmp_dict["U"])
  U_s_h_05.extend(tmp_dict["Us"])
  cv_h_05.extend(tmp_dict["Cv"])
  eta_h_05.extend(tmp_dict["eta"])
#endregion
tmp_dict={}
#region - h=1.0
M_h_10=[]
M_s_h_10=[]
U_h_10=[]
U_s_h_10=[]
cv_h_10=[]
eta_h_10=[]
for i in [0.30,0.55,0.75]:
  with open(data_path+'h= '+str(h[3])+'_field_eta_vs_data_until_'+ "{:.2f}".format(i)+'.txt','rt') as file:
     tmp_dict=json.load(file)
  M_h_10.extend(tmp_dict["M"])
  M_s_h_10.extend(tmp_dict["Ms"])
  U_h_10.extend(tmp_dict["U"])
  U_s_h_10.extend(tmp_dict["Us"])
  cv_h_10.extend(tmp_dict["Cv"])
  eta_h_10.extend(tmp_dict["eta"])
#endregion
tmp_dict={}
#region - h=0
M_h_0=[]
M_s_h_0=[]
U_h_0=[]
U_s_h_0=[]
cv_h_0=[]
eta_h_0=[]
for i in [0.30,0.55,0.80]:
  with open(data_path+'h= '+str(h[0])+'_field_eta_vs_data_until_'+ "{:.2f}".format(i)+'.txt','rt') as file:
     tmp_dict=json.load(file)
  M_h_0.extend(tmp_dict["M"])
  M_s_h_0.extend(tmp_dict["Ms"])
  U_h_0.extend(tmp_dict["U"])
  U_s_h_0.extend(tmp_dict["Us"])
  cv_h_0.extend(tmp_dict["Cv"])
  eta_h_0.extend(tmp_dict["eta"])
#endregion
tmp_dict={}
#region - h=0,Around critiacl
M_h_0_HR=[]
M_s_h_0_HR=[]
U_h_0_HR=[]
U_s_h_0_HR=[]
cv_h_0_HR=[]
eta_h_0_HR=[]
for i in [0.44,0.47]:
  with open(data_path+'h= '+str(h[0])+'_field_eta_vs_data_until_HR_'+ "{:.2f}".format(i)+'.txt','rt') as file:
     tmp_dict=json.load(file)
  M_h_0_HR.extend(tmp_dict["M"])
  M_s_h_0_HR.extend(tmp_dict["Ms"])
  U_h_0_HR.extend(tmp_dict["U"])
  U_s_h_0_HR.extend(tmp_dict["Us"])
  cv_h_0_HR.extend(tmp_dict["Cv"])
  eta_h_0_HR.extend(tmp_dict["eta"])
#endregion
#endregion
#region - Other important variables
n=32
N=n**2
#endregion
#region -Plotting functions
def Plot_Data(x,y,save_dir=' ',save_fig=0,new_fig=0,log=0):
    '''
    :param x: The data of X axis.
    :param y: The data of Y axis.
    :param save_dir: ' ' by defult, The path for the saving folder.
    :param save_fig: 0 by defult ,Type 1 if save intended.
    :param new_fig:  0 by defult ,Type 1 if you wish to create a nex figure.
    :param log: 0 by defult , Type 1 if yo wish to create semi-log ,Type -1 to create log
    :return: Plots the wanted figures
    '''

    if new_fig==0:
       plt.figure()
    print('The lagand of the ploted data: ')
    plot_name=input()
    plt.plot(x,y,'*',label=plot_name)
    if log==1:
        plt.yscale('log')
        plt.xscale('log')
        plt.grid()
    elif log==-1:
        plt.yscale('log')
        plt.grid()
    else:
        plt.grid()
    print('The x Axis label is: ')
    x_name=input()
    print('The y Axis label is: ')
    y_name=input()
    plt.xlabel(x_name,fontsize=15)
    plt.ylabel(y_name, fontsize=15)
    if save_fig==0:
       print('The plot file name :')
       name = input()
       plt.savefig(save_dir+name+'.png')
    return
def Onsager(eta_min,eta_max,num_res):
   '''
    :param eta_min: Minimum of eta list
    :param eta_max: Maximum of eta list
    :param num_res: number of elements in eta list
    :return m:|M/N| claculated by Onsager analytical solution
   '''
   eta=np.linspace(eta_min,eta_max,num_res,endpoint=True)
   z=np.exp(-2*eta)
   m=((1+z**2)**(1/4))*((1-6*(z**2)+z**4)**(1/8))/((1-z**2)**(1/2))
   return m

#endregion
 #region - Creating the final results
# x=np.arange(1,10,0.1)
# y=x**2
eta_h_0.extend(eta_h_0_HR)
M_h_0.extend(M_h_0_HR)
Plot_Data(eta_h_0,abs(np.divide(M_h_0,N)),save_dir='C:\\Users\\user\\Desktop\\',save_fig=1,new_fig=0,log=0)
plt.show()
#endregion