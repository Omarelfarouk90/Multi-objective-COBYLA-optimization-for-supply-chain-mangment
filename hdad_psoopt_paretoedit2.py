# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 15:46:28 2019

@author: omar.elfarouk
"""

import numpy as np
import math
from scipy.optimize import minimize
from scipy.optimize import NonlinearConstraint
import matplotlib.pyplot as plt
#Factors for the third party in the supply chain%%
U_Demand = 12000#normrnd(mn,std)%Monte carlo simulation mean 12000 unit%
alpha = 0.1#Percentage of rate of return of products from third party%
lamda = alpha * U_Demand
miu =0.1#Probability of havieng a returned product in a good working condition% 
gamma=0.7#Probability of having a rerurned part after disassembly with good working condition%
 
Q_TP=(lamda*(miu))+(lamda*(1-miu)*gamma)#Quantity of the third party%
#pd = makedist('Poisson','lamda',12000)
#U_Demand = random(pd)

std  = 200 # 10000 runslamda and %
#U_Demand = round(random('Normal',12000,109));%400000;
#JDO1=[265;30;30;65;105;40;190;100;120;40;85;50;202];% demand curve fitting
#for daikin australia
pd_1 = 600#round(random('Normal',600,24.49));q
pd_2 = 60#round(random('Normal',60,7.749));
Z_var = U_Demand - pd_1-pd_2 
#Objective equation start#
#l=l+1;%Integrate i value in the loop%
                                #wa(l)=w1;%Store weight 1%
                               # wb(l)=w2;%Store weight 2%
#Transportation Costs#
TS_s = 5000  #Transportation cost for the supplier(From alexandria road to downtown cairo)%
TS_m = 5000  #Transportation cost for the manufacturer(Assumed to be almost fixed due to practicallity)%
TS_d = 5000  #Transportation cost for the distributer%
TS_rt = 5000 #Transportation cost for the retailer%
TS_tp = 5000 #Transportation cost for the third party%
#collection Costs%%
C_tp = 5.1 #collection cost of recovering used product from the customer%
#facility opening Costs%%
F_rt = 10000000 #facility opening cost for the recovery center(Assumed to be 10 million  Egyptian pound)%
#Ordering Costs%%
S_s = 5.1
S_ms = 58.956
S_m1 = 700
S_m2 = 800
S_m3 = 850
S_d = 173.4
S_r = 204
S_tp = 42.5
#Holding Costs%%
H_s = 50.126
H_ms = 589.56
H_m = 1473.9
H_dr = 1734
H_rt = 2040
H_tp = 425.9571
#Production Rates%%
P_m1=200 #Production Rates assumed to be 200 unit per day for the power plant %%
P_m2=210
P_m3=220
#U_Demand = 400000 #Demand rate is asumed to be 400,000 unit per month%
P_m = P_m1+P_m2+P_m3 # Production rate of the manufacuter
#i_m #conunting of manufacturer%
#i_mp    
#i_d   #Counting of Distributer
 
##Factors for the third party in the supply chain##
alpha = 0.1 #Percentage of rate of return of products from third party%
lamda =(alpha*U_Demand) 
miu =0.1  #Probability of havieng a returned product in a good working condition% 
gamma=0.7 #Probability of having a rerurned part after disassembly with good working condition%
 
Q_TP =(lamda*(miu))+(lamda*(1-miu)*gamma)     #Quantity of the third party%
#Values of supplied chain quantities
n_s = 1                                           
n_m = 1                                           #1:2
n_d = 1 
#input of optimized models data 
(lambda x : x[0])
def Q_rt1(x):
    return x[0] # quantity of the retailer in the forward cycle
def Q_rt2(x):
    return x[1] # quantity of the retailer in the forward cycle
def Q_rt3(x):
    return x[2] # quantity of the retailer in the forward cycle
def Q_d1(x):
    return x[3] # Quantity of the distributer
def Q_d2(x):
    return x[4] # Quantity of the distributer
def Q_d3(x):
    return x[5] # Quantity of the distributer
def Q_m1(x):
    return x[6] # Quantity of the Manufacturer
def Q_m2(x):
    return x[7] # Quantity of the Manufacturer
def Q_m3(x):
    return x[8] # Quantity of the Manufacturer
def Q_s1(x):
    return x[9] #Quantity of Supplied Parts
def Q_s2(x):
    return x[10] #Quantity of Supplied Parts
def Q_s3(x):
    return x[11] #Quantity of Supplied Parts                                          #1:2
#Cycle time of the supply chain#
def t_r(x):
    return (U_Demand)/(x[0])   #cycle time of the retailer
def t_d(x):
    return n_d * t_r(x)            #cycle time of the Distribiter
def t_m(x):
    return (n_m * n_d* t_r(x))    #cycle time of the Manufacturer
def t_s(x):
    return  n_s *n_m *n_d *t_r(x)  #cycle time of the supplier
def t_tp(x):
    return t_s(x)                  #cycle time of the third party
#
S_jfS=30   #Job Index factor number of fixed jobs at the supplier assumed to be 30 fixed employees %
S_jfM=30   #Job index for the number of fixed jobs by Mamufacturer assumed to be 30 fixed employees %
S_jfD=30   #Job index for the number of fixed jobs by distributer assumed to be 30 fixed employees%
S_jfRT=30 #Job index for the number of fixed jobs by retialer assumed to be 30 fixed employees%
S_jfTP=20 #Job index for the number of fixed jobs by third party recovery assumed to be 20 fixed employees%
S_jvS=270 #Job Index factor number of variable jobs at the supplier assumed to be 270 workers per facility%
S_jvM=270 #Job index for the number of variable jobs by Mamufacturer  270 workers per facility%
S_jvD=270 #Job index for the number of variable jobs by distributer  270 workers per facility%
S_jvRT=270#Job index for the number of variable jobs by retialer  270 workers per facility%
S_jvTP=100#Job index for the number of variable jobs by third party recovery  100 workers per facility%
S_u=20    #Employee satisfaction factor of the refurbrished parts for the third party disassembler%
S_rt=30   #Customer satisfaction factor of the refurbrished parts%
#Number of lost days at work%
S_ds=5  # Number of lost days from injuries or work damage at the suppliers / month%
S_dm=5  #Number of lost days from injuries or work damage at the manufactuer%
S_dd=5  #Number of lost days from injuries or work damage at the distributer%
S_drt=5 #Number of lost days from injuries or work damage at the retailer%
S_dtp=5 #Number of lost days from injuries or work damage at the third party%
#Enviromental Aspect of the supply chain (Emissions calculated from carbon footprint)%
E_q=10   #Emission factor from production line
E_tp=10  #Emission from wastes removal%
#Transportation emission cost%
E_ts=20   #Emission from Transportation made by the supplier%
E_tm=20   #Emission from Transportation made by the manufacturer%
E_td=20   #Emission from Transportation made by the distributer%
E_trt=20  #Emission from Transportation made by the retailer%
E_ttp=20  #Emission from Transportation made by the third party%
#for w1 in w11:
   #starting of the loop#   
#i_s = fun()
#if type(i_s) == int :# if the return value is an integer
    #do this
#    elif type(i_s) == str:# if the retrun value is string#
        #do this
i_s = 1 
i_ss=np.arange(i_s,n_s+1,1)
tc_s1= list(range(i_s,n_s+1))  
for i_s in i_ss:
    def tc_s1(x):
        tc_s1 = np.sum(((i_ss)/n_s)*(Q_s1(x))*(t_s(x)))
        return tc_s1
    i_s=i_s + 1  # Adding value of Supplier integer#
def tc_s4(x):
        tc_s4 = (tc_s1(x))
        return tc_s4 
def TC_s1(x):
        TC_s1 =  (S_s*(1/(n_s*t_s(x))))+(((H_s+TS_s)/(n_s*(t_s(x))))*tc_s4(x)) #cost of the supplier for component 1%
        return TC_s1
i_s= 1    #starting of the loop#         
i_ss=np.arange(i_s,n_s+1,1)
#for w1 in w11:
tc_s2= list(range(i_s,n_s+1))  
for i_s in i_ss:
    def tc_s2(x):
        tc_s2=np.sum((i_ss/n_s)*Q_s2(x)*t_s(x)) #((x(11)) +Q_TP#
        return tc_s2
    i_s = i_s + 1   #Adding value of Supplier integer
def tc_s5(x):    
    tc_s5 = (tc_s2(x))
    return tc_s5
def TC_s2(x):
    TC_s2 = (S_s*(1/(n_s*t_s(x))))+(((H_s+TS_s)/(n_s*(t_s(x))))*tc_s5(x))
    return TC_s2
i_s=1    #starting of the loop#         
tc_s3= list(range(i_s,n_s+1))  
for i_s in i_ss:
    def tc_s3(x):
        tc_s3=np.sum((i_ss/n_s)*Q_s3(x)*t_s(x))  #((x(12)+ Q_TP))%
        return  tc_s3
    i_s = i_s + 1   #Adding value of Supplier integer (No addition for Q_TP )%
def tc_s6(x):
    tc_s6 = tc_s3(x)
    return tc_s6
def TC_s3(x):
    TC_s3 = (S_s*(1/(n_s*t_s(x))))+(((H_s+TS_s)/(n_s*(t_s(x))))*tc_s6(x))  #cost of the supplier for component 3%
    return TC_s3
i_m = 1    #starting of the loop#   
i_mm=np.arange(i_m,n_m+1,1)
#for w1 in w11:
tc_m2= list(range(i_m,n_m+1))  
for i_m in i_mm:
    tc_m1 =np.arange(1,n_m,1) #Defining range with starting and ending point
    def tc_m2(x):
        tc_m2 = np.sum((1-((i_mm)/(n_m)))*((Q_m1(x))+Q_TP)) #Defining range with start & ending point#
        return tc_m2
    i_m=i_m + 1  # Adding value of manufacturer integer#
def tc_m3(x):
    tc_m3=(tc_m2(x))        
    return tc_m3
tc_s7 =np.arange(1,n_s,1) 
#Total cost of manufacturer#
tc_m = sum(tc_m1)
tc_s8 = sum(tc_s7)
def TC_m(x):
     TC_m =(H_m*((0.5*(Q_m1(x)**2)*(1/P_m1))\
              +(tc_m*(Q_m1(x)*t_m(x)*(1/(n_m**2))))))\
              +((S_m1+TS_m)*(1/t_m(x)))+((S_ms+TS_tp)*(1/t_s(x)))\
              +(H_ms*(1/t_s(x))*(((((Q_s1(x)+Q_TP)*Q_m1(x))/P_m1))\
                       +(tc_s8*(((Q_s1(x))+Q_TP)/n_s)*(t_m(x)-(Q_m1(x)/P_m1)))))              
     return TC_m 
def TC_m2(x):
    TC_m2 =(H_m*((0.5*(Q_m2(x)**2)*(1/P_m2))\
              +(tc_m*(Q_m2(x)*t_m(x)*(1/(n_m**2))))))\
              +((S_m2+TS_m)*(1/t_m(x)))+((S_ms+TS_tp)*(1/t_s(x)))\
              +(H_ms*(1/t_s(x))*(((((Q_s2(x)+Q_TP)*Q_m2(x))/P_m2))\
                       +(tc_s8*(((Q_s2(x))+Q_TP)/n_s)*(t_m(x)-(Q_m2(x)/P_m2)))))
    return TC_m2
def TC_m3(x):
    TC_m3 =(H_m*((0.5*(Q_m3(x)**2)*(1/P_m3))\
              +(tc_m*(Q_m3(x)*t_m(x)*(1/(n_m**2))))))\
              +((S_m3+TS_m)*(1/t_m(x)))+((S_ms+TS_tp)*(1/t_s(x)))\
              +(H_ms*(1/t_s(x))*(((((Q_s3(x)+Q_TP)*Q_m3(x))/P_m3))\
                       +(tc_s8*(((Q_s3(x))+Q_TP)/n_s)*(t_m(x)-(Q_m3(x)/P_m2)))))
                   #Cost of the manufacturer for product 3
    return TC_m3
i_d=1    #starting of the loop#         
i_dd=np.arange(i_d,n_d+1,1)
#for w1 in w11:
tc_d1= list(range(i_d,n_d+1))
tc_d2= list(range(i_d,n_d+1))
tc_d3= list(range(i_d,n_d+1))  
for i_d in i_dd:
    def tc_d1(x):
        tc_d1=np.sum(((i_dd)/(n_d))*(Q_d1(x)))    #Cost of the Distributer for Product 1%%
        return tc_d1
    def tc_d2(x):
        tc_d2=np.sum(((i_dd)/(n_d))*(Q_d2(x)))   #Cost of the Distributer for Product 2%%
        return tc_d2
    def tc_d3(x):
        tc_d3=np.sum(((i_d)/(n_d))*(Q_d3(x)))  #Cost of the Distributer for Product 3%%
        return tc_d3
    i_d = i_d + 1
def tc_d_f(x):
    tc_d_f = (tc_d1(x))+(tc_d2(x))+(tc_d3(x))
    return tc_d_f
def TC_d(x):
    TC_d = (H_dr*(tc_d_f(x)/n_d))+((S_d+TS_d)*(1/t_d(x)))  #Total cost of the distributer of the supply chain%
    return TC_d 
#Total cost of retailer
def TC_rt(x):     
    TC_rt = (H_rt*((Q_rt1(x))/2))+ ((S_r+TS_rt)*(1/t_r(x)))  #Cost of the retailer%%
    return TC_rt
def TC_rt2(x):
    TC_rt2=(H_rt*((Q_rt2(x))/2))+ ((S_r+TS_rt)*(1/t_r(x))) #Cost of the retailer for product 2%%
    return TC_rt2
def TC_rt3(x):
    TC_rt3 = (H_rt*((Q_rt3(x))/2))+ ((S_r+TS_rt)*(1/t_r(x))) #Cost of the retailer for product 3%%
    return TC_rt3
#Total cost of third party recovery 
def TC_tp(x):
    TC_tp = ((H_tp/2)*Q_TP)+((S_tp+TS_tp)*(1/t_tp(x)))
    return TC_tp
S_jfS=30   #Job Index factor number of fixed jobs at the supplier assumed to be 30 fixed employees %
S_jfM=30   #Job index for the number of fixed jobs by Mamufacturer assumed to be 30 fixed employees %
S_jfD=30   #Job index for the number of fixed jobs by distributer assumed to be 30 fixed employees%
S_jfRT=30 #Job index for the number of fixed jobs by retialer assumed to be 30 fixed employees%
S_jfTP=20 #Job index for the number of fixed jobs by third party recovery assumed to be 20 fixed employees%
S_jvS=270 #Job Index factor number of variable jobs at the supplier assumed to be 270 workers per facility%
S_jvM=270 #Job index for the number of variable jobs by Mamufacturer  270 workers per facility%
S_jvD=270 #Job index for the number of variable jobs by distributer  270 workers per facility%
S_jvRT=270#Job index for the number of variable jobs by retialer  270 workers per facility%
S_jvTP=100#Job index for the number of variable jobs by third party recovery  100 workers per facility%
S_u=20    #Employee satisfaction factor of the refurbrished parts for the third party disassembler%
S_rt=30   #Customer satisfaction factor of the refurbrished parts%
#Number of lost days at work%
S_ds=5  # Number of lost days from injuries or work damage at the suppliers / month%
S_dm=5  #Number of lost days from injuries or work damage at the manufactuer%
S_dd=5  #Number of lost days from injuries or work damage at the distributer%
S_drt=5 #Number of lost days from injuries or work damage at the retailer%
S_dtp=5 #Number of lost days from injuries or work damage at the third party%
#Enviromental Aspect of the supply chain (Emissions calculated from carbon footprint)%
E_q=10   #Emission factor from production line
E_tp=10  #Emission from wastes removal%
#Transportation emission cost%
E_ts=20   #Emission from Transportation made by the supplier%
E_tm=20   #Emission from Transportation made by the manufacturer%
E_td=20   #Emission from Transportation made by the distributer%
E_trt=20  #Emission from Transportation made by the retailer%
E_ttp=20  #Emission from Transportation made by the third party%
#Cycle time%
t_r  #cycle time of the retailer
t_d  #cycle time of the Distribiter
t_m #cycle time of the Manufacturer
t_s #cycle time of the supplier
t_tp #cycle time of the third party

#def objective(x):
   # return  (w1*EQO)-(w2*LSC)+(w3*ESC)
def EQO(x):
    EQO = TC_s1(x)+TC_s2(x)+TC_s3(x)+TC_m(x)+TC_m2(x)+TC_m3(x)+TC_d(x)+TC_rt(x)\
            +TC_rt2(x)+TC_rt3(x)+TC_tp(x)
    return EQO
    #Economical aspect#
def LSC(x):
    LSC =(S_jfS+S_jfM+S_jfD+S_jfRT+S_jfTP)\
     +((S_jvS*Q_s1(x))+(S_jvD*Q_d1(x))+(S_jvM*Q_m1(x))\
       +(S_jvRT*Q_rt1(x))+(S_jvTP*Q_TP))\
     +(S_u*(U_Demand))+(S_rt*Q_rt1(x))-(S_ds*Q_s1(x))\
       +(S_dd*Q_d1(x))+(S_dm*Q_m1(x))+(S_drt*Q_rt1(x))\
       +(S_dtp*Q_TP)#Social aspect equation%
    return LSC
def ESC(x):
    ESC=(E_q*(Q_s1(x)+Q_d1(x)+Q_m1(x)+Q_rt1(x)))\
    +(E_ts*(1/t_s(x)))+(E_td*(1/t_d(x)))\
    +(E_tm*(1/t_m(x)))+(E_trt*(1/t_r(x)))\
    +(E_ts*(1/t_tp(x)))+(E_tp*Q_TP)  #Enviromental aspect
    return ESC
#def EQO(X):
 #   return EQO==TC_s1(x)+TC_s2(x)+TC_s3(x)+TC_m(x)+\
#TC_m2(x)+TC_m3(x)+TC_d(x)+TC_rt(x)+TC_rt2(x)+TC_rt3(x)+TC_tp(x)
    #Economical aspect#
#def LSC(x):
   # return LSC==((S_jfS+S_jfM+S_jfD+S_jfRT+S_jfTP)\
   #  +((S_jvS*Q_s1)+(S_jvD*Q_d1)+(S_jvM*Q_m1)\
    #   +(S_jvRT*Q_rt1)+(S_jvTP*Q_TP)))\
    # +(S_u*(U_Demand))+(S_rt*Q_rt1)-((S_ds*Q_s1)\
      # +(S_dd*Q_d1)+(S_dm*Q_m1)+(S_drt*Q_rt1)\
     # +(S_dtp*Q_TP))#Social aspect equation%
#def ESC(x):                                   #Enivromental aspect#
  #  return ESC==(E_q*(Q_s1+Q_d1+Q_m1+Q_rt1))\
  #  +(E_ts*(1/t_s))+(E_td*(1/t_d))\
  #  +(E_tm*(1/t_m))+(E_trt*(1/t_r))\
   # +(E_ts*(1/t_tp))+(E_tp*Q_TP)
f=1
j=1 
w11=np.arange(0.1,1.1,0.1)
solx=np.arange(1,11,1)
soly=np.arange(1,11,1)
solz=np.arange(1,11,1)
a0 = np.zeros(111)
a1 = np.zeros(111)
obj1 = np.zeros(111)
obj2 = np.zeros(111)
obj3 = np.zeros(111)
for w1 in w11:
#l=0;%Start indexing l value%
  w22=np.arange(0.1,1.1-w1,0.1)	
  for w2 in w22:
       w3=1.1-w1-w2
       def objective(x):
    #return  -(w2*LSC)
           return  (EQO(x)*w1)+(LSC(x)*w2)+(ESC(x)*w3)      
       U_Demand = 120000#normrnd(mn,std)%Monte carlo simulation mean 12000 unit%
       #dem = np.array([U_Demand])       #Storing demand #                 
       #Economic data%
       #solx[j]= EQO #Economical objective% 
       #Enviromental data%
       #soly[j]= ESC#Enviromental object
       #social data%
       #solz[j]= LSC#Economical objective%
            
       #j=j+1#Inceremental iteration
       #f=f+1 
       #Define constraint#
       def constraint1(x):
           return x[6]-U_Demand 
       def constraint2(x):
           return x[7]-U_Demand 
       def constraint3(x):
           return x[8]-U_Demand 
       #def constraint4(x):
           #return  U_Demand - (U_Demand-(3*std)) #Monte carlos Inequality#
       #def constraint5(x):
           #return (U_Demand+(3*std))-U_Demand
       def constraint6(x):
           return objective(x)-0
                    #C(1) = (n_s*(x[1]))-pd_1- x[0];
                    #C(2) = (n_s*(x[1]))-pd_1- x[0];
                    #C(3) = (n_mp*(x[2])) - x[1];
                    #C(4) = (n_mp*(x[2])) - x[1];
                    #C(5) = (n_d*x[3])- x[2];
                    #C(6) = (n_d*x[3])- x[2];
                    #C(7) = U_Demand - x[3];
                    #C(8) = U_Demand-x[3];
       def constraint7(x):
            return x[6]-(n_m*(x[3]))
       def constraint8(x):
            return x[7]-(n_m*(x[4]))
       def constraint9(x):
           return  x[8]-(n_m*(x[5]))
       def constraint13(x):
           return x[3]-(n_d*x[0])
       def constraint14(x):
           return x[4]- (n_d*x[1])
       def constraint15(x):
           return x[5]- (n_d*x[2])
       def constraint16(x):
            return ((x[9])+Q_TP)-(n_s*x[6])
       def constraint17(x):
            return ((x[10])+Q_TP)-(n_s*x[7])
       def constraint18(x):
            return ((x[11])+Q_TP)-(n_s*x[8])
       #def constraint10(x):
            #return x[9]-(n_s*x[6])
       #def constraint11(x):
            #return x[10]-(n_s*x[7])
       #def constraint12(x):
            #return x[11]-(n_s*x[8])


#Ceq(4) = (n_mp*(x[2])) - x[1];

# initial guesses
       n = 12
       x0 = np.zeros(n)
       x0[0] = 20000
       x0[1] = 20000
       x0[2] = 20000
       x0[3] = 20000
       x0[4] = 20000
       x0[5] = 20000
       x0[6] = 20000
       x0[7] = 20000
       x0[8] = 20000
       x0[9] = 20000
       x0[10] = 20000
       x0[11] = 20000
       # show initial objective
       print('Initial Objective: ' + str(objective(x0)))
       # optimize
       b = (1,100000)
       bnds = (b, b, b, b, b, b, b, b, b, b, b, b)
       #initial values
       initial_point=[20000.0,20000.0,20000.0,20000.0,
                      20000.0,20000.0,20000.0,20000.0,
                      20000.0,20000.0,20000.0,20000.0]    
       #lower and upper bound for variables
       bounds=[ [1,100000],[1,100000],[1,100000],[1,100000],
               [1,100000],[1,100000],[1,100000],[1,100000],
              [1,100000],[1,100000],[1,100000],[1,100000], ]

       #construct the bounds in the form of constraints
       cons = []
       for factor in range(len(bounds)):
           lower, upper = bounds[factor]
           l = {'type': 'ineq',
           'fun': lambda x, lb=lower, i=factor: x[i] - lb}
           u = {'type': 'ineq',
             'fun': lambda x, ub=upper, i=factor: ub - x[i]}
           cons.append(l)
           cons.append(u)
       con1 = {'type': 'ineq', 'fun': constraint1}
       cons.append(con1)
       con2 = {'type': 'ineq', 'fun': constraint2}
       cons.append(con2)
       con3 = {'type': 'ineq', 'fun': constraint3}
       cons.append(con3)
       #con4 = {'type': 'ineq', 'fun': constraint4}
       #con5 = {'type': 'ineq', 'fun': constraint5}
       con6 = {'type': 'ineq', 'fun': constraint6}
       cons.append(con6)
       con7 = {'type': 'ineq', 'fun': constraint7}
       cons.append(con7)
       con8 = {'type': 'ineq', 'fun': constraint8}
       cons.append(con8)
       con9 = {'type': 'ineq', 'fun': constraint9}
       cons.append(con9)
       #con10 = {'type': 'eq', 'fun': constraint10}
       #con11 = {'type': 'eq', 'fun': constraint11}
       #con12 = {'type': 'eq', 'fun': constraint12}
       con13 = {'type': 'ineq', 'fun': constraint13}
       cons.append(con13)
       con14 = {'type': 'ineq', 'fun': constraint14}
       cons.append(con14)
       con15 = {'type': 'ineq', 'fun': constraint15}
       cons.append(con15)
       con16 = {'type': 'ineq', 'fun': constraint16}
       cons.append(con16)
       con17 = {'type': 'ineq', 'fun': constraint17}
       cons.append(con17)
       con18 = {'type': 'ineq', 'fun': constraint18}
       cons.append(con18)
       
       #construct the bounds in the form of constraints


#similarly aditional constrains can be added
       solution = minimize(objective,x0,method='COBYLA',\
                    bounds=bnds,constraints=cons)
       #trust-constr #SLSQP options={'ftol': 1e-9, 'disp': True}
       x = solution.x
       # show final objective
       #print('Final Objective: ' + str(objective(x)))

# show final objective
       print('Final Objective: ' + str(objective(x)))
# print solution
       print('Solution')
       print('Qrt1 = ' + str(x[0]))
       print('Qrt2 = ' + str(x[1]))
       print('Qrt3 = ' + str(x[2]))
       print('Qd1 = ' + str(x[3]))
       print('Qd2 = ' + str(x[4]))
       print('Qd3 = ' + str(x[5]))
       print('Qm1 = ' + str(x[6]))
       print('Qm2 = ' + str(x[7]))
       print('Qm3 = ' + str(x[8]))
       print('Qs1 = ' + str(x[9]))
       print('Qs2 = ' + str(x[10]))
       print('Qs3 = ' + str(x[11]))
       print('w1 = ' +str(w1))
       print('w2 = ' +str(w2))
       print('w3 = ' +str(w3))
       a0[j]= str(objective(x))
       obj1[j]= str(EQO(x)*w1)
       obj2[j]=str(LSC(x)*w2) 
       obj3[j]=str(ESC(x)*w3  )
       a1[j]= str(x[0])
       j=j+1
plt.plot(a0,a1,label='Retailer 1 vs objective')
plt.xlabel('Retailer')
plt.ylabel('Objective fun')
plt.title('Retailer 1 vs objective function')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')




ax.scatter(obj1, obj2, obj3, c='r', marker='o')

ax.set_xlabel('Economic objective')
ax.set_ylabel('social objective')
ax.set_zlabel('enviromental objective')

plt.show()