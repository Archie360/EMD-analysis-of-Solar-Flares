import numpy as np
import pandas as pd


name = '21_0160363001/rgs12o12_b50_21_0160363001_decay_PF.dat'
L=name.split('/')[1].split('.')[0]+'variation_removed'
datContent = [i.strip().split() for i in 
              open(name).readlines()]

datContent = datContent[5:]

t= np.array(datContent)[:,0]
x1= np.array(datContent)[:,1]
e= np.array(datContent)[:,2]
t_,x_,e_=[],[],[]
for i in range(len(t)):
    t_.append(float(t[i]))
    x_.append(float(x1[i]))    
    e_.append(float(e[i]))    
    
t_ = np.array(t_)
x_ = np.array(x_)
e_ = np.array(e_)

def sine(A,w,t):
    return A*np.sin(w*t)

orbital_variation= sine(5,100,t_)

x=x_-orbital_variation

dict= {}
dict['time']=t_
dict['counts']=x
dict['error']=e



df=pd.DataFrame(dict)
df.to_csv(L+'.csv')
