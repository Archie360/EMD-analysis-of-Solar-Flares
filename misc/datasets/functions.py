import numpy as np
import math
import matplotlib.pyplot as plt
from statistics import mean


def exp(t,c1,m1,c2,m2,c3,m3,y0):
    return c1*np.exp(-t/m1) + c2*np.exp(-t/m2) + c3*np.exp(-t/m3)+y0

def scaling(t,c,m):
    return c*t**m

def ds(t,A,P,tau,phi,t0):
        return A*exp(-(t-t0)/tau)*np.sin(2*np.pi*(t-t0)/P -phi)
    
def bm(List):
    mins=[]
    maxs=[]
    for n, i in enumerate(List[1:-1]):
            if List[n] < List[n-1] and List[n] < List[n+1]:
                mins.append(List[n])
            elif List[n] > List[n-1] and i> List[n+1]:
                maxs.append(List[n])
    return len(mins)+len(maxs)

def log_fit(xi,yi,x,deg):  
    coefficients = np.polyfit(np.log10(xi), np.log10(yi), deg)
    polynomial = np.poly1d(coefficients)
    log10_y_fit = polynomial(np.log10(x))  # <-- Changed

    L=[]
    L.append(x)
    L.append(10**log10_y_fit)
    return L

def fit(xi,yi,x,deg,col='green'):  
    coefficients = np.polyfit((xi), (yi), deg)
    polynomial = np.poly1d(coefficients)
    log10_y_fit = polynomial((x))  # <-- Changed

    L=[]
    L.append(x)
    L.append(log10_y_fit)
    return L

def spread(x,k,N):
    return x+(k*(math.sqrt(2/N))*(np.e**(x/2))),x-(k*(math.sqrt(2/N))*(np.e**(x/2)))

def bypart(x,xi,yi,s,e,k) : 
    lf= log_fit(xi[s:e],yi[s:e],np.linspace(xi[s],xi[-1],len(x)),1)
    slope= (lf[1][1]-lf[1][0])/(lf[0][1]-lf[0][0])
    U,L=[],[]
    
    for i in range(len(lf[0])):
        u,l=spread(lf[1][i],k,len(lf[1])*10)
        L.append((l))
        U.append((u))
    return lf, slope, L, U

def scargle(p0, fxx, pxx):
    y0= -np.log(1-(1-p0)**(1/len(fxx)))*mean(pxx)
    return y0






if __name__ == '__main__':
    print('Hello')