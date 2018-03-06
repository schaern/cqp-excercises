# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 18:01:07 2018


"""

from numpy import *
from matplotlib.pyplot import *

m=1.
E=0.7
hbar=1.

amount_of_steps=105

def V(x,a):
    
    if ((0<x) and (x<=a)):
        return 0.
        #return 1.
        #return 4.*(x/a-x**2./a**2)
        #return 4*((x-a/2.)**2./(a/2.)-(x-a/2.)**4./(a/2.))
    else:
        return 0.
    


def k(E,x,a):
    v=V(x,a)
    return 2.*m*(E-v)/(hbar**2.)

def numerov_step(k0,k1,k2,psi0,psi1,dx):
   psi2=(2.*(1.-5.*dx**2*k1/12.)*psi1-(1.+dx**2*k0/12.)*psi0)/(1.+dx**2.*k2/12.)
   return psi2

#initial values
def init_psi(dx):
    #psi(a)
    psi1_init=1.
    #psi(a+dx)
    psi0_init=exp((0+1j)*dx*(2*m*E)/hbar**2.)
    return psi0_init, psi1_init



def numerov_multsteps(a,amount_of_steps):
    x=linspace(0,a,amount_of_steps)
    x=x[::-1]
    #set initial values at a and a+dx, need to delete a form list so that its not done twice
    x=x[1:]
    
    psi=[]
    dx=x[1]-x[0]
    #initial values
    psi0,psi1=init_psi(dx)
    k0=k(E,a+dx,a)
    k1=k(E,a,a)
    
    
    for i in range(len(x)):
        k2=k(E,x[i],a)
        psi2=numerov_step(k0,k1,k2,psi0,psi1,dx)
        psi.append(psi2)
        
        k0=k1
        k1=k2
        
        psi0=psi1
        psi1=psi2
    return psi

def prob(psi):
    psi_squared=[]
    for i in psi:
        #abs doesnt do anything but removing the warning about removing 0j i.e. casting complex nr to real etc...
        
        psi_squared.append(abs(conj(i)*(i)))
        
    #some kind of normalization...??
    maxx=max(psi_squared)
    for i in range(len(psi)):
        psi_squared[i]=psi_squared[i]/maxx
        
    return psi_squared
        
def make_plot(a,amount_of_steps):
    x=linspace(0,a,amount_of_steps)
    x=x[::-1]
    x=x[1:]
   
    
    psi=numerov_multsteps(5,amount_of_steps)
    psi_squared=prob(psi)
    #plot(x,psi_squared) 
    """
    def f(x):
        return 0.3*x
    

    """
    plot(x,real(psi))
    print(real(psi[0]))
    plot(x,imag(psi))
    
    VV=[]
    for n in x:
        VV.append(V(n,a))
    plot(x,VV)  
    plot(x,ones(len(x)))    
make_plot(1.,amount_of_steps)