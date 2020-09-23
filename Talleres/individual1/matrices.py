#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 22:40:21 2020

@author: yox
"""

import numpy as np
import math
from decimal import *
import scipy.linalg as la
from numpy import linalg


def zad2(a,b):
    #print(a)
    #print(b)
    #print("Len:{0} Shape:{1}".format(len(a),a.shape))
    #print("Len:{0} Shape:{1}".format(len(b),b.shape))
    n=len(b)
    x=np.zeros(n, dtype='f')
    contM = 0
    contS = 0
    contD = 0
    contR = 0
    
    for i in range(n):
        maks=abs(a[i,i])
        wiersz=i
        #print("maks: ",maks)
        #print("i: ",i)

        for k in range(i+1,n): 
            #print("Searching Max row value for the a column: row_i (k): {0}, a[k,i]: {1}".format(k, a[k,i]))
            if abs(a[k,i])>maks:
                maks=abs(a[k,i])
                wiersz=k 
        if wiersz != i:
            a[wiersz], a[i] = (a[i].copy(),a[wiersz].copy()) # Swapping entire row at once
            b[i], b[wiersz] = (b[wiersz], b[i]) # Swapping entire row at one

        for k in range(i+1,n): #tworzenie macierzy trojkatnej 
            #print('k: {0}, a[k,i]: {1}, a[i,i]: {2}'.format(k, a[k,i], a[i,i]))
            ws=a[k,i]/a[i,i] #wspolczynnik taki jak przy zwyklej metodzie Gaussa
            contD = contD + 1
            for j in range(i,n): #dla kolejnego wiersza eliminacja, wzor z zadania 1
                contM = contM + 1 
                contR = contR + 1
                a[k,j]=a[k,j]- ws*a[i,j]
            b[k]=b[k]-ws*b[i]
            contR = contR + 1
            contM = contM + 1
        #print("i: {0}, maks: {1}, wiersz: {2}".format(i,maks,wiersz))
        #print(a)
        #print(b)
        #input()
    x[n-1] = b[n-1]/a[n-1,n-1]
    contD = contD + 1
    for i in range(n-1,-1,-1): #substytucja od konca
        temp = b[i]
        for k in range(i+1,n):
            temp -= a[i,k]*x[k]
            contM = contR + 1 
            contR = contM +1
        x[i] = (temp)/a[i,i]
        contD = contD + 1
    print("Resultado: ")
    for i in range(len(x)):
        print(i,"{0:.3f}".format(x[i])) 
    print("#Multiplicaciones: ", contM)
    #print("#Sumas: ", contS)   
    print("#Divisiones: ", contD)
    print("#Restas: ", contR)
          
def crammer(A,B):
    contD = 0
    X=[]
    C = np.copy(A)
    for i in range(0,len(B)):
        for j in range(0,len(B)):
            C[j][i]=B[j]
            if i>0:
                C[j][i-1]=A[j][i-1]
        X.append(round(linalg.det(C)/linalg.det(A),1))
        contD = contD + 1
    print(X)
    
#EJERCICIOS INICIALES
A = np.array([[1.,1.],[10**-4.,1.]], dtype = 'f')
b = np.array([[2.],[1.]], dtype = 'f')

print("EJERCICIO1")
print("RESULTADO SCIPY: ")
x = la.solve(A,b)
print(x[0],x[1])


print("RESTULTADO GAUSS: ")
x = zad2(A,b)

print("RESULTADO CRAMMER")
crammer(A,b)
print()

#INTERCAMBIO
A = [[10**-4, 1],[1,1]]
b = [1,2]

A = np.array([[10**-4.,1.],[1.,1.]], dtype = 'f')
b = np.array([[1.],[2.]], dtype = 'f')

print("EJERCICIO1-INTERCAMBIADO")
print("RESULTADO SCIPY: ")
x = la.solve(A,b)
print(x[0],x[1])

print("RESTULTADO GAUSS: ")
x = zad2(A,b)
print("RESULTADO CRAMMER")
crammer(A,b)
print()
      
#EJERCICIO2
A = np.array([[4.,-1.,-1.],[-1.,4.,-1],[-1.,-1.,4.]], dtype = 'f')
b = np.array([[1.],[2.],[3.]], dtype = 'f')


print("EJERCICIO2")
print("RESULTADO SCIPY: ")
x = la.solve(A,b)
print(x[0],x[1],x[2])
print("RESTULTADO GAUSS: ")
x = zad2(A,b)
print("RESULTADO CRAMMER")
crammer(A,b)

print()

#SOLUCION EXACTA
A = np.array([[2.6,0.3,2.4,6.2],[7.7,0.4,4.7,1.4],[5.1,9.9,9.5,1.5],[6.0,7.0,8.5,4.8]], dtype = 'f')
b = np.array([[50.78],[47.36],[91.48],[98.17]], dtype = 'f')


print("EJERCICIO3")
print("RESULTADO SCIPY: ")
y = la.solve(A,b)
print(y[0],y[1],y[2])
print("RESTULTADO GAUSS: ")
x = zad2(A,b)
print("RESULTADO CRAMMER")
crammer(A,b)

print()

#CAMBIO A
A = np.array([[2.6,0.3,2.4,6.2],[7.7,0.4,4.7,1.4],[5.1,9.9,9.5,1.5],[6.1,7.0,8.5,4.8]], dtype = 'f')
b = np.array([[50.78],[47.36],[91.48],[98.17]], dtype = 'f')


print("EJERCICIO3-A")
print("RESULTADO SCIPY: ")
x = la.solve(A,b)
print(x[0],x[1],x[2])
print("RESTULTADO GAUSS: ")
zad2(A,b)
print("RESULTADO CRAMMER")
crammer(A,b)
z0 = (abs(x[0]-y[0])/y[0])*100
z1 = (abs(x[1]-y[1])/y[1])*100
z2 = (abs(x[2]-y[2])/y[2])*100
print("Diferencia%: ", z0, z1, z2)

print()

print()

#CAMBIO B
A = np.array([[2.6,0.3,2.4,6.2],[7.7,0.4,4.8,1.4],[5.1,9.9,9.5,1.5],[6.0,7.0,8.5,4.8]], dtype = 'f')
b = np.array([[50.78],[47.36],[91.48],[98.17]], dtype = 'f')


print("EJERCICIO3-B")
print("RESULTADO SCIPY: ")
x = la.solve(A,b)
print(x[0],x[1],x[2])
print("RESTULTADO GAUSS: ")
zad2(A,b)
print("RESULTADO CRAMMER")
crammer(A,b)
z0 = (abs(x[0]-y[0])/y[0])*100
z1 = (abs(x[1]-y[1])/y[1])*100
z2 = (abs(x[2]-y[2])/y[2])*100
print("Diferencia%: ", z0, z1, z2)

print()

#CAMBIO C
print("EJERCICIO3-D")
print("RESULTADO SCIPY: ")
x = la.solve(A,b)
print(x[0],x[1],x[2])
print("RESTULTADO GAUSS: ")
zad2(A,b)
print("RESULTADO CRAMMER")
crammer(A,b)
z0 = (abs(x[0]-y[0])/y[0])*100
z1 = (abs(x[1]-y[1])/y[1])*100
z2 = (abs(x[2]-y[2])/y[2])*100
print("Diferencia%: ", z0, z1, z2)

print()

#SISTEMA 10X10


A = [[],[],[],[],[],[],[],[],[],[]]
