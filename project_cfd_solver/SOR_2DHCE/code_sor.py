import numpy as np
import time as time
import matplotlib.pyplot as plt
lx=0.010
ly=0.015
nx=41
ny=41
dx=lx/(nx-1)
dy=ly/(ny-1)

x=np.zeros((nx,ny))
y=np.zeros((nx,ny))

a_s=1/(dy*dy)
a_w=1/(dx*dx)
a_e=1/(dx*dx)
a_n=1/(dy*dy)
a_p=1/(a_e + a_n)

#grid generation
for i in range(nx):
    for j in range(ny):
        x[i][j]=i*dx
        y[i][j]=j*dx

#initialization of T_old at the iterior nodes only.... boundary nodes will be taken care of at the boundary conditions step
T_old=np.zeros((nx,ny))

#iterations

start_time=time.time()
residuals=[]
interations=[]
final=10**(-6)
res=1
count=1
while res>=final:
    T_new=np.zeros((nx,ny))
    sum1=0
    for i in range(nx):
        x[i][ny-1]=100.00
        
    for i in range(nx):
        x[i][0]=0.00
        
    for j in range(1,)
        
    
