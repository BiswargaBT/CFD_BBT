import numpy as np
import time as time
#import matplotlib.pyplt as plt

lx=1.0
ly=1.5
nx=41
ny=41

dx=lx/(nx-1)
dy=ly/(ny-1)
dt=0.5
alpha=1.8

a_s=alpha*dt/(dy*dy*2)
a_w=alpha*dt/(dx*dx*2)
a_e=alpha*dt/(dx*dx*2)
a_n=alpha*dt/(dy*dy*2)
a_p=-2*alpha*dt*(1/(2*dx*dx) + 1/(2*dy*dy))

x=np.zeros((nx,ny))
y=np.zeros((nx,ny))

#grid generation

for i in range(nx):
    for j in range(ny):
        x[i][j]=dx*i
        y[i][j]=dy*j
        
#initializations
T_old=np.zeros((nx,ny))
T_oldtime=np.zeros((nx,ny))
T_newtime=np.zeros((nx,ny))
T_new=np.zeros((nx,ny))
b_old=np.zeros((nx,ny))
count=0
sum1=0
sum2=0
res1=1.0
res2=1.0
final_s=10**(-6)#final solver residual
final_t=10**(-1)#final time marching residual
iterations=[]
residuals=[]
tl=1 # at t=0, time level tl is 1
start_time=time.time()
#Initial guess for the solution at all the interior node
for i in range(1,nx-1):
    for j in range(1,ny-1):
        T_oldtime[i][j]=0

while res1>=final_t:     
    #Top surface boundary condition
    for i in range(nx):
        T_old[i][ny-1]=100.00
        T_oldtime[i][ny-1]=100.00
    #bottom boundary condition
    for i in range(nx):
        T_old[i][0]=0.00
        T_oldtime[i][0]=0.00
    #left boundary condition
    for j in range(1,ny-1):
        T_old[0][j]=0.00
        T_oldtime[0][j]=0.00
    #right boundary condition
    for j in range(1,ny-1):
        T_old[nx-1][j]=0.00
        T_oldtime[nx-1][j]=0.00
    
    
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            T_old[i][j]=T_oldtime[i][j]#the initial condition(guess) that we have specified, we are assigning that to T_old.... so that interations for solver can be started at time level tl=1
    
    #computation of b_old in outer loop because it will remain constant during the entire second loop, the second loop will be based on what value is obtained for B_old in the first loop
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            b_old[i][j]=T_oldtime[i][j] + a_s*T_oldtime[i][j-1] + a_w*T_oldtime[i-1][j] + (a_p-1)*T_oldtime[i][j] + a_e*T_oldtime[i+1][j] + a_n*T_oldtime[i][j+1]
            
    #solver
    sum2=0.0
    res2=1.0
    while res2>=final_s:
        for i in range(1,nx-1):
            for j in range(1,ny-1):
                T_new[i][j]=T_old[i][j] + (b_old[i][j]-(a_s*T_old[i][j-1] + a_w*T_old[i-1][j] + (a_p-1)*T_old[i][j] + a_e*T_old[i+1][j] + a_n*T_old[i][j+1]))/a_p
                sum2=sum2 + abs(T_new[i][j] - T_old[i][j])
                T_old[i][j]=T_new[i][j]
                
        res2=sum2
        count=count + 1
    
    #loop to assign the T_new to T_newtime
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            T_newtime[i][j]=T_new[i][j]
    #loop for setting outer loop residual and assigning newtime to oldtime
    for i in range(1,nx-1):
        for j in range(1,ny-1):       
            sum1=sum1+abs(T_newtime[i][j] - T_oldtime[i][j])
            res1=sum1
            T_oldtime[i][j]=T_newtime[i][j]
            
    tl=tl+1
    sum1=0.0#resetting it because not doing so will keep on increasing the value. we need to interate and iterate and go on untill a point of time comes when their difference
            #....is very small.
    
    
    
    
    
            
            
            
            
            
            
            
            
        
        
    
    