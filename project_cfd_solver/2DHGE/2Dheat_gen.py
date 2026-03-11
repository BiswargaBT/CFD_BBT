import numpy as np
lx=0.01 #[m]
ly=0.015 #[mm]
nx=41
ny=41
dx=lx/(nx-1)
dy=ly/(ny-1)
k=167#for aluminium 6061
q=5*(10**5)
alpha=-(q/k)

a_s=1/(dy*dy)
a_w=1/(dx*dx)
a_e=1/(dx*dx)
a_n=1/(dy*dy)
a_p=-2*(a_e + a_n)

x=np.zeros((nx,ny))
y=np.zeros((nx,ny))

#grid generation
for i in range(nx):
    for j in range(ny):
        x[i][j]=i*dx
        y[i][j]=j*dy
        
#initialization
T_old=np.zeros((nx,ny))

#iterations
res=1.0
count=0
final=(10)**(-3)
while res>final:
    T_new=np.zeros((nx,ny))
    sum1=0
    #for top boundary
    for i in range(nx):
        T_old[i][ny-1]=100.0
    
    #for bottom boundary
    for i in range(nx):
        T_old[i][0]=0.0
        
    #for left boundary
    for j in range(1,ny-1):
        T_old[0][j]=0.0
    
    #for right boundary
    for j in range(1,ny-1):
        T_old[nx-1][j]=0.0
        
    #point jacobi iteration updation algorithm for the inner nodes only
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            T_new[i][j]=(a_s*T_old[i][j-1] + a_w*T_old[i-1][j] + a_e*T_old[i+1][j] + a_n*T_old[i][j+1] + alpha)/a_p
            sum1=sum1+abs(T_new[i][j] - T_old[i][j])
            
    #point jacobi updation
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            T_old[i][j]=T_new[i][j]
            
    res=sum1
    count=count + 1
    
#for displaying values
print(f"{'x':<20}{'y':<20}{'T_old':<25}")
for i in range(nx):
    for j in range(ny):
        print(f"{x[i][j]:<20}{y[i][j]:<20}{T_old[i][j]:<25}")

