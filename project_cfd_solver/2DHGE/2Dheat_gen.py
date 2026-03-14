#solved with Gauss_seidal
import numpy as np
import time as time
import matplotlib.pyplot as plt #we aare importing the pyplot module of the matplotlib library
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
start_time=time.time()
res=1.0
count=0
final=(10)**(-3)
residuals=[]
iterations=[]
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
    #we are importing all the values of residuals and iteratin counts into a list so that plot can be made using matplotlib
    residuals.append(res)
    iterations.append(count)
    count=count + 1
    
end_time=time.time()
    
#for displaying values
#print(f"{'x':<20}{'y':<20}{'T_old':<25}")
#for i in range(nx):
    #for j in range(ny):
        #print(f"{x[i][j]:<20}{y[i][j]:<20}{T_old[i][j]:<25}")
        
#for making a .dat file for post processing work
with open("temp_field2.dat","w") as f:
    f.write('VARIABLES: "X", "Y", "T"\n')
    f.write('ZONE: I = {nx}, J = {ny}')
    for i in range(nx):
        for j in range(ny):
            f.write(f"{x[i][j]} {y[i][j]} {T_old[i][j]}\n")
show=True#in order to show the graph we are using this variable and since we have to show the graph, thus we are assigning it a TRUE value            
plt.figure()#starting of a matlab figure
plt.plot(iterations,residuals)#what to plot, what values on x axis and what values on the y axis
plt.xlabel("Iterations")#title of the x axis
plt.ylabel("Residuals")#title of the y axis
plt.title("Residuals v/s Iterations (Gauss-Seidal)")#title of the entire plot
plt.yscale("log")#we are using log scale because the residual decreases exponentially, so we are using the exponential scale
plt.grid(True)#this is used to state whether grid is neccessary or not... if yes than true if not then false
plt.show()#end of the matlab figure


print(end_time-start_time, "seconds to converge")
print("The final residual is ", res)
print("The number of iterations required to converge is ", count)
      
      
      
      
      
      
      
      
  
