import numpy as np
import time as time
import matplotlib.pyplot as plt
lx=1.0
ly=1.5
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
a_p=-2*(1/(dx*dx) + 1/(dy*dy))

#grid generation
for i in range(nx):
    for j in range(ny):
        x[i][j]=i*dx
        y[i][j]=j*dy

#initialization of T_old at the iterior nodes only.... boundary nodes will be taken care of at the boundary conditions step
T_old=np.zeros((nx,ny))

#iterations

start_time=time.time()
residuals=[]
iterations=[]
final=10**(-6)
res=1
count=1
omega=1.80
while res>=final:
    T_new=np.zeros((nx,ny))
    sum1=0
    for i in range(nx):
        T_old[i][ny-1]=100.00
        
    for i in range(nx):
        T_old[i][0]=0.00
        
    for j in range(1,ny-1):
        T_old[0][j]=0.00
        
    for j in range(1,ny-1):
        T_old[nx-1][j]=0.00
        
    #SOR Iterations
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            T_new[i][j] = (1.00 - omega)*T_old[i][j] + omega*(a_s * T_old[i][j-1] + a_w * T_old[i-1][j] + a_e * T_old[i+1][j] + a_n * T_old[i][j+1])/a_p
            sum1 = sum1 + abs(T_new[i][j] - T_old[i][j])
            T_old[i][j] = T_new[i][j]
    res=sum1       
    iterations.append(count)
    residuals.append(res)
    count=count + 1
    
end_time=time.time()

with open("temp_field3.dat", "w") as f:
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
plt.title("Residuals v/s Iterations (Successive Over Relaxation)")#title of the entire plot
plt.yscale("log")#we are using log scale because the residual decreases exponentially, so we are using the exponential scale
plt.grid(True)#this is used to state whether grid is neccessary or not... if yes than true if not then false
plt.show()#end of the matlab figure

print(end_time-start_time, "seconds to converge")
print("The final residual is ", res)
print("The number of iterations required to converge is ", count)
            
            
    

    
                                                    
    
    
        
        
    
