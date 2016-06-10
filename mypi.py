import numpy as np
import math
import matplotlib.pyplot as plt
import timeit
points=[100,1000,10000,100000]
a=[0,points[0]-1]
b=[points[0],points[1]-1]
c=[points[1],points[2]-1]
d=[points[2],points[3]-1]
bins=[a,b,c,d]
print(bins)
sqrt=math.sqrt
inside=0
for npoints in points:
    for i in range(0,npoints):
        circ=np.random.rand(1,2)
        x=circ[:,0]
        y=circ[:,1]
        #print('x is', x)
        r=sqrt(x**2+y**2)
        #print('circ is', circ)
        #print('this was',i)
        if r<=1:
            inside+=1
            #print('good radius is ',r)
    pi=4*inside/npoints
    print('pi is', pi)
times=np.empty((4,1))
for i in range(0,4):
    time=timeit.timeit(number=points[i])
    print('calculates in',time, ' seconds')
    times[i]=time
print(times)

N=4 #four separate bars
ind = np.arange(N) #x location for bar groups
width=10 #width of each bar graph
rects = plt.bar(points, times, width=1)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Number of Points')
plt.ylabel('Time (s)')
plt.suptitle("The amount of time it takes Python to calculate pi according to number of sources")
#plt.bar(points, times)
plt.show()


