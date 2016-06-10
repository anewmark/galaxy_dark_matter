## still working on this
import numpy as np
import math
import matplotlib.pyplot as plt
import timeit
points=[10, 100,1000,10000,100000]
sqrt=math.sqrt
inside=0
repeat=1
Nruns=1000
if repeat:
	points=[10, 100,1000,10000,100000]
	length=len(points)
	pis=[]
	pisst=[]
	for i in range(0,length):
		pisj=[]
		for j in range(0,Nruns):
			N=points[i]
			x=-1+2*np.random.rand(N)
			y=-1+2*np.random.rand(N)
			inside=(x**2+y**2<1).sum()
			pi=4*inside/N
			#print('pi is', pi)
			pisj.append(pi)
		#print(pism)
		pisj=np.array(pisj)
		print(pisj.mean(), pisj.std())
		#print(pisj)
		pis.append(pi)
		pisst.append(pisj)
	print(pisst)
	print(len(pisst))
	pisst=np.array(pisst)
	#print(pisst.mean(), pisst.std())
	



darts=0
if darts:
	plt.plot(pis, "bD", xdata=pis, ydata=0)
	plt.xlabel('Value of Pi')
	pix=math.pi
	piy=0
	plt.plot(pix, piy, 'ro', label='pi')
	plt.show()
    
times=np.empty((5,1))
for i in range(0,5):
    time=timeit.timeit(number=points[i])
    print('calculates in',time, ' seconds')
    times[i]=time
#print(times)

timetable=0
if timetable:
	N=5 #four separate bars
	width=10 #width of each bar graph
	rects = plt.bar(points, times, width=10, color='b')
	plt.xscale('log')
	#plt.yscale('log')
	plt.xlabel('Number of Points (Log10)')
	plt.ylabel('Time (s)')
	plt.suptitle("Time Calculation based on Number of Points", fontsize=12)
	plt.show()


##another way to calculate pi
