## first day of learning to use Python
#convert degrees to radians and arrange them into columns
degint=10
print('IS THIS WORKING?!')
print('  Degrees  Radians')
for x in range(1,19):
	degrees = x * degint
	pi= 3.14159
	radians=degrees*pi/180
	degs=str(degrees)
	if len(degs)==2:
		degs= "   "+degs
	elif len(degs)==3:
		degs="  "+degs
	rads=str(radians)
	deg= 'degrees= ' + degs
	rad= 'radians= ' + rads
#	print(deg, rad)
	print(degs+"      "+rads[:7])	#this only shows the first 7 places of rads
#	data1=[['degrees', 'radians'],[degs, rads]]
#	data=[degs, rads]
#	print(data)
#	print([[header],[data]])

print('Now the import math exercise')
import math
pi=math.pi
print('  Degrees  Radians    Cosine')
for x in range(1,19):
	degrees = x * degint
	radians=degrees*pi/180
	degs=str(degrees)
	if len(degs)==2:
		degs= "   "+degs
	elif len(degs)==3:
		degs="  "+degs
	rads=str(radians)
	cosrads=math.cos(radians)
	cosrad=str(cosrads)
#	print(deg, rad)
	print(degs+"      "+rads[:7]+"    "+cosrad)

print('Working with lists')
import math
pi=math.pi
listhead=['Degrees', 'Radians', 'Cosine']
radlist=[]	#start with empty lists
coslist=[]
for x in range(1,19):
	degrees = x * degint
	radians=degrees*pi/180
	radlist.append(radians)	#adds each radian to the list
	cosrads=math.cos(radians)
	coslist.append(cosrads)
print(radlist)
print(coslist)

print('creating an array')
import numpy
a=numpy.zeros((6,3))
for i in range(6):
	for j in range (3):
		a[i][j]=i+j
print(a)

print('creating a rotation matrix of 30º') #useful when creating array of non-integers
rot=numpy.empty((2,2))
anglerad=30*pi/180	
rot[0][0]=math.cos(anglerad)
rot[0][1]=math.sin(anglerad)
rot[1][0]=-math.sin(anglerad)
rot[1][1]=math.cos(anglerad)
print(rot)

print('Creating a circle')
circlearr=numpy.empty((18,2))
b=numpy.linspace(0,360,18)
for i in range(18):
	angle=0
	incr=2*pi/18	#already in radians, no need to convert
	angle+=incr	#original angle plus the increments of incr. it is additive
	#can also do:
	#	ang=0
	#	incr=2*pi/18
	#	angle=ang+i*incr
	x=math.cos(angle)
	y=math.sin(angle)
circlearr[i][0]=x
circlearr[i][1]=y
print(circlearr)


print('Now rotating this circle and scaling by 2')
rotcirc=numpy.dot(circlearr, rot)
scalerotcirc=2.*rotcirc
print(scalerotcirc)


print('Working with matplotlib')
import numpy
import matplotlib.pyplot as plt
x=numpy.array(range(10))
y=3*x
plt.plot(x,y)
#plt.show()


