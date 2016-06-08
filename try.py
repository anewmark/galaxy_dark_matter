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

#print('Defining Functions that formats radians')	#****************
#import math
#pi=math.pi
#for x in range(1,19):
#	degrees = x * degint
#	radians=degrees*pi/180
#	degs=str(degrees)
#	rads=str(radians)
#	def fixformat(rads)
#	while len(rads)<7:
#		rads=rads+"  "
#		return s
#	r=fixformat(rads)
#	print(r)