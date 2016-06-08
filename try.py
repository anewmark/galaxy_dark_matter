degint=10
print('IS THIS WORKING?!')
for x in range(0,18):
	degrees = x * degint
	pi= 3.14159
	radians=degrees*pi/180
	degs=str(degrees)
	rads=str(radians)
	deg= 'degrees= ' + degs
	rad= 'radians= ' + rads
	print(deg, rad)
print('THE END')