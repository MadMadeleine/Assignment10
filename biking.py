import numpy as np
import matplotlib.pyplot as pp
import gpxpy as gp

# Define the different classes used

class class_coordinate:
	def __init__(self,lat,lon,ele):
		self.lat = lat
		self.lon = lon
		self.ele = ele

class class_point:
	def __init__(self,x,y,z):
		self.x = x
		self.y = y
		self.z = z

class class_marker:
	def __init__(self,index,distance,distance_along,grade):
		self.i = index
		self.d = distance
		self.da = distance_along
		self.g = grade

def GeocentricLatitude(latitude):
	e2 = 0.00669437999014
	return np.atan((1-e2) * np.tan(latitude))

def Location_to_Point(c):
	lat = c.lat*np.pi/180
	lon = c.lon*np.pi/180
	radius = 6467 * 1000
	clat = GeocentricLatitude(lat)

	cosLon = np.cos(lon)
	sinLon = np.sin(lon)
	cosLat = np.cos(clat)
	sinLat = np.sin(clat)

	x = radius * cosLon * cosLat
	y = radius * sinLon * cosLat
	z = radius * sinLat

	cosGlat = np.cos(lat)
	sinGlat = np.sin(lat)

	nx = cosGlat * cosLon
	ny = cosGlat * sinLon
	nz = sinGlat

	x += c.ele * nx
	y += c.ele * ny
	z += c.ele * nz

	point = class_point(x,y/100000,z)

	return point

def Marker_Generator(index, p1, p2, da_previous):
	dx = p2.x - p1.x
	dy = p2.y - p1.y
	dz = p2.z - p1.z

	distance = np.sqrt((dx**2) + (dy**2) + (dz**2))
	da = distance + da_previous
	grade = dz / (np.sqrt((dx**2) + (dy**2)) * 100)

	marker = class_marker(index, distance, da, grade)

	return(marker)

# Define the functions used in solving for the motion of the bicycle

def grade_given_pos(position,mlist):
	iterator = 0
	while(mlist[iterator].da < position):
		iterator +=1
	return(mlist[iterator].g)

def gravity(mass,grade):
	return -1*mass*9.8*np.sin(np.atan(grade))

def drag_aerodynamic(v,cd,den,a):
	return ((-1*0.00002*a*(v/2))+(-0.5*cd*den*a*v*v))

def ode(power,mass,velocity,Cd,density,A,grade,roll_resist):
	return(power/(mass*velocity)+((drag_aerodynamic(velocity,Cd,density,A))/mass)+(gravity(mass,grade)/mass)-(roll_resist/mass))
	
def euler(vn,ts,p,m,cd,den,a,g,rr):
	return vn + (ts*ode(p,m,vn,cd,den,a,g,rr))


#open and parse the gpx file
gpx_file = open('SS-Getaria.gpx','r')
gpx = gp.parse(gpx_file)

grade_list = []
distpoints = []


clist = []

lat_temp = 0
lon_temp = 0
ele_temp = 0 

# Parse the gpx file for the longitudes, longitudes, and elevations

for track in gpx.tracks:
	for segment in track.segments:
		for point in segment.points:
			lat_temp = point.latitude
			lon_temp = point.longitude
			ele_temp = point.elevation

			c = class_coordinate(lat_temp,lon_temp,ele_temp)

			clist.append(c)

# Convert the list of coordinates into a list of points and markers

plist = []
mlist = []

for i in range(len(clist)):
	plist.append(Location_to_Point(clist[i]))
	if(i == 0):
		mlist.append(Marker_Generator(i, plist[i], plist[i-1], 0))
	else:
		mlist.append(Marker_Generator(i, plist[i], plist[i-1], mlist[i-1].da))

# Remove the duplicate points

removal_list = []

for i in range(len(mlist)):
	if(mlist[i].d == 0):
		removal_list.append(i)

plist = np.delete(plist,removal_list)
clist = np.delete(clist,removal_list)
mlist = np.delete(mlist,removal_list)

# Initialize the figure

figure = pp.figure(figsize=(7, 5), layout="constrained")
spec = figure.add_gridspec(ncols=2, nrows=2)

# Define some constants

Cd = 0.9 #Drag Coeffcicient for aerodynamic drag equation
density = 1.29 #Air Density
A = 0.33 #Value for A in aerodynamic drag equation

initial_velocity = 1
new_velocity = initial_velocity
velocity = [initial_velocity] #list for all the velocities
max_velocity = initial_velocity

initial_position = 0
new_position = initial_position
position = [initial_position]

mass = 70
power = 1162#Power output of the rider (watts)
tire_resistance = 25 #rolling resistance for the two tires combined (watts)

timestep = 0.1
time = [0] #list for all the times
	
# Start creating the velocity and position lists

iterator = 0
index_final = len(mlist) - 1

glist = [0]

sum_of_distances = 0

for i in range(len(mlist)):
	sum_of_distances += mlist[i].d
	if(mlist[i].g > 3):
		mlist[i].g = 3
	elif(mlist[i].g < -3):
		mlist[i].g = -3
	
while(position[iterator] < sum_of_distances):
	time.append(iterator * timestep)
	
	grade = grade_given_pos(position[iterator],mlist)
	glist.append(grade_given_pos(position[iterator],mlist))

	velocity.append(euler(new_velocity,timestep,power,mass,Cd,density,A,grade,tire_resistance))
	new_velocity = euler(new_velocity,timestep,power,mass,Cd,density,A,grade_given_pos(position[iterator-1],mlist),tire_resistance)

	position.append(position[iterator] + (new_velocity * timestep))
	iterator += 1

position_graph = figure.add_subplot(spec[0,0])
position_graph.set_xlabel('Time (seconds)')
position_graph.set_ylabel('Position (meters)')

position_graph.plot(time,position,'r-')

velocity_graph = figure.add_subplot(spec[0,1])
velocity_graph.set_xlabel('Time (seconds)')
velocity_graph.set_ylabel('Velocity (m/s)')

velocity_graph.plot(time,velocity,'r-')

grade_graph = figure.add_subplot(spec[1,0])
grade_graph.set_xlabel('Time (seconds)')
grade_graph.set_ylabel('Grade (degrees)')

grade_graph.plot(time,glist,'r-')

pp.show()
