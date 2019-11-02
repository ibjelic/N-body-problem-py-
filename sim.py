import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
import math, time, matplotlib, random


#Some important variables
dt = 1e-2 #second
#assume that all distance units are in meters
T=10#time of simulation in seconds
animation_speed=3 #how much times faster should animation be
e0 = 8.85e-12; #epsilon 0 
gravity_constant = 6.67e-11; #gravity constant

#For animation
fig = plt.figure() #plot
ax1 = fig.add_subplot(1,2,1) #graph for bodies movment
ax2 = fig.add_subplot(1,2,2) #graph velocities=f(time)

class body:
	def __init__(self, x, y, m, q, r, vx, vy):
		self.x = x; #x coordinate
		self.y = y; #y coordinate
		self.m = m; #mass
		self.q = q; #charge
		self.r = r; #radius
		self.v = [vx,vy]; #velocity
		self.f = [0.0, 0.0]; #force acting on body
		self.a = [0.0, 0.0]; #acceleration
	
	def angle(self, body): #angle between this and other body in radians
		x_dist = self.x-body.x; #one side of right angle triangle
		y_dist = self.y-body.y; #second side of same triangle
		#prevent division by 0
		if(x_dist==0 and y_dist<0): 
			angle=3.14/2
		elif(x_dist==0 and y_dist>0):
			angle=-3.14/2
		else:
			angle = np.arctan(y_dist/x_dist); #division of sides gives tan(angle)
		return angle
    
	def dist(self,body): #distance between two points (self and other body)
		dis = math.sqrt((self.x-body.x)**2+(self.y-body.y)**2) 
		return dis
    
	def force(self, body): #calculate force on body
		for i in range(len(body)):
			distance = self.dist(body[i]) 
			if(distance!=0): #prevent division by 0 (no force for one body on that same body)
				self.f[0]=self.f[0]+gravity_constant*self.m*body[i].m/distance**2*np.sin(self.angle(body[i])) #gravity force x
				self.f[0]=self.f[0]+1/(4*3.14*e0)*self.q*body[i].q/distance**2*np.sin(self.angle(body[i])) #electric force x
				self.f[1]=self.f[1]+gravity_constant*self.m*body[i].m/distance**2*np.cos(self.angle(body[i])) #gravity force y
				self.f[1]=self.f[1]+1/(4*3.14*e0)*self.q*body[i].q/distance**2*np.cos(self.angle(body[i])) #electric force y
		return self.f
	
	def acceleration(self):
		self.a[0]=self.f[0]/self.m; #x acceleration from F=m*a
		self.a[1]=self.f[1]/self.m; #y 
		self.v[0]=self.v[0]+self.a[0]*dt #x calculate new velocity
		self.v[1]=self.v[1]+self.a[1]*dt #y
		return self.a; 
    
	def energy(self,body):
		for i in range(len(body)):
			distance = self.dist(body[i])
			if(distance!=0):
				ep=-1/(4*3.14*e0)*self.q*body[i].q/distance-gravity_constant*self.m*body[i].m/distance
		ek=0.5*self.m*(self.v[0]**2+self.v[1]**2)
		return ek,ep

	def collision(self,body): #detect collision
		for i in range(len(body)):
			distance = self.dist(body[i]) 	
			if(distance<=(self.r/2+body[i].r/2) and distance>0): #we are working with circles with radius so we dont get division by zero when bodies collide
				#print(distance)
				v1=(self.m-body[i].m)/(self.m+body[i].m)*math.sqrt(self.v[0]**2+self.v[1]**2)+2*body[i].m/(self.m+body[i].m)*math.sqrt(body[i].v[0]**2+body[i].v[1]**2)
				v2=-(self.m-body[i].m)/(self.m+body[i].m)*math.sqrt(body[i].v[0]**2+body[i].v[1]**2)+2*self.m/(self.m+body[i].m)*math.sqrt(self.v[0]**2+self.v[1]**2)
				self.v[0]=v1*np.cos(self.angle(body[i]))
				self.v[1]=v1*np.sin(self.angle(body[i]))
				body[i].v[0]=v2*np.cos(body[i].angle(self))
				body[i].v[1]=v2*np.sin(body[i].angle(self))
				while(self.dist(body[i])<=(self.r/2+body[i].r/2)): #prevent another body to detect collision
					self.move()
					body[i].move()
	def move(self):
		self.x=self.x+self.v[0]*dt #move for time step x
		self.y=self.y+self.v[1]*dt #move for time step y
		return [self.x, self.y] #return position for debugging



def simulation(bodies, t):
	start = time.time()
	N=int(t/dt)
	print("%.0f iterations" % N)
	simulation_matrix = np.zeros([N,len(bodies),5])
	for k in range(N):
		for i in range(len(bodies)):
			bodies[i].force(bodies) #calculate force on each body
			bodies[i].collision(bodies) #check for collisions
			bodies[i].acceleration() #calculate new velocity for each body
			bodies[i].move() #move all bodies for time step
			#add positions and velocities into one matrix 
			simulation_matrix[k][i][0], simulation_matrix[k][i][1] = bodies[i].x, bodies[i].y
			simulation_matrix[k][i][2], simulation_matrix[k][i][3] = bodies[i].v[0], bodies[i].v[1]
			simulation_matrix[k][i][4] = bodies[i].r
		
		#proogress info
		if(k==int(N/4)):  
			end_ = time.time()
			print("25%s done, ETA = %.3f s, Time in simulation: %.3f" %("%",(end_-start)*3,k*dt)) 
		if(k==int(N/2)): 
			end_ = time.time()
			print("50%s done, ETA = %.3f s, Time in simulation: %.3f" %("%",(end_-start),k*dt))
		if(k==int(N*3/4)):
			end_ = time.time()
			print("75%s done, ETA = %.3f s, Time in simulation: %.3f" %("%",(end_-start)/4,k*dt))
		if(k==0): 
			end_ = time.time()
			print('Time for one iteration: %.3f seconds, ETA: %.3f seconds' %(end_-start,(end_-start)*N))
		
	end = time.time()
	print("100%s done, time taken to simulate: %.3f seconds" % ("%",end-start))
	return simulation_matrix #get back the data

#generate random coordinates
test = []
number_of_bodies = 10
xcor, ycor = random.sample(range(100), number_of_bodies), random.sample(range(100), number_of_bodies)
mass, charge = random.sample(range(50),number_of_bodies), random.sample(range(-50,50), number_of_bodies)
for i in range(number_of_bodies):
	test.append(body(xcor[i],ycor[i],mass[i],charge[i]*1e-1,1,0,0))

#simulation start
print("Starting simulation with %.0f bodies for time period of %.2f seconds" %(len(test),T))
simulation_data = simulation(test,T) 


#Create animation
colors = cm.rainbow(np.linspace(0, 1, len(test))) #create colors list 

def update(frame):	#animation frames update
	frame=frame*animation_speed #less frames -> faster animation
	ax1.clear() #clear all dots and plot them again, imitating movment
	ax2.set_xlabel('T [s]')
	ax2.set_ylabel('V [pixel/s]')
	ax1.set_xlabel('x [pixel]')
	ax1.set_ylabel('y [pixel]')
	for i in range(len(test)):
		ax1.plot(simulation_data[frame][i][0],simulation_data[frame][i][1],c=colors[i],marker='o',markersize=5)
	for i in range(len(test)):
		ax2.plot(frame*dt,math.sqrt(simulation_data[frame][i][2]**2+simulation_data[frame][i][3]**2), c=colors[i],marker='o', markersize=5)



animation = FuncAnimation(fig, update, interval=44, frames=int(len(simulation_data)/animation_speed), repeat=False)
plt.show()

