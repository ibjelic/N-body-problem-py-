import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
import math, time

dt = 1e-3 #second
#assume that all distance units are in meters
T=2#time of simulation in seconds
animation_speed=10 #how much times faster should animation be
e0 = 8.85e-12; #epsilon 0 
gravity_constant = 6.67e-11; #gravity constant
fig = plt.figure() #plot
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

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
		if(x_dist==0 and y_dist<0):
			angle=3.14/2
		elif(x_dist==0 and y_dist>0):
			angle=-3.14/2
		else:
			angle = np.arctan(y_dist/x_dist); #division of sides gives tan(angle)
		return angle
    
	def dist(self,body):
		dis = math.sqrt((self.x-body.x)**2+(self.y-body.y)**2) #distance between two points
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
				while(self.dist(body[i])<=(self.r/2+body[i].r/2)):
					self.move()
					body[i].move()
				print(distance)
				print(self.x,body[i].x)
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
		if(k==int(N/4)):  #debugging info 
			end_ = time.time()
			print("25%s done, ETA = %.3f s, Time in simulation: %.3f" %("%",(end_-start)*3,k*dt)) 
		if(k==int(N/2)): 
			end_ = time.time()
			print("50%s done, ETA = %.3f s, Time in simulation: %.3f" %("%",(end_-start),k*dt))
		if(k==int(N*3/4)):
			end_ = time.time()
			print("75%s done, ETA = %.3f s, Time in simulation: %.3f" %("%",(end_-start)/4,k*dt))
		for i in range(len(bodies)):
			bodies[i].force(bodies) #calculate force on each body
			bodies[i].collision(bodies) #check for collisions
			bodies[i].acceleration() #calculate new velocity for each body
			bodies[i].move() #move all bodies for time step
			simulation_matrix[k][i][0], simulation_matrix[k][i][1] = bodies[i].x, bodies[i].y
			simulation_matrix[k][i][2], simulation_matrix[k][i][3] = bodies[i].v[0], bodies[i].v[1]
			simulation_matrix[k][i][4] = bodies[i].r
	end = time.time()
	print("100%s done, time taken to simulate: %.3f seconds" % ("%",end-start))
	return simulation_matrix

test = [body(0,10.001, 20,1e-4,5,0,-10),body(0,0,1,1e-5,5,0,0),body(50,100,1,-1e-2,5,-10,-20) ]
print("Starting simulation with %.0f bodies for time period of %.2f seconds" %(len(nis),T))

simulation_data = simulation(nis,T) #simulation start
colors = cm.rainbow(np.linspace(0, 1, len(nis))) #create colors list 

def update(frame):	#animation frames update
	frame=frame*animation_speed #less frames -> faster animation
	ax1.clear() #clear all dots and plot them again, imitating movment
	for i in range(len(nis)):
		ax1.plot(simulation_data[frame][i][0],simulation_data[frame][i][1],c=colors[i],marker='o',markersize=5)
	for i in range(len(nis)):
		ax2.plot(frame*dt,math.sqrt(simulation_data[frame][i][2]**2+simulation_data[frame][i][3]**2), c=colors[i],marker='o', markersize=5)



animation = FuncAnimation(fig, update, interval=1, frames=int(len(simulation_data)/animation_speed), repeat=False)
plt.show()

