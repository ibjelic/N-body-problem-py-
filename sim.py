import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.animation import FuncAnimation

dt = 0.1 #second
#assume that all distance units are in meters

e0 = 8.85e-12; #epsilon 0 
gravity_constant = 6.67e-11; #gravity constant
fig = plt.figure() #plot
ax1 = fig.add_subplot(1,1,1)

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
		if(y_dist<1e-5):
			angle=3.14/2*np.sign(x_dist)
		else:
			angle = np.arctan(x_dist/y_dist); #division of sides gives tan(angle)
		return angle;
    
	def dist(self,body):
		dis = math.sqrt((self.x-body.x)**2+(self.y-body.y)**2) #distance between two points
		return dis;
    
	def force(self, body): #calculate force on body
		for i in range(len(body)):
			distance = self.dist(body[i]) 
			if(distance!=0): #prevent division by 0 (no force for one body on that same body)
				self.f[0]=self.f[0]+gravity_constant*self.m*body[i].m/distance**2*np.sin(self.angle(body[i])) #gravity force x
				self.f[0]=self.f[0]+1/(4*3.14*e0)*self.q*body[i].q/distance**2*np.sin(self.angle(body[i])) #electric force x
				self.f[1]=self.f[1]+gravity_constant*self.m*body[i].m/distance**2*np.cos(self.angle(body[i])) #gravity force y
				self.f[1]=self.f[1]+1/(4*3.14*e0)*self.q*body[i].q/distance**2*np.cos(self.angle(body[i])) #electric force y
		
		return self.f;
	
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
				self.v[0]=v1*np.sin(self.angle(body[i]))
				self.v[1]=v1*np.cos(self.angle(body[i]))
				body[i].v[0]=v2*np.sin(body[i].angle(self))
				body[i].v[1]=v2*np.cos(body[i].angle(self))
				while(self.dist(body[i])<=(self.r/2+body[i].r/2)):
					self.move()
					body[i].move()
					print(distance)
	def move(self):
		self.x=self.x+self.v[0]*dt #move for time step x
		self.y=self.y+self.v[1]*dt #move for time step y
		return [self.x, self.y] #return position for debugging




def simulation(bodies, time):
	N=int(time/dt)
	bodies_next = bodies.copy();
	xcor1, ycor1 = np.zeros(N), np.zeros(N)
	xcor2, ycor2 = np.zeros(N), np.zeros(N)
	xcor3, ycor3 = np.zeros(N), np.zeros(N)
	for k in range(N):
		for i in range(len(bodies)):
			bodies[i].force(bodies) #calculate force on each body
			bodies[i].collision(bodies) #check for collisions
			bodies[i].acceleration() #calculate new velocity for each body
			bodies[i].move() #move all bodies for time step
		xcor1[k] = bodies[0].x
		ycor1[k] = bodies[0].y
		xcor2[k] = bodies[1].x
		ycor2[k] = bodies[1].y
		xcor3[k] = bodies[2].x
		ycor3[k] = bodies[2].y
	return xcor1, ycor1, xcor2, ycor2, xcor3,ycor3
nis = [body(10,100, 1,0,1,-1,0),body(0,0,1,0,1,0,0),body(0,100,1,0,1,1,0) ]
x1,y1,x2,y2,x3,y3= simulation(nis,1000)

def update(frame):	
	ax1.clear()
	ax1.plot(x1[frame],y1[frame], c='r',marker='o',markersize=10)
	ax1.plot(x2[frame],y2[frame], c='g',marker='o',markersize=10)
	ax1.plot(x3[frame],y3[frame], c='b',marker='o',markersize=10)
animation = FuncAnimation(fig, update, interval=50)
plt.show()

