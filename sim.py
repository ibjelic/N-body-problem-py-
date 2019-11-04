import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
import math, time, matplotlib, random

matplotlib.use('Qt5Agg')
#Some important variables
dt = 1e-3#second
#assume that all distance units are in meters
T=1#time of simulation in seconds
animation_speed=1 #how much times faster should animation be
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

	def movement_angle(self): #angle of vector in our coordinate system
		if(math.isnan(self.v[0])):
			self.v[0]=0.0

		if(self.v[0]==0.0 and self.v[1]>0):
			angle = np.pi/2
		elif(self.v[0]==0.0  and self.v[1]<0):
			angle = -np.pi/2
		elif(self.v[1]==0):
			angle=0
		else:
			angle = np.arctan(self.v[1]/self.v[0])
		return angle

	def angle(self, body): #angle between this and other body in radians, also named contact angle
		x_dist = self.x-body.x; #one side of right angle triangle
		y_dist = self.y-body.y; #second side of same triangle
		#prevent division by 0
		if(x_dist==0 and y_dist<0): 
			angle=-np.pi/2
		elif(x_dist==0 and y_dist>0):
			angle=np.pi/2
		#elif(x_dist==0 and y_dist==0):
		#	angle=0
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
				if(self.angle(body[i])==np.pi/2 or self.angle(body[i])==-np.pi/2):
					pass
				else:
					self.f[0]=self.f[0]-gravity_constant*self.m*body[i].m/distance**2*(np.cos(self.angle(body[i]))) #gravity force x
				self.f[1]=self.f[1]-gravity_constant*self.m*body[i].m/distance**2*(np.sin(self.angle(body[i]))) #gravity force y

				if(self.angle(body[i])==np.pi/2 or self.angle(body[i])==-np.pi/2):
					pass
				else:
					self.f[0]=self.f[0]+1/(4*3.14*e0)*self.q*body[i].q/distance**2*(np.cos(self.angle(body[i])))#*(self.x-body[i].x)/np.abs(self.x-body[i].x) #electric force x
				self.f[1]=self.f[1]+1/(4*3.14*e0)*self.q*body[i].q/distance**2*np.sin(self.angle(body[i]))#*(self.y-body[i].y)/np.abs(self.y-body[i].y) #electric force y

		return self.f

	def rk4(self,f,dx): #dy/dx=f()
		k1=dx*f
		k2=dx*(k1*0.5+f)
		k3=dx*(k2*0.5+f)
		k4=dx*(k3+f)
		return (k1+2*k2+2*k3+k4)/6

	def acceleration(self):
		self.a[0]=self.f[0]/self.m; #x acceleration from F=m*a
		self.a[1]=self.f[1]/self.m; #y 

		self.v[0]=self.v[0]+self.rk4(self.a[0],dt) #x calculate new velocity
		self.v[1]=self.v[1]+self.rk4(self.a[1],dt) #y
		return self.a; 
    
	def energy(self,body):
		for i in range(len(body)):
			distance = self.dist(body[i])
			if(distance!=0):
				ep=-1/(4*3.14*e0)*self.q*body[i].q/distance-gravity_constant*self.m*body[i].m/distance
		ek=0.5*self.m*(self.v[0]**2+self.v[1]**2)
		return ek,ep

	def box_collision(self,box): #collision with box walls
		if(self.x>=box.right):
			self.v[0]=self.v[0]*(-1)
		if(self.x<=box.left):
			self.v[0]=self.v[0]*(-1)
		if(self.y>=box.up):
			self.v[1]=self.v[1]*(-1)
		if(self.y<=box.down):
			self.v[1]=self.v[1]*(-1)

	def collision(self,body): #detect collision
		for i in range(len(body)):
			distance = self.dist(body[i]) 	
			if(distance<=(self.r/2+body[i].r/2) and distance>0): #we are working with circles with radius so we dont get division by zero when bodies collide
				#print(distance)
				v1, v2=math.sqrt(self.v[0]**2+self.v[1]**2),math.sqrt(body[i].v[0]**2+body[i].v[1]**2)
				m1, m2=self.m, body[i].m
				phi = self.angle(body[i])
				theta1, theta2 = self.movement_angle(), body[i].movement_angle()
				v1xs, v1ys, v2xs, v2ys = np.sign(self.v[0]), np.sign(self.v[1]), np.sign(body[i].v[0]), np.sign(body[i].v[1])

				if(phi!=np.pi/2 and phi!=-np.pi/2):
					self.v[0]   =(v1*np.cos(theta1-phi)*(m1-m2)+2*v2*m2*np.cos(theta2-phi))/(m1+m2)*np.cos(phi)+v1*np.sin(theta1-phi)*np.cos(phi+np.pi/2)
					body[i].v[0]=(v2*np.cos(theta2-phi)*(m2-m1)+2*v1*m1*np.cos(theta1-phi))/(m2+m1)*np.cos(phi)+v2*np.sin(theta2-phi)*np.cos(phi+np.pi/2)
				self.v[1]   =(v1*np.cos(theta1-phi)*(m1-m2)+2*v2*m2*np.cos(theta2-phi))/(m1+m2)*np.sin(phi)+v1*np.sin(theta1-phi)*np.sin(phi+np.pi/2)
				body[i].v[1]=(v2*np.cos(theta2-phi)*(m2-m1)+2*v1*m1*np.cos(theta1-phi))/(m2+m1)*np.sin(phi)+v2*np.sin(theta2-phi)*np.sin(phi+np.pi/2)
				
				
				if(v1xs>v2xs and v1ys>v2ys):
					self.v[0]=self.v[0]*(-1)	
					self.v[1]=self.v[1]*(-1)
				if(v2xs>v1xs and v2ys>v1ys):
					body[i].v[0]=body[i].v[0]*(-1)
					body[i].v[1]=body[i].v[1]*(-1)
				if(v1xs>v2xs):
					self.v[0]=self.v[0]*(-1)
				if(v2xs>v1xs):
					body[i].v[0]=body[i].v[0]*(-1)
				if(v1xs<v2xs and v1ys>v2ys):
					body[i].v[1]=body[i].v[1]*(-1)
				if(v2xs<v1xs and v2ys>v1ys):
					self.v[1]=self.v[1]*(-1)
				
				while(self.dist(body[i])<=(self.r/2+body[i].r/2)): #prevent another body to detect collision
					self.move()


	def move(self):
		self.x=self.x+self.rk4(self.v[0],dt) #move for time step x
		self.y=self.y+self.rk4(self.v[1],dt) #move for time step y
		return [self.x, self.y] #return position for debugging

class box():
	def __init__(self, left,right,up,down):
		self.left = left
		self.right = right
		self.up = up
		self.down = down


def simulation(bodies, t, box):
	start = time.time()
	N=int(t/dt)
	print("%.0f iterations" % N)
	simulation_matrix = np.zeros([N,len(bodies),5])
	for k in range(N):
		for i in range(len(bodies)):
			bodies[i].force(bodies) #calculate force on each body
			bodies[i].box_collision(box)
			bodies[i].acceleration() #calculate new velocity for each body
			bodies[i].collision(bodies) #check for collisions
			bodies[i].move() #move all bodies for time step
			#add positions and velocities into one matrix 
			simulation_matrix[k][i][0], simulation_matrix[k][i][1] = bodies[i].x, bodies[i].y
			simulation_matrix[k][i][2], simulation_matrix[k][i][3] = bodies[i].v[0], bodies[i].v[1]
			simulation_matrix[k][i][4] = np.sign(bodies[i].q)
		
		#proogress info
		if(k==int(N/4)):  
			end_ = time.time()
			print("25%s done, ETA = %.3f seconds" %("%",(end_-start)*3)) 
		if(k==int(N/2)): 
			end_ = time.time()
			print("50%s done, ETA = %.3f seconds" %("%",(end_-start)))
		if(k==int(N*3/4)):
			end_ = time.time()
			print("75%s done, ETA = %.3f seconds" %("%",(end_-start)/4))
		
	end = time.time()
	print("100%s done, time taken to simulate: %.3f seconds" % ("%",end-start))
	return simulation_matrix #get back the data

#create box:
box = box(-100,100,100,-100)

#generate random coordinates
test = []

number_of_bodies = 20

xcor, ycor = random.sample(range(box.left+5,box.right-5), number_of_bodies), random.sample(range(box.down+5,box.up-5), number_of_bodies)
mass, charge = random.sample(range(100),number_of_bodies), random.sample(range(-100,100), number_of_bodies)
for i in range(number_of_bodies):
	test.append(body(xcor[i],ycor[i],mass[i],charge[i]*1e-4,1,0,0))

#test = [body(-10,10,1,1e-2,1,5,-5),body(10,-10,1,1e-2,1,-5,5)]

#simulation start
print("Starting simulation with %.0f bodies for time period of %.2f seconds" %(len(test),T))
simulation_data = simulation(test,T, box) 


#Create animation
colors = cm.rainbow(np.linspace(0, 1, len(test))) #create colors list 

def update(frame):	#animation frames update
	frame=frame*animation_speed #less frames -> faster animation
	ax1.clear() #clear all dots and plot them again, imitating movment
	ax1.set_xlim(box.left,box.right)
	ax1.set_ylim(box.down,box.up)
	ax2.set_xlabel('T [s]')
	ax2.set_ylabel('V [pixel/s]')
	ax1.set_xlabel('x [pixel]')
	ax1.set_ylabel('y [pixel]')
	for i in range(len(test)):
		ax1.plot(simulation_data[frame][i][0],simulation_data[frame][i][1],c=colors[i],marker='o',markersize=5)
		if(simulation_data[frame][i][4]==-1):
			ax1.annotate('-',(simulation_data[frame][i][0],simulation_data[frame][i][1]))
		elif(simulation_data[frame][i][4]==+1):
			ax1.annotate('+',(simulation_data[frame][i][0],simulation_data[frame][i][1]))
	for i in range(len(test)):
		ax2.plot(frame*dt,math.sqrt(simulation_data[frame][i][2]**2+simulation_data[frame][i][3]**2), c=colors[i],marker='*', markersize=2)


try:
	animation = FuncAnimation(fig, update, interval=50, frames=int(len(simulation_data)/animation_speed), repeat=False)
except Exception as e: print(e)
plt.show()

