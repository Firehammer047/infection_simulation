'''
Infectious disease simulation in enclosed space.

Dynamic time series simulation inside a box, with particles
entering and exiting. 

- Only contagious people can infect others.
- People who are infected inside are not contagious
- Transmissibility decays as 1/r or 1/r^2
- People take the shortest route possible to their destination

Version 1.7

Copyright Tony Cabrera 2020

'''
import sys
import numpy as np
import matplotlib.pyplot as plt
import time
t1 = time.time()

MAKE_PLOT = 0
PRINT_DATA = 0

R0 = 1.0	# infection is exp(-20.0x/R0) for 1 < R0 < 6.5

####### INIT ##########
MAX_X = 5	# X dimensions of dining area
MAX_Y = 5	# Y dimensions of dining area
X_IN_DOOR = 0.5*MAX_X 
Y_IN_DOOR = 0.0
X_OUT_DOOR = 0.5*MAX_X
Y_OUT_DOOR = 0.0
#Y_OUT_DOOR = MAX_Y
REAL_MAX_CAP = 100
CAP_FACTOR = 0.25
MAX_CAP = REAL_MAX_CAP * CAP_FACTOR
TABLE_COORD = [ [1,1],[1,2],[1,3],[1,4],[2,1],[2,2],[2,3],[2,4],[3,1],[3,2],[3,3],[3,4],[4,1],[4,2],[4,3],[4,4] ]
X_TABLE_COORD = []
#TABLE_COORD = [ [1,1],[1,3],[2,2],[2,4],[3,1],[3,3],[4,2],[4,4] ]
#X_TABLE_COORD = [ [1,2], [1,4], [2,1], [2,3], [2,4], [3,2], [3,4], [4,1], [4,3] ] 
ALL_TABLE_COORD = [ [1,1],[1,2],[1,3],[1,4],[2,1],[2,2],[2,3],[2,4],[3,1],[3,2],[3,3],[3,4],[4,1],[4,2],[4,3],[4,4] ]
TABLES = len(TABLE_COORD)
TABLE_RADIUS = 0.2
AVG_V = 0.11	# walking speed in meters/sec
EPSILON = 0.15
AIR_EXCHANGES = 3
VOLUME = MAX_X*MAX_Y*3

DT = 1
END_TIME = 3600*8	# in seconds
#END_TIME = 300	# in seconds
TIME_STEP = int(1.0*END_TIME / 40)

AVG_VISIT_TIME = 30	# in minutes
SD_VISIT_TIME = 10	# in minutes
	
AVG_CUST_PER_DAY =  100
A1 = AVG_CUST_PER_DAY / 53476.0	# AVG_CUST/455 for 5 min day ; AVG_CUST/53476 for 8 hour day
#A1 = AVG_CUST_PER_DAY / 455.0
A2 = A1
C1 = 3600*1.0	# in seconds (1 hour  after open)
#C1 = 60*1.0
C2 = C1*7.0		# in seconds (7 hours  after open)
S1 = 3600.0		# SD of gaussian 1
#S1 = 60.0
S2 = S1			# SD of gaussian 2

PCT_INFECTED = 0.01

AVG_PEOPLE_IN_GROUP = 4
SD_PEOPLE_IN_GROUP = 1

#######################

class Person:

	def __init__(self,group,infected=0,contagious=0,table=-1,dest_x=0,dest_y=0,vec_x=0,vec_y=0,vx=0,vy=0,visit_time=0,action="WAITING",x=X_IN_DOOR,y=Y_IN_DOOR):
		self.group = group
		self.infected = infected
		self.contagious = contagious
		self.table = table
		self.dest_x = dest_x
		self.dest_y = dest_y
		self.vec_x = vec_x
		self.vec_y = vec_y
		self.vx = vx
		self.vy = vy
		self.visit_time = visit_time
		self.action = action
		self.x = x
		self.y = y
		self.active = 1
		self.elapsed_time = 0
	
	def print_stats(self):
		print("group: "),
		print(self.group)
		#print "infected: ",
		#print(self.infected)
		#print "contagious: ",
		#print(self.contagious)
		print "table: ",
		print(self.table)
		#print "vec_x: ",
		#print(self.vec_x)
		#print "vec_y: ",
		#print(self.vec_y)
		#print "vx: ",
		#print(self.vx)
		#print "vy: ",
		#print(self.vy)
		print "visit_time: ",
		print(self.visit_time)
		print "action: ",
		print(self.action)
		#print "x: ",
		#print(self.x)
		#print "y: ",
		#print(self.y)
		print "active: ",
		print(self.active)
		print "elapsed_time: ",
		print(self.elapsed_time)

	def d(self):
		return np.sqrt((self.dest_x - self.x)**2 + (self.dest_y - self.y)**2)
	
	def walk(self):
		self.x += self.vx * DT
		self.y += self.vy * DT

	def set_new_vectors(self):
		self.vec_x = self.dest_x - self.x
		self.vec_y = self.dest_y - self.y
		self.vx,self.vy = get_vx_vy(self.vec_x,self.vec_y)
		

def get_a_group(avg,sd):
	n = np.random.normal(avg,sd)
	return int(n)

def get_group_visit_time(avg,sd):
	t = np.random.normal(avg*60.0,sd*60.0)
	if t < 1:
		t = 1
	return t

def is_infected():
	r = np.random.rand()
	if(r < PCT_INFECTED):
		return 1
	return 0

def get_table():
	if(customers_inside >= MAX_CAP):
		return -1
	for i in range(TABLES):
		if(tables[i] == 0):
			tables[i] = 1
			return i
	return -1

def get_vx_vy(vec_x,vec_y):
	tangent = np.exp(100)
	if(vec_x != 0.0):
		tangent = 1.0*vec_y/vec_x
	sx = 1.0
	sy = 1.0
	if(vec_x < 0):
		sx = -1.0
	if(vec_y < 0):
		sy = -1.0
	vx = sx*AVG_V / np.sqrt(1 + tangent**2)
	vy = sy*np.sqrt(AVG_V**2 - vx**2)
	return (vx,vy)

def get_ortho_vx_vy(vx,vy,x,table_x):
	s = 1.0
	if(table_x > x):
		s = -1.0
	v = np.sqrt(vx**2 + vy**2)
	ortho_vx = s*vy/v
	ortho_vy = -(s*vx)/v
	return (ortho_vx,ortho_vy)

def get_infection(x1,y1,x2,y2):
	r = np.sqrt((y2-y1)**2 + (x2-x1)**2)
	p = np.exp(-20.0*r/R0)	
	if(np.random.rand() < p):
		return 1
	
	return 0
	
	#C = S*(1 - np.exp(-(K*I*t)/(AIR_EXCHANGES*VOLUME)))

def get_group_prob(t):
	g1 = A1*np.exp(-0.5*(1.0*(t-C1)/S1)**2) 
	g2 = A2*np.exp(-0.5*(1.0*(t-C2)/S2)**2)
	pdf = g1 + g2
	#if(pdf>0.4):
		#print 'PDF: {}'.format(pdf)
	
	p = np.random.rand()
	if(p<pdf):
		return 1

	return 0

########################

people = []
tables = [0]*TABLES
sim_time = 0
group = 0
next_group = 1
seating_group = 0
table = -1
total_customers = 0
total_infected = 0
new_infected = 0
customers_inside = 0

########################

scale = "0%" + " "*(40-4) + "100%"
print scale
sys.stdout.write("[")
sys.stdout.flush()

########################
while(sim_time <= END_TIME):
	
	x_healthy = []
	y_healthy = []
	x_infected = []
	y_infected = []
	x_contagious = []
	y_contagious = []
	
	if(sim_time in range(TIME_STEP,END_TIME+TIME_STEP,TIME_STEP)):
		#pct_done = sim_time*100.0 / END_TIME
		#print '{} % done'.format(pct_done)
		#print ("%0.2f" % pct_done) + " % done"
		sys.stdout.write("#")
		sys.stdout.flush()
	
	group_arrived = get_group_prob(sim_time)
		
	if(group_arrived == 1):
		# group of people arrive	
		num_people = get_a_group(AVG_PEOPLE_IN_GROUP,SD_PEOPLE_IN_GROUP)
		group += 1
		infected = 0
		contagious = 0
		table = -1
		dest_x = 0
		dest_y = 0
		vec_x = 0
		vec_y = 0
		vx = 0
		vy = 0
		visit_time = get_group_visit_time(AVG_VISIT_TIME,SD_VISIT_TIME)
		action = "WAITING"
		num_contagious = 0
		for _ in range(num_people):
			infected = is_infected()
			contagious = infected
			num_contagious += contagious
			p = Person(group,infected,contagious,table,dest_x,dest_y,vec_x,vec_y,vx,vy,visit_time,action)
			people.append(p)
			#p.print_stats()
			#print ""
		if(group == 6):
			if(num_contagious == 0):
				people[-1].contagious = 1
				people[-1].infected = 1
	
	for p in people:
		
		if(p.active == 0):
			continue

		if(p.action == "WAITING"):
			if(seating_group == 0):
				table = get_table() # only do this once!
				#print ">>>>>>>>>>>>>>    GOT A TABLE: ",
				#print table
			if(table != -1):
				if(p.group == next_group): # seat these fuckers
					p.action = "ARRIVING"
					p.table = table
					p.dest_x = TABLE_COORD[table][0]
					p.dest_y = TABLE_COORD[table][1]
					p.x += EPSILON*np.random.normal(0,1)
					p.y += EPSILON*np.random.normal(0,1)
					p.vec_x = p.dest_x - p.x
					p.vec_y = p.dest_y - p.y
					#string = "Vector: (" + str(p.vec_x) + "," + str(p.vec_y) + ")"
					#print string
					p.vx,p.vy = get_vx_vy(p.vec_x,p.vec_y)
					total_infected += p.infected
					customers_inside += 1
					total_customers += 1
					seating_group = 1

		if(p.elapsed_time >= p.visit_time):
			if(p.action != "LEAVING"):
				p.action = "LEAVING"
				p.dest_x = X_OUT_DOOR
				p.dest_y = Y_OUT_DOOR
				p.vec_x = p.dest_x - p.x
				p.vec_y = p.dest_y - p.y
				#string = "Vector: (" + str(p.vec_x) + "," + str(p.vec_y) + ")"
				#print string
				p.vx,p.vy = get_vx_vy(p.vec_x,p.vec_y)
				tables[p.table] = 0

		# people arrive
		if(p.action == "ARRIVING"):
			if(p.d() > EPSILON):
				p.walk()
			else:
				angle = 2.0*np.pi*np.random.rand()
				p.x = p.dest_x + TABLE_RADIUS*np.cos(angle)
				p.y = p.dest_y + TABLE_RADIUS*np.sin(angle)
				p.action = "SEATED"
		
		# people leave		
		
		if(p.action == "LEAVING"):
			if(p.d() > EPSILON):
				p.walk()
			else:
				p.x = X_OUT_DOOR
				p.y = Y_OUT_DOOR
				p.active = 0
				p.action = "LEFT"
				customers_inside -= 1
		

		#print("(%0.2f , %0.2f)" % (p.x,p.y))
		
		if(p.active == 1):
			if(p.action != "WAITING"):
				p.elapsed_time += DT
		
			if(p.infected == 0):
				x_healthy.append(p.x)
				y_healthy.append(p.y)
			if(p.infected == 1):
				x_infected.append(p.x)
				y_infected.append(p.y)
			if(p.contagious == 1):
				x_contagious.append(p.x)
				y_contagious.append(p.y)
	
	########### END PEOPLES LOOP ###############
	for p in people:
		if(p.active == 0):
			continue
		# Calculate new infections
		if(p.contagious == 1):
			if(p.action != "WAITING"):	# only contagious people inside can infect others
				for p2 in people:
					if(p.active == 0):
						continue
					if(p2.infected == 0):
						if(p2.action != "WAITING"):	# can't infect outside waiting
							infected = get_infection(p.x,p.y,p2.x,p2.y)
							p2.infected = infected
							total_infected += infected
							new_infected += infected
							#if(infected == 1):
								#print ">>>>>>> New infected!"
	########### DONE NEW INFECTIONS #################
	
	if(seating_group == 1):	# get next group for seating
		seating_group = 0
		next_group += 1
	
	if(PRINT_DATA):	
		print x_healthy
		print y_healthy
		print x_infected
		print y_infected
		print x_contagious
		print y_contagious
	
	if(MAKE_PLOT):
		for c in TABLE_COORD:
			x=c[0]
			y=c[1]
			plt.plot(x,y,'ko',ms=15)

		for c in X_TABLE_COORD:
			x=c[0]
			y=c[1]
			plt.plot(x,y,'ko',ms=15,mfc='none')

		#plt.plot(x_healthy,y_healthy,'bo',ms=4, label='Sano')
		#plt.plot(x_infected,y_infected,'ro',ms=4, label='Infectado')
		plt.plot(x_healthy,y_healthy,'bo',ms=4, label='Healthy')
		plt.plot(x_infected,y_infected,'ro',ms=4, label='Infected')
		plt.plot(x_contagious,y_contagious,'o',color='coral',ms=6, label='Contagious')
		plt.xlim(0,MAX_X)
		plt.ylim(0,MAX_Y)
		plt.legend(loc='upper right')
		#title = "Total Contagiosos: " + str(total_infected - new_infected) + "    Total Infectados: " + str(total_infected)
		title = "Total Contagious: " + str(total_infected - new_infected) + "    Total Infected: " + str(total_infected)
		plt.title(title)
		#xlabel = "Tiempo simulado: " + str(sim_time) + "s"
		xlabel = "Simulated Time: " + str(sim_time) + "s"
		plt.xlabel(xlabel)
		filename = "data/plot-" + ("%03d" % sim_time) + ".jpg"
		plt.savefig(filename)
		plt.close()
		#plt.show()


	sim_time += DT

###### END WHILE LOOP #########
sys.stdout.write("]\n")
sys.stdout.flush()
print "Total Arrivals: ",
print len(people)
print "Total Customers: ",
print total_customers
print "Contagious: ",
print str(total_infected - new_infected)
print "Newly Infected: ",
print new_infected
print "Total Infected: ",
print total_infected
print "% Infected: ",
print str(round(100.0*new_infected/total_customers,1))
t2 = time.time()
print "Exec Time: ",
print str(t2-t1)

