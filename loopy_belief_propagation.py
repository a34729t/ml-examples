from __future__ import division # fix python division
import sys
from math import log
import re
# Implementation of Loopy Belief Propagation

def uniquifyList(seq, idfun=None): # TODO: move this to utilities file
	# From http://www.peterbe.com/plog/uniqifiers-benchmark
    # order preserving 
    if idfun is None: 
        def idfun(x): return x 
    seen = {} 
    result = [] 
    for item in seq: 
        marker = idfun(item) 
        if marker in seen: continue 
        seen[marker] = 1 
        result.append(item) 
    return result

class distribution:
	# Implement a CPD as a list (see pg 359 Kohler and Friedman)
	# get the variable names by doing factor.values.keys() or factor.stride.keys()
	def __init__(self,variableList,cardinality,cpdList):
		self.variables = variableList # variables in the CPD (ie A,B,C)
		self.cardinality = cardinality # the number of values for a specific variable (basically array offset since we are flattening out the CPD)
		self.cpdList = cpdList
		self._stride = {} # numer of steps to see next assignment to a variable
		
		lastStride = 1;
		for var in self.variables:
			self._stride[var] = lastStride
			lastStride *= self.cardinality[var]
	
	def __len__(self):
		return len(self.cpdList)
	
	def stride(self,index):
		try:
			x = self._stride[index]
			return x
		except:
			return 0
			
	def __getitem__(self,index):
		return self.cpdList[index]
		
	def __setitem__(self,index, value):
		self.cpdList[index] = value
		
	def printOut(self):
		print "todo"
		
	def marginal(self,marginalVar): # marginalize over a given variable
		
		stride = self.stride(marginalVar) # add error catcher?
		card = self.cardinality[marginalVar] # add error catcher?
		newCpd = [0] * card # cpd is initialized as a list of zeroes of the proper length
		
		strideCount = 0
		cardCount = 0
		for i in range(len(self)):
			
			if strideCount < stride:
				strideCount += 1
			else:
				strideCount = 1
				cardCount += 1
				
			if cardCount >= card:
				cardCount = 0
				
			newCpd[cardCount] += self[i]	
					
		# Construct a new factor that only goes over the variable we marginalized across
		cardinality = {}
		cardinality[marginalVar] = card
		newDist = distribution([marginalVar],cardinality,newCpd)
		newDist.normalize()
		return newDist
		
	def __mul__(self,phi2):
		# Multiply two factors and return a new one
		# Algorithm 10.A.1 from Kohler and Friedman 2009
		# phi1 is current factor
		# phi2 is multiplicand factor
		phi1 = self
		
		# Check distributions are compatible (ie same cardinality for same variables) and get cardinality and the variables in the new factor
		card = {}
		variables = []
		for var1 in phi1.variables:
			for var2 in phi2.variables:
				if var1 == var2:					
					c1 = phi1.cardinality[var1]
					c2 = phi1.cardinality[var2]

					if c1 != c2:
						print "ERROR: factor multiplication, cardinality mismatch:"," ",var1,"(",c1,") != ",var2,"(",c2,")"
						
					variables.extend(var1) # only need to do 1
					card[var1] = phi1.cardinality[var1]
					
				else:
					variables.extend(var1)
					variables.extend(var2)
					card[var1] = phi1.cardinality[var1]
					card[var2] = phi2.cardinality[var2]
		variables = uniquifyList(variables) # get rid of duplicates in the variable list
		
		lenNewCpd = 1
		for c in card.keys():
			lenNewCpd *= card[c]
		newCpd = [0] * lenNewCpd # cpd is initialized as a list of zeroes of the proper length

		# start of Algorithm 10.A.1 #################
		j = 0
		k = 0
		
		assignment = {}
		for l in card.keys():
			assignment[l] = 0
		
		for i in range(lenNewCpd):
			#print "<i,j,k>=",i,",",j,",",k,",   ",assignment
			newCpd[i]= phi1[j] * phi2[k]
			for l in variables: # Advances to the next assignment to the variables in psi and calculates the indexes into each factor array
				assignment[l] = assignment[l]+1
				if assignment[l] == card[l] or assignment[l] % card[l] == 0:
					j = j - (card[l]-1)*phi1.stride(l)
					k = k - (card[l]-1)*phi2.stride(l)
				else:
					j = j + phi1.stride(l)
					k = k + phi2.stride(l)
					break
		# end of Algorithm 10.A.1 #################
		
		newDist = distribution(variables,card,newCpd)
		newDist.normalize()
		return newDist
		
	def normalize(self):
		# sum across all elements in array
		div = sum(self)
		for i in range(len(self)):
				self[i] /= div

	def setUniform(self):
		self = [1] * len(self)
	
	def setZeros(self):
		self = [0] * len(self)
		
	def kld(self,other): # K-L Divergence Score
		if other == None: # starting case
			return 100
		else:
			ret = 0.0
			self.normalize()
			other.normalize();
			for i in range(len(self)):
				if self[i] > 0.0:
					ret += self[i] * log(self[i]/other[i])
			return ret;
		
class factorGraph:
	def __init__(self):
		x = 1

class Message:
	def __init__(self,fromId,toId,dist):
		self.fromId = fromId
		self.toId = toId
		self.dist = dist

class Factor:
	def __init__(self, id, name, dist):
		self.id = id
		self.name = name
		self.dist = dist
		self.links = []
		self.incoming = []

	def link(self,links):
		self.links.extend(links)	

	def send(self):
		for link in self.links:
			#print "Factor: Generating outgoing messages:",self.name,"->",link.name
			dist = self.dist
			for inMsg in self.incoming:
				if inMsg.fromId != link.id:
					#print "\t\tmessage from",inMsg.fromId
					#print "\t\tmsg:",inMsg.dist.cpdList
					#print "\t\tdist:",dist.cpdList
					dist *= inMsg.dist
					#print "\t\tdist*msg:",dist.cpdList
			m = dist.marginal(link.name)
			#print "\tmarginal("+link.name+"):",m.cpdList
			message = Message(self.id,link.id,dist.marginal(link.name))
			link.incoming.append(message)
		self.incoming = [] # clean up

class Variable:
	def __init__(self, id, name):
		self.id = id
		self.name = name
		self.links = []
		self.incoming = []

	def link(self,links):
		self.links.extend(links)	

	def initialize(self):
		dist = distribution([self.name],{self.name:2}, [1,1])
		for link in self.links:
			message = Message(self.id,link.id,dist)
			link.incoming.append(message)
	
	def send(self):
		flag = True
		dist = None
		for link in self.links:
			#print "Variable: Generating outgoing messages:",self.name,"->",link.name
			curr = distribution([self.name],{self.name:2}, [1,1])
			for inMsg in self.incoming:
				if inMsg.fromId != link.id:
					curr *= inMsg.dist
					#print "\t\tmessage from",inMsg.fromId
					#print "\t\tmsg:",inMsg.dist.cpdList
					#print "\t\tcurr*msg:",curr.cpdList

			#print "\t\tcurr:",curr.cpdList
			message = Message(self.id,link.id,curr)
			link.incoming.append(message)
			
			if flag: # capture the first outgoing message so we can return it (see below, need for k-l divergence)
				flag = False
				dist = curr

		dist *= self.incoming[0].dist
		self.incoming = [] # clean up
		return dist
	
class Network:
	def __init__(self):
		self.name = ""
		self.variables = []
		self.factors = []
	
def loadBIF(filename):

	# for sprinkler.bif:

	# C = Variable(1, "c")
	# S = Variable(2, "s")
	# R = Variable(3, "r")
	# W = Variable(4, "w")

	# c = Factor(5, "c", distribution(["c"],{"c":2}, [0.5,0.5]))
	# s = Factor(6, "s", distribution(["s","c"],{"s":2, "c":2}, [0.5,0.9,0.5,0.1]))
	# r = Factor(7, "r", distribution(["r","c"],{"r":2, "c":2}, [0.8,0.2,0.2,0.8]))
	# w = Factor(8, "w", distribution(["w","s","r"],{"w":2,"r":2, "s":2}, [0.99,0.01,0.1,0.9,0.1,0.9,0.01,0.99]))
	
	# C.link([c,s,r])
	# S.link([s,w])
	# R.link([r,w])
	# W.link([w])

	# c.link([C])
	# s.link([C,S])
	# r.link([C,R])
	# w.link([S,R,W])
	
	# variables = [C,S,R,W]
	# factors = [c,s,r,w]

	f=open(filename, 'r')
	
	network = Network()
	
	
	p_var = re.compile('^variable\s+(\S*)\s+\{')
	p_net = re.compile('^network\s+(\S*)\s+\{')
	
	file = f.readlines()
	for i in range(len(file)):
		line = file[i]
		var = p_var.match(line)
		net = p_net.match(line)
		
		if var:
			varName = var.group(1)
			print "adding Variable:",varName
		
		if net:
			netName = net.group(1)
			network.name = netName
			print "Network name:",netName
	
def run():
	print "Starting belief propagation ..."

	# the way my arrays are set out is first column changes first:
	# for example for [a1,a2,a3,b1,b2,c1,c2]
	# a1 b1 c1
	# a2 b1 c1
	# a3 b1 c1
	# a1 b2 c1
	# a2 b2 c1
	# a3 b2 c1
	# a1 b1 c2
	# .....
	#f1 = distribution(["a","b"],{"a":3,"b":2},[0.5,0.1,0.3,0.8,0,0.9])
	#f2 = distribution(["b","c"],{"b":2,"c":2},[0.5,0.1,0.7,0.2])
	#f3 = distribution(["a","b","c"],{"a":3,"b":2,"c":2},[0.25,0.05,0.15, 0.08,0,0.09, 0.35,0.07,0.21, 0.16,0,0.18])
	
	# Hand coded example, for now (read BIF files later on)
	
	#varA = Variable(1, "a")
	#varB = Variable(2, "b")
	#facA = Factor(3, "a", distribution(["a"],{"a":2}, [0.5,0.5]))
	#facB = Factor(4, "b", distribution(["b","a"],{"b":2, "a":2}, [0.5,0.4,0.3,0.7]))
	
	#varA.links.append(facA)
	#varA.links.append(facB)
	#varB.links.append(facB)
	#facA.links.append(varA)
	#facB.links.append(varA)
	#facB.links.append(varB)

	#variables = [varA, varB]
	#factors = [facA, facB]	

	#bn = loadBIF('./sprinkler.bif')
	#sys.exit();
	
	C = Variable(1, "c")
	S = Variable(2, "s")
	R = Variable(3, "r")
	W = Variable(4, "w")

	c = Factor(5, "c", distribution(["c"],{"c":2}, [0.5,0.5]))
	s = Factor(6, "s", distribution(["s","c"],{"s":2, "c":2}, [0.5,0.9,0.5,0.1]))
	r = Factor(7, "r", distribution(["r","c"],{"r":2, "c":2}, [0.8,0.2,0.2,0.8]))
	w = Factor(8, "w", distribution(["w","s","r"],{"w":2,"r":2, "s":2}, [0.99,0.01,0.1,0.9,0.1,0.9,0.01,0.99]))
	
	C.link([c,s,r])
	S.link([s,w])
	R.link([r,w])
	W.link([w])

	c.link([C])
	s.link([C,S])
	r.link([C,R])
	w.link([S,R,W])
	
	variables = [C,S,R,W]
	factors = [c,s,r,w]
	
	# Initialize all messages to [1,1]
	for var in variables:
		var.initialize()
	
	iter = 0
	MAX_ITER = 100
	MIN_CHANGE = 0.000001

	currDists = [None] * len(variables)	

	while True:
		biggestChange = 0.0
		print "i=",iter

		# Print messages in each factor
		print "\t","Factors- incoming messages"
		for factor in factors:
			for msg in factor.incoming:
				print "\t\t",factor.name+": ",msg.fromId,"->",msg.toId," = ",msg.dist.cpdList

		# For each factor, receive and send messages
		for factor in factors:
			factor.send()
		
		# Print messages in each variable
		print "\t","Variables- incoming messages"
		for var in variables:
			for msg in var.incoming:
				print "\t\t",var.name+": ",msg.fromId,"->",msg.toId," = ",msg.dist.cpdList

		# For each variable, receive and send messages
		for i in range(len(variables)):
			newDist = variables[i].send()

			# Track change in distributions
			kld = newDist.kld(currDists[i])
			if kld > biggestChange:
				biggestChange = kld
			currDists[i] = newDist
	
		

		print "\t","Biggest change =",biggestChange
	
		iter += 1
		if iter == MAX_ITER or biggestChange < MIN_CHANGE:
			break


if __name__ == "__main__":
    run()
	
