import matplotlib.font_manager
import matplotlib.pyplot as plt
import matplotlib as mpl
from numpy import array
from numpy import pi
from numpy import sqrt
from numpy import meshgrid
from numpy import arange
from numpy import linspace
from numpy import eye
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

import string, sys, os, shutil, time
from numpy import *
from scipy.special import *
from scipy.linalg import *
from scipy.integrate import odeint
from scipy.io import savemat
from scipy.io import loadmat
from numpy.random import *

mpl.rcParams['font.size'] = 20
mpl.rcParams['legend.fontsize'] = 20
mpl.rcParams['axes.labelsize'] = 20

from copy import copy,deepcopy
import pickle
import sys

## SERVER COMMUNICATION #######################################################

def submit_ampl(model,identifier,description,watch_variables,objective='',solver_options='',solver_name='knitro'):
	import httplib, urllib
	
	params = urllib.urlencode({'model': model, 'identifier': identifier, 'description': description, 'solver_name': solver_name, 'solver_options': solver_options, 'watch_variables': watch_variables, 'email':'jruths@gmail.com', 'objective':objective})
	response = urllib.urlopen('http://amlmax.seas.wustl.edu/submit/', params)
	
	return response.read()

def result_ampl(job_id):
	import httplib, urllib

	params = urllib.urlencode({'job_id': job_id})
	response = urllib.urlopen('http://amlmax.seas.wustl.edu/result/', params)

	return response.read()
	
def job_description(job_id):
	import httplib, urllib

	#params = urllib.urlencode({'job_id': job_id})
	response = urllib.urlopen('http://amlmax.seas.wustl.edu/description/%s' % job_id)

	return response.read()
	
## AMPL #######################################################################

def array2ampl(pv,varname,a,endit=True):
	if a.ndim == 1: # one dimensional
		s= '%s %s :=\n' % (pv,varname)
		d, = a.shape
		for i in arange(d):
			s=s+ '  %i\t%f\n' % (i+1,a[i])
	elif a.ndim == 2: # two dimensional
		s = '%s %s' % (pv,varname)
		s=s+ array2D2ampl(a)
	elif a.ndim == 3: # three dimensional
		s= '%s %s ' % (pv,varname)
		d = a.shape
		d1 = d.index(min(d))
		stars = ['*','*','*']
		stars[d1] = '%i'
		stars = '[' + ','.join(stars) + ']'
		
		for i in arange(min(d)):
			s=s+ stars % (i+1)
			if d1 == 0:
				s=s+ array2D2ampl(a[i,:,:])
			elif d1 == 1:
				s=s+ array2D2ampl(a[:,i,:])
			else:
				s=s+ array2D2ampl(a[:,:,i])
	if endit: s=s+';'
	s=s+ '\n'
	return s

def array2D2ampl(a):
	d = a.shape
	d1 = d.index(min(d))
	d2 = d.index(max(d))
	
	# Decide if need transpose
	if d2 > d1: 
		s = ' (tr)'
		a = a.transpose()
	else: 
		s = ''
	s=s+ ' :\n'
	
	# Print first row of indexes for shorter dim
	for i in arange(min(d)):
		s=s+ '\t%i' % (i+1,)
	s=s+ ' :=\n'
	
	# Print out 2D array
	for i in arange(max(d)):
		s=s+ '\t%i' % (i+1,)
		for j in arange(min(d)):
			s=s+ '\t%f' % (a[i,j])
		s=s+ '\n'
		
	return s
	
def _read_one_dim(lines,c,contents):
	while True:
		line = lines[c]
		line = line.strip()
		c += 1

		if line == ';':
			break
		else:
			data = line.split()
			for i in range(0,len(data),2):
				idx = int(data[i])
				if len(contents) < idx:
					contents.extend([0]*(idx - len(contents)))
				contents[idx-1] = float(data[i+1])
	return c

def _read_two_dim(lines,c,contents):
	while True:
		line = lines[c]
		line = line.strip()
		c += 1

		if line == ';':
			break
		elif len(line) == 0 or line.startswith(':'):
			continue
		else:
			data = line.split()
			idx = int(data[0])
			if len(contents) < idx:
				contents.append([])

			contents[idx-1].extend([float(d) for d in data[1:]])
	return c

def _read_two_dim_idx(lines,c,contents):
	num_rows = 1
	num_cols = 1
	while True:
		line = lines[c]
		line = line.strip()
		c += 1

		if line == ';':
			break
		elif len(line) == 0:
			continue
		else:
			data = line.split()
			if(int(data[0]) > num_rows): num_rows = int(data[0])
			if(int(data[1]) > num_cols): num_cols = int(data[1])
			contents.append(float(data[2]))
	contents = array(contents).reshape((num_rows,num_cols))
	return (c,contents)

def ampl2dictionary(s):

	variables = {}

	lines = s.split('\n')
	c = 0
	while c < len(lines):
		line = lines[c]
		line = line.strip()

		# skip empty and comment lines
		if line.startswith('#') or len(line) == 0:
			c += 1
			continue
		else:
			data = line.split()
			varname = data[0]
			dim_str = data[1]
			contents = []
			#variables[varname] = contents

			if dim_str == '[*]':
				c = _read_one_dim(lines,c+1,contents)
			elif dim_str == '[*,*]':
				c = _read_two_dim(lines,c+1,contents)
			elif dim_str == ':=':
				(c,contents) = _read_two_dim_idx(lines,c+1,contents)
			elif dim_str == '=':
				contents.append(float(data[2]))
				c += 1
			else:
				print line
				raise Exception, 'Unrecognized variable dimension'
			variables[varname] = contents

	return variables
	
## PSEUDOSPECTAL ##############################################################

def Dmatrix(N):
	t, w = LGLnodes(N)
	L_t = legendre(N)(t)
	L_t[0]=(-1)**N
	L_t[-1]=1
	D = zeros((N+1,N+1))
	for i in arange(N+1):
		for j in arange(N+1):
			if i != j:
				D[i,j] = L_t[i]/L_t[j]/( t[i]-t[j] )
			else:
				if i == 0:
					D[i,j] = -N*(N+1.)/4
				elif i == N:
					D[i,j] = N*(N+1.)/4
				else:
					D[i,j] = 0
	return D
	
def LGLnodes(N):
	# Start by guessing the Chebyshev-Gauss-Lobbatoo nodes
	x = cos(pi*arange(N+1)/N)
	# Legendre Vandermonde matrix
	P = zeros((N+1,N+1))
	# Compute P_(N) using the recursion relation
	# Compute its first and second derivatives and 
	# update x using the Newton-Raphson method.
	xold = 2

	while max(abs(x-xold)) > 10**(-10):
		xold = x
		P[:,0] = ones(N+1)
		P[:,1] = x

		for k in arange(2,N+1):
			P[:,k] = ( (2*k-1)*x*P[:,k-1] - (k-1)*P[:,k-2] )/k;

		x = xold - ( x*P[:,N]-P[:,N-1] )/( (N+1)*P[:,N] )
		w = 2/( (N+1)*N*(P[:,N]**2) )

	x.sort()
	return x,w

## USEFUL TOOLS ###############################################################

class OptSolution:
	def __init__(self):
		self.param = {}
	
	def getDict(self):
		ret = self.param
		ret['t']=self.t
		ret['controls']=self.controls
		ret['states']=self.states
		return ret
		
	def __str__(self):
		s = ''
		for p,v in self.param.iteritems():
			s=s+ '%s:\t\t%s\n' % (p,str(v))
		return s

def clip_vector(v,maximum):
	for i in range(len(v)):
		if v[i] > maximum:
			v[i] = maximum
		if v[i] < -maximum:
			v[i] = -maximum
	return v

def clip_amplitude(u,A):
	a = sqrt(u[0,:]**2+u[1,:]**2)
	phi = arctan2(u[1,:],u[0,:])
	a = clip_vector(a,A)
	u[0,:] = a*cos(phi)
	u[1,:] = a*sin(phi)
	return u

def write2file(fname,s):
	f = open(fname,'w')
	f.write(s)
	f.close()
	
def save_queue(q,filename):
	f = open(filename,'w')
	pickle.dump(q, f)
	f.close()
	
def load_queue(filename):
	f = open(filename,'r')
	queue = pickle.load(f)
	f.close()
	if type(queue) is list:
		return queue
	else:
		print 'creating queue'
		return []