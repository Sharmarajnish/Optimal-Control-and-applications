from numpy import *
from scipy.special import *

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
	
D = Dmatrix(10)
print array2ampl('param','D',D)