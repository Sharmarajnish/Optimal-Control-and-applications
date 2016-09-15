from tools import *

## AMPL #######################################################################


## PROCESS ####################################################################

def bloch_process_ampl(ampl):
	from scipy.interpolate import lagrange
	from scipy.interpolate import InterpolatedUnivariateSpline
	
	ampl_dict = ampl2dictionary(ampl)
	opt = OptSolution()
	
	A = ampl_dict['A'][0]
	N = int(ampl_dict['N'][0])
	T = ampl_dict['T'][0]
	
	opt.param['x0'] = array([ampl_dict['x0'][0],ampl_dict['y0'][0],ampl_dict['z0'][0]])
	u = ampl_dict['u']
	v = ampl_dict['v']
	
	x = ampl_dict['x']
	y = ampl_dict['y']
	z = ampl_dict['z']
	opt.states_nodes = vstack([x,y,z])
	
	opt.controls_nodes = vstack([u,v])
	opt.nodes = ((LGLnodes(N))[0]+1)*T/2
	
	opt.t = linspace(0,T,1001)
	if N <= 20:
		U = lagrange(opt.nodes,u)(opt.t)
		V = lagrange(opt.nodes,v)(opt.t)
	else:
		U = InterpolatedUnivariateSpline(opt.nodes,u)(opt.t)
		V = InterpolatedUnivariateSpline(opt.nodes,v)(opt.t)
	opt.controls = clip_amplitude(vstack([U,V]),A)
	
	opt.states = bloch_evolution2(opt.t,U,V,opt.param['x0'])
	opt.performance = bloch_plot(opt)
	return opt

def bloch_evolution(t,U,V,x0):
	X = zeros((3,len(t)))
	X[:,0] = x0
	
	for k in arange(1,len(t)):
		dt = t[k]-t[k-1]
		u = U[k]
		v = V[k]
		
		A = array([[0,0,u], [0,0,-v], [-u,v,0]])
		X[:,k] = dot(expm(A*dt),X[:,k-1])
	return X

def bloch_evolution2(t,U,V,x0):
	def bloch_dynamics(x,t,t_in,u_in,v_in):
		u = interp(t,t_in,u_in)
		v = interp(t,t_in,v_in)
		dx = zeros(x.shape)
		
		dx[0] = u*x[2]
		dx[1] = -v*x[2]
		dx[2] = -u*x[0] + v*x[1]
		return dx
	
	X = transpose( odeint(bloch_dynamics,x0,t, args = (t,U,V)) )
	return X

def bloch_plot(opt):
	
	plt.figure(figsize=(8,4))
	plt.plot(opt.t/2/pi,opt.controls[0,:],'k',label='$u(t)$',linewidth=2)
	plt.hold(True)
	plt.plot(opt.t/2/pi,opt.controls[1,:],'k--',label='$v(t)$',linewidth=2)
	control_amp = sqrt(opt.controls[0,:]**2+opt.controls[1,:]**2)
	plt.hold(False)
	plt.xlabel('$t$')
	#plt.ylabel('controls')
	plt.axis([0,opt.t[-1,]/2/pi,opt.controls.min(),control_amp.max()])
	plt.grid(True)
	plt.legend()
	plt.show()
	
	plt.figure()
	plt.plot(opt.t/2/pi,opt.states[0,:],'r')
	plt.hold(True)
	plt.plot(opt.t/2/pi,opt.states[1,:],'g')
	plt.plot(opt.t/2/pi,opt.states[2,:],'b')
	plt.plot(opt.nodes/2/pi,opt.states_nodes[0,:],'ro')
	plt.plot(opt.nodes/2/pi,opt.states_nodes[1,:],'go')
	plt.plot(opt.nodes/2/pi,opt.states_nodes[2,:],'bo')
	plt.xlabel('t')
	#plt.ylabel('states')
	plt.grid(True)
	plt.hold(False)
	plt.axis([0,opt.t[-1]/2/pi,opt.states.min(),opt.states.max()])
	plt.show()
	
	return average

## RUN ########################################################################
f = open('amploutput','r')
amploutput = f.read()
f.close()
bloch_process_ampl(amploutput)