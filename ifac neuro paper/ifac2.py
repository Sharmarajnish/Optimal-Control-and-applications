from tools import *

def ifac_neuro_ampl(ampl):
	from scipy.interpolate import InterpolatedUnivariateSpline
	
	ampl_dict = ampl2dictionary(ampl)
	opt = optsolution()
	
	N = int(ampl_dict['N'][0])
	Tmax = ampl_dict['T'][0]
	eps = int(ampl_dict['eps'][0]
	
	x1 = ampl_dict['x1']
	x2 = ampl_dict['x2']
	x3 = ampl_dict['x3']
	x4 = ampl_dict['x4']

	
	
	
	opt.controls_nodes = vstack([u])
	opt.nodes = ((LGLnodes(N))[0]+1)*T/2
	opt.t = linspace(0,T,1001)
	U = InterpolatedUnivariateSpline(opt.nodes,u)(opt.t)
	
	opt.states = bloch_evolution2(opt.t,U,V,opt.param['x0'])
	opt.performance = bloch_plot(opt)
	return opt
	
def neuro_evolution(t,u,x0):
	def neuro_dynmaic(x,t,i_in,u_in)
		u = interp(t,t_in,u_in)
		dx = zeros(x.shape)
		
		dx[0] = 
		dx[1] = 
		dx[2] =
		dx[3] =
		
		return dx
	
	X = transpose( odeint(neuro_dynamics,x0,t, args = (t,U)))
	return X
	
def neuro_plot(opt)

	plt.figure(figsize=(8,4))
	