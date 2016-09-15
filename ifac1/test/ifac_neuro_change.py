from tools import *

def ifac_neuro_ampl(ampl):
	from scipy.interpolate import InterpolatedUnivariateSpline
	
	ampl_dict = ampl2dictionary(ampl)
	opt = optsolution()
	
	N = int(ampl_dict['N'][0])
	Tmax = ampl_dict['T'][0]
	eps = int(ampl_dict['eps'][0]
	
	u = ampl_dict['u']
	x1 = ampl_dict['x'][0]
	x2 = ampl_dict['x[nIN]']
	x3 = ampl_dict['x[nTC]']
	x4 = ampl_dict['x[nRE]']
	
	opt.states_nodes= vstack[(x1,x2,x3,x4)]
	opt.nodes = ((LGLnodes(int(opt.p['N'])))[0]+1)*opt.p['T']/2
	opt.controls_nodes = opt.p['u'].T
	opt.t = linspace(opt.nodes[0],opt.nodes[-1],1001)

	U = InterpolatedUnivariateSpline(opt.nodes,u)(opt.t)
	opt.states = neuro_evolution(opt.t,opt.controls)
	plot_result(opt,dosave)


def neuro_evolution(t,u,x0):
		def neuro_dynmaic(x,t,i_in,u_in)
			u = interp(t,t_in,u_in)
			dx = zeros(x.shape)
		
			fPY = 1/(1+p['eps']**(-x[xPY]))
			fIN = 1/(1+p['eps']**(-x[xIN]))
			fTC = 1/(1+p['eps']**(-x[xTC]))
			fRE = 1/(1+p['eps']**(-x[xRE]))
			
			dx[nPY] = tau[nPY]*(h[nPY] - x[nPY] + fPY*C[1] - fIN*C[3] + fTC*C[9] + u)
			dx[nIN] = tau[nIN]*(h[nIN] - x[nIN] + fPY*C[2] + u)
			dx[nTC] = tau[nTC]*(h[nTC] - fTC*C[6] + fPY*C[7]
			dx[nRE] = tau[nRE]*(h[nRE] - fRE*C[4] - fTC*C[5] + fPY*C[8])
			
			return dx
	
		X = transpose( odeint(neuro_dynamics,x0,t, args = (t,U)))
		return X

f = open('ifac_neuro_output.txt','r')
amploutput = f.read()
f.close()
ifac_neuro_ampl(ifac_neuro_output)