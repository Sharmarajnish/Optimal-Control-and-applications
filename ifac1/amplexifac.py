from tools import *

def neuron_ampl(N,T,x0,options):
	C = array([1.8,4.0,1.5,0.2,10.0,1.5,3.0,3.0,1.0])
	tau = array([1.0,1.25,0.1,0.1])
	h = array([-0.35,-3.4,-2.0,-5.0])
	eps = 250000

	xT = array([0.1695,0.1639,-0.0914,0.0031])

	# MODEL ----------
	ampl = 'param N = %i;\n' % N
	ampl += 'param Tmax = %1.4f;\n' % T
	ampl += 'param eps = %1.4f;\n\n'% eps
	
	ampl += 'param nPY = 1;\n'
	ampl += 'param nIN = 2;\n'
	ampl += 'param nTC = 3;\n'
	ampl += 'param nRE = 4;\n'
	
	ampl += 'set nodes := 1..(N+1);\n'
	ampl += 'set states := 1..4;\n\n'
	
	ampl += 'param D {nodes,nodes};\n'
	ampl += 'param times {nodes};\n'
	ampl += 'param wts {nodes};\n\n'
	
	ampl += 'param x0 {states};\n'
	ampl += 'param xT {states};\n'
	ampl += 'param tau {states};\n'
	ampl += 'param h {states};\n'
	ampl += 'param C {1..9};\n\n'

	ampl += 'var T >=0 <= Tmax;\n'
	ampl += 'var x {states,nodes} >= -2 <= 2;\n'
	
	ampl += 'var u {nodes} >= -2 <= 2;\n'
	
	ampl += 'maximize cost: '
	if 'minT' in options:
		ampl += ' - %f*(T/Tmax)' % (options['minT'],)
	if 'minE' in options:
		ampl += ' - %f*(sum{t in nodes} (u[t]^2)*wts[t])' % (options['minE'],)
	if 'minErr' in options:
		ampl += ' - (sum{j in states} (x[j,N+1] - xT[j])^2)'
	ampl += ';\n\n'
	
	fPY = '(1/(1+eps^(-x[nPY,t])))'
	fIN = '(1/(1+eps^(-x[nIN,t])))'
	fTC = '(1/(1+eps^(-x[nTC,t])))'
	fRE = '(1/(1+eps^(-x[nRE,t])))'
	
	ampl += 'subject to dynamics_PY {t in nodes}:\n'
	ampl += '\ttau[nPY]*(h[nPY] - x[nPY,t] + C[1]*%s - C[3]*%s + C[9]*%s) + u[t] = (2/T)*(sum{k in nodes} D[t,k]*x[nPY,k]);\n' % (fPY,fIN,fTC)
	ampl += 'subject to dynamics_IN {t in nodes}:\n'
	ampl += '\ttau[nIN]*(h[nIN] - x[nIN,t] + C[2]*%s) + u[t] = (2/T)*(sum{k in nodes} D[t,k]*x[nIN,k]);\n' % fPY
	ampl += 'subject to dynamics_TC {t in nodes}:\n'
	ampl += '\ttau[nTC]*(h[nTC] - x[nTC,t] - C[6]*%s + C[7]*%s) = (2/T)*(sum{k in nodes} D[t,k]*x[nTC,k]);\n' % (fRE,fPY)
	ampl += 'subject to dynamics_RE {t in nodes}:\n'
	ampl += '\ttau[nRE]*(h[nRE] - x[nRE,t] - C[4]*%s + C[5]*%s + C[8]*%s) = (2/T)*(sum{k in nodes} D[t,k]*x[nRE,k]);\n' % (fRE,fTC,fPY)

	ampl += 'subject to initialCondition {j in states}: x[j,1] = x0[j];\n\n'
	
	if 'fixedT' in options:
		ampl += 'subject to fixedTime: T=Tmax;\n\n'
	
	#ampl += 'subject to terminalCondition {t in N-10..N+1}: (sum{j in states} (x[j,t] - xT[j])^2) <= 0.0001;\n\n'
	ampl += 'subject to terminalCondition: (sum{j in 1..4} (x[j,N+1] - xT[j])^2) <= 0.0001;\n\n'
	ampl += 'subject to zerou0: u[1] = 0;\n\n'
	ampl += 'subject to zerouT: u[N+1] = 0;\n\n'

	# DATA ---------
	ampl += 'data;\n\n'
	
	ampl += array2ampl('param','C',C)
	ampl += array2ampl('param','h',h)
	ampl += array2ampl('param','tau',tau)
	ampl += array2ampl('param','x0',x0)
	ampl += array2ampl('param','xT',xT)
	
	wts = (LGLnodes(N))[1]
	ampl += array2ampl('param','wts',wts)
	
	times = ((LGLnodes(N))[0]+1.0)*float(T)/2
	ampl += array2ampl('param','times',times)
	
	D = Dmatrix(N)
	ampl += array2ampl('param','D',D)

	ampl += 'solve;\n'
	for v in 'N;Tmax;x0;xT;T;eps;C;h;tau;x;u;cost'.split(';'):
		ampl += 'display %s;\n' % v

	return ampl
	
opts = dict()
opts['minT'] = 0.1
opts['minE'] = 1
opts['fixedT'] = True
x0 = array([0,0,0,0])
amplstr = neuron_ampl(24,4,x0,opts)

f=open('ifac_ampl.txt','w')
f.write(amplstr)
f.close()