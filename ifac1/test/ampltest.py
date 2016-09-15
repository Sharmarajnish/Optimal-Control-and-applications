
from tools import *

## AMPL #######################################################################
def generate_ampl(N,T,options):
	mat = loadmat('data/patient28.mat')
	eps = 250000
	Ne = mat['Ne']
	
	# MODEL ----------
	ampl = 'param N = %i;\n' % N
	ampl += 'param Ne = %i;\n' % Ne
	ampl += 'param Tmax = %1.4f;\n' % T
	ampl += 'param eps = %1.4f;\n\n'% eps
	
	ampl += 'set nodes := 1..(N+1);\n'
	ampl += 'set electrodes := 1..Ne;\n\n'
	
	ampl += 'param D {nodes,nodes};\n'
	ampl += 'param times {nodes};\n'
	ampl += 'param wts {nodes};\n\n'
	
	ampl += 'param PY_0 {electrodes};\n'
	ampl += 'param IN_0 {electrodes};\n'
	ampl += 'param TC_0;\n'
	ampl += 'param RE_0;\n'
	ampl += 'param PY_T {electrodes};\n'
	ampl += 'param IN_T {electrodes};\n'
	ampl += 'param TC_T;\n'
	ampl += 'param RE_T;\n'
	
	ampl += 'param tau {1..4};\n'
	
	ampl += 'param h_PY {electrodes};\n'
	ampl += 'param h_IN {electrodes};\n'
	ampl += 'param h_TC;\n'
	ampl += 'param h_RE;\n'
	
	ampl += 'param C_PYPY {electrodes,electrodes};\n'
	ampl += 'param C_PYIN {electrodes,electrodes};\n'
	ampl += 'param C_INPY {electrodes,electrodes};\n'
	
	ampl += 'param C_RERE;\n'
	ampl += 'param C_TCRE;\n'
	ampl += 'param C_RETC;\n'
	
	ampl += 'param C_PYTC {electrodes};\n'
	ampl += 'param C_PYRE {electrodes};\n'
	ampl += 'param C_TCPY {electrodes};\n\n'
		
	ampl += 'var T >=0 <= Tmax;\n'
	ampl += 'var PY {electrodes,nodes} >= -1 <= 1;\n'
	ampl += 'var IH {electrodes,nodes} >= -1 <= 1;\n'
	ampl += 'var TC {nodes} >= -1 <= 1;\n'
	ampl += 'var RE {nodes} >= -1 <= 1;\n'
	ampl += 'var u {electrodes,nodes} >= -2 <= 2;\n\n'
	
	if 'b:binary' in options:
		ampl += 'var b {electrodes} binary;\n\n'	
	elif 'b:fixed' in options:
		ampl += 'param b{electrodes} = 1;\n\n'
	else:
		ampl += 'var b {electrodes} >= 0 <= 1;\n\n'
	
	ampl += 'maximize cost: '
	if 'minT' in options:
		ampl += ' - %f*(T/Tmax)' % (options['minT'],)
	if 'minE' in options:
		ampl += ' - %f*(sum{e in electrodes} (b[e]*(sum{t in nodes} u[e,t]^2*wts[t])) )' % (options['minE'],)
	if 'minB' in options:
		ampl += ' - %f*(sum{e in electrodes} (b[e]) )' % (options['minB'],)
	ampl += ';\n\n'
	
	fPY = '(1/(1+eps^(-PY[j,t])))'
	fIN = '(1/(1+eps^(-IH[j,t])))'
	fTC = '(1/(1+eps^(-TC[t])))'
	fRE = '(1/(1+eps^(-RE[t])))'
	
	ampl += 'subject to dynamics_PY {t in nodes, e in electrodes}:\n'
	ampl += '\ttau[1]*(h_PY[e] - PY[e,t] + (sum{j in electrodes} C_PYPY[e,j]*%s) - (sum{j in electrodes} C_INPY[e,j]*%s) + C_TCPY[e]*%s ) + u[e,t] = (2/T)*(sum{k in nodes} D[t,k]*PY[e,k]);\n' % (fPY,fIN,fTC)
	ampl += 'subject to dynamics_IN {t in nodes, e in electrodes}:\n'
	ampl += '\ttau[2]*(h_IN[e] - IH[e,t] + (sum{j in electrodes} C_PYIN[e,j]*%s)) + u[e,t] = (2/T)*(sum{k in nodes} D[t,k]*IH[e,k]);\n' % fPY
	ampl += 'subject to dynamics_TC {t in nodes}:\n'
	ampl += '\ttau[3]*(h_TC - TC[t] - C_RETC*%s + (sum{j in electrodes} C_PYTC[j]*%s) ) = (2/T)*(sum{k in nodes} D[t,k]*TC[k]);\n' % (fRE,fPY)
	ampl += 'subject to dynamics_RE {t in nodes}:\n'
	ampl += '\ttau[4]*(h_RE - RE[t] - C_RERE*%s + C_TCRE*%s + (sum{j in electrodes} C_PYRE[j]*%s) ) = (2/T)*(sum{k in nodes} D[t,k]*RE[k]);\n\n' % (fRE,fTC,fPY)
	
	ampl += 'subject to initialConditionPY {e in electrodes}: PY[e,1] = PY_0[e];\n'
	ampl += 'subject to initialConditionIN {e in electrodes}: IH[e,1] = IN_0[e];\n'
	ampl += 'subject to initialConditionTC: TC[1] = TC_0;\n'
	ampl += 'subject to initialConditionRE: RE[1] = RE_0;\n\n'
	
	if 'fixedT' in options:
		ampl += 'subject to fixedTime: T=Tmax;\n\n'
	
	ampl += 'subject to finalConditionPY {e in electrodes}: PY[e,N+1] = PY_T[e];\n'
	ampl += 'subject to finalConditionIN {e in electrodes}: IH[e,N+1] = IN_T[e];\n'
	ampl += 'subject to finalConditionTC: TC[N+1] = TC_T;\n'
	ampl += 'subject to finalConditionRE: RE[N+1] = RE_T;\n\n'
	
	# DATA ---------
	ampl += 'data;\n\n'
	
	for nm in ['PY_0','IN_0','TC_0','RE_0','PY_T','IN_T','TC_T','RE_T','tau','h_PY','h_IN','h_TC','h_RE','C_PYPY','C_INPY','C_PYIN','C_PYRE','C_PYTC','C_RERE','C_RETC','C_TCPY','C_TCRE']:
		ampl += array2ampl('param',nm,mat[nm])
		
	wts = (LGLnodes(N))[1]
	ampl += array2ampl('param','wts',wts)
	
	times = ((LGLnodes(N))[0]+1.0)*float(T)/2
	ampl += array2ampl('param','times',times)
	
	D = Dmatrix(N)
	ampl += array2ampl('param','D',D)
	
	ampl += '\nsolve;\n'
	for v in 'N;Tmax;T;eps;Ne;u;cost;PY_0;IN_0;TC_0;RE_0;PY_T;IN_T;TC_T;RE_T;tau;h_PY;h_IN;h_TC;h_RE;C_PYPY;C_INPY;C_PYIN;C_PYRE;C_PYTC;C_RERE;C_RETC;C_TCPY;C_TCRE;PY;IH;TC;RE'.split(';'):
		ampl += 'display %s;\n' % v
	
	return ampl
