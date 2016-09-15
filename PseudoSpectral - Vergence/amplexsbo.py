from tools import *

def sbo_ampl(N,T):
	#N=10
	A=1;
	#T=1.5708
	x0=0
	y0=0
	z0=1

	# MODEL ----------
	ampl = 'param N = %i;\n' % N
	ampl += 'param T = %1.4f;\n' % T
	ampl += 'param A = %1.4f;\n\n'% A
	
	
	ampl += 'param x0 = %1.4f;\n' % x0
	ampl += 'param y0 = %1.4f;\n' % y0
	ampl += 'param z0 = %1.4f;\n' % z0
	
	ampl += 'set nodes := 1..(N+1);\n\n'
	
	ampl += 'param D {nodes,nodes};\n\n'
	
	ampl += 'var x {nodes} >= -1, <= 1;\n'
	ampl += 'var y {nodes} >= -1, <= 1;\n'
	ampl += 'var z {nodes} >= -1, <= 1;\n'
	ampl += 'var u {nodes} >= -A, <= A;\n'
	ampl += 'var v {nodes} >= -A, <= A;\n\n'

	ampl += 'maximize cost: x[N+1];\n\n'
	
	
	ampl += 'subject to dynamics_x {t in nodes}:\n'
	ampl += '\tu[t]*z[t] = (2/T)*(sum{k in nodes} D[t,k]*x[k]);\n\n'
	
	ampl += 'subject to dynamics_y {t in nodes}:\n'
	ampl += '\t-v[t]*z[t]= (2/T)*(sum{k in nodes} D[t,k]*y[k]);\n\n'
	
	ampl += 'subject to dynamics_z {t in nodes}:\n'
	ampl += '\tv[t]*y[t]-u[t]*x[t]= (2/T)*(sum{k in nodes} D[t,k]*z[k]);\n\n'
	
	ampl += 'subject to initialCondition_x {t in nodes}: x[1] = x0;\n'
	ampl += 'subject to initialCondition_y {t in nodes}: y[1] = y0;\n'
	ampl += 'subject to initialCondition_z {t in nodes}: z[1] = z0;\n\n'
	
	ampl += 'subject to amplitudeBound {t in nodes}:\n'
	ampl += '\tu[t]^2+v[t]^2 <= A^2;\n\n'

	# DATA ---------
	ampl += 'data;\n\n'

	D = Dmatrix(N)
	ampl += array2ampl('param','D',D)
	
	ampl += 'solve;\n'
	ampl +='display A;\n'
	ampl +='display T;\n'
	ampl +='display N;\n'
	ampl +='display u;\n'
	ampl +='display v;\n'
	ampl +='display cost;\n'
	ampl +='display x0;\n'
	ampl +='display y0;\n'
	ampl +='display z0;\n'
	ampl +='display x;\n'
	ampl +='display y;\n'
	ampl +='display z;\n'
	return ampl
	
amplstr = sbo_ampl(10,3.1428)

f = open('sbo_ampl.txt','w')
f.write(amplstr)
f.close()