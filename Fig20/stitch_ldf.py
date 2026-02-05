import numpy as np
from scipy.optimize import minimize

files_ldf = [ './Data/Batch1/ldf_3walker.txt', './Data/Batch2/ldf_3walker.txt', './Data/Batch3/ldf_3walker.txt' ]
####################
### COLLECT DATA ###
####################
nldf = len(files_ldf)
ldf_list = [ ]
q_list = [ ]
for filename in files_ldf:
   data = np.genfromtxt(filename)
   q_list += [ data[:,0] ]
   ldf_list += [ data[:,1] ]
######################################
### CREATE A SINGLE ABSCISSA ARRAY ###
######################################
dq = q_list[0][1] - q_list[0][0]
qmin = np.inf
qmax = - np.inf
for i in range(nldf):
    qmin = min(qmin, np.amin(q_list[i]))
    qmax = max(qmax, np.amax(q_list[i]))
q = np.arange(qmin, qmax + dq, dq)
nq = len(q)
#####################################################
### COLLAPSE THE DIFFERENT LDFS BY VERTICAL SHIFT ###
#####################################################
ldf = np.zeros(nq)
computed = np.zeros(nq)
for i in range(1, nldf):
    jmin = np.argmin(np.absolute(q_list[i - 1] - q_list[i][0]))
    jmax = np.argmin(np.absolute(q_list[i - 1] - q_list[i][-1]))
    def f(a):       
        return np.sum(( ldf_list[i][:jmax - jmin] + a[0] - ldf_list[i - 1][jmin:jmax] ) ** 2)
    res = minimize(f, np.array([0.]))
    ldf_list[i] += res.x
##################################################################
### AVERAGE THE DIFFERENT LDF FOR THE MULTI-VALUATED ABSCISSAE ###
##################################################################
for j in range(nq):
    ldf_values = [ ]
    for i in range(nldf):
        jmin = np.argmin(np.absolute(q[j] - q_list[i]))
        if np.absolute(q[j] - q_list[i][jmin]) < 0.99 * dq:
            ldf_values += [ ldf_list[i][jmin] ]
    if len(ldf_values) > 0:
        computed[j] = 1.
        ldf[j] = np.mean(ldf_values)
np.savetxt('./Data/ldf_3walker_full.txt', np.vstack((q[computed > 0.], ldf[computed > 0.])).transpose())