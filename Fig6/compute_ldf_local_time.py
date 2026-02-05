import numpy as np
from scipy.optimize import root_scalar, root

pc = 2.
qc = 1.5

def tanhc(u):
    if u != 0:
        return np.tanh(u) / u
    else:
        return 1.
        
def dres2(u,p):
    return tanhc(u) - 2. / p
    
def uc(p):
    if p > pc:
        tol = 1e-10
        umin = tol
        umax = 0.5 * p
        up = root_scalar(lambda u: dres2(u,p), bracket = (umin, umax), x0 = umax, xtol = tol)
        return up.root
    else:
        return 0.
x = np.linspace(0., 3.1, 20000)
x_below_qc = x[x < qc]
ldf_below_qc = pc * x_below_qc
x_above_qc = x[x >= qc]
ldf_above_qc = np.zeros(len(x_above_qc))
warning = np.zeros(len(x_above_qc))
p = np.linspace(0., 100., 20000)
mu = np.array([uc(pp) ** 2 for pp in p])
p = p[np.isnan(mu) == False]
mu = mu[np.isnan(mu) == False]
np.savetxt('./Analytics/scgf_local_time_theory.txt', np.vstack((p, mu)).transpose())
for i, xx in enumerate(x_above_qc):
    j = np.argmax(p * xx - mu)
    ldf_above_qc[i] = p[j] * xx - mu[j]
    if j == len(p) - 1:
        warning[i] = 1.
with open('./Analytics/ldf_local_time_theory.txt', 'w') as f:
    np.savetxt(f, np.vstack((x_below_qc, ldf_below_qc)).T, fmt='%.10f')
    f.write('\n\n')
    np.savetxt(f, np.vstack((x_above_qc[warning < 1], ldf_above_qc[warning < 1])).T, fmt='%.10f')