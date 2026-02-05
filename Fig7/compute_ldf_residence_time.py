import numpy as np
from scipy.optimize import root_scalar, root

xstar = root_scalar(lambda x : np.tan(x) * x - 1., bracket = (0, np.pi / 3.), x0 = 0.)
pc = 9. * xstar.root ** 2
qc = 3. * ( 2. + xstar.root ** 2 ) / ( 12. + 11. * xstar.root ** 2 )

def tanhc(u):
    if u != 0:
        return np.tanh(u) / u
    else:
        return 1.

def tanc(u):
    if u != 0:
        return np.tan(u) / u
    else:
        return 1.

def K(u):
    if u != 0:
        return ( 1. - np.tanh(u) ** 2 - tanhc(u) ) / u ** 2
    else:
        return - 2. / 3.

def J(u):
    if u != 0:
        return ( 1. + np.tan(u) ** 2 - tanc(u) ) / np.tan(u) ** 2
    else:
        return 2. / 3.

def dres(u,p):
    return ( p - 2. * u * u ) * tanhc(2. * u / 3.) * tanc(np.sqrt(p - u * u) / 3.) - 9.

def dres2(y,p):
    ym2 = p / 9.
    return tanc(y) - 1. / ( 2. * y * y - ym2 ) / tanhc(2. * np.sqrt(ym2 - y * y))

def uc(p):
    if p > pc:
        tol = 1e-10
        if np.sqrt(p) < 3. * np.pi / 2.:
            ymin = np.sqrt(p / 18.) + tol
            ymax = np.sqrt(p) / 3. - tol
        elif np.sqrt(p) < 3. * np.pi / np.sqrt(2.):
            ymin = np.sqrt(p / 18.) + tol
            ymax = np.pi / 2. - tol
        elif np.sqrt(p) < 3 * np.sqrt(2.) * np.pi:
            ymin = np.pi / 2. + tol
            ymax = np.sqrt(p / 18.) - tol
        else:
            ymin = np.pi / 2. + tol
            ymax = np.pi
        yp = root_scalar(lambda y: dres2(y,p), bracket = (ymin, ymax), x0 = ymax, xtol = tol)
        return np.sqrt(p - 9. * yp.root ** 2)
    else:
        return 0.
x = np.linspace(0., 0.88, 10000)
x = x[:-1]
x_below_qc = x[x < qc]
ldf_below_qc = pc * x_below_qc
x_above_qc = x[x >= qc]
ldf_above_qc = np.zeros(len(x_above_qc))
warning = np.zeros(len(x_above_qc))
p = np.linspace(0., 100., 20000)
mu = np.array([uc(pp) ** 2 for pp in p])
p = p[np.isnan(mu) == False]
mu = mu[np.isnan(mu) == False]
np.savetxt('./Analytics/scgf_residence_time_theory.txt', np.vstack((p, mu)).transpose())
for i, xx in enumerate(x_above_qc):
    j = np.argmax(p * xx - mu)
    ldf_above_qc[i] = p[j] * xx - mu[j]
    if j == len(p) - 1:
        warning[i] = 1.
with open('./Analytics/ldf_residence_time_theory.txt', 'w') as f:
    np.savetxt(f, np.vstack((x_below_qc, ldf_below_qc)).T, fmt='%.10f')
    f.write('\n\n')
    np.savetxt(f, np.vstack((x_above_qc[warning < 1], ldf_above_qc[warning < 1])).T, fmt='%.10f')
