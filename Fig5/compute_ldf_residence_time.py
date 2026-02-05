import numpy as np
from scipy.optimize import root_scalar

filename = './Analytics/ldf_residence_time_theory.txt'

xstar = root_scalar(lambda x: np.tan(x) * x - 1., bracket=(0, np.pi / 3.), x0=0.)
pc = 1. + xstar.root ** 2
qc = 0.5 * ( 1. + 1. / pc )

def s(p):
    xf = root_scalar(lambda x: np.tan(x) * x - p, bracket=(0, np.pi / 2. - 1e-10), x0=0.)
    yp = root_scalar(lambda x: np.tan(x) * x - np.sqrt(p - x * x), bracket=(0, xf.root), x0=0.)
    return p - yp.root ** 2 - 1.
    
def ldf(q):
    if q < qc:
        return pc * q
    elif q < 1.:
        yp = root_scalar(lambda x: np.tan(x) * (np.tan(x) + x * (1. + np.tan(x) ** 2)) - q / (1. - q), bracket=(1e-10, np.pi / 2. - 1e-10), x0=1.).root
        pstar = yp ** 2 * (1. + np.tan(yp) ** 2)
        return pstar * q - s(pstar)
        
x = np.linspace(0., 1., 10000)
x = x[:-1]
x_below_qc = x[x < qc]
ldf_below_qc = np.array([ldf(xx) for xx in x_below_qc])
x_above_qc = x[x >= qc]
ldf_above_qc = np.array([ldf(xx) for xx in x_above_qc])
with open(filename, 'w') as f:
    np.savetxt(f, np.vstack((x_below_qc, ldf_below_qc)).T, fmt='%.10f')
    f.write('\n\n')
    np.savetxt(f, np.vstack((x_above_qc, ldf_above_qc)).T, fmt='%.10f')