import numpy as np
from scipy.integrate import quad
from scipy.special import gammainc

#################
### FUNCTIONS ###
#################
def gt(a, t):
    """
    Compute the function g_t(x) for the calculation of the exact pdf.
    """
    return np.exp(- 3. * a * a * 0.25 / t ** 3) / np.sqrt(4. * np.pi * t ** 3 / 3.)

def proba(a, t):
    """
    Compute the exact pdf of the area under the trajectory of the walker before diing.
    """
    g = gt(a, t)
    return np.exp(- t) * g + quad(lambda tt: np.exp(- tt) * gt(a, tt), 0., t)[0]

##################################
### DEFINE THE HISTOGRAM RANGE ###
##################################
Amin = -1.45
Amax = 1.45
T = 100.
dt = 0.01
dA = np.sqrt(2. * dt) / T ** 2 * 10
As = np.arange(Amin, Amax + dA, dA)

#############################################
### COMPUTE AND STORE THE THEORETICAL LDF ###
#############################################
ldf = np.array([ - np.log(proba(AAA * T ** 2, T) * T ** 2) / T for AAA in As ])
np.savetxt('./Analytics/ldf_area_theory.txt', np.vstack((As, ldf)).transpose())
