import numpy as np
import matplotlib.pyplot as plt


from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)



dt = 1e-2;
mass = 1.;
charge = 1.;
vAc = 3e-4;

duration = 1000;

v = np.array([0., 1., 0.]);
x = np.array([-1., 0., 0.]);

B = np.array([0., 0., 1.]);
E = np.array([0., 0., 0.]);

X = np.zeros((duration,3)) 
V = np.zeros((duration,3)) 

for time in range(duration):
    t = charge / mass * B * 0.5 * dt;
    s = 2. * t / (1. + t*t);
    v_minus = v + charge / (mass * vAc) * E * 0.5 * dt;
    v_prime = v_minus + np.cross(v_minus,t);
    v_plus = v_minus + np.cross(v_prime,s);
    v = v_plus + charge / (mass * vAc) * E * 0.5 * dt;
    x += v * dt;
    X[time,:] = x;
    V[time,:] = v; 



plt.plot(X[:,0],X[:,1],'k',linewidth=2.0); 
plt.xlabel(r'$x/d_{\rm p}$',fontsize=16)
plt.ylabel(r'$y/d_{\rm p}$',fontsize=16)

plt.show()