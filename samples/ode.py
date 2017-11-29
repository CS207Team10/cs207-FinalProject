import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


#2NO + Br2 <-> 2NOBr
def rhs(y, t, kf, kb):
    rf = kf * y[0]**2 * y[1] #forward reaction rate
    rb = kb * y[2]**2 #backward reaction rate
    return [2*(rb - rf), rb - rf, 2*(rf - rb)]


#time interval
tout = np.linspace(0, 20)

#reaction coefficients
k_vals = 0.42, 0.17  # arbitrary in this case
#Initial concentration for each specie
y0 = [1, 1, 0]
#concentration
yout = odeint(rhs, y0, tout, k_vals)  
print(yout)
#concentration plot
plt.plot(tout, yout)
plt.show()

