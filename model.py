# Python Implementation

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

# System of equations
def f(s, t):
	pi = 2
	u_s = 0.0
	u_i = 0.01
	a = 1
	b = 0.5
	d = 1
	e = 0.3
	Ns = s[0]
	Ni = s[1]
	Nss, Nsi, Nii = s[2], s[3], s[4]
	dNsdt = pi - (u_s*Ns) - (2*(a/e)*(Ns**2)) - ((a/e)*Ns*Ni) + 2*(u_s+(d/e))*Nss + (u_i+(d/e))*Nsi
	dNidt = -(u_i*Ni) - (2*(a/e)*(Ni**2)) - ((a/e)*Ns*Ni) + 2*(u_i+(d/e))*Nii + (u_s+(d/e))*Nsi
	dNssdt = -(2*u_s+(d/e))*Nss + (a/e)*(Ns**2)
	dNsidt = -(u_i+u_s+(d/e)+b)*Nsi + (a/e)*Ns*Ni
	dNiidt = -(2*u_i+(d/e))*Nii + (b*Nsi) + ((a/e)*(Ni**2))
	return [dNsdt, dNidt, dNssdt, dNsidt, dNiidt]

# Solver
t = np.linspace(1,1000, num=500)
s0 = [100, 10, 0, 0, 0]
s = odeint(f, s0, t)
Ns, Ni, Nss, Nsi, Nii = s[:,0], s[:,1], s[:,2], s[:,3], s[:,4]

# Initial plot
plot_axes = plt.axes([0.1, 0.2, 0.8, 0.70])
plt.axes(plot_axes)
susceptible, = plt.plot(np.log(t), Ns+Nsi+2*Nss, 'b-', linewidth=2.0, label='Susceptible')
infected, = plt.plot(np.log(t), Ni+Nsi+2*Nii, 'r-', linewidth=2.0, label='Infected')
plt.xlabel('ln(t)')
plt.ylabel('Ns(t), Ni(t)')
plt.ylim(0, 250)
plt.legend(loc='best')

# Create sliders
sus_init = plt.axes([0.20, 0.02, 0.65, 0.03])
inf_init = plt.axes([0.20, 0.07, 0.65, 0.03])
initial_sus = Slider(sus_init, 'Initial Susceptible', 0, 200, valinit=s0[0])
initial_inf = Slider(inf_init, 'Initial Infected', 0, 200, valinit=s0[1])

# Recalculate and update plot action
def update(val):
	s0[0] = initial_sus.val
	s0[1] = initial_inf.val
	s = odeint(f, s0, t)
	Ns, Ni, Nss, Nsi, Nii = s[:,0], s[:,1], s[:,2], s[:,3], s[:,4]
	susceptible.set_ydata(Ns+Nsi+2*Nss)
	infected.set_ydata(Ni+Nsi+2*Nii)
	fig.canvas.draw_idle()

# Set trigger action
initial_sus.on_changed(update)
initial_inf.on_changed(update)
	
plt.show()


