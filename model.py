# Python Implementation

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# System of equations
def f(s, t, a0, a1, a2, a3, a4, a5, a6):	
	Ns = s[0] # susceptible population
	Ni = s[1] # infected population
	Nr = s[2] # recovered population
	Nss, Nsi, Nsr, Nii, Nir, Nrr = s[3], s[4], s[5], s[6], s[7], s[8]	
	pi = a0 # recruitment rate
	u_i = a1 # death rate
	a = a2 # association rate
	b = a3 # transmission rate
	d = a4 # dissociation rate
	e = a5 # association/dissociation scaling
	y = a6 # recovery rate
	dNsdt = pi - (2*(a/e)*(Ns**2)) - ((a/e)*Ns*Ni) - ((a/e)*Ns*Nr) + (2*(d/e)*Nss) + (u_i+(d/e))*Nsi + ((d/e)*Nsr)
	dNidt = -((u_i+y)*Ni) - (2*(a/e)*(Ni**2)) - ((a/e)*Ns*Ni) - ((a/e)*Ni*Nr) + ((d/e)*Nir) + 2*(u_i+(d/e))*Nii + ((d/e)*Nsi)
	dNrdt = -(2*(a/e)*(Nr**2)) - ((a/e)*Ni*Nr) - ((a/e)*Ns*Nr) + (y*Ni) + ((d/e)*Nsr) + (u_i+(d/e)*Nir) + (2*(d/e)*Nrr)
	dNssdt = -((d/e)*Nss) + (a/e)*(Ns**2)
	dNsidt = -(u_i+(d/e)+b+y)*Nsi + ((a/e)*Ns*Ni)
	dNsrdt = -((d/e)*Nsr) + ((a/e)*Ns*Nr) + (y*Nsi)
	dNiidt = -(2*u_i+(d/e)+2*y)*Nii + (b*Nsi) + ((a/e)*(Ni**2))
	dNirdt = -((u_i+(d/e)+y)*Nir) + ((a/e)*Ni*Nr) + (2*y*Nii)
	dNrrdt = -((d/e)*Nrr) + ((a/e)*(Nr**2)) + (y*Nir)
	return [dNsdt, dNidt, dNrdt, dNssdt, dNsidt, dNsrdt, dNiidt, dNirdt, dNrrdt]

# Solver
t = np.linspace(1,1000, num=500)
s0 = [100, 10, 0, 0, 0, 0, 0, 0, 0]
args = (2, 0.01, 1, 0.5, 1, 1, 0.05)
s = odeint(f, s0, t, args)
Ns, Ni, Nr, Nss, Nsi, Nsr, Nii, Nir, Nrr = s[:,0], s[:,1], s[:,2], s[:,3], s[:,4], s[:,5], s[:,6], s[:,7], s[:,8]

# Initial plot
plot_axes = plt.axes([0.1, 0.19, 0.8, 0.68])
plt.axes(plot_axes)
susceptible, = plt.plot(np.log(t), Ns+Nsi+Nsr+2*Nss, 'b-', linewidth=2.0, label='Susceptible')
infected, = plt.plot(np.log(t), Ni+Nsi+Nir+2*Nii, 'r-', linewidth=2.0, label='Infected')
recovered, = plt.plot(np.log(t), Nr+Nsr+Nir+2*Nrr, 'g-', linewidth=2.0, label='Recovered')
plt.xlabel('ln(t)')
plt.ylabel('Ns(t), Ni(t), Nr(t)')
plt.ylim(0, 250)
plt.legend(loc='best')

# Create sliders
sus_init = plt.axes([0.175, 0.02, 0.25, 0.03])
inf_init = plt.axes([0.175, 0.07, 0.25, 0.03])
recov_rate = plt.axes([0.67, 0.02, 0.25, 0.03])
rec_rate = plt.axes([0.67, 0.07, 0.25, 0.03])
trans_rate = plt.axes([0.2, 0.9, 0.25, 0.03])
assoc_rate = plt.axes([0.67, 0.9, 0.25, 0.03])
dissoc_rate = plt.axes([0.2, 0.95, 0.25, 0.03])
dth_rate = plt.axes([0.67, 0.95, 0.25, 0.03])
initial_sus = Slider(sus_init, 'Initial Susceptible', 0, 200, valinit=s0[0])
initial_inf = Slider(inf_init, 'Initial Infected', 0, 200, valinit=s0[1])
recover_rate = Slider(recov_rate, 'Recovery Rate', 0, 1, valinit=0.05)
recruitment_rate = Slider(rec_rate, 'Recruitment Rate', 0, 15, valinit=2)
transmission_rate = Slider(trans_rate, 'Transmission Rate', 0, 1, valinit=0.5)
association_rate = Slider(assoc_rate, 'Association Rate', 0, 1, valinit=1)
dissociation_rate = Slider(dissoc_rate, 'Dissociation Rate', 0, 1, valinit=1)
death_rate = Slider(dth_rate, 'Death Rate', 0, 1, valinit=0.01)

# Recalculate and update plot action
def update(val):
	s0[0] = initial_sus.val
	s0[1] = initial_inf.val
	args =(recruitment_rate.val, death_rate.val, association_rate.val, transmission_rate.val, dissociation_rate.val, 1, recover_rate.val)	
	s = odeint(f, s0, t, args)
	Ns, Ni, Nr, Nss, Nsi, Nsr, Nii, Nir, Nrr = s[:,0], s[:,1], s[:,2], s[:,3], s[:,4], s[:,5], s[:,6], s[:,7], s[:,8]
	susceptible.set_ydata(Ns+Nsi+Nsr+2*Nss)
	infected.set_ydata(Ni+Nsi+Nir+2*Nii)
	recovered.set_ydata(Nr+Nsr+Nir+2*Nrr)
	fig.canvas.draw_idle()

# Set trigger action
initial_sus.on_changed(update)
initial_inf.on_changed(update)
recover_rate.on_changed(update)
recruitment_rate.on_changed(update)
transmission_rate.on_changed(update)
association_rate.on_changed(update)
dissociation_rate.on_changed(update)
death_rate.on_changed(update)

plt.show()


