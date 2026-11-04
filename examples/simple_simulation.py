#
# Simulation: un bloc attaché à un ressort sur un sol plat
#
# Charles Pecheur 2018-2025
#

import math
import matplotlib.pyplot as plt
import numpy as np

### Constantes

g = 9.81         # gravitation [m/s**2]

### Paramètres du système

m = 0.1          # masse du bloc [kg]
mu = 0.03        # coefficient de frottement visqueux [N*s/m]
k = 0.5          # coefficient du ressort [N/m]

### Paramètres de la simulation

step = 0.001     # pas (dt) [s]
end = 20.0       # durée [s]
x_0 = -5.0       # position initiale [m]
v_0 = 0.0        # vitesse initiale [m/s]

t = np.arange(0, end, step)

### Solution analytique (sans frottement)

omega = math.sqrt(k/m)
x_t = x_0*np.cos(omega*t) + v_0/omega*np.sin(omega*t)

### Simulation

x = np.empty_like(t)          
v = np.empty_like(t)
a = np.empty_like(t)
    
# conditions initiales
x[0] = x_0
v[0] = v_0

for i in range(len(t)-1):

    dt = step
    
    # calcul de la force totale
    f_ressort = -k * x[i]
    f_frott = -mu * v[i]
    f_tot = f_ressort + f_frott
    
    # calcul accélération, vitesse, position
    a[i] = f_tot / m
    v[i+1] = v[i] + a[i] * dt
    x[i+1] = x[i] + v[i] * dt
    a[i+1] = a[i]

### Graphique trajectoire

plt.figure(1)
plt.plot(t, x, label="simulation")
plt.plot(t, x_t, label="analytique")
plt.xlabel("t [s]")
plt.ylabel("x [m]")
plt.legend()
plt.show()

### Graphique énergies

e_ressort = k * x**2 / 2
e_cin = m * v**2 / 2
e_tot = e_ressort + e_cin

plt.figure(2)
plt.plot(t, e_ressort, label="ressort")
plt.plot(t, e_cin, label="cinétique")
plt.plot(t, e_tot, label="total")
plt.legend()
plt.xlabel("t [s]")
plt.ylabel("E [J]")
plt.show()

### Graphique de phase

plt.figure(3)
plt.plot(x, v)
plt.xlabel("x [m]")
plt.ylabel("v [m/s]")
plt.show()

