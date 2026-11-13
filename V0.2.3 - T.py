import numpy as np
import matplotlib.pyplot as plt

# Paramètres physiques
g = 9.81
m = 1.0
fd = 0.4
H = -2.75
d = 15
alpha = np.arctan(H/d)
f = np.sqrt(H**2 + d**2)
L = 15.5

# Définition de la forme de la piste v(u) et de ses dérivées
a = L/2
b = np.sqrt(L**2 - f**2)/2

def track_v(u):
    u_clipped = np.clip(u, -a, a)
    return -b * np.sqrt(1.0 - (u_clipped**2 / a**2))

def track_v_prime(u):
    v_val = track_v(u)
    if v_val == 0:
        return 0.0
    return -(u * b**2) / (v_val * a**2)

def track_v_double_prime(u):
    v_val = track_v(u)
    vp_val = track_v_prime(u)
    if v_val == 0:
        return 0.0
    return -(b**2 / a**2) * (v_val - u * vp_val) / (v_val**2)

# Paramètres de simulation
h = 0.001
t_max = 30.0
n_steps = int(t_max / h)

# Initialisation des tableaux
t = np.arange(n_steps) * h
u = np.zeros(n_steps)
w = np.zeros(n_steps)

# Conditions initiales
u[0] = -a + 0.001
w[0] = 0.0

# Fonction d'accélération
def acceleration_u(u_current, w_current):
    vp = track_v_prime(u_current)
    vpp = track_v_double_prime(u_current)

    terme_courbure = (vp * vpp / (1.0 + vp**2)) * w_current**2
    terme_frottement = (fd / m) * np.sqrt(1.0 + vp**2) * w_current
    terme_gravite = g * (np.sin(alpha) + vp * np.cos(alpha)) / (1.0 + vp**2)

    u_ddot = -terme_courbure - terme_frottement - terme_gravite
    return u_ddot

# Boucle d'intégration numérique (Méthode d'Euler)
for i in range(n_steps - 1):
    a_u = acceleration_u(u[i], w[i])

    # Mise à jour de la position
    u[i + 1] = u[i] + w[i] * h

    # Mise à jour de la vitesse
    w[i + 1] = w[i] + a_u * h

    # Arrêt à l'extrémité de l'ellipse
    if u[i + 1] >= a - 0.001:
        u = u[:i+2]
        w = w[:i+2]
        t = t[:i+2]
        break

# Calcul de la trajectoire
v = track_v(u)

# Conversion du repère (u,v) vers (x,y)
x = u * np.cos(alpha) - v * np.sin(alpha)
y = u * np.sin(alpha) + v * np.cos(alpha)

# Calcul des énergies
e_p = m * g * (y - np.min(y))
e_c = m * w**2 / 2
e_t = e_p + e_c

# Affichage de la trajectoire
fig, axd = plt.subplot_mosaic([['speed','speed'],
                               ['pos','pos'],
                               ['posv','posv'],
                               ['energy','energy%']], layout='constrained', figsize=(8, 9))

axd['speed'].plot(t, w, label="Vitesse du téléphérique")
axd['speed'].set_xlabel("Temps [s]")
axd['speed'].set_ylabel("Vitesse [m/s]")
axd['speed'].set_title("Vitesse du téléphérique")
axd['speed'].grid()

axd['pos'].plot(x, y, label="Positions du téléphérique")
axd['pos'].plot(x[::100], y[::100], "+", label="Positions du téléphérique")
axd['pos'].set_xlabel("Position horizontale x [m]")
axd['pos'].set_ylabel("Position verticale y [m]")
axd['pos'].set_title("Positions du téléphérique")
axd['pos'].grid()

axd['posv'].plot(t, v, label="Hauteur du téléphérique (repère v u)")
axd['posv'].set_xlabel("Temps [s]")
axd['posv'].set_ylabel("Hateur [m]")
axd['posv'].set_title("Hauteur du téléphérique (repère v u)")
axd['posv'].grid()

axd['energy'].plot(t, e_c, label="Energie cinétique")
axd['energy'].plot(t, e_p, label="Energie potentielle")
axd['energy'].plot(t, e_t, label="Energie totale")

axd['energy'].set_xlabel("Temps [s]")
axd['energy'].set_ylabel("Energie [J]")
axd['energy'].set_title("Energie du téléphérique")
axd['energy'].legend()
axd['energy'].grid()


axd['energy%'].plot(t, e_t / e_t[0] * 100, label="Rendement")

axd['energy%'].set_xlabel("Temps [s]")
axd['energy%'].set_ylabel("Rendement (%)")
axd['energy%'].set_title("Rendement du téléphérique")
axd['energy%'].grid()

plt.show()
