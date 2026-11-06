import numpy as np
import matplotlib.pyplot as plt
 
# Paramètres physiques
g = 9.81       # Accélération de la pesanteur
m = 1.0        # Masse
fd = 0.4       # Coefficient de frottement
alpha = 0.0    # Angle d'inclinaison (tan(alpha) = h / d)
d = 15         # distance horizontale entre les 2 points d'accroches
h = -2.5         # différence de hauteur entre les 2 points d'accroches

# Définition de la forme de la piste v(u) et de ses dérivées
# Equation de l'éllipse avec les foyers a et b
f=sqrt(d**2+h**2)      # Distance entre les 2 points d'accroches     
L=1.01*f               # Longueur du câble 
a= L/2
b= sqrt(L**2-f**2)/2
 
def track_v(u):
    u_clipped = np.clip(u, -a, a)
    return b * np.sqrt(1.0 - (u_clipped**2))
 
def track_v_prime(u):
    return 
 
def track_v_double_prime(u):
    return -0.5 * np.cos(u)
 
# Paramètres de simulation
h = 0.001       # Pas de temps (dt)
t_max = 30.0   # Temps total
n_steps = int(t_max / h)
 
# Initialisation des tableaux
t = np.linspace(0, t_max, n_steps)
u = np.zeros(n_steps)  # Coordonnée u
w = np.zeros(n_steps)  # Vitesse u_dot (dérivée de u)
 
# Conditions initiales
u[0] = 0.0     # Position u initiale
w[0] = 0.0     # Vitesse initiale
 
# Fonction d'accélération (u_ddot) basée sur l'équation de l'annexe
def acceleration_u(u_current, w_current):
    vp = track_v_prime(u_current)
    vpp = track_v_double_prime(u_current)
    
        # Termes de l'équation
    terme_courbure = (vp * vpp / (1.0 + vp**2)) * (w_current**2)
    terme_frottement = (fd / m) * np.sqrt(1.0 + vp**2) * w_current
    terme_gravite = g * (np.sin(alpha) + vp * np.cos(alpha)) / np.sqrt(1.0 + vp**2)
    
    # Équation isolant u_ddot (ou ä)
    u_ddot = -terme_courbure - terme_frottement - terme_gravite
    return u_ddot
 
# Boucle d'intégration numérique (Méthode d'Euler)
for i in range(n_steps - 1):
    # Calcul de l'accélération à l'instant i
    a_u = acceleration_u(u[i], w[i])
    
    # Mise à jour de la position u : u_{i+1} = u_i + w_i * h
    u[i + 1] = u[i] + w[i] * h
    
    # Mise à jour de la vitesse w : w_{i+1} = w_i + a_u * h
    w[i + 1] = w[i] + a_u * h
 
# Si vous avez besoin de la position verticale y ou de la position réelle s le long de la courbe :
y = track_v(u)
 
# Affichage des résultats
plt.plot(t, u, label="Position u(t)")
plt.plot(t, w, label="Vitesse u_dot(t)")
plt.xlabel("Temps (s)")
plt.legend()
plt.show()
