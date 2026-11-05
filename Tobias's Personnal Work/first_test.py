import numpy as np
import matplotlib.pyplot as plt

### CONSTS GENERAL
logs = True
g = 9.81         # gravitation [m/s**2]

### CONSTS SYSTEM (START) - TO BE MODIFIED
ETAT_0 = np.array([[
    0.0,  # t
    0.0,  # x
    9.0,  # y
    0.0,  # vx
    -2.0,  # vy
    1.0,  # ax
    1.0   # ay
]])

### SIMU PARAMS
DT = 0.001     # pas (dt) [s]
START_T = 0.0 # durée [s]
END_T = 20.0       # durée [s]

MIN_X = 0.0
MAX_X = 10.0
MIN_Y = 0.0
MAX_Y = 10.0

### EQUATIONS
equation_xy = lambda x, y, vx, vy, ax, ay, dt: (x + vx * dt + 0.5 * ax * dt**2, y + vy * dt + 0.5 * ay * dt**2)
equation_v = lambda x, y, vx, vy, ax, ay, dt: (vx + ax * dt, vy + ay * dt)
equation_a = lambda x, y, vx, vy, ax, ay, dt: (ax, ay)

### FUNCTIONS
def log(text):
    if logs:
        print(text)

def check_limits(state):
    x, y = state[1], state[2]
    if x < MIN_X or x > MAX_X or y < MIN_Y or y > MAX_Y:
        log(f"Warning: State out of bounds at t={state[0]}: x={x}, y={y}")
        return False
    return True

### SETUP arrays
state = np.empty((int((END_T - START_T) / DT), 7))  # t, x, y, vx, vy, ax, ay
state[0] = ETAT_0[0]
state[:, 0] = np.arange(START_T, END_T, DT)  # Fill the first column with time values

### RUN SIMULATION
log("Simulation start")

for i in range(len(state) - 1):
    x, y, vx, vy, ax, ay = state[i, 1], state[i, 2], state[i, 3], state[i, 4], state[i, 5], state[i, 6]
    
    # Update position
    new_x, new_y = equation_xy(x, y, vx, vy, ax, ay, DT)
    
    # Update velocity
    new_vx, new_vy = equation_v(x, y, vx, vy, ax, ay, DT)
    
    # Update acceleration (for this example we keep it constant)
    new_ax, new_ay = equation_a(x, y, vx, vy, ax, ay, DT)
    
    # Store the new state
    state[i+1] = [state[i+1][0], new_x, new_y, new_vx, new_vy, new_ax, new_ay]
    if not check_limits(state[i+1]):
        log(f"Simulation stopped at t={state[i+1][0]} due to out-of-bounds state.")
        state = state[:i+2]  # Keep only the valid states
        break

log("Simulation end")

### PLOT RESULTS
t = state[:, 0]
x, y = state[:, 1], state[:, 2]
v = np.sqrt(state[:, 3]**2 + state[:, 4]**2)
a = np.sqrt(state[:, 5]**2 + state[:, 6]**2)


fig, ax = plt.subplot_mosaic([['up', 'up'],
                               ['lowleft', 'lowright']], figsize=(10, 8))
ax['up'].plot(x, y)
ax['up'].set_title("Trajectory")
ax['up'].set_xlabel("X Position [m]")
ax['up'].set_ylabel("Y Position [m]")

ax['lowleft'].plot(t, v)
ax['lowleft'].set_title("Velocity over Time")
ax['lowleft'].set_xlabel("Time [s]")
ax['lowleft'].set_ylabel("Velocity [m/s]")

ax['lowright'].plot(t, a)
ax['lowright'].set_title("Acceleration over Time")
ax['lowright'].set_xlabel("Time [s]")
ax['lowright'].set_ylabel("Acceleration [m/s^2]")

plt.show()
