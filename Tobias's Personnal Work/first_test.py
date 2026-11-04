import numpy as np
import matplotlib.pyplot as plt

### CONSTS GENERAL
logs = True
g = 9.81         # gravitation [m/s**2]

### CONSTS SYSTEM (START)


### SIMU PARAMS
dt = 0.001     # pas (dt) [s]
START_T = 0.0 # durée [s]
END_T = 20.0       # durée [s]


### EQUATIONS
equation_x = lambda x, v, a, dt: x + v * dt + 0.5 * a * dt**2
equation_v = lambda x, v, a, dt: v + a * dt
equation_a = lambda x, v, a, dt: a

### FUNCTIONS
def log(text):
    if logs:
        print(text)

### RUN SIMULATION
log("Simulation start")

