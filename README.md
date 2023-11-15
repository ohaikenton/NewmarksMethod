# Newmark's Method

Sample implementation in Python of Newmark's Method for solving the equation of motion of a single degree of freedom system.

## Sample Results
Parameters for the sample results:

```angular2html
Phi = np.matrix([[0.3338, -0.8953, 1.1726, 1.0776, -0.6405],
                 [0.6405, -1.1726, 0.3338, -0.8953, 1.0776],
                 [0.8953, -0.6405, -1.0776, -0.3338, -1.1726],
                 [1.0776, 0.3338, -0.6405, 1.1726, 0.8953],
                 [1.1726, 1.0776, 0.8953, -0.6405, -0.3338]])
gamma = 0.5
beta = 1/6
dt = 0.1
K = np.matrix([[31.3,0,0,0,0],
              [0,266.4,0,0,0],
              [0,0,0,662,0],
              [0,0,0,1092.6,0],
              [0,0,0,0,1421.3]])
C = np.matrix([[0.5592,0,0,0,0],
              [0,1.6322,0,0,0],
              [0,0,2.573,0,0],
              [0,0,0,3.3054,0],
              [0,0,0,0,3.77]])
M = np.identity(5)
p = np.zeros([5, 1])
u0 = np.zeros((5, 1))
v0 = np.zeros([5, 1])
timesteps = 22
u_g0 = 193
```

Usage:

```
time, displacements, velocities, accelerations = newmarks_method(M, C, K, Phi, u0, v0, dt, gamma, beta, timesteps, u_g0)
```