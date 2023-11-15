import numpy as np


def newmarks_method(M, C, K, Phi, u0, v0, dt, gamma, beta, timesteps, u_g0):
    # Step 1.0: Initial calculations
    q = u0
    q_dot = v0
    q_double_dot = np.linalg.solve(M, np.zeros([5,1]) - np.dot(C, q_dot) - np.dot(K, q))

    # Step 1.5: Define coefficients
    a1 = (1 / (beta * dt ** 2)) * M + (gamma / (beta * dt)) * C
    a2 = (1 / (beta * dt)) * M + ((gamma / beta) - 1) * C
    a3 = ((1 / (2 * beta)) - 1) * M + dt * ((gamma / (2 * beta)) - 1) * C

    # Step 1.6: Formulate modified stiffness matrix
    H_K = K + a1

    # Initialize arrays to store results
    num_dofs = M.shape[0]
    time = np.zeros(timesteps)
    displacements = np.zeros((timesteps, num_dofs))
    velocities = np.zeros((timesteps, num_dofs))
    accelerations = np.zeros((timesteps, num_dofs))

    # Store initial values
    time[0] = 0
    displacements[0] = q.flatten()
    velocities[0] = q_dot.flatten()
    accelerations[0] = q_double_dot.flatten()

    # Step 2.0: Calculations for each time step
    for i in range(1, timesteps):
        # Step 2.1: Calculate modified equivalent force vector
        u_g_double_dot = u_g0 * np.sin(2 * np.pi * (i-1) * dt)
        p = -0.2591 * np.array([1,1,1,1,1]).reshape(-1,1) * u_g_double_dot

        print(np.dot(Phi.T, p))
        H_P = np.dot(Phi.T, p) + np.dot(a1, q) + np.dot(a2, q_dot) + np.dot(a3, q_double_dot)

        # Step 2.2: Solve the equation H_K * q(i+1) = H_P to determine displacements at the next time step
        q = np.linalg.solve(H_K, H_P)

        # Step 2.3: Calculate velocities at the next time step
        q_dot = (gamma / (beta * dt)) * (q.T - displacements[i - 1]) + (1 - (gamma / beta)) * velocities[i - 1] + dt * (
                1 - (gamma / (2 * beta))) * accelerations[i - 1]
        q_dot = q_dot.T

        # print (f"q_dot: {q_dot}")
        # print (f"q: {q}")
        # print (f"displacements: {displacements[i - 1]}")
        # print (f"velocities: {velocities[i - 1]}")
        # print (f"accelerations: {accelerations[i - 1]}")
        # q_dot = (q.T - displacements[i - 1]) + velocities[i - 1] + accelerations[i - 1]
        # q_dot = q_dot.T

        # Reshape q_dot to match the shape of velocities[i]
        # q_dot = q_dot.reshape(velocities[i - 1].shape)

        # Step 2.4: Calculate accelerations at the next time step
        q_double_dot = (1 / (beta * dt ** 2)) * (q.T - displacements[i - 1]) - (1 / (beta * dt)) * velocities[i - 1] - (
                (1 / (2 * beta)) - 1) * accelerations[i - 1]
        q_double_dot = q_double_dot.T

        # Step 2.5: Calculate response vector
        u = np.dot(Phi, q)

        # Store results
        time[i] = i * dt
        displacements[i] = q.flatten()
        velocities[i] = q_dot.flatten()
        accelerations[i] = q_double_dot.flatten()

    return time, displacements, velocities, accelerations


# Example usage
# Define system matrices: M, C, K
# Define mode shape matrix: Phi
# Define applied force vector: p
# Define initial displacements and velocities: u0, v0
# Define time step size: dt
# Define gamma and beta coefficients
# Define number of time steps: timesteps

# time, displacements, velocities, accelerations = newmarks_method(M, C, K, Phi, p, u0, v0, dt, gamma, beta, timesteps)


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

# print(f"C:\n{C}")
# print(f"q_dot:\n{v0}")
# print(f"dot product:\n{np.dot(C, v0)}")


time, displacements, velocities, accelerations = newmarks_method(M, C, K, Phi, u0, v0, dt, gamma, beta, timesteps, u_g0)
print(time)
print(displacements)
print(velocities)
print(accelerations)
