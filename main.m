function [time, displacements, velocities, accelerations] = newmarks_method()
    % Step 1.0: Initial calculations

    Phi = [0.3338, -0.8953, 1.1726, 1.0776, -0.6405;
       0.6405, -1.1726, 0.3338, -0.8953, 1.0776;
       0.8953, -0.6405, -1.0776, -0.3338, -1.1726;
       1.0776, 0.3338, -0.6405, 1.1726, 0.8953;
       1.1726, 1.0776, 0.8953, -0.6405, -0.3338];
    gamma = 0.5;
    beta = 1/6;
    dt = 0.1;
    K = [31.3, 0, 0, 0, 0;
         0, 266.4, 0, 0, 0;
         0, 0, 0, 662, 0;
         0, 0, 0, 1092.6, 0;
         0, 0, 0, 0, 1421.3];
    C = [0.5592, 0, 0, 0, 0;
         0, 1.6322, 0, 0, 0;
         0, 0, 2.573, 0, 0;
         0, 0, 0, 3.3054, 0;
         0, 0, 0, 0, 3.77];
    M = eye(5);
    p = zeros(5, 1);
    timesteps = 22;
    u_g0 = 193;

    q = zeros(5, 1);
    q_dot = zeros(5, 1);
    q_double_dot = (M \ zeros(5,1)) - C*q_dot - K*q;

    % Step 1.5: Define coefficients
    a1 = (1 / (beta * dt^2)) * M + (gamma / (beta * dt)) * C;
    a2 = (1 / (beta * dt)) * M + ((gamma / beta) - 1) * C;
    a3 = ((1 / (2 * beta)) - 1) * M + dt * ((gamma / (2 * beta)) - 1) * C;

    % Step 1.6: Formulate modified stiffness matrix
    H_K = K + a1;

    % Initialize arrays to store results
    num_dofs = size(M, 1);
    time = zeros(timesteps, 1);
    displacements = zeros(timesteps, num_dofs);
    velocities = zeros(timesteps, num_dofs);
    accelerations = zeros(timesteps, num_dofs);

    % Store initial values
    time(1) = 0;
    displacements(1, :) = q';
    velocities(1, :) = q_dot';
    accelerations(1, :) = q_double_dot';

    % Step 2.0: Calculations for each time step
    for i = 2:timesteps
        % Step 2.1: Calculate modified equivalent force vector
        u_g_double_dot = u_g0 * sin(2 * pi * (i-1) * dt);
        p = -0.2591 * ones(5, 1) * u_g_double_dot;

        H_P = Phi' * p + a1 * q + a2 * q_dot + a3 * q_double_dot;

        % Step 2.2: Solve the equation H_K * q(i+1) = H_P to determine displacements at the next time step
        q = H_K \ H_P;

        % Step 2.3: Calculate velocities at the next time step
        q_dot = (gamma / (beta * dt)) * (q' - displacements(i - 1, :)) + (1 - (gamma / beta)) * velocities(i - 1, :) + dt * (1 - (gamma / (2 * beta))) * accelerations(i - 1, :);

        % Step 2.4: Calculate accelerations at the next time step
        q_double_dot = (1 / (beta * dt^2)) * (q' - displacements(i - 1, :)) - (1 / (beta * dt)) * velocities(i - 1, :) - ((1 / (2 * beta)) - 1) * accelerations(i - 1, :);

        % Step 2.5: Calculate response vector
        u = Phi * q;

        % Store results
        time(i) = (i - 1) * dt;
        displacements(i, :) = q';
        velocities(i, :) = q_dot';
        accelerations(i, :) = q_double_dot';
    end
end