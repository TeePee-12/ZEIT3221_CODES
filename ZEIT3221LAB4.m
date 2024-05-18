%% Section 1 - making a state space model
%declaring the matixes
E = [0 1 0 0 0 0 ;0 -0.25 10.32 0 -0.42 0;0 0 0 1 0 0 ;0 0.16 -43.38 0 35.96 0 ;0 0 0 0 0 1;0 -0.01 190.92 0 -188.73 0];
f = transpose([0 0.125 0 -0.0803 0 0.0105]);
G = [0 0.002 -35.2238 0 26.6023 0; 1 0 0 0 0 0 ;0 1 0 0 0 0 ];
h = transpose([-0.002003 0 0 ]);
Ts = 0.02;

% CT State-Space Model
sysc = ss(E,f,G,h);
% Sampled  Data Model
sysd = c2d(sysc,Ts,'zoh');

% extact the ABCD matrixes from the sampled data model to suit the form
% q[n + 1] = Aq[n] + Bx[n]
% y[n] = Cq[n] + Dx[n]
[A, B, C, D] = ssdata(sysd)

%Determine the eigenvalues and eigenvectors
[V, Lambda] = eig(A)

%Controllability
syscont = ctrb(A,B)
ctrlrank = rank(syscont)

%Observability
sysobs = obsv(A,C)
obsvrank = rank(sysobs)
%% Zero Response

% Define the initial condition q0
q0 = [0.1; 0; 0.05; 0; 0.1; 0];
n=100000;
qn = A^(n)*q0;
An = A^(n)
An_2 = real(V*(Lambda^n)*inv(V))

% Final time
Tf = 10;

% Generate zero-input response using initial function
[yt, t, qt] = initial(sysd, q0, 0:Ts:Tf);

% Extract state variables q and ˙θ1 from the state trajectory matrix qt
qdot = qt(:, 1);  % State variables q
theta_dot_1 = qt(:, 4);  % ˙θ1

% Plot the outputs
figure;
subplot(2, 1, 1);
stairs(t, yt);
xlabel('Time (s)');
ylabel('Outputs');
title('Output Plots');

% Plot state variables q and ˙θ1
subplot(2, 1, 2);
stairs(t, qdot);
hold on;
stairs(t, theta_dot_1);
xlabel('Time (s)');
ylabel('State Variables');
legend('q','˙θ1');
title('State Variable Plots');

%% General resposne


A0 = 2;     
Omega0 = 3; 
n = 100000;
q = [0.1; 0; 0.05; 0; 0.1; 0];

% Calculate the state at the nth time step using the general response equation
for k = 0:n
    x_k = A0 * sin(Omega0 * k * Ts);
    % Update state using the general response equation
    q = A * q + B * x_k;
end
disp(q);

Tf = 10;
t = 0:Ts:Tf;
x = A0 * sin(Omega0 * t);
% Simulate the system response using lsim
[yt, t, qt] = lsim(sysd, x, t, q0);

% Extract state variables q and ˙θ1 from the state trajectory matrix qt
qstate = qt(:, 1);  % State variables q
theta_dot_1 = qt(:, 4);  % ˙θ1

% Plot the outputs and input
figure;
subplot(2, 1, 1);
plot(t, x, 'b', 'LineWidth', 1.5); % Input signal
hold on;
plot(t, yt, 'r', 'LineWidth', 1.5); % Output signal
xlabel('Time (s)');
ylabel('Input and Output');
title('Input and Output Plots');
legend('x[n]', 'y[n]');
grid on;

% Plot state variables q and ˙θ1
subplot(2, 1, 2);
stairs(t, qstate);
hold on;
stairs(t, theta_dot_1);
xlabel('Time (s)');
ylabel('State Variables');
legend('q', '˙θ1');
title('State Variable Plots');
grid on;

%% Section 2a - Transfer Functions
% Generate the transfer functions from the state-space model

% Display the transfer functions
fprintf('Transfer Functions:\n');
tf(sysd)

% Express transfer functions in zero-pole-gain form
G_zpk = zpk(sysd);

% Display transfer functions in zero-pole-gain form
fprintf('Transfer Functions in Zero-Pole-Gain Form:\n');
disp(G_zpk);

%% Section 2b - Pole-Zero Plots
% Generate pole-zero plots for each input-output transfer function
iopzmap(sysd);

% Comment on pole and zero locations

%% Section 2c - Poles
% Obtain poles of the system
poles = pole(sysd);

% Display poles
fprintf('Poles of the System:\n');
disp(poles);

% Comment on pole locations and their relationship with eigenvalues

%% Section 2d - Frequency Response
% Generate magnitude and phase plots of the frequency response
bode(sysd);

%% Section 3 - Stochastic State-Space Model and Simulation
% Define the standard deviations for process and measurement noises
process_noise_stdev = 0.006;
measurement_noise_stdev = 0.02;

% Number of time steps for the simulation
numSteps = 1000;  % Choose an appropriate number

% Generate process noise and measurement noise signals
w = process_noise_stdev * randn(1, numSteps);
zeta = measurement_noise_stdev * randn(1, numSteps);

% Create the matrix of input signals, including control input, process noise, and measurement noise
input_matrix = zeros(numSteps, 10);  % 10 inputs in total
input_matrix(:, 1) = A0 * sin(Omega0 * (0:Ts:(numSteps-1) * Ts));  % Sinusoidal input
input_matrix(:, 7) = w;
input_matrix(:, 8) = zeta;

% Create the state-space model with stochastic inputs
A_stochastic = A;
B_stochastic = [B, eye(6), zeros(6, 3)];
C_stochastic = C;
D_stochastic = [D, zeros(3, 6), eye(3)];

sysdstochastic = ss(A_stochastic, B_stochastic, C_stochastic, D_stochastic, Ts);

% Simulate the system dynamics with stochastic inputs
[yt_stochastic, t_stochastic, qt_stochastic] = lsim(sysdstochastic, input_matrix);

% Extract state variables q and ˙θ1 from the output matrix yt_stochastic
q_stochastic = qt_stochastic(:, 1);
theta_dot_1_stochastic = qt_stochastic(:, 4);



% Plot the outputs
figure;
subplot(1, 2, 1);
plot(0:Ts:(numSteps-1) * Ts, yt_stochastic);
xlabel('Time (s)');
ylabel('Outputs');
title('Stochastic System Output');
grid on;

subplot(1, 2, 2);
stairs(0:Ts:(numSteps-1) * Ts, q_stochastic);
hold on
stairs(0:Ts:(numSteps-1) * Ts, theta_dot_1_stochastic);

xlabel('Time (s)');
ylabel('State Variables');
title('State Variable Plots');
legend('q','θ1');
grid on;

%% Observer Matrix
p_noise = 0.006*randn(6,501);
m_noise = 0.02*randn(3,501);
p = [0.7 0.78 0.7246+0.2294i 0.7246-0.2294i 0.0860+0.7562i 0.0860-0.7562i];
L = -place(A',C',p)';
A_LC = A + L*C;
eig(A_LC);
B_o = [B+L*D, -L];
ss(A_LC,B_o,C,[D zeros(3,3)]);
q_bar = [0 0 0 0 0 0]';
[observer_output, t_stochastic, observer_states] = lsim(sysdstochastic,[transpose(x) transpose(p_noise) transpose(m_noise)],t,q_bar)
% Extract state variables q_hat and θ1_hat from the observer_states matrix
q_hat = observer_states(:, 1);
theta1_hat = observer_states(:, 4);

% Plot the observer's estimated state variables
figure;
subplot(2, 1, 1);
stairs(t_stochastic, q_hat);
xlabel('Time (s)');
ylabel('Estimated State Variables (q)');
title('Observer State Estimates for q');

subplot(2, 1, 2);
stairs(t_stochastic, theta1_hat);
xlabel('Time (s)');
ylabel('Estimated State Variables (θ1)');
title('Observer State Estimates for θ1');

%% 4c
% Define the initial estimate of the state for the observer
q_hat_initial = [0; 0; 0; 0; 0; 0];

% Simulate the observer using outputs from the stochastic system
[observer_output, t_stochastic, observer_states] = lsim(observer_model, [yt_stochastic, q_stochastic], t_stochastic, q_hat_initial);

% Extract state variables q_hat and θ1_hat from the observer_states matrix
q_hat = observer_states(:, 1:6);
theta1_hat = observer_states(:, 7:9);

% Plot the observer's estimated state variables against the 'true' state variables
figure;
subplot(2, 1, 1);
stairs(t_stochastic, q_hat, 'b', 'LineWidth', 1.5); % Estimated state variables
hold on;
stairs(t_stochastic, q_stochastic, 'r', 'LineWidth', 1.5); % True state variables
xlabel('Time (s)');
ylabel('State Variables');
title('State Observer Estimates vs. True State Variables (q)');

legend('Estimated q','True q');
grid on;

subplot(2, 1, 2);
stairs(t_stochastic, theta1_hat, 'b', 'LineWidth', 1.5); % Estimated state variables
hold on;
stairs(t_stochastic, theta_dot_1_stochastic, 'r', 'LineWidth', 1.5); % True state variables
xlabel('Time (s)');
ylabel('State Variables');
title('State Observer Estimates vs. True State Variables (θ1)');

legend('Estimated θ1', 'True θ1');
grid on;

%% 5 State Feedback Control
p = [0.9,0.95,0.9057+0.2867i,0.9057-0.2867i,0.095+0.9452i,0.095-0.9452i];
%g is the feedback gain vector
g = -place(A, B, p)

A = A + (B*g);

%5b
sysdfb = ss(A, B, C, D, Ts)
poles = pole(sysdfb)
iopzmap(sysdfb)

% Define the initial condition q0
q0 = [0.1; 0; 0.05; 0; 0.1; 0];
n=100000;
qn = A^(n)*q0;
An = A^(n)
An_2 = real(V*(Lambda^n)*inv(V))

% Final time
Tf = 3;

% Generate zero-input response using initial function
[yt, t, qt] = initial(sysdfb, q0, 0:Ts:Tf);

% Extract state variables q and ˙θ1 from the state trajectory matrix qt
qdot = qt(:, 1);  % State variables q
theta_dot_1 = qt(:, 4);  % ˙θ1

% Plot the outputs
figure;
subplot(2, 1, 1);
stairs(t, yt);
xlabel('Time (s)');
ylabel('Outputs');
title('Output Plots');

% Plot state variables q and ˙θ1
subplot(2, 1, 2);
stairs(t, qdot);
hold on;
stairs(t, theta_dot_1);
xlabel('Time (s)');
ylabel('State Variables');
legend('q','˙θ1');
title('State Variable Plots');
%%
tf(sysdfb)
bode(sysdfb)
