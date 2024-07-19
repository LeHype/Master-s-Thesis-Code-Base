function [f, f_ipopt,cost_function,constraint_function_left,constraint_function_g,constraint_function_right] = ConstructOCP_X()
%% Make sure the correct models and scripts are in the path
import casadi.*
addpath(['src'])
currentDir = pwd;
parentDir = fileparts(currentDir);
modelsDir = fullfile(parentDir, 'Models');
addpath(modelsDir);

%% Define main properties
nx              = 4;           % Size of the state
nu              = 1;           % Size of input u
foh             = 0;           % First order hold reconstruction is on. (1)nSteps = 100;                    % Number of steps in MPC
timestep = 0.1;               % Discretization
nSteps   = 20;
uMax            = 0.85;        % input voltage is limited !When using dry friction compensation this needs to be reduced accordingly!

% Hardcode l_bar
l_bar = 0.4;
% Weighing vector for each component w(1) = u, w(2) = stabilization, w(3)
% tracking cost 
w= [0.02,0,1];
% wQ are the weights for all 4 states, i.e. wQ(1) = xc/yc 
wQ = [1,0.3,0.01,0.01];

% Construct solver
ocp = casadi.Opti();
qpsol_options = struct;
qpsol_options.print_iter = false;
qpsol_options.print_info = false;
qpsol_options.print_lincomb = false;
qpsol_options.print_header= false;
s_opts = struct("qpsol", 'qrqp', "max_iter", 2000, 'qpsol_options',qpsol_options,'print_time',false,'print_header',false,'verbose',0,'print_status',false);
ocp.solver('sqpmethod',s_opts);
ocp.subject_to( )  

% Load parametric trajectory for trajectory tracking. 
% Here, scenario 3 is loaded
% load('Scenario_1_trajectory.mat')
% load('Scenario_2_trajectory.mat')
load('Scenario_3_trajectory.mat')


% The trajectory only depends on the time t_p that changes every mpc
% increment. 
t_p = ocp.parameter(1,1);
x_ref = casadi.MX(2*nx,nSteps);
u_ref = casadi.MX(2*nu,nSteps);

for i = 1:nSteps
x_ref(:,i) = evaluateSines(sine_params_trajectory,(t_p+timestep*(i-1)));
u_ref(:,i) = evaluateSines(sine_params_trajectory_u,(t_p+timestep*(i-1)));
end

%% split for subsystems xX and xY
x_X_ref = x_ref([1,3,5,7],:);
u_X_ref = u_ref([1],:);
x_Y_ref = x_ref([2,4,6,8],:);
u_Y_ref = u_ref([2],:);
%%
% Create ode with integrator This ODE is in form x_dot = Ax + B 
% But form x_plus = Ax+ Bu could also be used. Makes no difference for
% runtime here

[A_noL B_noL] = linearize_system([0,0,l_bar,0,0,0,0,0,0,0]',[0,0,0]');
sys_cont = ss(A_noL, B_noL,zeros(2,size(A_noL,1)),zeros(2, size(B_noL, 2)));
sys_disc = c2d(sys_cont, timestep, 'zoh');
A_d = sys_disc.A;
B_d = sys_disc.B;
C_d = sys_disc.C;
D_d = sys_disc.D;
% Construct 2 decoupled systems
indices1 = [1, 2, 3, 4];
indices2 = [5, 6, 7, 8];

AX = A_d(indices1, indices1);
AY = A_d(indices2, indices2);

BX = B_d(indices1, 1);
BY = B_d(indices2, 2);
ode_fun_X = @(x,u) [AX*x + BX*u];
ode_fun_Y = @(x,u) [AY*x + BY*u];
%%

% initial state as parameter
x0_p = ocp.parameter(nx,1);

% control input as decision variable
u_dec = ocp.variable(nu,  nSteps+foh);
u = u_dec;

% apply dynamic constraint
idx = 1;

% Here, full discretization is used. This simplifies a specific step 
% later on the algorithm. The slightly increased computation time 
% is negligable since all optimal solutions are only calcuatled once and
% then stored in the hashmap
x_dec = casadi.MX(nx,nSteps-1);
x=[x0_p, x_dec]; 
for Node = 1:nSteps-1
    x(:,Node+1) = ode_fun_X(x(:,Node), u(:, Node));
end

%%
% Define all components of cost function
ref_tracking_error = 0;
ref_tracking_u_error = 0;

Q_w = wQ.*eye(4);
R_w = eye(1);

ref_tracking_error = wQ(1)*sum(sum((x(1,:)-x_X_ref(1,:)).^2));
ref_tracking_error =ref_tracking_error + wQ(2)*sum(sum((x(2,:)-x_X_ref(2,:)).^2));
ref_tracking_error =ref_tracking_error + wQ(3)*sum(sum((x(3,:)-x_X_ref(3,:)).^2));
ref_tracking_error =ref_tracking_error + wQ(4)*sum(sum((x(4,:)-x_X_ref(4,:)).^2));
ref_tracking_u_error =sum(sum((u(1,:)-u_X_ref(1,:)).^2));

% Additional cost to keep the movement of the pendulum relativ to the cart
% minimal (penalize any angle/angle derivative)
cost_stabilize =  sum(sum(x([2,4],:).^2));

cost =   w(1)*ref_tracking_u_error + w(2)*cost_stabilize   + w(3)*ref_tracking_error;
ocp.minimize(cost);


ocp.subject_to(-uMax <= u <= uMax);

% Apply admissible set constraints
lower_bound = [0, -pi/5, -0.5, -pi/2.5];
upper_bound = [0.5,  pi/5, 0.5,  -pi/2.5];
ocp.subject_to(lower_bound'<= x < upper_bound');

% construct function for object handling
f = ocp.to_function('f1',{ocp.p},{full(full(u_dec)),full(cost)});
opts.ipopt.print_level = 0;
ocp.solver('ipopt',opts);


% ipopt as solver for backup of qrqp fails
f_ipopt = ocp.to_function('f1',{ocp.p},{full(full(u_dec)),full(cost)});
cost_function = Function('f',{ocp.p,u_dec},{full(ocp.f)});
constraint_function_left= Function('leq',{ocp.p,u_dec},{ocp.lbg});
constraint_function_g = Function('g',{ocp.p,u_dec},{ocp.g});
constraint_function_right= Function('laq',{ocp.p,u_dec},{ocp.ubg});
end


function evaluated_sines = evaluateSines(params, t)
    % This function evaluates sine waves given sine parameters and a time vector t
    % params - a matrix where each row contains the parameters [Amplitude, Frequency, Phase, Vertical Shift] for a sine wave
    % t - a scalar or vector of time points at which to evaluate the sine waves

    % Initialize an array to hold the evaluated sines
    num_trajectories = size(params, 1); % Number of trajectories (rows of params)
    evaluated_sines = casadi.MX(num_trajectories, numel(t));

    % Loop over each set of parameters (each trajectory)
    for j = 1:num_trajectories
        % Extract the parameters for this trajectory
        amplitude = params(j, 1);
        frequency = params(j, 2);
        phase = params(j, 3);
        vertical_shift = params(j, 4);

        % Evaluate the sine wave at each time point in t
        evaluated_sines(j, :) = amplitude * sin(frequency * t + phase) + vertical_shift;
    end

end
