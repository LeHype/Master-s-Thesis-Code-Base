function [A_noL B_noL] = linearize_system_A(x0,u0)
% x0 = [0,0,pi/4,pi/4,0,0,-2,-2]';
% u0 = [0,0,0]';
persistent VF Jx Ju x u
syms xT yT  l theta_x theta_y
syms dxT dyT dzT dl dtheta_x dtheta_y
syms  un1 un2 un3
% params = struct('Tx', 39.108, 'Tsx_u', 2.7741, 'Y_Gain',9.5952, 'static_friction_force_y', 0.0018091, ...
%                 'Ty', 19.549, 'Tsy_u', 1.0234, 'X_Gain', 59.065, 'static_friction_force_x', 0.0016816, ...
%                 'Tz', 539.37, 'Tsz_u', 8.3732, 'Z_Gain', 41.809, 'static_friction_force_z', 13.176, ...
%                 'FaX', 0.013824, 'FaY', 0.0071167, 'l0', 0.18606 + 0.547, 'mc', 0.468, 'mY', 1.155, 'mX', 3.3550, 'g', 9.81);

params = struct('Tx', 39.108, 'Tsx_u', 2.7741, 'Y_Gain',9.5952, 'static_friction_force_y', 0.0018091, ...
                'Ty', 19.549, 'Tsy_u', 1.0234, 'X_Gain', 59.065, 'static_friction_force_x', 0.0016816, ...
                'Tz', 539.37, 'Tsz_u', 8.3732, 'Z_Gain', 41.809, 'static_friction_force_z', 13.176, ...
                'FaX', 0.013824, 'FaY', 0.0071167, 'l0', 0.18606, 'mc', 0.468, 'mY', 1.155, 'mX', 3.3550, 'g', 9.81);

mL = params.mc;
mT = params.mX; % X
mY = params.mY; % Y
l0 = params.l0;
Fax = params.FaX;
Fay = params.FaY;
g = params.g;
Tx = params.Tx;
Ty = params.Ty;
Tz = params.Tz;
Txs_u = params.Tsx_u;
Tys_u = params.Tsy_u;
Tzs_u = params.Tsz_u;
X_Gain = params.X_Gain;
Y_Gain = params.Y_Gain;
Z_Gain = params.Z_Gain;
static_friction_force_z = params.static_friction_force_z;
static_friction_force_y = params.static_friction_force_y;
static_friction_force_x = params.static_friction_force_x;

if isempty(VF)
    disp('linearizeing ...')
%Rotation x first
xL = xT +l*sin(theta_x)*cos(theta_y);
yL = yT + l*sin(theta_y);
zL = -l*cos(theta_x)*cos(theta_y);
%Rotation y first
% xL = xT + l*sin(theta_x);
% yL = yT + l*sin(theta_y);
% zL = -l*cos(theta_x)*cos(theta_y);
r_L = [xL; yL; zL];

% Step 1: Define time-dependent temporary variables
syms xT_temp(t) yT_temp(t) l_temp(t) theta_x_temp(t) theta_y_temp(t)

% Step 2: Construct r_L using the temporary variables
% xL_temp = xT_temp + l_temp*sin(theta_x_temp)*cos(theta_y_temp);
% yL_temp = yT_temp + l_temp*sin(theta_y_temp)*cos(theta_x_temp);
% zL_temp = -l_temp*cos(theta_x_temp)*cos(theta_y_temp);
xL_temp = xT_temp + l_temp*sin(theta_x_temp)*cos(theta_y_temp);
yL_temp = yT_temp + l_temp*sin(theta_y_temp);
zL_temp = -l_temp*cos(theta_x_temp)*cos(theta_y_temp);
r_L_temp = [xL_temp; yL_temp; zL_temp];

% Step 3: Differentiate r_L_temp with respect to t
dr_L_temp = diff(r_L_temp, t);

% Step 4: Replace the temporary variables with the respective symbols

% Perform the replacements
dr_L = subs(dr_L_temp, [diff(xT_temp(t), t), diff(yT_temp(t), t), diff(l_temp(t), t), diff(theta_x_temp(t), t), diff(theta_y_temp(t), t)], [dxT, dyT, dl, dtheta_x, dtheta_y]);
dr_L = subs(dr_L, [xT_temp(t), yT_temp(t), l_temp(t), theta_x_temp(t), theta_y_temp(t)], [xT, yT, l, theta_x, theta_y]);



% Now use the new derivative symbols in the kinetic energy T
T = (1/2)*mT*dxT^2 + (1/2)*mY*dyT^2 + (1/2)*mL*dr_L.'*dr_L;

% Potential energy V remains the same as zL is not differentiated
V = mL*g*zL;


L   = T - V;
X = {xT dxT yT dyT l dl theta_x dtheta_x theta_y dtheta_y};
Q_i = {-dxT*Tx - 0*tanh(300*dxT)*Txs_u,-dyT*Ty - 0*tanh(300*dyT)*Tys_u,-dl*Tz - 0*tanh(300*dl)*Tzs_u,-dtheta_x*Fax,-dtheta_y*Fay};
Q_e = {un1,un2,un3,0,0};
R   = 0;
par = {mT mY mL g  Tx Txs_u Ty Tys_u Tz Tzs_u Fax Fay};
% VF  = EulerLagrange(L,X,Q_i,Q_e,R,par,'m','s','GantryCraneODE');
VF  = EulerLagrange(L,X,Q_i,Q_e,R,par);
x = transpose([xT dxT yT dyT l dl theta_x dtheta_x theta_y dtheta_y]);

u = transpose([un1 un2 un3]);

    % Calculate the Jacobians
    Jx = jacobian(VF, x);
    Ju = jacobian(VF, u);
end


   

 X_correct_order = transpose([xT yT l theta_x theta_y dxT dyT dl dtheta_x dtheta_y]) ;
    A_sub = double(subs(Jx, [X_correct_order; u], [x0; u0]));
    B_sub = double(subs(Ju, [X_correct_order; u], [x0; u0]));
x = transpose([xT dxT yT dyT l dl theta_x dtheta_x theta_y dtheta_y]);
 
newOrder = [1, 7, 2, 8, 3, 9, 4,10,5,6];
n = length(newOrder); % Assuming 'n' is the size of your state vector
% Create the permutation matrix P
P = eye(n); % Start with an identity matrix of size 'n'
P = P(newOrder, :); % Reorder rows according to 'newOrder'
% Apply the permutation
A_new = P * A_sub * inv(P); % Reorder A
B_new = P * B_sub; % Reorder B
    
B_new(:,1) = B_new(:,1)*(X_Gain/mT);
B_new(:,2) = B_new(:,2)*(Y_Gain/mY);
% Substitute the linearization point and default input
% A_noL = A_new(1:8,1:8);
% B_noL = B_new(1:8,1:2); 
A_noL = A_new;
B_noL = B_new;
% X_correct_order = transpose([xT yT theta_x theta_y dxT dyT dtheta_x dtheta_y])
%     A_lin = double(subs(A_new, [X_correct_order; u], [x0; u0]));
%     B_lin = double(subs(B_new, [X_correct_order; u], [x0; u0]));

end
