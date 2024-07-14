function x_dot = Linearized_Crane_Model(x,u)
params = struct('Tx', 39.108, 'Tsx_u', 2.7741, 'Y_Gain',9.5952, 'static_friction_force_y', 0.0018091, ...
                'Ty', 19.549, 'Tsy_u', 1.0234, 'X_Gain', 59.065, 'static_friction_force_x', 0.0016816, ...
                'Tz', 539.37, 'Tsz_u', 8.3732, 'Z_Gain', 41.809, 'static_friction_force_z', 13.176, ...
                'FaX', 0.013824, 'FaY', 0.0071167, 'l0', 0.18606, 'mc', 0.468, 'mY', 1.155, 'mX', 3.3550, 'g', 9.81);

% params = struct('Tx', 0, 'Tsx_u', 2.7741, 'Y_Gain',9.5952, 'static_friction_force_y', 0.0018091, ...
%                 'Ty', 19.549, 'Tsy_u', 1.0234, 'X_Gain', 59.065, 'static_friction_force_x', 0.0016816, ...
%                 'Tz', 539.37, 'Tsz_u', 8.3732, 'Z_Gain', 41.809, 'static_friction_force_z', 13.176, ...
%                 'FaX', 0, 'FaY', 0, 'l0', 0.18606, 'mc', 0.468, 'mY', 1.155, 'mX', 3.3550, 'g', 9.81);

% Now use params.mc, params.mw, etc., in your function logic
mL = params.mc;
mT = params.mX; % X
mY = params.mY; % Y
l0 = params.l0 +  0.5478;
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
un = [(u(1)/mT)*X_Gain, (u(2)/(mY))*Y_Gain];
% xT = x(1); 
yT = x(2);  
% l = x(3); 
% theta_x = x(3); 
theta_y = x(4);
% dxT = x(5);
dyT = x(6);
% dl = x(8);
% dtheta_x = x(7);
dtheta_y = x(8);
% un1 = un(1);
un2 = un(2);
% un3 = un(3);
% t2 = 1.0./mL;
% t3 = 1.0./mT;
% t4 = 1.0./mY;
% x_dot = [dxT;dyT;dl;dtheta_x;dtheta_y;t3.*un1+Fax.*dtheta_x.*t3.*(5.0e+1./9.0)-Tx.*dxT.*t3;t4.*un2+Fay.*dtheta_y.*t4.*(5.0e+1./9.0)-Ty.*dyT.*t4;t2.*un3-Tz.*dl.*t2;g.*theta_x.*(-5.0e+1./9.0)-t3.*un1.*(5.0e+1./9.0)-dtheta_x.*(Fax.*t2.*3.08641975308642e+1+Fax.*t3.*3.08641975308642e+1)+Tx.*dxT.*t3.*(5.0e+1./9.0);g.*theta_y.*(-5.0e+1./9.0)-t4.*un2.*(5.0e+1./9.0)-dtheta_y.*(Fay.*t2.*3.08641975308642e+1+Fay.*t4.*3.08641975308642e+1)+Ty.*dyT.*t4.*(5.0e+1./9.0)];
t2 = 1.0./l0;
t4 = 1.0./mL;
t5 = 1.0./mT;
t6 = 1.0./mY;
t3 = t2.^2;
static_friction_force_y = params.static_friction_force_y;
static_friction_force_x = params.static_friction_force_x;
x_dot = [dyT;dtheta_y;t6.*un2-Ty.*dyT.*t6+Fay.*dtheta_y.*t2.*t6;-g.*t2.*theta_y-t2.*t6.*un2+Ty.*dyT.*t2.*t6-Fay.*dtheta_y.*t3.*t4.*t6.*(mL+mY)];
% epsilon = 0.01;
% if (abs(x_dot(5)) <= static_friction_force_x && abs(dxT) <= epsilon)
%     x_dot(1) = 0;
%     x_dot(5)=0;
% end
% 
% if (abs(x_dot(6)) <= static_friction_force_y && abs(dyT) <= epsilon)
%     x_dot(6) = 0;
%     x_dot(2)=0;
% end


end
