function [derivatives] = crane3DModel_derivation(x, u, varargin)

% Define default values is no parameters are provided. 

persistent params
if isempty(params)
defaultParams = struct('Tx', 39.108, 'Tsx_u', 2.7741, 'Y_Gain',9.5952, 'static_friction_force_y', 0.0018091, ...
                'Ty', 19.549, 'Tsy_u', 1.0234, 'X_Gain', 59.065, 'static_friction_force_x', 0.0016816, ...
                'Tz', 539.37, 'Tsz_u', 8.3732, 'Z_Gain', 41.809, 'static_friction_force_z', 13.176, ...
                'FaX', 0.013824, 'FaY', 0.0071167, 'l0', 0.18606, 'mc', 0.468, 'mY', 1.155, 'mX', 3.3550, 'g', 9.81);
if ~isempty(varargin)
    % Override default values with passed arguments
    params = varargin{1};
    fields = fieldnames(params);
    for i = 1:length(fields)
        defaultParams.(fields{i}) = params.(fields{i});

    end
    params = defaultParams;
else
    params = defaultParams;
end
end

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
Tx_d = params.Tsx_u;
Ty_d = params.Tsy_u;
Tz_d = params.Tsz_u;
X_Gain = params.X_Gain;
Y_Gain = params.Y_Gain;
Z_Gain = params.Z_Gain;
static_friction_force_z = params.static_friction_force_z;
static_friction_force_y = params.static_friction_force_y;
static_friction_force_x = params.static_friction_force_x;



% Input normalization
un = [(u(1)/mT)*X_Gain, (u(2)/(mY))*Y_Gain, (u(3)/mL)*Z_Gain];
un1 = un(1);
un2 = un(2);
un3 = un(3);
t2 = x(6); %xdot
t3 = x(7); %ydot
t4 = x(3)+params.l0; %r
t5 = x(8); %r_d
t6 = x(4); %xa
t7 = mL.^2;
t8 = mL.^3;
t11 = x(9); %xad
t12 = x(5); %ya
t13 = x(10);%yad
t14 = -mL;
t15 = 1.0./mL;
t16 = -mT;
t9 = cos(t6);
t10 = sin(t6);
t17 = t4.^2;
t18 = cos(t12);
t19 = sin(t12);
t20 = 1.0./t4;
t22 = t11.^2;
t23 = t13.^2;
t28 = t2.*3.0e+2;
t29 = t3.*3.0e+2;
t30 = t5.*3.0e+2;
t21 = 1.0./t17;
t24 = t9.^2;
t25 = t9.^3;
t26 = t10.^2;
t27 = t10.^3;
t31 = t18.^2;
t32 = t18.^3;
t34 = t19.^2;
t35 = t19.^3;
t38 = 1.0./t18;
t39 = tanh(t28);
t40 = tanh(t29);
t41 = tanh(t30);
t33 = t31.^2;
t36 = t34.^2;
t37 = mL.*t24;
t42 = mL.*t31;
t43 = mL.*t34;
t44 = t14.*t31;
t45 = t14.*t34;
t46 = t26.*t42;
t47 = t26.*t43;
t48 = mL+mY+t44+t45;
t49 = 1.0./t48;
t50 = t14+t16+t37+t46+t47;
t51 = 1.0./t50;
et1 = t8.*t17.*t23+t4.*t7.*un3-Tz.*t4.*t5.*t7-Tz_d.*t4.*t7.*t41+mL.*mT.*t4.*un3+mL.*mY.*t4.*un3+mT.*mY.*t4.*un3+mT.*t7.*t17.*t23+mY.*t7.*t17.*t23+mT.*t4.*t44.*un3-t8.*t17.*t23.*t24+t8.*t17.*t22.*t31-t8.*t17.*t23.*t31-t8.*t17.*t22.*t33-t8.*t17.*t23.*t34-t4.*t7.*t19.*un2-t4.*t7.*t24.*un3-t4.*t7.*t31.*un3-Fax.*t7.*t9.*t10.*t11-Fay.*t7.*t13.*t18.*t19+Tz.*mT.*t4.*t5.*t14+Tz.*mT.*t4.*t5.*t42+Tz_d.*mT.*t4.*t14.*t41+Tz_d.*mT.*t4.*t41.*t42+Tz.*mY.*t4.*t5.*t14+Tz.*mY.*t4.*t5.*t16+Tz.*mY.*t4.*t5.*t37+Tz.*mY.*t4.*t5.*t47+Tz_d.*mY.*t4.*t14.*t41+Tz_d.*mY.*t4.*t16.*t41+Tz_d.*mY.*t4.*t37.*t41+Tz_d.*mY.*t4.*t41.*t47+Ty.*t3.*t4.*t7.*t19;
et2 = Ty_d.*t4.*t7.*t19.*t40+Tz.*t4.*t5.*t7.*t24+Tz.*t4.*t5.*t7.*t31+Tz_d.*t4.*t7.*t24.*t41+Tz_d.*t4.*t7.*t31.*t41+mL.*mT.*mY.*t17.*t23+mT.*mY.*t17.*t22.*t42+g.*t4.*t8.*t9.*t18-g.*t4.*t8.*t9.*t32-g.*t4.*t8.*t18.*t25+g.*t4.*t8.*t25.*t32+mT.*t7.*t17.*t22.*t31-mY.*t7.*t17.*t23.*t24+mY.*t7.*t17.*t22.*t31+mT.*t4.*t14.*t19.*un2+mY.*t4.*t14.*t24.*un3+mY.*t4.*t26.*t45.*un3+t7.*t16.*t17.*t23.*t31+t7.*t16.*t17.*t22.*t33+t7.*t16.*t17.*t23.*t34-t8.*t17.*t22.*t24.*t31+t8.*t17.*t23.*t24.*t31+t8.*t17.*t22.*t24.*t33-t8.*t17.*t22.*t26.*t31-t8.*t17.*t23.*t26.*t31+t8.*t17.*t22.*t26.*t33+t8.*t17.*t23.*t24.*t34+t8.*t17.*t23.*t26.*t33-t8.*t17.*t23.*t26.*t34+t8.*t17.*t23.*t26.*t36;
et3 = -t8.*t17.*t22.*t31.*t34-t4.*t7.*t10.*t18.*un1+t4.*t7.*t10.*t32.*un1+t4.*t7.*t19.*t24.*un2+t4.*t7.*t24.*t31.*un3-t4.*t7.*t26.*t34.*un3+t4.*t7.*t26.*t35.*un2+g.*mT.*t4.*t7.*t9.*t18+g.*mY.*t4.*t7.*t9.*t18-g.*mY.*t4.*t7.*t18.*t25-g.*t4.*t8.*t9.*t18.*t26+g.*t4.*t7.*t9.*t16.*t32-g.*t4.*t8.*t9.*t18.*t34+g.*t4.*t8.*t9.*t26.*t32+g.*t4.*t8.*t18.*t25.*t34-mY.*t7.*t17.*t22.*t24.*t31-mY.*t7.*t17.*t22.*t26.*t31-mY.*t7.*t17.*t23.*t26.*t31-mY.*t7.*t17.*t23.*t26.*t34+mY.*t4.*t10.*t14.*t18.*un1+t7.*t16.*t17.*t22.*t31.*t34+t8.*t17.*t22.*t24.*t31.*t34+t8.*t17.*t22.*t26.*t31.*t34+t8.*t17.*t23.*t26.*t31.*t34.*2.0+t4.*t7.*t10.*t18.*t34.*un1+t4.*t7.*t19.*t26.*t31.*un2;
et4 = Fay.*mT.*t13.*t14.*t18.*t19+Fax.*mY.*t9.*t10.*t11.*t14+Fax.*t7.*t9.*t10.*t11.*t31+Fax.*t7.*t9.*t10.*t11.*t34+Fay.*t7.*t13.*t18.*t19.*t24+Fay.*t7.*t13.*t18.*t19.*t26+Ty.*mL.*mT.*t3.*t4.*t19+Ty_d.*mL.*mT.*t4.*t19.*t40+Tx.*t2.*t4.*t7.*t10.*t18-Tx.*t2.*t4.*t7.*t10.*t32+Tx_d.*t4.*t7.*t10.*t18.*t39-Tx_d.*t4.*t7.*t10.*t32.*t39-Ty.*t3.*t4.*t7.*t19.*t24-Ty.*t3.*t4.*t7.*t26.*t35-Ty_d.*t4.*t7.*t19.*t24.*t40-Ty_d.*t4.*t7.*t26.*t35.*t40-Tz.*t4.*t5.*t7.*t24.*t31+Tz.*t4.*t5.*t7.*t26.*t34-Tz_d.*t4.*t7.*t24.*t31.*t41+Tz_d.*t4.*t7.*t26.*t34.*t41+Fay.*mL.*mY.*t13.*t18.*t19.*t26+Tx.*mL.*mY.*t2.*t4.*t10.*t18+Tx_d.*mL.*mY.*t4.*t10.*t18.*t39+g.*mL.*mT.*mY.*t4.*t9.*t18-Tx.*t2.*t4.*t7.*t10.*t18.*t34;
et5 = -Tx_d.*t4.*t7.*t10.*t18.*t34.*t39-Ty.*t3.*t4.*t7.*t19.*t26.*t31-Ty_d.*t4.*t7.*t19.*t26.*t31.*t40-g.*mY.*t4.*t7.*t9.*t18.*t26+g.*t4.*t7.*t9.*t16.*t18.*t34+g.*t4.*t8.*t9.*t18.*t26.*t34;
et6 = Fax.*t11.*t14+Fax.*t11.*t16+Fax.*t11.*t46+Fax.*t11.*t47-g.*t4.*t7.*t10.*t18+g.*t4.*t7.*t27.*t32-mT.*t4.*t5.*t11.*t42.*2.0-t4.*t5.*t7.*t11.*t31.*2.0+t4.*t9.*t14.*t18.*un1+t4.*t9.*t10.*t42.*un3+g.*mT.*t4.*t10.*t14.*t18+g.*t4.*t7.*t10.*t24.*t32+g.*t4.*t7.*t18.*t27.*t34+t4.*t5.*t7.*t11.*t24.*t31.*2.0+t7.*t11.*t13.*t17.*t18.*t19.*2.0+t4.*t5.*t7.*t11.*t26.*t33.*2.0-t7.*t9.*t10.*t17.*t22.*t31+t7.*t9.*t10.*t17.*t22.*t33+Tx.*mL.*t2.*t4.*t9.*t18+Tx_d.*mL.*t4.*t9.*t18.*t39+Tz.*t4.*t5.*t9.*t10.*t44+Tz_d.*t4.*t9.*t10.*t41.*t44+Fay.*mL.*t9.*t10.*t13.*t18.*t19+mL.*mT.*t11.*t13.*t17.*t18.*t19.*2.0+g.*t4.*t7.*t10.*t18.*t24.*t34-t7.*t11.*t13.*t17.*t18.*t19.*t24.*2.0+t4.*t5.*t7.*t11.*t26.*t31.*t34.*2.0;
et7 = t7.*t11.*t13.*t17.*t19.*t26.*t32.*-2.0-t7.*t11.*t13.*t17.*t18.*t26.*t35.*2.0+t7.*t9.*t10.*t17.*t22.*t31.*t34;
et8 = -Fay.*t7.*t13.*t18+mT.*t4.*t44.*un2-t4.*t7.*t31.*un2+Fay.*mT.*t13.*t14.*t18+Fay.*mT.*t13.*t18.*t43+Fay.*mY.*t13.*t14.*t18+Fay.*mY.*t13.*t16.*t18+Fay.*mY.*t13.*t18.*t37+Fay.*t7.*t13.*t18.*t24+Fay.*t7.*t13.*t18.*t34+Fay.*t7.*t13.*t26.*t32+Ty.*mT.*t3.*t4.*t42+Ty_d.*mT.*t4.*t40.*t42+Ty.*t3.*t4.*t7.*t31+Ty_d.*t4.*t7.*t31.*t40+mT.*t4.*t19.*t42.*un3-t4.*t5.*t8.*t13.*t18.*2.0+t4.*t5.*t8.*t13.*t32.*2.0-t8.*t17.*t19.*t22.*t31+t8.*t17.*t19.*t22.*t33+t8.*t17.*t22.*t31.*t35+t4.*t7.*t19.*t31.*un3+t4.*t7.*t24.*t31.*un2+t4.*t7.*t26.*t33.*un2+mT.*mY.*t17.*t19.*t22.*t44-g.*t4.*t8.*t9.*t18.*t19+g.*t4.*t8.*t9.*t19.*t32+g.*t4.*t8.*t9.*t18.*t35+g.*t4.*t8.*t18.*t19.*t25-g.*t4.*t8.*t19.*t25.*t32;
et9 = -g.*t4.*t8.*t18.*t25.*t35-mT.*t4.*t5.*t7.*t13.*t18.*2.0+mT.*t4.*t5.*t7.*t13.*t32.*2.0+mT.*t7.*t17.*t19.*t22.*t33+mT.*t7.*t17.*t22.*t31.*t35-mY.*t4.*t5.*t7.*t13.*t18.*2.0-mY.*t7.*t17.*t19.*t22.*t31+mY.*t4.*t19.*t26.*t44.*un3+t4.*t5.*t8.*t13.*t18.*t24.*2.0+t4.*t5.*t8.*t13.*t18.*t34.*2.0-t4.*t5.*t8.*t13.*t24.*t32.*2.0+t4.*t5.*t8.*t13.*t26.*t32.*2.0+t7.*t16.*t17.*t19.*t22.*t31+t8.*t17.*t19.*t22.*t24.*t31-t8.*t17.*t19.*t22.*t24.*t33+t8.*t17.*t19.*t22.*t26.*t31-t8.*t17.*t19.*t22.*t26.*t33-t8.*t17.*t22.*t24.*t31.*t35-t8.*t17.*t22.*t26.*t31.*t35+t4.*t7.*t10.*t18.*t19.*un1-t4.*t7.*t10.*t19.*t32.*un1-t4.*t7.*t10.*t18.*t35.*un1-t4.*t7.*t19.*t24.*t31.*un3-t4.*t7.*t19.*t26.*t31.*un3;
et10 = t4.*t7.*t26.*t31.*t34.*un2-t4.*t5.*t8.*t13.*t26.*1.0./t38.^5.*2.0+Fay.*mL.*mY.*t13.*t26.*t32+Fax.*t7.*t9.*t10.*t11.*t19-Fax.*t7.*t9.*t10.*t11.*t35-Fay.*t7.*t13.*t18.*t24.*t34+Tz.*mT.*t4.*t5.*t19.*t44+Tz_d.*mT.*t4.*t19.*t41.*t44+Tz.*mY.*t4.*t5.*t19.*t46+Tz_d.*mY.*t4.*t19.*t41.*t46-Ty.*t3.*t4.*t7.*t24.*t31-Ty.*t3.*t4.*t7.*t26.*t33-Ty_d.*t4.*t7.*t24.*t31.*t40-Ty_d.*t4.*t7.*t26.*t33.*t40-Tz.*t4.*t5.*t7.*t19.*t31-Tz_d.*t4.*t7.*t19.*t31.*t41+Fax.*mL.*mY.*t9.*t10.*t11.*t19-Fax.*t7.*t9.*t10.*t11.*t19.*t31-Tx.*t2.*t4.*t7.*t10.*t18.*t19+Tx.*t2.*t4.*t7.*t10.*t19.*t32+Tx.*t2.*t4.*t7.*t10.*t18.*t35-Tx_d.*t4.*t7.*t10.*t18.*t19.*t39+Tx_d.*t4.*t7.*t10.*t19.*t32.*t39;
et11 = Tx_d.*t4.*t7.*t10.*t18.*t35.*t39-Ty.*t3.*t4.*t7.*t26.*t31.*t34-Ty_d.*t4.*t7.*t26.*t31.*t34.*t40+Tz.*t4.*t5.*t7.*t19.*t24.*t31+Tz.*t4.*t5.*t7.*t19.*t26.*t31+Tz_d.*t4.*t7.*t19.*t24.*t31.*t41+Tz_d.*t4.*t7.*t19.*t26.*t31.*t41-mL.*mT.*mY.*t4.*t5.*t13.*t18.*2.0+g.*mT.*t4.*t7.*t9.*t19.*t32+g.*mT.*t4.*t7.*t9.*t18.*t35-g.*mY.*t4.*t7.*t9.*t18.*t19+g.*mY.*t4.*t7.*t18.*t19.*t25+g.*t4.*t7.*t9.*t16.*t18.*t19+g.*t4.*t8.*t9.*t18.*t19.*t26-g.*t4.*t8.*t9.*t19.*t26.*t32-g.*t4.*t8.*t9.*t18.*t26.*t35+mL.*mY.*t4.*t10.*t18.*t19.*un1+mT.*t4.*t5.*t7.*t13.*t18.*t34.*2.0+mY.*t4.*t5.*t7.*t13.*t18.*t24.*2.0+mY.*t4.*t5.*t7.*t13.*t26.*t32.*2.0+mY.*t7.*t17.*t19.*t22.*t24.*t31+mY.*t7.*t17.*t19.*t22.*t26.*t31;
et12 = t4.*t5.*t8.*t13.*t18.*t24.*t34.*-2.0+t4.*t5.*t8.*t13.*t18.*t26.*t34.*2.0-t4.*t5.*t8.*t13.*t18.*t26.*t36.*2.0-t4.*t5.*t8.*t13.*t26.*t32.*t34.*4.0+Tx.*mY.*t2.*t4.*t10.*t14.*t18.*t19+Tx_d.*mY.*t4.*t10.*t14.*t18.*t19.*t39+g.*mT.*mY.*t4.*t9.*t14.*t18.*t19+g.*mY.*t4.*t7.*t9.*t18.*t19.*t26+mY.*t4.*t5.*t7.*t13.*t18.*t26.*t34.*2.0;
derivatives = [t2;t3;t5;t11;t13;t20.*t38.*t51.*(-Fax.*t9.*t11-t4.*t18.*un1+Tx.*t2.*t4.*t18+Tx_d.*t4.*t18.*t39+t10.*t17.*t22.*t44+t4.*t10.*t31.*un3+Fay.*t10.*t13.*t18.*t19-Tz.*t4.*t5.*t10.*t31-Tz_d.*t4.*t10.*t31.*t41+mL.*t10.*t17.*t22.*t33+t10.*t17.*t22.*t34.*t42+g.*mL.*t4.*t9.*t10.*t32+g.*t4.*t9.*t10.*t14.*t18+g.*t4.*t9.*t10.*t18.*t43);t20.*t49.*(t4.*un2+Fay.*t13.*t18-Ty.*t3.*t4-Ty_d.*t4.*t40-t4.*t19.*un3+Tz.*t4.*t5.*t19+Tz_d.*t4.*t19.*t41);-t15.*t20.*t49.*t51.*(et1+et2+et3+et4+et5);-(t15.*t21.*t51.*(et6+et7))./t31;-t15.*t21.*t38.*t49.*t51.*(et8+et9+et10+et11+et12)]';

epsilon = 0.005; % Numeric treshhold for static dry friction
if (abs(derivatives(8)) <= static_friction_force_z && abs(t5) <= epsilon)
    derivatives(3) = 0;
    derivatives(8)=0;
end

if (abs(derivatives(7)) <= static_friction_force_y && abs(t3) <= epsilon)
    derivatives(2) = 0;
    derivatives(7)=0;
end

if (abs(derivatives(6)) <= static_friction_force_x && abs(t2) <= epsilon)
    derivatives(1) = 0;
    derivatives(6)=0;

end


end

