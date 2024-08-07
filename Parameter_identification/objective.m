function [cost] = objective(params, training,plotting)
arguments
    params struct = struct() % Hand over candidate parameters
    training logical = true  % Specify whether to use training (true) or test data (false)
    plotting logical = false % Plot progress of simulation
end


load('Parameter_identification_meneuvers.mat');

Ts_plant = 0.01; % Constant due to Hardware
RMSEIDX = 1:10;  % Minimize error on these states
rng(3);          % Set the seed for reproducibility

% Draw 60 unique indices from 1 to 200 for test/training split
indices = randperm(size(Meassurement, 1), floor(size(Meassurement, 1) * 0.3));

X_Estimate = [];
X_Val = [];
Val_PWM = [];
% Test if models are in path
try
    crane_Model_non_linear(zeros(10,1),zeros(3,1));
catch
    error('Please add models to path')
end
for kk = 1:size(Meassurement, 1)
    % Exclude problematic measurements and manage training/testing data split
    if kk == 60 || kk == 152 || (~ismember(kk, indices) && ~training) || (ismember(kk, indices) && training)
        continue
    end

    % Extract current scenario
    rt_OnlyPlant = Meassurement{kk, 1}{1};
    rt_TrueGaussNoise = Meassurement{kk, 2}{1};
    X_True = rt_OnlyPlant.signals.values(:,:);
    X_Val = [X_Val; X_True];

    % Input Signals
    PWM_Plant = [rt_TrueGaussNoise.signals.values];

    % Time-step for integration
    Ts = 0.0002;
    T_begin = rt_OnlyPlant.time(1);
    T_end = rt_OnlyPlant.time(end);

    try
        F = griddedInterpolant(rt_TrueGaussNoise.time, PWM_Plant);
    catch
        continue
    end

    Val_PWM = [Val_PWM; PWM_Plant];

    X_0 = rt_OnlyPlant.signals.values(1, :);
    clear crane3DModel_derivation
    options= odeset('OutputFcn',@odeplot,'RelTol',1e-4);
    if plotting == false
        options.OutputFcn = '';
    end

    [time_ode23t, X_Estimate_ode23t] = ode45(@(t, y) crane3DModel_derivation(y, [F(t)]', params)', T_begin:Ts_plant:T_end, X_0,options);
    X_Estimate = [X_Estimate; X_Estimate_ode23t];
end

% Standardize the data
X_Val_std = (X_Val - mean(X_Val)) ./ std(X_Val);
X_Estimate_std = (X_Estimate - mean(X_Val)) ./ std(X_Val);

% Additional lowpass filter because of mesurement artifacts (spikes)
% These artifacts were only noticed after the measurements were taken

FC = 60;
X_Val_std(:,6) = lowpass(X_Val_std(:,6),0.01,FC);
X_Val_std(:,7) = lowpass(X_Val_std(:,7),0.01,FC);
X_Val_std(:,9) = lowpass(X_Val_std(:,9),0.01,FC);
X_Val_std(:,10) = lowpass(X_Val_std(:,10),0.01,FC);

% Compute RMSE
Error = rmse(X_Val_std(1:length(X_Estimate_std), RMSEIDX), X_Estimate_std(:, RMSEIDX), 'omitnan');
cost = sum(Error) / length(Error);


end
