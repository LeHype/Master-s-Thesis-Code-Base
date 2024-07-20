function [testLoss] = Experiment1_smallnet(params, layerSize, numLayers)
saving_true = 1;  
X_switch = true;
%% Load Data
    % Load the structured data containing Training and Test sets
    if X_switch
    data = load('S3XX_Training_Data.mat');
    filename = ['S3XX_NL_' num2str(numLayers) 'LS_' num2str(layerSize) '.mat'];
    else
    filename = ['S3YY_NL_' num2str(numLayers) 'LS_' num2str(layerSize) '.mat'];
    
    data = load('S3YY_Training_Data.mat');
    end
    Training = data.Training;
    Test = data.Test;

    %% Define Network Architecture
    % Create the input layer
    layers = [
        featureInputLayer(size(Training.Input, 2), 'Normalization', 'none', 'Name', 'input')
    ];

    % Add fully connected layers according to numLayers and layerSize
    for i = 1:numLayers
        layers = [% Display the best hyperparameters for each configuration
disp('Best Hyperparameters for each configuration:');
disp(bestResults);

            layers,
            fullyConnectedLayer(layerSize, 'Name', ['fc' num2str(i)]),
            reluLayer('Name', ['relu' num2str(i)])
        ];
    end

    % Append output layer
    layers = [
        layers,
        fullyConnectedLayer(size(Training.Target, 2), 'Name', 'output'), % Output layer for regression
        regressionLayer('Name', 'regressionoutput')
    ];


    %% Define Training Options
options = trainingOptions('adam', ...
    'MaxEpochs', 50, ...
    'MiniBatchSize', params.MiniBatchSize, ...
    'InitialLearnRate', params.InitialLearnRate, ...
    'ValidationFrequency', 30, ...
    'Verbose', false, ...
    'Shuffle', 'every-epoch', ...
    'ValidationData', {Test.Input, Test.Target}, ...
    'OutputFcn', @customStopTraining, ... % Add the custom stop function
    'ExecutionEnvironment', 'gpu', ...
    'Plots', 'none');
if saving_true
options = trainingOptions('adam', ...
    'MaxEpochs', 50, ...
    'MiniBatchSize', params.MiniBatchSize, ...
    'InitialLearnRate', params.InitialLearnRate, ...
    'ValidationFrequency', 30, ...
    'Verbose', false, ...
    'Shuffle', 'every-epoch', ...
    'ValidationData', {Test.Input, Test.Target}, ...
    'OutputFcn', @customStopTraining, ... % Add the custom stop function
    'ExecutionEnvironment', 'gpu', ...
    'Plots', 'training-progress');
end
    %% Train Network
    [net, info] = trainNetwork(Training.Input, Training.Target, layers, options);
    
    % Calculate loss on Test Data
    testPred = predict(net, Test.Input);
    testLoss = mean((testPred - Test.Target).^2); % Assuming MSE for regression
    
    % Output the final test loss
    output = struct('FinalTestLoss', testLoss, 'TrainedNet', net);
    testLoss = output.FinalTestLoss;
    
    if saving_true
    
    save(filename,'net',"numLayers","layerSize","options");
    end
end
function stop = customStopTraining(info)
    persistent minLoss numIterationsWithoutImprovement
    
    % Initialize persistent variables during the first iteration
    if info.State == "start"
        minLoss = inf;
        numIterationsWithoutImprovement = 0;
    end

    % Check the validation loss and update the minimum loss found so far
    if ~isempty(info.ValidationRMSE)
        if info.ValidationRMSE < minLoss
            % Update minimum loss and reset the counter
            if minLoss - info.ValidationRMSE > 0.01
                minLoss = info.ValidationRMSE;
                numIterationsWithoutImprovement = 0;
            else
                numIterationsWithoutImprovement = numIterationsWithoutImprovement + 1;
            end
        else
            % Increment the counter if no significant improvement
            numIterationsWithoutImprovement = numIterationsWithoutImprovement + 1;
        end
        
        % Stop condition
        if numIterationsWithoutImprovement >= 300
            disp(['Training will stop. No significant improvement in the last ' num2str(numIterationsWithoutImprovement) ' iterations.']);
            stop = true;
            return;
        end
    end
    
    stop = false;
end