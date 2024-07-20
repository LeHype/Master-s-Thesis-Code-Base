% Define discrete layer sizes and number of layers
layerSizes = [8 16 32 64 128];
numLayersOptions = [1 2 3];

% Initialize a table to store the best hyperparameters for each configuration
bestResults = table();

% Loop over number of layers
for numLayers = numLayersOptions
    % Loop over each layer size
    for layerSize = layerSizes
        % % Define the optimization variables for learning rate and batch size
        vars = [
            optimizableVariable('InitialLearnRate', [1e-4, 1e-2], 'Transform', 'log'), ...
            optimizableVariable('MiniBatchSize', [64, 512], 'Type', 'integer')
        ];
        initialGuess = table(0.00012559, 99, 'VariableNames', {'InitialLearnRate', 'MiniBatchSize'});

        % Objective function for Bayesian optimization
        % Pass layer size and number of layers as additional parameters
        objFcn = @(params) Experiment1_smallnet(params, layerSize, numLayers);

        % Perform Bayesian optimization
        results = bayesopt(objFcn, vars, ...
            'AcquisitionFunctionName', 'expected-improvement-plus', ...
            'MaxObjectiveEvaluations', 10, ...
            'Verbose', 0, ...
            'UseParallel', false,InitialX=initialGuess);

        % Extract the best hyperparameters

        bestParams = bestPoint(results);
         [val_worst worstResultsIndex] = max(results.ObjectiveTrace);
         [val_best best_idx] = min(results.ObjectiveTrace);
          
         worstParams = results.XTrace(worstResultsIndex, :);

      
        % Add results to the table
        bestResults = [bestResults; table(numLayers, layerSize, bestParams.InitialLearnRate, bestParams.MiniBatchSize, results.MinObjective, 'VariableNames', {'NumLayers', 'LayerSize', 'LearnRate', 'BatchSize', 'MinLoss'})];
    % Print the summary of this iteration
        fprintf('Test with layer size %d, and %d number of layers complete.\n', layerSize, numLayers);
        fprintf('Best result, %.4f,  with learning rate %.4f, and Batch size %d\n',val_best, bestParams.InitialLearnRate, bestParams.MiniBatchSize);
        fprintf('Worst result, %.4f, with learning rate %.4f, and Batch size %d\n',val_worst, worstParams.InitialLearnRate, worstParams.MiniBatchSize);

    end
end
