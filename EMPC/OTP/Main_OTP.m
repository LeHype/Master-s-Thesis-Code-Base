%% Implementation of the OTP algorithm, configured for subsystem xX

% Create a results folder with the timestamp
dateTimeString = datestr(now, 'yyyy-mm-dd_HH-MM');
resultsFolder = ['Results/Results_' dateTimeString];
mkdir(resultsFolder);
% Add paths
addpath('Results')
addpath('src')
addpath(['src'])
currentDir = pwd;
parentDir = fileparts(currentDir);
parentDir = fileparts(parentDir);
modelsDir = fullfile(parentDir, 'Construct_OCPs/');
addpath(modelsDir);
% Initialize parallel pool
numCores = 16;
parpool_current =gcp('nocreate');
if isempty(parpool_current)
parpool_current = parpool(numCores);
end
%% Construct f-OCP and f-constraints to be used thoughout. 

import casadi.*
global f f_ipopt cost_function constraint_function_left constraint_function_g constraint_function_right
[f, f_ipopt,cost_function,constraint_function_left,constraint_function_g,constraint_function_right] = Construct_OCP_X();
%%
terminalLeaves = cell({}); % initialize terminal leaves to be empty
lower_bound = [0,0.05, -pi/6, -0.45, -pi/2.5];
upper_bound = [2*pi/3.5,0.55,  pi/6, 0.45, pi/2.5];

paramSpace = [lower_bound;upper_bound]';

% sample for tree leave structure. I am avare this is an ineffcient
% datastructure. Due to the long time needed for calculation this was not
% refactored and re-run. It was also confirmed, that 95 percent of the time
% to calcualte came from the NLP to find K0 and h.
treeRoot = struct(...
    'paramSpace', paramSpace, ...
    'SplitDimension', 3, ...
    'midpoint', NaN, ...
    'K0', [], ...
    'h0', [],...
    'Message',[],...
    'epsilon',[],...
    'terminationFlag', false, ...
    'feasibility',0,...
    'leftChild', [], ...
    'rightChild', []);
leafNodes = {'treeRoot'};
%% convert paramspaces from leaf nodes to array for parpool
iteration = 0;
while(~isempty(leafNodes) || iteration >= 19) % Stop after 20 splits
    iteration = iteration +1;
    num_terminatedNodes = 0;
    num_splitNodes = 0;
    lower_bound_array = zeros(length(leafNodes),length(paramSpace));
    upper_bound_array = zeros(length(leafNodes),length(paramSpace));
    K0_array = cell(length(leafNodes),1);
    h0_array = cell(length(leafNodes),1);
    epsilon_array = cell(length(leafNodes),1);
    Message_array = cell(length(leafNodes),1);
    feasibility_array = cell(length(leafNodes),1);
    nextSplitArray = zeros(length(leafNodes),1);
    for i = 1:length(leafNodes)
        leafpath = leafNodes{i};
        node = eval([leafpath]);
        lower_bound_array(i,:) = node.paramSpace(:,1)';
        upper_bound_array(i,:) = node.paramSpace(:,2)';
    end
    % run the investigation of the hyperrecangles in paralell
    parfor i =  1:length(leafNodes) 
        [K0_sol, h0_sol, epsilon, feasibility,Message,nextSplit] = CheckHyperRectangle(lower_bound_array(i,:),upper_bound_array(i,:));
        K0_array{i} = K0_sol;
        h0_array{i} = h0_sol;
        epsilon_array{i} = epsilon;
        Message_array{i} = Message;
        feasibility_array{i} = feasibility;
        nextSplitArray(i,1) = nextSplit;
    end
    %
    % Update Nodes and get new leaf node paths
    new_leaves = cell({});
    statusReports = cell(length(leafNodes), 1);
    for i = 1:length(leafNodes)
        leafpath = leafNodes{i};
        node = eval([leafpath]);
        eval([ leafpath '.K0 = K0_array{i};']);
        eval([ leafpath '.h0 = h0_array{i};']);
        eval([ leafpath '.epsilon = epsilon_array{i};']);
        eval([ leafpath '.feasibility = feasibility_array{i};']);
        eval([ leafpath '.Message = Message_array{i};']);
        eval([ leafpath '.SplitDimension = nextSplitArray(i);']);
        statusReports{i} = sprintf('Node %d: %s, with epsilon %f', i, Message_array{i}, epsilon_array{i});
        node = eval([leafpath]);
        if iteration < 16
            if  (node.feasibility && node.epsilon < 0.4) || strcmp(node.Message,'Infeasible_region_found')
                eval([ leafpath '.terminationFlag = 1;']);
                terminalLeaves{end+1} = leafpath;
                num_terminatedNodes = num_terminatedNodes+1;
                              
            else
               num_splitNodes = num_splitNodes +1;
                [midpoint, leftChild, rightChild, rootpath_right ,rootpath_left] = partition_Node(leafpath, node);

                eval([ leafpath '.midpoint = midpoint;']);
                eval([ leafpath '.leftChild = leftChild;']);
                eval([ leafpath '.rightChild = rightChild;']);
                new_leaves{end+1} = rootpath_right;
                new_leaves{end+1} = rootpath_left; 
            end
        else
            if ~node.feasibility || node.epsilon <= 0.5 || strcmp(node.Message,'Infeasible_region_found')
                eval([ leafpath '.terminationFlag = 1;']);
                terminalLeaves{end+1} = leafpath;
                num_terminatedNodes = num_terminatedNodes+1;
            else 
                num_splitNodes = num_splitNodes +1;
                [midpoint, leftChild, rightChild, rootpath_right ,rootpath_left] = partition_Node(leafpath, node);

                eval([ leafpath '.midpoint = midpoint;']);
                eval([ leafpath '.leftChild = leftChild;']);
                eval([ leafpath '.rightChild = rightChild;']);
                new_leaves{end+1} = rootpath_right;
                new_leaves{end+1} = rootpath_left;

            end
        end

    end
    dateTimeString = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    logFilename = [resultsFolder filesep 'StatusReport_' dateTimeString '.txt'];
    fileID = fopen(logFilename, 'w');
    header = sprintf('Finished iteration %d with %d nodes. %d were split and %d were terminated. Total number of new nodes to investigate: %d.\n\n', ...
        iteration, numel(leafNodes), num_splitNodes, num_terminatedNodes, numel(new_leaves));
    fprintf(fileID, '%s', header);

    for i = 1:length(statusReports)
        fprintf(fileID, '%s\n', statusReports{i});
    end
    fclose(fileID);
    old_leaves = leafNodes;
    leafNodes = new_leaves;
    dateTimeString = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    savename = [resultsFolder filesep 'Preliminary_Resutls' dateTimeString '.mat'];
    save(savename,'treeRoot',"leafNodes","terminalLeaves","old_leaves");
end
% eval([ leafpath '.K0 = K0_sol;']);
% eval([ leafpath '.h0 = h0_sol;']);
% eval([ leafpath '.epsilon = epsilon;']);
% eval([ leafpath '.feasibility = feasibility;']);
% eval([ leafpath '.Message = Message;']);
%  fprintf(fileID, '%d, %f, %d, %s\n', i, epsilon, feasibility, Message);

%%
% dateTimeString = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
% logFilename = ['logFile_' dateTimeString '.txt']; % Creates a filename like 'logFile_2024-03-29_15-45-30.txt'
% fileID = fopen(logFilename, 'w'); % Opens a new file for writing
%
% for i = 1:length(leafNodes)
% leafpath = leafNodes{i};
% node = eval([leafpath]);
% [K0_sol, h0_sol, epsilon, feasibility,Message] = CheckPolygon(node.paramSpace(:,1)', node.paramSpace(:,2)');
% eval([ leafpath '.K0 = K0_sol;']);
% eval([ leafpath '.h0 = h0_sol;']);
% eval([ leafpath '.epsilon = epsilon;']);
% eval([ leafpath '.feasibility = feasibility;']);
% eval([ leafpath '.Message = Message;']);
%  fprintf(fileID, '%d, %f, %d, %s\n', i, epsilon, feasibility, Message);
% end
%
% fclose(fileID); % Close the log file
% %%
% % Example parameter vector theta
% theta = rand(10,1);
%
% paramSpace = traverseAndPrintLeaf(treeRoot, theta)