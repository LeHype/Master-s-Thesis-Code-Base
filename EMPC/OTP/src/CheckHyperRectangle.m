function [K0_sol, h0_sol, epsilon, feasibility,Message,idx_to_besplit] = CheckHyperRectangle(lower_bound, upper_bound)
% Initialize feasibility to true
% addpath('C:\Users\Lennart Heib\Downloads\casadi-3.6.5-windows64-matlab2018b');
feasibility = 1;
Message = 'sucess';
K0_sol = [];
h0_sol = [];
epsilon = NaN;
idx_to_besplit = -1;
import casadi.*
persistent f  f_ipopt  cost_function  constraint_function_left  constraint_function_g  constraint_function_right
%due to the paralell computing the function has to be re-constructed for some threads
if isempty(f)
    [f, f_ipopt, cost_function, constraint_function_left, constraint_function_g, constraint_function_right] = YConstruct_OCP_X();

end
n = length(lower_bound);
num_vertices = 2^n;
size_u_Star = [1,20];   % Size of optimal solution

% Preallocate storage for results
u_star_all = zeros(num_vertices, size_u_Star(1)*size_u_Star(2))*NaN;
V_star_all = zeros(num_vertices, 1)*NaN;
vertices_all = zeros(num_vertices, n)*NaN;

% Iterate over each vertex
for i = 0:(num_vertices - 1)
    binaryString = dec2bin(i, n);
    scaled_vertex = lower_bound + (upper_bound - lower_bound) .* (binaryString - '0');
    vertices_all(i+1,:) = scaled_vertex;
    % Get OPC solution. CachedF calcualtes it and stores the solution in a
    % hashmap so it does not have to be recomputed
    [u_star, V_star] = cachedF(scaled_vertex);

    % Store the results. These are the vertex points
    u_star_all(i + 1, :) = full(reshape(u_star,size_u_Star(1)*size_u_Star(2),1));
    V_star_all(i + 1) = full(V_star);
end

if any(isnan(V_star_all)) % Any vertex infeasible
    feasibility = 0;
    Message = 'Infeasibly Vertex found';
end
if ~any(~isnan(u_star_all),'all')
    % All vertices infeasible, possible infeasible region. If midpoint also
    % infeasible --> region is deemed infeasible
    middle_point = lower_bound +((upper_bound -lower_bound) / 2); %% Check midpoint if infeasible
    try
        f(middle_point);
    catch
        Message = 'Infeasible_region_found';
        return;
    end

end
%% Start eplison calculation for region
% Calculate inner rectangle, i.e. test points
middle_point = lower_bound +((upper_bound -lower_bound) / 2);
third_step = (upper_bound-lower_bound)/3;
lower_bound_ = lower_bound+third_step;
upper_bound_ = upper_bound-third_step;
u_star_inner = zeros(num_vertices+1, size_u_Star(1)*size_u_Star(2))*NaN; % Adjust 'size_of_u_star' as necessary
V_star_inner  = zeros(num_vertices+1, 1)*NaN;
vertices_inner = zeros(num_vertices+1, n)*NaN;
% Iterate over each additional vertex theta_test
for i = 0:(num_vertices )
    if i == num_vertices
        scaled_vertex = middle_point;
    else
        binaryString = dec2bin(i, n);
        scaled_vertex = lower_bound_ + (upper_bound_ - lower_bound_) .* (binaryString - '0');
    end
    vertices_inner(i+1,:) = scaled_vertex;
    % Call the function 'cachedF' with the scaled vertex
    [u_star, V_star] = cachedF(scaled_vertex);
    % Store the results these are the Test points
    u_star_inner(i + 1, :) = full(reshape(u_star,size_u_Star(1)*size_u_Star(2),1));
    V_star_inner(i + 1) = full(V_star);
end

% Calcualte Epsilon
epsilon_g = 0.01;
import casadi.*
opti = casadi.Opti();
opti.subject_to();
K0 = opti.variable(size_u_Star(1)*size_u_Star(2),n);
h0 = opti.variable(size_u_Star(1)*size_u_Star(2),1);
cost = 0;
u_cost = 0;
for i = 1:length(V_star_all)
    u_Temp_ = K0*vertices_all(i,:)' + h0;
    u_Temp = reshape(u_Temp_,size_u_Star(1),size_u_Star(2));
    theta_temp = vertices_all(i,:)';
    temp_cost = cost_function(theta_temp,u_Temp) - V_star_all(i) ;
    temp_u_cost = 0.8*norm(u_Temp_-u_star_all(i,:)');
    u_cost = u_cost + temp_u_cost+ u_cost;
    cost = cost+ temp_cost;
    opti.subject_to(constraint_function_left(theta_temp,u_Temp) - constraint_function_g(theta_temp,u_Temp)<= epsilon_g)
    opti.subject_to(constraint_function_g(theta_temp,u_Temp) - constraint_function_right(theta_temp,u_Temp)<=epsilon_g)
    opti.minimize(cost)
end
opts.ipopt.print_level = 0;
opti.solver('ipopt',opts);
% Find K0 and h
try
    opti.solve()
catch
    Message = 'Linear approximation infeasible';
    feasibility = 0;
    epsilon = NaN;

end
estimate_epsilon = zeros(num_vertices*2 +1 ,1)*NaN;
% If approximation can be found do 
if feasibility
    estimate_epsilon = zeros(num_vertices*2 +1 ,1);
    K0_sol = opti.value(K0);
    h0_sol = opti.value(h0);
    %%
    u_hat = @(x) reshape(K0_sol * x+ h0_sol, size_u_Star(1), size_u_Star(2));
    %% Estimate error
    for i = 1:length(vertices_all)
        u_est = u_hat(vertices_all(i,:)');
        estimate_epsilon(i,1) = (u_est(1)-u_star_all(i,1))^2;
    end
    idx = 1;
    for i = length(vertices_all)+1:length(estimate_epsilon)

        u_est = u_hat(vertices_inner(idx,:)');
        estimate_epsilon(i,1) = (u_est(1)-u_star_inner(idx,1))^2;
        idx = idx+1;
    end
    epsilon = sqrt(sum(estimate_epsilon,'all'));

end


%% Next Split accoring to heuristic splitting rule
vertices_combined = [vertices_all ; vertices_inner];
midpoints = lower_bound + (upper_bound-lower_bound)/2;
error_after_split = zeros(n,2);
for i = 1:n
    error_after_split(i,1) = sum(estimate_epsilon((vertices_combined(:,i) < midpoints(i))));
    error_after_split(i,2) = sum(estimate_epsilon((vertices_combined(:,i) >= midpoints(i))));
end
[val,idx_to_besplit] = max(abs(error_after_split(:,1)-error_after_split(:,2)));
% If approximation is infeasible choose random from P without t
if isnan(val)
    feasibility_relevant_idx = [2,3,4,5];
    idx_to_besplit = feasibility_relevant_idx(randi(length(feasibility_relevant_idx))); %% infeasibility comes from x
end


end


function [u_star, V_star] = cachedF(scaled_vertex)
% Persistent cache for storing computation results

persistent cache
if isempty(cache)
    cache = dictionary('init',{[0 0 0 0 0]});
end

key = mat2str(scaled_vertex);
if isKey(cache, key)
    result = cache{key};  % Access using curly braces for cell contents
    u_star = result{1};
    V_star = result{2};
else
    [u_star, V_star] = getF(scaled_vertex);
    cache{key} = {u_star, V_star};  % Store cell array in dictionary
end
end

function [u_star, V_star] = getF(scaled_vertex)
% Calculate OPC solution. If qrqp fails sometime ipopt finds a feasible
% solution.
persistent f  f_ipopt
if isempty(f)
    [f, f_ipopt, ~, ~, ~, ~] = Construct_OCP_X();

end
try
    try
        [u_star, V_star] = f(scaled_vertex);
    catch
        [u_star, V_star] = f_ipopt(scaled_vertex);
        u_star = NaN*u_star;
        V_star = V_star*NaN;
    end
catch
    u_star = NaN;
    V_star = NaN;
end
end

