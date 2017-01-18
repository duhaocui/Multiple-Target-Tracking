clc;clear
dbstop if error

% Generate model
model = gen_model;

lambdau = model.lambdab;
xu = model.xb;
Pu = model.Pb;

% Monte Carlo simulations
numTrial = 1;
K = 100; % simulation time steps

% GOSPA parameters
gospa_p= 1;
gospa_c= 100;
gospa_alpha= 2;
gospa_vals= zeros(100,4,numTrial);

for trial = 1:numTrial
    
    % Generate ground truth
    truth = gen_truth(model);
    
    % Generate measurements
    meas =  gen_meas(model,truth);
    
    % Multi-Bernoulli representation
    dim = model.x_dim;
    n = 0;
    r = zeros(0,1);
    x = zeros(dim,n);
    P = zeros(dim,dim,n);
    
    % Loop through time
    for t = 1:K
        t
        % Predict step
        [r,x,P,lambdau,xu,Pu] = predict(r,x,P,lambdau,xu,Pu,model);
        
        % Update step
        [lambdau,xu,Pu,r,x,P,x_est] = updating(lambdau,xu,Pu,r,x,P,meas.Z{t},model);
        
        % Performance evaluation using GOSPA metric
        [gospa_vals(t,1,trial), gospa_vals(t,2,trial),...
            gospa_vals(t,3,trial), gospa_vals(t,4,trial)] = ...
            gospa_dist(get_comps(truth.X{t},[1 3]),get_comps(x_est,[1 3])...
            ,gospa_c,gospa_p,gospa_alpha);
    end
    
end

averGospa = mean(gospa_vals,3);
mean(averGospa)

