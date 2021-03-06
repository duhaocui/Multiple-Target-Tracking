clc;clear
dbstop if error
% Generate model
model= gen_model(0.98,10);

% Monte Carlo simulations
numTrial = 50;
K = 100; % time steps

% GOSPA parameters
gospa_p= 1;
gospa_c= 100;
gospa_alpha= 2;
gospa_vals= zeros(K,4,numTrial);

parfor trial = 1:numTrial
    %Generate ground truth
    truth= gen_truth(model);
    
    % Generate measurements
    meas=  gen_meas(model,truth);
    
    % State dimension
    dim = model.x_dim;
    % Multi-Bernoulli representation
    r = cell(0,1);
    r{1} = zeros(0,1);
    x = cell(0,1);
    x{1} = zeros(dim,0);
    P = cell(0,1);
    P{1} = zeros(dim,dim,0);
    w_update = 1;
    
    % Unknown target PPP parameters
    lambdau{1} = model.lambdab;
    xu{1} = model.xb;
    Pu{1} = model.Pb;
    
    % Loop through time
    for t = 1:K
        % Predict
        [r,x,P,lambdau,xu,Pu] = predict(r,x,P,lambdau,xu,Pu,model);
        
        % Update
        [lambdau,xu,Pu,r,x,P,x_est,w_update] = ...
            updating(lambdau,xu,Pu,r,x,P,meas.Z{t},model,w_update);
        
        % Performance evaluation using GOSPA metric
        [gospa_vals(t,:,trial)] = gospa_dist(get_comps(truth.X{t},[1 3]),...
            get_comps(x_est,[1 3]) ,gospa_c,gospa_p,gospa_alpha);
    end 
end

averGospa = mean(gospa_vals,3);



