clc;clear
dbstop if error
% Generate model
model= gen_model(0.75,10);
% Monte Carlo simulations
numTrial = 100;
% GOSPA parameters
gospa_p= 1;
gospa_c= 100;
gospa_alpha= 2;

gospa_vals= zeros(100,4,numTrial);

parfor trial = 1:numTrial
    % Generate ground truth
    truth= gen_truth(model);
    
    % Generate measurements
    meas=  gen_meas(model,truth);
    
    % State dimension
    dim = model.x_dim;
    % Multi-Bernoulli representation
    n = 0;
    r = zeros(0,1);
    x = zeros(dim,n);
    P = zeros(dim,dim,n);
    
    % Unknown target PPP parameters
    lambdau = model.lambdab;
%     lambdau = model.lambdau;
    xu = model.xb;
    Pu = model.Pb;
    
    % Loop through time
    for t = 1:100
        
        % Predict
        [r,x,P,lambdau,xu,Pu] = predict(r,x,P,lambdau,xu,Pu,model);
        
        % Update
        [lambdau,xu,Pu,r,x,P,x_est] = updating(lambdau,xu,Pu,r,x,P,meas.Z{t},model);
        
        % Performance evaluation using GOSPA metric
        [gospa_vals(t,:,trial)] = gospa_dist(get_comps(truth.X{t},[1 3]),...
            get_comps(x_est,[1 3]),gospa_c,gospa_p,gospa_alpha);
    end

end

averGospa = mean(gospa_vals,3);
save pmb_10_75 averGospa



