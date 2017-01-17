clc;clear
dbstop if error
% Generate model
model= gen_model2;
% Monte Carlo simulations
numTrial = 1;
K = 101;
% GOSPA parameters
gospa_p= 1;
gospa_c= 100;
gospa_alpha= 2;
gospa_vals= zeros(K,4,numTrial);

for trial = 1:numTrial
    
    % Generate ground truth
    truth= gen_truth2(model);
    
    % Generate measurements
    meas=  gen_meas(model,truth);
    
    % State dimension
    dim = model.x_dim;
    % Multi-Bernoulli representation
    n = 0;
    r = cell(0,1);
    r{1} = zeros(n,1);
    x = cell(0,1);
    x{1} = zeros(dim,n);
    P = cell(0,1);
    P{1} = zeros(dim,dim,n);
    w_update = 1;
    
    % Unknown target PPP parameters
    lambdau = model.lambdau;
    xu = model.xb;
    Pu = model.Pb;
    
    for t = 1:K
        t
        % Predict
        [r,x,P,lambdau,xu,Pu] = predict(r,x,P,lambdau,xu,Pu,model);
        
        % Update
        [lambdau,xu,Pu,r,x,P,x_est,w_update] = ...
            updating(lambdau,xu,Pu,r,x,P,meas.Z{t},model,w_update);
        
        % GOSPA
        [gospa_vals(t,:,trial)] = ...
            gospa_dist(get_comps(truth.X{t},[1 3]),get_comps(x_est,[1 3])...
            ,gospa_c,gospa_p,gospa_alpha);
    end
end

averGospa = mean(gospa_vals,3);
mean(averGospa)


