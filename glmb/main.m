clc;clear
dbstop if error
% averGospa = zeros(100,4,4);

% parfor i = 1:4
%     if i == 1
%         pd = 0.98;lambda = 10;
%     elseif i == 2
%         pd = 0.75;lambda = 10;
%     elseif i == 3
%         pd = 0.98;lambda = 30;
%     elseif i == 4
%         pd = 0.75;lambda = 30;
%     end
    
    numTrial = 20;
    % GOSPA parameters
    gospa_p= 1;
    gospa_c= 100;
    gospa_alpha= 2;
    gospa_vals= zeros(100,4,numTrial);
    model= gen_model(0.75,30);
    
    
    parfor trial = 1:numTrial
        truth= gen_truth(model);
        meas=  gen_meas(model,truth);
        est=   run_filter(model,meas);
        
        for k=1:100
            [gospa_vals(k,:,trial)] = ...
                gospa_dist(get_comps(truth.X{k},[1 3]),get_comps(est.X{k},[1 3])...
                ,gospa_c,gospa_p,gospa_alpha);
        end
    end
    averGospa = mean(gospa_vals,3);
    
    % averGospa = mean(gospa_vals,3);
    % save glmb_30_75_3000 averGospa
% end