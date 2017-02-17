clc;clear
dbstop if error
averGospa = zeros(100,4,4);

parfor i = 1:4
    if i == 1
        pd = 0.98;lambda = 10;
    elseif i == 2
        pd = 0.75;lambda = 10;
    elseif i == 3
        pd = 0.98;lambda = 30;
    elseif i == 4
        pd = 0.75;lambda = 30;
    end
    
    numTrial = 100;
    model = gen_model(pd,lambda);
    % GOSPA parameters
    gospa_p = 1;
    gospa_c = 100;
    gospa_alpha = 2;
    gospa_vals = zeros(100,4,numTrial);
    
    for trial = 1:numTrial
        truth = gen_truth(model);
        meas = gen_meas(model,truth);
        est = run_filter(model,meas);
        
        for k=1:100
            [gospa_vals(k,:,trial)] = gospa_dist(get_comps(truth.X{k},[1 3]),...
                get_comps(est.X{k},[1 3]),gospa_c,gospa_p,gospa_alpha);
        end
    end
    
    averGospa(:,:,i) = mean(gospa_vals,3);
    
%     if i == 1
%         save lmb_merge_coal_10_75 averGospa
%     elseif i == 2
%         save lmb_merge_coal_10_75 averGospa
%     elseif i == 3
%         save lmb_merge_coal_30_98 averGospa
%     elseif i == 4
%         save lmb_merge_coal_30_75 averGospa
%     end
    
end