numTrial = 1;
model= gen_model;
truth= gen_truth(model);
% GOSPA parameters
gospa_p= 1;
gospa_c= 100;
gospa_alpha= 2;
gospa_vals= zeros(truth.K,4,numTrial);

for trial = 1:numTrial
    meas=  gen_meas(model,truth);
    est=   run_filter(model,meas);
    trial
    
    for k=1:meas.K
        [gospa_vals(k,1,trial), gospa_vals(k,2,trial),...
            gospa_vals(k,3,trial), gospa_vals(k,4,trial)] = ...
            gospa_dist(get_comps(truth.X{k},[1 3]),get_comps(est.X{k},[1 3])...
            ,gospa_c,gospa_p,gospa_alpha);
    end
end

averGospa = mean(gospa_vals,3);
mean(averGospa)