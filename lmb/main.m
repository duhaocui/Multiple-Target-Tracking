dbstop if error
numTrial = 100;
model= gen_model2;
% GOSPA parameters
gospa_p= 1;
gospa_c= 100;
gospa_alpha= 2;
gospa_vals= zeros(101,4,numTrial);

parfor trial = 1:numTrial
    truth= gen_truth2(model);
    meas=  gen_meas(model,truth);
    est=   run_filter(model,meas);

    for k=1:101
        [gospa_vals(k,:,trial)] = gospa_dist(get_comps(truth.X{k},[1 3]),...
            get_comps(est.X{k},[1 3]),gospa_c,gospa_p,gospa_alpha);
    end
end

averGospa = mean(gospa_vals,3);
mean(averGospa)