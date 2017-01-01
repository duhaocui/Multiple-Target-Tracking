function truth= gen_truth(model)

%variables
truth.K= 100;                   %length of data/number of scans
truth.X= cell(truth.K,1);             %ground truth for states of targets  
truth.N= zeros(truth.K,1);            %ground truth for number of targets
truth.L= cell(truth.K,1);             %ground truth for labels of targets (k,i)
truth.track_list= cell(truth.K,1);    %absolute index target identities (plotting)
truth.total_tracks= 0;          %total number of appearing tracks

%target initial states and birth/death times
nbirths= 12;
xstart = zeros(4,nbirths);
for i = 1:nbirths
    xstart(:,i) = mvnrnd(zeros(4,1),eye(4));
end

tbirth = [1 1 1 20 20 20 40 40 60 60 80 80];
tdeath = (truth.K+1)*ones(1,nbirths);
tdeath(1) = 70; tdeath(3) = 70;

%generate the tracks
for targetnum=1:nbirths
    targetstate = xstart(:,targetnum);
    for k=tbirth(targetnum):min(tdeath(targetnum),truth.K)
        targetstate = gen_newstate_fn(model,targetstate,'noise');
        truth.X{k}= [truth.X{k} targetstate];
        truth.track_list{k} = [truth.track_list{k} targetnum];
        truth.N(k) = truth.N(k) + 1;
     end
end
truth.total_tracks= nbirths;
