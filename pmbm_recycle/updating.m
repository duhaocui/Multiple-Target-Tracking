function [lambdau,xu,Pu,r_update,x_update,p_update,x_est,w_update] = ...
    updating(lambdau,xu,Pu,r,x,P,z,model,w_update)
%UPDATE: CONSTRUCT COMPONENTS OF DISTRIBUTION UPDATED WITH MEASUREMENTS

% Extract number of global hypothesis
len = length(r);

% Allocate memory for existing tracks
r_update = cell(0,1);
x_update = cell(0,1);
p_update = cell(0,1);

rRecycle = cell(0,1);
xRecycle = cell(0,1);
pRecycle = cell(0,1);

if isempty(r{1})
    [~,rnew,xnew,Pnew] = ppp_update(lambdau{1},xu{1},Pu{1},z,model);
    lambdau{1} = (1-model.Pd)*lambdau{1};
    [r_update{1},x_update{1},p_update{1},r_recycle,x_recycle,P_recycle] = ...
        track_pruning(rnew,xnew,Pnew,model);
    lambdau{1} = [lambdau{1};r_recycle];
    xu{1} = [xu{1} x_recycle];
    Pu{1} = cat(3,Pu{1},P_recycle);
    x_est = state_extract(w_update,r_update,x_update);
    return;
end

% Loop through the global hypothesis
w = cell(len,1);
for l = 1:len
    % Gating
    [gatingGroup,gating_mb,idx_out] = group_gating(z,r{l},x{l},P{l},model);
    numGroups = length(gatingGroup);
    [wout,rout,xout,Pout] = ppp_update(lambdau{l},xu{l},Pu{l},z(:,idx_out),model);
    
    group_r = cell(numGroups,1);
    group_x = cell(numGroups,1);
    group_P = cell(numGroups,1);
    group_w = cell(numGroups,1);
    sets = cell(1,numGroups);
    for g = 1:numGroups
        % Update unknown targets
        [wnew,rnew,xnew,Pnew] = ppp_update(lambdau{l},xu{l},Pu{l},z(:,gatingGroup{g}),model);
        
        % Extract number of tracks and measurements
        n = length(gating_mb{g}.r);
        m = size(z(:,gatingGroup{g}),2);
        
        % Update existing tracks
        [wupd,rupd,xupd,Pupd] = mbm_update(z(:,gatingGroup{g}),n,gating_mb{g},model);
        
        % Calculate cost function
        cost = -log(wupd(:,2:end)./repmat(wnew',n,1));
        
        % Find M-best assignment
        Mt = ceil(model.H_max*sqrt(w_update(l)/sum(sqrt(w_update))));
        
        [bestAssign, nCost] = mbestwrap_updt_custom(cost,Mt,wupd(:,1));
        
        % Update single target hypothesis
        [group_r{g},group_x{g},group_P{g}] = make_assign(bestAssign,rupd,xupd,Pupd,...
            rnew,xnew,Pnew,n,m,model);
        
        Wp = exp(-nCost)'*prod(wnew);
        Wp2 = Wp/sum(Wp);
        [Wp3,order] = sort(Wp2,'descend');
        Y = cumsum(Wp3);
        pos = find(Y>=0.999,1);
        group_w{g} = Wp(order(1:pos));
        group_r{g} = group_r{g}(:,order(1:pos));
        group_x{g} = group_x{g}(:,:,order(1:pos));
        group_P{g} = group_P{g}(:,:,:,order(1:pos));
        sets{g} = 1:length(group_w{g});
    end
    
    lambdau{l} = (1-model.Pd)*lambdau{l};
    
    c = allcombs(sets{:});
    [N,M] = size(c);
    wMurty = zeros(N,1);
    rMurty = cell(N,1);
    xMurty = cell(N,1);
    pMurty = cell(N,1);
    
    r_re = cell(N,1);
    x_re = cell(N,1);
    P_re = cell(N,1);
    
    for i = 1:N
        rr = zeros(1,0);
        xx = zeros(4,0);
        PP = zeros(4,4,0);
        wc = 1;
        
        for k = 1:M
            wc = wc*group_w{k}(c(i,k));
            rr = [rr;group_r{k}(:,c(i,k))];
            xx = [xx group_x{k}(:,:,c(i,k))];
            PP = cat(3,PP,group_P{k}(:,:,:,c(i,k)));
        end
        
        wMurty(i) = wc;
        [rMurty{i},xMurty{i},pMurty{i},r_re{i},x_re{i},P_re{i}] = ...
            track_pruning([rr;rout],[xx xout],cat(3,PP,Pout),model);
        r_re{i} = [r_re{i};lambdau{l}];
        x_re{i} = [x_re{i} xu{l}];
        P_re{i} = cat(3,P_re{i},Pu{l});
    end
    
    w{l} = wMurty*w_update(l)*prod(wout);
    r_update = cat(1,r_update,rMurty);
    x_update = cat(1,x_update,xMurty);
    p_update = cat(1,p_update,pMurty);
    
    rRecycle = cat(1,rRecycle,r_re);
    xRecycle = cat(1,xRecycle,x_re);
    pRecycle = cat(1,pRecycle,P_re);
    
end

lambdau = rRecycle;
xu = xRecycle;
Pu = pRecycle;

% Normalisation
w_update = normalize(cell2mat(w));

% Pruning
[w_update,r_update,x_update,p_update,lambdau,xu,Pu] = pruning(w_update,r_update,...
    x_update,p_update,lambdau,xu,Pu,model.H_threshold);

% Capping
[w_update,r_update,x_update,p_update,lambdau,xu,Pu] = capping(w_update,r_update,...
    x_update,p_update,lambdau,xu,Pu,model.H_max);

% Best state extraction
x_est = state_extract(w_update,r_update,x_update);