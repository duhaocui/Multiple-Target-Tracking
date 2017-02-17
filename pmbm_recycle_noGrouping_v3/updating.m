function [r_recycle,x_recycle,p_recycle,r_update,x_update,p_update,x_est,w_update] = ...
    updating(lambdau,xu,Pu,r,x,P,z,model,w_update)
%UPDATE: CONSTRUCT COMPONENTS OF DISTRIBUTION UPDATED WITH MEASUREMENTS

% Extract parameters from model
Pd = model.Pd;
H_max = model.H_max;
H_threshold = model.H_threshold;

% Only PPP update at the first time step
if isempty(r{1})
    [~,rr,xx,PP] = ppp_update(lambdau{1},xu{1},Pu{1},z,model);
    idx = rr > model.threshold;
    r_update{1} = rr(idx);
    x_update{1} = xx(:,idx);
    p_update{1} = PP(:,:,idx);
    x_est = x_update{1}(:,r_update{1}>0.5);
    
    idx = (rr < model.threshold) & (rr > H_threshold);
    r_recycle{1} = [(1-Pd)*lambdau{1};rr(idx)];
    x_recycle{1} = [xu{1} xx(:,idx)];
    p_recycle{1} = cat(3,Pu{1},PP(:,:,idx));
    return;
end

% Extract number of global hypothesis
len = length(r);

% Allocate memory for existing tracks
r_update = cell(0,1);
x_update = cell(0,1);
p_update = cell(0,1);
r_recycle = cell(0,1);
x_recycle = cell(0,1);
p_recycle = cell(0,1);

% Loop through the global hypothesis
w = cell(len,1);
for l = 1:len
    
    % Update unknown targets
    [wnew,rnew,xnew,Pnew] = ppp_update(lambdau{l},xu{l},Pu{l},z,model);
    lambdau{l} = (1-Pd)*lambdau{l};
    
    % Gating
    [valid_idx,idx_out] = gate_meas_gms(z,model,x{l},P{l});
    
    % Update existing tracks
    [wupd,rupd,xupd,Pupd] = mbm_update(z(:,valid_idx),r{l},x{l},P{l},model);
    
    % Calculate cost function
    cost = -log(wupd(:,2:end)./repmat(wnew(valid_idx)',length(r{l}),1));
    
    % Find M-best assignment
    Mt = ceil(H_max*sqrt(w_update(l)/sum(sqrt(w_update))));
    [bestAssign, nCost] = mbestwrap_updt_custom(cost,Mt,wupd(:,1));
    
    % Update single target hypothesis
    [rMurty,xMurty,pMurty,lambdauMurty,xuMurty,PuMurty] = ...
        hypo_update(bestAssign,rupd,xupd,Pupd,rnew,xnew,Pnew,...
        nCost,valid_idx,idx_out,lambdau{l},xu{l},Pu{l},model);
    
    w{l} = exp(-nCost)'*prod(wnew)*w_update(l);
    r_update = cat(1,r_update,rMurty);
    x_update = cat(1,x_update,xMurty);
    p_update = cat(1,p_update,pMurty);
    r_recycle = cat(1,r_recycle,lambdauMurty);
    x_recycle = cat(1,x_recycle,xuMurty);
    p_recycle = cat(1,p_recycle,PuMurty);
    
end

% Normalisation
w_update = normalize(cell2mat(w));

% Pruning
[w_update,r_update,x_update,p_update,r_recycle,x_recycle,p_recycle] = pruning(w_update,r_update,...
    x_update,p_update,r_recycle,x_recycle,p_recycle,H_threshold);

% Capping
[w_update,r_update,x_update,p_update,r_recycle,x_recycle,p_recycle] = capping(w_update,r_update,...
    x_update,p_update,r_recycle,x_recycle,p_recycle,H_max);

% Best state extraction
x_est = state_extract( w_update, r_update, x_update);
