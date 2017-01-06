function [lambdau,xu,Pu,r_update,x_update,p_update,x_est,w_update] = ...
    updating(lambdau,xu,Pu,r,x,P,z,model,w_update)
%UPDATE: CONSTRUCT COMPONENTS OF DISTRIBUTION UPDATED WITH MEASUREMENTS

% Extract parameters from model
Pd = model.Pd;
H_max = model.H_max;
H_threshold = model.H_threshold;
 
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
    [z_gate,valid_idx,idx_out] = gate_meas_gms(z,model,x{l},P{l});
    
    % Extract number of tracks and measurements
    n = length(r{l});
    m = size(z_gate,2);
    
    % Update existing tracks
    [wupd,rupd,xupd,Pupd] = mbm_update(z_gate,n,r{l},x{l},P{l},model);
    
    % Calculate cost function
    %cost = -log(wupd(:,2:end)./repmat(wupd(:,1),1,m) + eps);
    cost = -log(wupd(:,2:end)./repmat(wnew(valid_idx)',n,1));
    
    % Find M-best assignment
    Mt = ceil(H_max*sqrt(w_update(l)/sum(sqrt(w_update))));
    
    [bestAssign, nCost] = mbestwrap_updt_custom(cost,Mt,wupd(:,1));
    if length(nCost)<Mt
        M = length(nCost);
    else
        M = Mt;
    end
    
    % Update single target hypothesis
    [rMurty,xMurty,pMurty,lambdau_upd,xu_upd,pu_upd] = ...
        hypo_update(bestAssign,rupd,xupd,Pupd,rnew,xnew,...
        Pnew,n,m,nCost,M,model,lambdau{l},xu{l},Pu{l},valid_idx,idx_out);
    
    % Normalisation
    wMurty = exp(-nCost)'*prod(wnew);
    
    if size(wMurty,1) ~= 0
        w{l} = wMurty*w_update(l);
    end
    
    r_update = cat(1,r_update,rMurty);
    x_update = cat(1,x_update,xMurty);
    p_update = cat(1,p_update,pMurty);
    
    r_recycle = cat(1,r_recycle,lambdau_upd);
    x_recycle = cat(1,x_recycle,xu_upd);
    p_recycle = cat(1,p_recycle,pu_upd);
    
end

lambdau = r_recycle;
xu = x_recycle;
Pu = p_recycle;

if size(wMurty,1) ~= 0
    w_update = cell2mat(w);
else
    w_update = 1;
end

% Normalisation
w_update = normalize(w_update);

% Pruning
[w_update,r_update,x_update,p_update,lambdau,xu,Pu] = pruning(w_update,r_update,...
    x_update,p_update,lambdau,xu,Pu,H_threshold);

% Capping
[w_update,r_update,x_update,p_update,lambdau,xu,Pu] = capping(w_update,r_update,...
    x_update,p_update,lambdau,xu,Pu,H_max);

% Best state extraction
x_est = state_extract( w_update, r_update, x_update);
