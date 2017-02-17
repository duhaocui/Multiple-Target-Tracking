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

% Update unknown targets
[w_new,r_new,x_new,P_new] = ppp_update(lambdau,xu,Pu,z,model);
lambdau = (1-Pd)*lambdau;

% Loop through the global hypothesis
w = cell(len,1);
for l = 1:len
    
    % Gating
    [valid_idx,idx_out] = gate_meas_gms(z,model,x{l},P{l});
    
    % Update unknown targets
%     [wnew,rnew,xnew,Pnew] = ppp_update(lambdau,xu,Pu,z(:,valid_idx),model);
%     [wout,rout,xout,Pout] = ppp_update(lambdau,xu,Pu,z(:,idx_out),model);
%     lambdau = (1-Pd)*lambdau;
    wnew = w_new(valid_idx);rnew = r_new(valid_idx);
    xnew = x_new(:,valid_idx);Pnew = P_new(:,:,valid_idx);
    wout = w_new(idx_out);rout = r_new(idx_out);
    xout = x_new(:,idx_out);Pout = P_new(:,:,idx_out);
    
    % Extract number of tracks and measurements
    n = length(r{l});
    m = size(z(:,valid_idx),2);
    
    % Update existing tracks
    [wupd,rupd,xupd,Pupd] = mbm_update(z(:,valid_idx),n,r{l},x{l},P{l},model);
    
    % Calculate cost function
    cost = -log(wupd(:,2:end)./repmat(wnew',n,1));
    
    % Find M-best assignment
    Mt = ceil(H_max*sqrt(w_update(l)/sum(sqrt(w_update))));
    
    [bestAssign, nCost] = mbestwrap_updt_custom(cost,Mt,wupd(:,1));

    if length(nCost)<Mt
        M = length(nCost);
    else
        M = Mt;
    end
    
    % Update single target hypothesis
    [rMurty,xMurty,pMurty] = hypo_update(bestAssign,rupd,xupd,Pupd,...
        rnew,xnew,Pnew,rout,xout,Pout,n,m,nCost,M,model);
    
    % Normalisation
    wMurty = exp(-nCost)'*prod(wnew);
    
    if size(wMurty,1) ~= 0
        w{l} = wMurty*w_update(l)*prod(wout);
    end
    
    r_update = cat(1,r_update,rMurty);
    x_update = cat(1,x_update,xMurty);
    p_update = cat(1,p_update,pMurty);
    
end

if size(wMurty,1) ~= 0
    w_update = cell2mat(w);
else
    w_update = 1;
end

% Normalisation
w_update = normalize(w_update);

% Pruning
[w_update,r_update,x_update,p_update] = pruning(w_update,r_update,...
    x_update,p_update,H_threshold);

% Capping
[w_update,r_update,x_update,p_update] = capping(w_update,r_update,...
    x_update,p_update,H_max);

% Best state extraction
x_est = state_extract( w_update, r_update, x_update);