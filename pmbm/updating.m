function [lambdau,xu,Pu,r_update,x_update,p_update,x_est,w_update] = ...
    updating(lambdau,xu,Pu,r,x,P,z,model,w_update)
%UPDATE: CONSTRUCT COMPONENTS OF DISTRIBUTION UPDATED WITH MEASUREMENTS

% Extract number of global hypothesis
len = length(r);

% Allocate memory for existing tracks
r_update = cell(0,1);
x_update = cell(0,1);
p_update = cell(0,1);

if isempty(r{1})
    [~,rnew,xnew,Pnew] = ppp_update(lambdau,xu,Pu,z,model);
    lambdau = (1-model.Pd)*lambdau;
    [r_update{1},x_update{1},p_update{1}] = track_pruning(rnew,xnew,Pnew,model);
    x_est = state_extract(w_update,r_update,x_update);
    return;
end

% Loop through the global hypothesis
w = cell(len,1);
for l = 1:len
    
    % Gating
    [z_gate,idx_out] = gate_meas_gms(z,model,x{l},P{l});
    
    % Update unknown targets
    [wnew,rnew,xnew,Pnew] = ppp_update(lambdau,xu,Pu,z_gate,model);
    [~,rout,xout,Pout] = ppp_update(lambdau,xu,Pu,z(:,idx_out),model);
    lambdau = (1-model.Pd)*lambdau;
    
    % Extract number of tracks and measurements
    n = length(r{l});
    m = size(z_gate,2);
    
    % Update existing tracks
    [wupd,rupd,xupd,Pupd] = mbm_update(z_gate,n,r{l},x{l},P{l},model);
    
    % Calculate cost function
    cost = -log(wupd(:,2:end)./repmat(wnew',n,1));
    
    % Find M-best assignment
    Mt = ceil(model.H_max*sqrt(w_update(l)/sum(sqrt(w_update))));
    
    [bestAssign, nCost] = mbestwrap_updt_custom(cost,Mt,wupd(:,1));
    
    % Update single target hypothesis
    [rMurty,xMurty,pMurty] = hypo_update(bestAssign,rupd,xupd,Pupd,...
        rnew,xnew,Pnew,rout,xout,Pout,n,m,model);
    
    % Normalisation
    wMurty = exp(-nCost)'*prod(wnew);
    
    w{l} = exp(-wMurty)*w_update(l);
    
    r_update = cat(1,r_update,rMurty);
    x_update = cat(1,x_update,xMurty);
    p_update = cat(1,p_update,pMurty);
    
end


% Normalisation
w_update = normalize(cell2mat(w));

% Pruning
[w_update,r_update,x_update,p_update] = pruning(w_update,r_update,...
    x_update,p_update,model.H_threshold);

% Capping
[w_update,r_update,x_update,p_update] = capping(w_update,r_update,...
    x_update,p_update,model.H_max);

% Best state extraction
x_est = state_extract(w_update,r_update,x_update);