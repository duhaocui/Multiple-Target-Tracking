function [lambdau,xu,Pu,r_update,x_update,p_update,x_est] = ...
    updating(lambdau,xu,Pu,r,x,P,z,model)
%UPDATE: CONSTRUCT COMPONENTS OF DISTRIBUTION UPDATED WITH MEASUREMENTS

% Gating
[z_gate,idx_out] = gate_meas_gms(z,model,x,P);

% Update unknown targets
[wnew,rnew,xnew,Pnew] = ppp_update(lambdau,xu,Pu,z_gate,model);
[~,rout,xout,Pout] = ppp_update(lambdau,xu,Pu,z(:,idx_out),model);
lambdau = (1-model.Pd)*lambdau;

% Extract number of tracks and measurements
n = length(r);
m = size(z_gate,2);

% Update existing tracks
[wupd,rupd,xupd,Pupd] = mbm_update(z_gate,n,r,x,P,model);

% Calculate cost function
cost = -log(wupd(:,2:end)./repmat(wupd(:,1),1,m));

% Find M-best assignment
[bestAssign, nCost] = mbestwrap_updt_custom(cost,model.M,wnew);

% Update single target hypothesis
[r_update,x_update,p_update] = hypo_update(bestAssign,rupd,xupd,Pupd,...
    rnew,xnew,Pnew,rout,xout,Pout,n,m,model,nCost);

% Prune and recycle low track weights
idx = r_update <= model.recycleThreshold;
lambdau = [lambdau;r_update(idx)];
xu = [xu x_update(:,idx)];
Pu = cat(3,Pu,p_update(:,:,idx));

% Best state extraction
x_est = state_extract(r_update,x_update);

idx = r_update > model.recycleThreshold;
r_update = r_update(idx);
x_update = x_update(:,idx);
p_update = p_update(:,:,idx);





