function [lambdau,xu,Pu,r_update,x_update,p_update,x_est] = ...
    updating(lambdau,xu,Pu,r,x,P,z,model)
%UPDATE: CONSTRUCT COMPONENTS OF DISTRIBUTION UPDATED WITH MEASUREMENTS

% Gating
[gatingGroup,gating_mb,idx_out] = group_gating(z,r,x,P,model);
numGroups = length(gatingGroup);
if numGroups == 0
    [~,rnew,xnew,Pnew] = ppp_update(lambdau,xu,Pu,z,model);
    lambdau = (1-model.Pd)*lambdau;
    [r_update,x_update,p_update,r_recycle,x_recycle,p_recycle] = pruning(rnew,xnew,Pnew,model);
    lambdau = [lambdau;r_recycle];
    xu = [xu x_recycle];
    Pu = cat(3,Pu,p_recycle);
    x_est = state_extract(r_update,x_update);
    return;
end

group_r = cell(numGroups,1);
group_x = cell(numGroups,1);
group_P = cell(numGroups,1);
group_w = cell(numGroups,1);
sets = cell(1,numGroups);

for g = 1:numGroups
    % Extract number of tracks and measurements
    r = gating_mb{g}.r;
    x = gating_mb{g}.x;
    P = gating_mb{g}.P;
    n = length(r);
    m = length(gatingGroup{g});
    
    % Update unknown targets
    [wnew,rnew,xnew,Pnew] = ppp_update(lambdau,xu,Pu,z(:,gatingGroup{g}),model);
    
    % Update existing tracks
    [wupd,rupd,xupd,Pupd] = mbm_update(z(:,gatingGroup{g}),n,r,x,P,model);
    
    % Calculate cost function
    cost = -log(wupd(:,2:end)./repmat(wnew',n,1));
    
    % Find M-best assignment
    [bestAssign, nCost] = mbestwrap_updt_custom(cost,model.M,wupd(:,1));
    
    % Making assignments
    [group_r{g},group_x{g},group_P{g}] = make_assign(bestAssign,rupd,xupd,Pupd,...
        rnew,xnew,Pnew,n,m,model);
    
    Wp = exp(-nCost)'*prod(wnew);
    
    % Prune low-weight hypothesis in each gating group
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

% Try all combinations of different groups, select unique single target
% hypothesis and find the correspondence between Bernoulli component in MB
% mixture and approximated MB
[pn,ph,phi,h_r,h_x,h_p] = group_combs(sets,group_w,group_r,group_x,group_P);

% Variational approximation
[r_hat,x_hat,P_hat] = variational_approx(pn,ph,phi,h_r',h_x,h_p,model);
[~,rnew_out,xnew_out,Pnew_out] = ppp_update(lambdau,xu,Pu,z(:,idx_out),model);
rr = cat(1,r_hat,rnew_out);
xx = cat(2,x_hat,xnew_out);
PP = cat(3,P_hat,Pnew_out);

% Update unknown PPP intensity
lambdau = (1-model.Pd)*lambdau;

% Prune low-weight tracks
[r_update,x_update,p_update,r_recycle,x_recycle,p_recycle] = pruning(rr,xx,PP,model);
lambdau = [lambdau;r_recycle];
xu = [xu x_recycle];
Pu = cat(3,Pu,p_recycle);
% Best state extraction
x_est = state_extract(r_update,x_update);
