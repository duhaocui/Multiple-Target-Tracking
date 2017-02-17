function [lambdau,xu,Pu,r_update,x_update,p_update,x_est] = ...
    updating(lambdau,xu,Pu,r,x,P,z,model)
%UPDATE: CONSTRUCT COMPONENTS OF DISTRIBUTION UPDATED WITH MEASUREMENTS

% Gating
[gatingGroup,gating_mb,idx_out] = group_gating(z,r,x,P,model);
[~,rnew_out,xnew_out,Pnew_out] = ppp_update(lambdau,xu,Pu,z(:,idx_out),model);
numGroups = length(gatingGroup);
if numGroups == 0
    [~,rnew,xnew,Pnew] = ppp_update(lambdau,xu,Pu,z,model);
    lambdau = (1-model.Pd)*lambdau;
    [r_update,x_update,p_update] = pruning(rnew,xnew,Pnew,model);
    [lambdau,xu,Pu] = recycling(rnew,xnew,Pnew,lambdau,xu,Pu,model);
    x_est = state_extract(r_update,x_update);
    return;
end

group_r = cell(numGroups,1);
group_x = cell(numGroups,1);
group_P = cell(numGroups,1);
rr = rnew_out;
xx = xnew_out;
PP = Pnew_out;

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
    costs = -log(wupd(:,2:end)./repmat(wnew',n,1));
    
    % Find M-best assignment
    [bestAssign, nCost] = mbestwrap_updt_custom(costs,model.M,wupd(:,1));
    Wp = exp(-nCost-logsumexp(-nCost))';
    
    % Making assignments
    [group_r{g},group_x{g},group_P{g}] = make_assign(bestAssign,rupd,xupd,Pupd,...
        rnew,xnew,Pnew,n,m,model);
    
    % Variational mixture reduction
    [pn,ph,phi,h_r,h_x,h_p] = hypo_all(Wp,group_r{g},group_x{g},group_P{g},m+n,model);
    [C,r_hat,x_hat,P_hat] = cost(phi,h_r,h_x,h_p,model);
    [Cmin,phi] = LP_transport(C,pn,ph);
    temp = Cmin;
    maxIteration = 100;
    iteration = 0;
    while(1)
        iteration = iteration + 1;
        [C,r_temp,x_temp,P_temp] = cost(phi,h_r,h_x,h_p,model);
        [Cmin,phi] = LP_transport(C,pn,ph);
        if temp - Cmin < 1e-3 && temp >= Cmin || iteration == maxIteration
            r_hat = r_temp;
            x_hat = x_temp;
            P_hat = P_temp;
            break;
        else
            temp = Cmin;
        end
    end
    
    % Combine results from independent groups
    rr = [rr;r_hat];
    xx = [xx x_hat];
    PP = cat(3,PP,P_hat);
end

% Update unknown PPP intensity
lambdau = (1-model.Pd)*lambdau;
% Prune low-weight tracks
[r_update,x_update,p_update] = pruning(rr,xx,PP,model);
% Recycle low-weight tracks
[lambdau,xu,Pu] = recycling(rr,xx,PP,lambdau,xu,Pu,model);
% Best state extraction
x_est = state_extract(r_update,x_update);

