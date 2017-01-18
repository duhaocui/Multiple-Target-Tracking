function [ r_update,x_update,p_update] = hypo_update( bestAssign,rupd,xupd,Pupd,...
    rnew,xnew,Pnew,n,m,model,nCost,idx_in,idx_out )
%Update single target hypothesis according to the assignment

% Making assignments
[rr,xx,PP] = make_assign(bestAssign,rupd,xupd,Pupd,rnew(idx_in),...
    xnew(:,idx_in),Pnew(:,:,idx_in),n,m,model);
w = exp(-nCost-logsumexp(-nCost))';
if isempty(w)
    rr = cat(1,rr,rnew(idx_out));
    xx = cat(2,xx,xnew(:,idx_out));
    PP = cat(3,PP,Pnew(:,:,idx_out));
else
    % Generate q(h,j) r_h,x_h,P_h
    [pn,ph,phi,h_r,h_x,h_p] = hypo_all(w,rr,xx,PP,m,n,model);
    
    % Generate cost function
    [C,r_hat,x_hat,P_hat] = cost(phi,h_r,h_x,h_p,model);
    
    % Solve linear transport problem
    [Cmin,phi] = LP_transport(C,pn,ph);
    temp = Cmin;
    indicator = 0;
    while indicator == 0
        [C,r_temp,x_temp,P_temp] = cost(phi,h_r,h_x,h_p,model);
        [Cmin,phi] = LP_transport(C,pn,ph);
        if temp - Cmin < 1e-3 && temp >= Cmin
            indicator = 1;
            r_hat = r_temp;
            x_hat = x_temp;
            P_hat = P_temp;
        else
            temp = Cmin;
        end
    end
    
    rr = cat(1,r_hat,rnew(idx_out));
    xx = cat(2,x_hat,xnew(:,idx_out));
    PP = cat(3,P_hat,Pnew(:,:,idx_out));
end

% Prune low track weights
idx = rr > model.threshold;
r_update = rr(idx);
x_update = xx(:,idx);
p_update = PP(:,:,idx);


end

