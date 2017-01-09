function [ r_update,x_update,p_update] = ...
    hypo_update( bestAssign,rupd,xupd,Pupd,rnew,xnew,Pnew,rout,xout,Pout,n,m,model,nCost )
%Update single target hypothesis according to the assignment

% Making assignments
[rr,xx,PP] = make_assign(bestAssign,rupd,xupd,Pupd,rnew,...
    xnew,Pnew,n,m,model);
w = exp(-nCost-logsumexp(-nCost))';
if isempty(w)
    r_update = cat(1,rr,rout);
    x_update = cat(2,xx,xout);
    p_update = cat(3,PP,Pout);
else
    % Generate q(h,j) r_h,x_h,P_h
    [pn,ph,phi,h_r,h_x,h_p] = hypo_all(w,rr,xx,PP,m,n,model);
    
    % Generate cost function
    [C,r_hat,x_hat,P_hat] = cost(phi,h_r,h_x,h_p,model);
    
    % Solve linear transport problem
    [Cmin,phi] = LP_transport(C,pn,ph);
    temp = Cmin;
    indicator = 0;
    maxIterations = 1e2;
    iteration = 0;
    while indicator == 0
        [C,r_temp,x_temp,P_temp] = cost(phi,h_r,h_x,h_p,model);
        [Cmin,phi] = LP_transport(C,pn,ph);
        iteration = iteration + 1;
        if (abs(temp - Cmin) < 1e-3) || (iteration > maxIterations)
            indicator = 1;
            r_hat = r_temp;
            x_hat = x_temp;
            P_hat = P_temp;
        else
            temp = Cmin;
        end
    end
    
    r_update = cat(1,r_hat,rout);
    x_update = cat(2,x_hat,xout);
    p_update = cat(3,P_hat,Pout);
end


end

