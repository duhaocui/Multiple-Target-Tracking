function [ r_update,x_update,p_update] = hypo_update( bestAssign,rupd,xupd,Pupd,...
    rnew,xnew,Pnew,rout,xout,Pout,n,m,model,nCost )
%Update single target hypothesis according to the assignment

% Making assignments
[rr,xx,PP] = make_assign(bestAssign,rupd,xupd,Pupd,rnew,...
    xnew,Pnew,n,m,model);
w = exp(-nCost-logsumexp(-nCost))';
if isempty(w)
    rr = cat(1,rr,rout);
    xx = cat(2,xx,xout);
    PP = cat(3,PP,Pout);
else
    % Generate q(h,j) r_h,x_h,P_h
    [pj,phi,h_r,h_x,h_p] = hypo_all(w,rr,xx,PP,m,n,model);
    
    % Generate cost function
    [r_hat,x_hat,P_hat] = Mstep(phi,h_r,h_x,h_p,model);
    
    [w,Cmin] = Estep(r_hat,x_hat,P_hat,rr,xx,PP,pj,w);
    temp = Cmin;
    indicator = 0;
    maxIterations = 1e2;
    iteration = 0;
    while indicator == 0
        [pj,phi,h_r,h_x,h_p] = hypo_all(w,rr,xx,PP,m,n,model);
        [r_temp,x_temp,P_temp] = Mstep(phi,h_r,h_x,h_p,model);
        [w,Cmin] = Estep(r_temp,x_temp,P_temp,rr,xx,PP,pj,w);
        iteration = iteration + 1;
        if (abs((temp-Cmin)/temp) < 1e-2) || (iteration > maxIterations) || (Cmin>temp)
            indicator = 1;
            r_hat = r_temp;
            x_hat = x_temp;
            P_hat = P_temp;
        else
            temp = Cmin;
        end
    end
    
    
    % Solve linear transport problem
%     [Cmin,phi] = LP_transport(C,pn,ph);
%     temp = Cmin;
%     indicator = 0;
%     maxIterations = 1e2;
%     iteration = 0;
%     while indicator == 0
%         [C,r_temp,x_temp,P_temp] = cost(phi,h_r,h_x,h_p,model);
%         [Cmin,phi] = LP_transport(C,pn,ph);
%         iteration = iteration + 1;
%         if (abs(temp - Cmin) < 1e-3) || (iteration > maxIterations)
%             indicator = 1;
%             r_hat = r_temp;
%             x_hat = x_temp;
%             P_hat = P_temp;
%         else
%             temp = Cmin;
%         end
%     end
    
    rr = cat(1,r_hat,rout);
    xx = cat(2,x_hat,xout);
    PP = cat(3,P_hat,Pout);
end

% Prune low track weights
idx = rr > model.threshold;
r_update = rr(idx);
x_update = xx(:,idx);
p_update = PP(:,:,idx);


end

