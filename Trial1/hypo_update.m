function [ lambdau,xu,Pu,r_update,x_update,p_update] = hypo_update( lambdau,xu,Pu,bestAssign,rupd,xupd,Pupd,...
    rnew,xnew,Pnew,rout,xout,Pout,n,m,model,nCost )
%Update single target hypothesis according to the assignment

% Making assignments
[rr,xx,PP] = make_assign(bestAssign,rupd,xupd,Pupd,rnew,...
    xnew,Pnew,n,m,model);
w = exp(-nCost-logsumexp(-nCost))';
if isempty(w)
    r = cat(1,rr,rout);
    x = cat(2,xx,xout);
    P = cat(3,PP,Pout);
else
    [r_hat,x_hat,P_hat] = mixtureReduction(w,rr,xx,PP);
%     [~,~,phi,h_r,h_x,h_p] = hypo_all(w,rr,xx,PP,model);
%     [r_hat,x_hat,P_hat] = mixture_reduction(phi,h_r,h_x,h_p,model);
    [rr,xx,PP,temp] = expectation(w,rr,xx,PP,r_hat,x_hat,P_hat);
    maxIteration = 10;
    iteration = 0;
    while(1)
        iteration = iteration + 1;
        [r_hat,x_hat,P_hat] = mixtureReduction(w,rr,xx,PP);
        [rr,xx,PP,Cmin] = expectation(w,rr,xx,PP,r_hat,x_hat,P_hat);
        if Cmin < temp && temp - Cmin < 1 || iteration == maxIteration
            break;
        else
            temp = Cmin;
        end
    end
    
    r = cat(1,r_hat,rout);
    x = cat(2,x_hat,xout);
    P = cat(3,P_hat,Pout);
end

% idx = rr < 0.1;
% lambdau = [lambdau rr(idx)'];
% xu = [xu xx(:,idx)];
% Pu = cat(3,Pu,PP(:,:,idx));
% Prune low track weights
idx = r > model.threshold;
r_update = r(idx);
x_update = x(:,idx);
p_update = P(:,:,idx);


end

