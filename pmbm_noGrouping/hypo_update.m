function [ rMurty,xMurty,pMurty ] = hypo_update( bestAssign,rupd,xupd,Pupd,...
    rnew,xnew,Pnew,rout,xout,Pout,n,m,nCost,model )
%Update single target hypothesis according to the assignment

M = length(nCost);
rMurty = cell(M,1);
xMurty = cell(M,1);
pMurty = cell(M,1);


for assign = 1:M
    [rr,xx,PP] = make_assign(assign,bestAssign,rupd,xupd,Pupd,...
        rnew,xnew,Pnew,n,m);
    
    rr = cat(1,rr,rout);
    xx = cat(2,xx,xout);
    PP = cat(3,PP,Pout);
    
    % Prune low track weights
    idx = rr > model.threshold;
    rMurty{assign} = rr(idx);
    xMurty{assign} = xx(:,idx);
    pMurty{assign} = PP(:,:,idx);
    
end

end

