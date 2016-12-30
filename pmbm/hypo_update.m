function [ rMurty,xMurty,pMurty ] = hypo_update( bestAssign,rupd,xupd,Pupd,...
    rnew,xnew,Pnew,rout,xout,Pout,n,m,model )
%Update single target hypothesis according to the assignment

M = size(bestAssign,1);
rMurty = cell(M,1);
xMurty = cell(M,1);
pMurty = cell(M,1);

for k = 1:M
    [rr,xx,PP] = make_assign(k,bestAssign,rupd,xupd,Pupd,...
        rnew,xnew,Pnew,n,m);
    
    rr = cat(1,rr,rout);
    xx = cat(2,xx,xout);
    PP = cat(3,PP,Pout);
    
    % Prune low track weights
    [rMurty{k},xMurty{k},pMurty{k}] = track_pruning(rr,xx,PP,model);
    
end

end

