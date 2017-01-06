function [ rMurty,xMurty,pMurty ] = hypo_update( bestAssign,rupd,xupd,Pupd,...
    rnew,xnew,Pnew,n,m,nCost,M,model,valid_idx,idx_out )
%Update single target hypothesis according to the assignment
if isempty(nCost)
    M = 1;
    rMurty = cell(1,1);
    xMurty = cell(1,1);
    pMurty = cell(1,1);
else
    rMurty = cell(M,1);
    xMurty = cell(M,1);
    pMurty = cell(M,1);
end

for assign = 1:M
    [rr,xx,PP] = make_assign(assign,bestAssign,rupd,xupd,Pupd,...
        rnew(valid_idx),xnew(:,valid_idx),Pnew(:,:,valid_idx),n,m);
    
    rr = cat(1,rr,rnew(idx_out));
    xx = cat(2,xx,xnew(:,idx_out));
    PP = cat(3,PP,Pnew(:,:,idx_out));
    
    % Prune low track weights
    idx = rr > model.threshold;
    rMurty{assign} = rr(idx);
    xMurty{assign} = xx(:,idx);
    pMurty{assign} = PP(:,:,idx);
     
end

end

