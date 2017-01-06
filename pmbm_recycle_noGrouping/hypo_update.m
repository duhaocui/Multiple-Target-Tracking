function [ rMurty,xMurty,pMurty,lambdau_upd,xu_upd,pu_upd ] = hypo_update( bestAssign,rupd,xupd,Pupd,...
    rnew,xnew,Pnew,n,m,nCost,M,model,lambdau,xu,Pu,valid_idx,idx_out )
%Update single target hypothesis according to the assignment
if isempty(nCost)
    M = 1;
    rMurty = cell(1,1);
    xMurty = cell(1,1);
    pMurty = cell(1,1);
    lambdau_upd = cell(1,1);
    xu_upd = cell(1,1);
    pu_upd = cell(1,1);
else
    rMurty = cell(M,1);
    xMurty = cell(M,1);
    pMurty = cell(M,1);
    lambdau_upd = cell(M,1);
    xu_upd = cell(M,1);
    pu_upd = cell(M,1);
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
    
    % Recycle low weight track
    idx_recycle = rr <=model.threshold;
    
    % Append recycling MB component to PPP
    lambdau_upd{assign} = [lambdau;rr(idx_recycle)];
    xu_upd{assign} = [xu xx(:,idx_recycle)];
    pu_upd{assign} = cat(3,Pu,PP(:,:,idx_recycle));
    
end


end

