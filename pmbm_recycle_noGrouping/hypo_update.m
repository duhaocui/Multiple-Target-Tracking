function [ rMurty,xMurty,pMurty,lambdauMurty,xuMurty,PuMurty ] = hypo_update( bestAssign,rupd,xupd,Pupd,...
    rnew,xnew,Pnew,nCost,valid_idx,idx_out,lambdau,xu,Pu,model )
%Update single target hypothesis according to the assignment

M = length(nCost);
rMurty = cell(M,1);
xMurty = cell(M,1);
pMurty = cell(M,1);

lambdauMurty = cell(M,1);
xuMurty = cell(M,1);
PuMurty = cell(M,1);

for i = 1:M
    [rr,xx,PP] = make_assign(bestAssign(i,:),rupd,xupd,Pupd,...
        rnew(valid_idx),xnew(:,valid_idx),Pnew(:,:,valid_idx));
    
    rr = cat(1,rr,rnew(idx_out));
    xx = cat(2,xx,xnew(:,idx_out));
    PP = cat(3,PP,Pnew(:,:,idx_out));
    
    % Prune low track weights
    idx = rr > model.H_threshold;
    rMurty{i} = rr(idx);
    xMurty{i} = xx(:,idx);
    pMurty{i} = PP(:,:,idx);
    
    % Recycling
    idx = (rr < model.threshold) & idx;
    lambdauMurty{i} = [lambdau;rr(idx)];
    xuMurty{i} = [xu xx(:,idx)];
    PuMurty{i} = cat(3,Pu,PP(:,:,idx));
    
%     [lambdauMurty{i},xuMurty{i},PuMurty{i}] = ...
%         gaus_merge(lambdauMurty{i},xuMurty{i},PuMurty{i},1);
    
end


end

