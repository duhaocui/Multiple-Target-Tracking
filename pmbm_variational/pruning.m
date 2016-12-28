function [ r_update,x_update,p_update ] = pruning( r,x,P,model )

% Prune low track weights
idx1 = r > model.threshold;
idx2 = true(length(r),1);
for i = 1:length(r)
    [~,p] = chol(P(:,:,i));
    if p==1
        idx2(i) = false;
    end
end

idx = idx1&idx2;
r_update = r(idx);
x_update = x(:,idx);
p_update = P(:,:,idx);


end

