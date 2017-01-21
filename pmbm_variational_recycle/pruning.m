function [ r_update,x_update,p_update ] = pruning( r,x,P,model )

% Prune low track weights
idx = r > 0.1;
r_update = r(idx);
x_update = x(:,idx);
p_update = P(:,:,idx);


end

