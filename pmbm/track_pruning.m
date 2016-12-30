function [ r,x,P ] = track_pruning( r,x,P,model )

idx = r > model.threshold;
r = r(idx);
x = x(:,idx);
P = P(:,:,idx);

end

