function [ r,x,P,r_recycle,x_recycle,P_recycle ] = track_pruning( r,x,P,model )

idx = r <= model.recycleThreshold;
r_recycle = r(idx);
x_recycle = x(:,idx);
P_recycle = P(:,:,idx);

idx = r > model.recycleThreshold;
r = r(idx);
x = x(:,idx);
P = P(:,:,idx);

end

