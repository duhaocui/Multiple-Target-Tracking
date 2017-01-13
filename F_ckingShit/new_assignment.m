function [ w,rr_hat,xx_hat,PP_hat ] = new_assignment( ph,h_r,h_x,h_p,phi,C,N,model )

C = C./(phi+eps);
P0 = C';
P0 = P0 - min(P0(:));
[assignments, costs] = murty_custom(P0,model.M);
costs = costs + (min(P0(:)).*sum(assignments>0,2))';
w = exp(-costs-logsumexp(-costs))';
idx = w ~= 0;
w = w(idx);
assignments = assignments(idx,:);
len = length(w);
rr_hat = zeros(N,len);
xx_hat = zeros(4,N,len);
PP_hat = zeros(4,4,N,len);
for i = 1:len
    rr_hat(:,i) = h_r(assignments(i,:)).*ph(assignments(i,:))';
    xx_hat(:,:,i) = h_x(:,assignments(i,:)).*ph(assignments(i,:))';
    for j = 1:N
        PP_hat(:,:,j,i) = h_p(:,:,assignments(i,j))*ph(assignments(i,j));
    end
end


end

