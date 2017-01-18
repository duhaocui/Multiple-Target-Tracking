function [ pn,ph,phi,h_r,h_x,h_p ] = hypo_all( w,rr,xx,PP,N,model )

dim = model.x_dim;
M = length(w);
r = reshape(rr,numel(rr),1);
x = reshape(xx,dim,size(xx,2)*size(xx,3));
P = reshape(PP,dim,dim,size(PP,3)*size(PP,4));
idx = rr~=0;
r = r(idx);
x = x(:,idx);
P = P(:,:,idx);
[~,idx] = unique(x(1,:));
h_r = r(idx)';
h_x = x(:,idx);
h_p = P(:,:,idx);

H = length(h_r);
phi = zeros(H,N);

for i = 1:H
    for k = 1:N
        for j = 1:M
            if isequal(h_x(1,i),xx(1,k,j)) && isequal(h_r(i),rr(k,j))
                phi(i,k) = phi(i,k) + w(j);
            end
        end
    end
end

ph = sum(phi,2);
pn = sum(phi,1)';

end

