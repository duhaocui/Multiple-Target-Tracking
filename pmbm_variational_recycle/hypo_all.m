function [ pn,ph,phi,h_r,h_x,h_p ] = hypo_all( w,rr,xx,PP,m,n,model )

dim = model.x_dim;
M = length(w);
r = reshape(rr,numel(rr),1);
x = reshape(xx,dim,size(xx,2)*size(xx,3));
P = reshape(PP,dim,dim,size(PP,3)*size(PP,4));
idx = rr~=0;
r = r(idx);
x = x(:,idx);
P = P(:,:,idx);
[~,idx_unique,~] = unique(x(1,:));
h_r = r(idx_unique)';
h_x = x(:,idx_unique);
h_p = P(:,:,idx_unique);

H = length(h_r);
phi = zeros(H,n+m);

for i = 1:H
    for k = 1:n+m
        temp = 0;
        for j = 1:M
            if isequal(h_x(1,i),xx(1,k,j)) && isequal(h_r(i),rr(k,j))
                temp = temp + w(j);
            end
        end
        phi(i,k) = temp;
    end
end

phi(:,all(phi==0,1))=[];
ph = sum(phi,2);
pn = sum(phi,1)';

%-----
% ph(ph>=1) = 1;
% pn(pn>=1) = 1;
%-----


end

