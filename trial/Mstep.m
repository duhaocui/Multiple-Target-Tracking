function [ r_hat,x_hat,P_hat ] = Mstep(phi,h_r,h_x,h_p,model)

dim = model.x_dim;
[H,N] = size(phi);
x_hat = zeros(dim,N);
P_hat = zeros(dim,dim,N);
v = zeros(dim,H,N);

r_hat = (h_r*phi)';
r_hat(r_hat>=1) = 1-eps;

for j = 1:N
    x_hat(:,j) = sum(repmat(h_r.*phi(:,j)',dim,1).*h_x,2)/r_hat(j);
    for h = 1:H
        v(:,h,j) = h_x(:,h) - x_hat(:,j);
        P_hat(:,:,j) = P_hat(:,:,j) + phi(h,j)*h_r(h)*(h_p(:,:,h)+v(:,h,j)*v(:,h,j)');
    end
end
P_hat = abs(P_hat)./reshape(kron(r_hat,ones(dim))',dim,dim,N);


end

