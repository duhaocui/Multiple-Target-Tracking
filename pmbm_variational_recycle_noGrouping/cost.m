function [ C,r_hat,x_hat,P_hat ] = cost(phi,h_r,h_x,h_p,model)

dim = model.x_dim;
[H,N] = size(phi);
x_hat = zeros(dim,N);
P_hat = zeros(dim,dim,N);
v = zeros(dim,H,N);

r_hat = (h_r*phi)';

for j = 1:N
    x_hat(:,j) = sum(repmat(h_r.*phi(:,j)',dim,1).*h_x,2)/r_hat(j);
    for h = 1:H
        v(:,h,j) = h_x(:,h) - x_hat(:,j);
        P_hat(:,:,j) = P_hat(:,:,j) + phi(h,j)*h_r(h)*(h_p(:,:,h)+v(:,h,j)*v(:,h,j)');
    end
end
P_hat = P_hat./reshape(kron(r_hat,ones(dim))',dim,dim,N);

C = zeros(H,N);
for j = 1:N
    for h = 1:H
        try
            inv_P = pinv(P_hat(:,:,j));
            temp = trace(inv_P*h_p(:,:,h)) + v(:,h,j)'*inv_P*v(:,h,j) + log(det(2*pi*P_hat(:,:,j)));
        catch
            C(h,j) = 1/eps;
        end
        if r_hat(j)==1 && h_r(h)==1
            C(h,j) = h_r(h)/2*temp;
        elseif r_hat(j)==1 && h_r(h)~=1
            C(h,j) = 1/eps;
        elseif r_hat(j)==0
            C(h,j) = 1/eps;
        else
            C(h,j) = -(1-h_r(h))*log(1-r_hat(j)) - h_r(h)*log(r_hat(j)) + h_r(h)/2*temp;
        end
    end
end

end

