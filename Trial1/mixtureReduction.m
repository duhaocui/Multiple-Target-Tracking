function [ r_hat,x_hat,P_hat ] = mixtureReduction( w,rr,xx,PP )

[N,M] = size(rr);
r_hat = zeros(N,1);
x_hat = zeros(4,N);
P_hat = zeros(4,4,N);

for j = 1:N
    for i = 1:M
        r_hat(j) = r_hat(j) + w(i)*rr(j,i);
        x_hat(:,j) = x_hat(:,j) + w(i)*rr(j,i)*xx(:,j,i);
    end
    x_hat(:,j) = x_hat(:,j)/r_hat(j);
    v = xx(:,j,i) - x_hat(:,j);
    for i = 1:M
        P_hat(:,:,j) = P_hat(:,:,j) + w(i)*rr(j,i)*(PP(:,:,j,i) + v*v');
    end
    P_hat(:,:,j) = P_hat(:,:,j)/r_hat(j);
end


end

