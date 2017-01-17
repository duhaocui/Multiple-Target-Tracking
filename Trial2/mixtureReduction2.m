function [ r_hat,x_hat,P_hat ] = mixtureReduction2( w,qa,rr_hat,xx_hat,PP_hat )

[N,A,M] = size(rr_hat);
r_hat = zeros(N,1);
x_hat = zeros(4,N);
P_hat = zeros(4,4,N);

for j = 1:N
    for i = 1:M
        for a = 1:A
            r_hat(j) = r_hat(j) + w(i)*qa(a,i)*rr_hat(j,a,i);
            x_hat(:,j) = x_hat(:,j) + w(i)*qa(a,i)*rr_hat(j,a,i)*xx_hat(:,j,a,i);
        end
    end
    x_hat(:,j) = x_hat(:,j)/r_hat(j);
    for i = 1:M
        for a = 1:A
            v = xx_hat(:,j,a,i) - x_hat(:,j);
            P_hat(:,:,j) = P_hat(:,:,j) + w(i)*qa(a,i)*rr_hat(j,a,i)*(PP_hat(:,:,j,a,i) + v*v');
        end
    end
    P_hat(:,:,j) = P_hat(:,:,j)/r_hat(j);
end


end

