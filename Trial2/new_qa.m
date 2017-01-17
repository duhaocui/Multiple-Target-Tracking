function [ qa,Cmin ] = new_qa( w,r_hat,x_hat,P_hat,rr_hat,xx_hat,PP_hat )

[N,A,M] = size(rr_hat);
qa = zeros(A,M);
C = zeros(A,M);

for i = 1:M
    for a = 1:A
        for j = 1:N
            v = xx_hat(:,j,a,i) - x_hat(:,j);
            temp = trace(P_hat(:,:,j)\PP_hat(:,:,j,a,i)) + v'/P_hat(:,:,j)*v + log(abs(det(2*pi*P_hat(:,:,j))));
            C(a,i) = C(a,i) - (1-rr_hat(j,a,i))*log(1-r_hat(j)) - rr_hat(j,a,i)*log(r_hat(j)) + rr_hat(j,a,i)/2*temp;
        end
    end
    qa(:,i) = exp(-C(:,i)-logsumexp(-C(:,i)));
end

Cmin = 0;
for i = 1:M
    for a = 1:5
        Cmin = Cmin + w(i)*qa(a,i)*log(qa(a,i)) + w(i)*qa(a,i)*C(a,i);
    end
end


end

