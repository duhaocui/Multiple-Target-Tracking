function [ rr_new,xx_new,PP_new,qa,Cmin ] = expectation( w,rr,xx,PP,r_hat,x_hat,P_hat )
warning('off','all');

[N,M] = size(rr);
C = zeros(N);
rr_new = zeros(N,5,M);
xx_new = zeros(4,N,5,M);
PP_new = zeros(4,4,N,5,M);
qa = zeros(5,M);
nCost = zeros(5,M);

for i = 1:M
    for h = 1:N
        for j = 1:N
            if rr(h,i)==0
                C(h,j) = -log(1-r_hat(j));
            else
                v = xx(:,h,i) - x_hat(:,j);
                temp = trace(P_hat(:,:,j)\PP(:,:,h,i)) + v'/P_hat(:,:,j)*v + log(abs(det(2*pi*P_hat(:,:,j))));
                C(h,j) = -(1-rr(h,i))*log(1-r_hat(j)) - rr(h,i)*log(r_hat(j)) + rr(h,i)/2*temp;
            end
        end
    end
    C(isnan(C)) = 1/eps;
    x = min(min(C));
    C = C - x;
    [assignments, nCost(:,i)] = murty_custom([abs(C) inf(N,N)],5);
    nCost(:,i) = nCost(:,i) + (x.*sum(assignments>0,2));
    qa(:,i) = exp(-nCost(:,i)-logsumexp(-nCost(:,i)))';
    r = rr(:,i);x = xx(:,:,i);P = PP(:,:,:,i);
    for a = 1:5
        rr_new(:,a,i) = r(assignments(a,:));
        xx_new(:,:,a,i) = x(:,assignments(a,:));
        PP_new(:,:,:,a,i) = P(:,:,assignments(a,:));
    end
end

Cmin = 0;
for i = 1:M
    for a = 1:5
        Cmin = Cmin + w(i)*qa(a,i)*log(qa(a,i)) + w(i)*qa(a,i)*nCost(a,i);
    end
end

end

