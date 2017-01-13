function [ rr_new,xx_new,PP_new,Cmin ] = expectation( w,rr,xx,PP,r_hat,x_hat,P_hat )
warning('off','all');

[N,M] = size(rr);
C = zeros(N);
rr_new = zeros(N,M);
xx_new = zeros(4,N,M);
PP_new = zeros(4,4,N,M);
costs = zeros(M,1);
nCost = zeros(M,1);

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
    nCost(i) = sum(diag(C));
    C(isnan(C)) = 1/eps;
    x = min(min(C));
    C = C - x;
    [assignments, costs(i)] = murty_custom(C,1);
    costs(i) = costs(i) + (x.*sum(assignments>0,2))';
    r = rr(:,i);x = xx(:,:,i);P = PP(:,:,:,i);
    rr_new(:,i) = r(assignments);
    xx_new(:,:,i) = x(:,assignments);
    PP_new(:,:,:,i) = P(:,:,assignments);
end

Cmin = sum(w.*costs);

end

