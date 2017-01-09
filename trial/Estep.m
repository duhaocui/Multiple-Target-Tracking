function [ w_hat,Cmin ] = Estep( r_hat,x_hat,P_hat,rr,xx,PP,pj,w )

r_hat(r_hat>=1) = 1-eps;
r_hat(r_hat<=eps) = eps;
[N,M] = size(rr);
w_hat = zeros(M,1);
C = zeros(M,1);
Cmin = 0;
for i = 1:M
    for j = 1:N
        v = xx(:,j,i) - x_hat(:,j);
        temp = trace(P_hat(:,:,j)\PP(:,:,j,i)) + v'*inv(abs(P_hat(:,:,j)))*v + log(abs(det(2*pi*P_hat(:,:,j))));
        C(i) = C(i) + (1-rr(j,i))*log(1-r_hat(j)) + rr(j,i)*log(r_hat(j)) - rr(j,i)/2*temp;
    end
    w_hat(i) = exp(sum(C));
    Cmin = Cmin + sum(C)*w(i);
end
w_hat = w_hat/sum(w_hat);
Cmin = -Cmin;

end

