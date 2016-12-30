function [ rr, xx, PP ] = make_assign(bestAssign, rupd, xupd, Pupd, ...
    rnew, xnew, Pnew, n, m ,model)
%Make measurement to track assignment
dim = model.x_dim;
M = size(bestAssign,1);
rr = zeros(n+m,M);
xx = zeros(dim,n+m,M);
PP = zeros(dim,dim,n+m,M);

for k = 1:M
    r = zeros(n+m,1);
    x = zeros(dim,n+m);
    P = zeros(dim,dim,n+m);
    for i = 1:n
        if bestAssign(k,i) ~= 0
            r(i) = 1;
            x(:,i) = xupd(:,i,bestAssign(k,i)+1);
            P(:,:,i) = Pupd(:,:,i,bestAssign(k,i)+1);
        else
            r(i) = rupd(i,1);
            x(:,i) = xupd(:,i,1);
            P(:,:,i) = Pupd(:,:,i,1);
        end
    end
    
    for i = n+1:n+m
        if ~any(bestAssign(k,:)==i-n)
            r(i) = rnew(i-n);
            x(:,i) = xnew(:,i-n);
            P(:,:,i) = Pnew(:,:,i-n);
        end
    end
    
    rr(:,k) = r;
    xx(:,:,k) = x;
    PP(:,:,:,k) = P;
    
end

end

            