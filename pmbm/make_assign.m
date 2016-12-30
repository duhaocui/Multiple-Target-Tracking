function [ r, x, P ] = make_assign( assign, bestAssign, rupd, xupd, Pupd, ...
    rnew, xnew, Pnew, n, m )
%Make measurement to track assignment
dim = 4;
r = zeros(n+m,1);
x = zeros(dim,n+m);
P = zeros(dim,dim,n+m);

for i=1:n
    if bestAssign(assign,i) ~= 0
        r(i) = 1;
        x(:,i) = xupd(:,i,bestAssign(assign,i)+1);
        P(:,:,i) = Pupd(:,:,i,bestAssign(assign,i)+1);
    else
        r(i) = rupd(i,1);
        x(:,i) = xupd(:,i,1);
        P(:,:,i) = Pupd(:,:,i,1);
    end
end

for i=n+1:n+m
    if ~any(bestAssign(assign,:)==i-n)
        r(i) = rnew(i-n);
        x(:,i) = xnew(:,i-n);
        P(:,:,i) = Pnew(:,:,i-n);
    end
end

end

