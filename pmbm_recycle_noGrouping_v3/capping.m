function [ w, r, x, P ,lambdau,xu,Pu] = capping( w, r, x, P, lambdau,xu,Pu,threshold )

if length(w) > threshold
    [~,idxsort]= sort(w,'descend');
    idxkeep=idxsort(1:threshold);
    w= w(idxkeep);
    r = r(idxkeep);
    x = x(idxkeep);
    P = P(idxkeep);
    lambdau = lambdau(idxkeep);
    xu = xu(idxkeep);
    Pu = Pu(idxkeep);
end
w = normalize(w);

end

