function [ w, r, x, p, lambdau, xu, Pu ] = pruning( w, r, x, p, lambdau,xu,Pu,threshold )

if ~isempty(w)
    idxPrune = w > threshold;
    w = w(idxPrune);
    r = r(idxPrune);
    x = x(idxPrune);
    p = p(idxPrune);
    w = normalize(w);
    lambdau = lambdau(idxPrune);
    xu = xu(idxPrune);
    Pu = Pu(idxPrune);
end

end

