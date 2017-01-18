function [ pn,ph,phi,h_r,h_x,h_p ] = group_combs( sets,group_w,group_r,group_x,group_P )

c = allcombs(sets{:});
[N,M] = size(c);

Wupd = zeros(N,1);
num = zeros(N,1);
mbm_upd = cell(N,1);

h_r = zeros(1,0);
h_x = zeros(4,0);
h_p = zeros(4,4,0);
for i = 1:N
    r = zeros(1,0);
    x = zeros(4,0);
    P = zeros(4,4,0);
    wc = 1;
    
    for k = 1:M
        wc = wc*group_w{k}(c(i,k));
        r = [r;group_r{k}(:,c(i,k))];
        x = [x group_x{k}(:,:,c(i,k))];
        P = cat(3,P,group_P{k}(:,:,:,c(i,k)));
    end
    
    h_r = [h_r;r];
    h_x = [h_x x];
    h_p = cat(3,h_p,P);
    
    Wupd(i) = wc;
    num(i) = length(r);
    mbm_upd{i}.r = r;
    mbm_upd{i}.x = x;
    mbm_upd{i}.P = P;
end

maxN = max(num);
Wupd = Wupd/sum(Wupd);
% Limit the number of global association hypotheses
threshold = 100;
if N > threshold
    [~,idxsort]= sort(Wupd,'descend');
    idxkeep=idxsort(1:threshold);
    Wupd= Wupd(idxkeep);
    mbm_upd = mbm_upd(idxkeep);
end
Wupd = Wupd/sum(Wupd);
N = length(Wupd);

idx1 = find(h_r~=0);
[~,idx2] = unique(h_x(1,:));
idx = intersect(idx1,idx2);
h_r = h_r(idx);
h_x = h_x(:,idx);
h_p = h_p(:,:,idx);
% [h_r,h_x,h_p] = pruning(h_r,h_x,h_p,model);

H = length(h_r);
phi = zeros(H,maxN);
for i = 1:H
    for k = 1:maxN
        temp = 0;
        for j = 1:N
            if k<=length(mbm_upd{j}.r)
                if isequal(h_x(1,i),mbm_upd{j}.x(1,k))
                    temp = temp + Wupd(j);
                end
            end
        end
        phi(i,k) = temp;
    end
end
ph = sum(phi,2);
pn = sum(phi,1)';


end

