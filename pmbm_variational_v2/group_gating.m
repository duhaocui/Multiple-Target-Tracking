function [gatingGroup,gating_mb,idx_out] = group_gating(z,r,x,P,model)

idx_in = [];
zlength = size(z,2);
plength = size(x,2);

tentativeGroup = cell(plength,1);
for j=1:plength
    Sj = model.H*P(:,:,j)*model.H' + model.R;
    nu = z - model.H*repmat(x(:,j),[1 zlength]);
    dist = diag(nu'/Sj*nu);
    idx_in= union(idx_in,find(dist < model.gamma & dist > 0));
    tentativeGroup{j} = find(dist < model.gamma & dist > 0);
end
idx_out = setdiff(1:zlength,idx_in);

% Merging gating groups with common measurements
label = 1:plength;
for i = 1:plength-1
    for j = i+1:plength
        if intersect(tentativeGroup{i},tentativeGroup{j})
            label(j) = min(label(i),label(j));
        end
    end
end
uniqueLabel = unique(label);
n = length(uniqueLabel);

gatingGroup = cell(n,1);
gating_mb = cell(n,1);
for i = 1:n
    gating_mb{i}.r = zeros(0,1);
    gating_mb{i}.x = zeros(4,0);
    gating_mb{i}.P = zeros(4,4,0);
    for j = 1:plength
        if uniqueLabel(i) == label(j)
            gatingGroup{i} = union(gatingGroup{i},tentativeGroup{j});
            gating_mb{i}.r = [gating_mb{i}.r;r(j)];
            gating_mb{i}.x = [gating_mb{i}.x x(:,j)];
            gating_mb{i}.P = cat(3,gating_mb{i}.P,P(:,:,j));
        end
    end
end

end