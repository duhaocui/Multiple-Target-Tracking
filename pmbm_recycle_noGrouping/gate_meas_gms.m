function [z_gate,valid_idx,idx_out]= gate_meas_gms(z,model,m,P)

valid_idx = [];
zlength = size(z,2);
plength = size(m,2);

for j=1:plength
    Sj= model.R + model.H*P(:,:,j)*model.H';
    Vs= chol(Sj);
    inv_sqrt_Sj= inv(Vs);
    nu= z- model.H*repmat(m(:,j),[1 zlength]);
    dist= sum((inv_sqrt_Sj'*nu).^2);
    valid_idx= union(valid_idx,find( dist < model.gamma ));
end
z_gate = z(:,valid_idx);
[~,idx_out] = setdiff(z(1,:),z_gate(1,:),'stable');