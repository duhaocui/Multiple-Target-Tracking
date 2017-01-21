function [ lambdau,xu,Pu ] = recycling( rr,xx,PP,lambdau,xu,Pu,model )

idx = rr <= 0.1;
lambdau = [lambdau rr(idx)'];
xu = [xu xx(:,idx)];
Pu = cat(3,Pu,PP(:,:,idx));



end

