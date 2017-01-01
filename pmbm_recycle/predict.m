function [r,x,P,lambdau,xu,Pu] = predict(r,x,P,lambdau,xu,Pu,model)
%PREDICT: PREDICT MULTI-BERNOULLI AND POISSON COMPONENTS

% Get multi-Bernoulli prediction parameters from model
F = model.F;
Q = model.Q;
Ps = model.Ps;

% Get birth parameters from model
lambdab = model.lambdab;
xb = model.xb;
Pb = model.Pb;

% Interpret length of inputs
len = length(r);

% Implement prediction algorithm for each global hypothesis
for l = 1:len
    n = length(r{l});
    % Predict existing tracks
    r{l} = Ps*r{l};
    x{l} = F*x{l};
    for i = 1:n
        P{l}(:,:,i) = F*P{l}(:,:,i)*F' + Q;
    end
    
    % Predict existing PPP intensity
    lambdau{l} = Ps*lambdau{l};
    xu{l} = F*xu{l};
    nu = size(xu{l},2);
    for k = 1:nu
        Pu{l}(:,:,k) = F*Pu{l}(:,:,k)*F' + Q;
    end
    
    % Incorporate birth intensity into PPP
    lambdau{l} = [lambdau{l};lambdab];
    xu{l} = [xu{l} xb];
    Pu{l} = cat(3,Pu{l},Pb);
    
    % Truncate low weight components
    ss = lambdau{l} > model.H_threshold;
    lambdau{l} = lambdau{l}(ss);
    xu{l} = xu{l}(:,ss);
    Pu{l} = Pu{l}(:,:,ss);
    
end




