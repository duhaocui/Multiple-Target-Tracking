function [ r_hat,x_hat,P_hat ] = variational_approx( pn,ph,phi,h_r,h_x,h_p,model )

% Generate cost function
[C,r_hat,x_hat,P_hat] = cost(phi,h_r,h_x,h_p,model);

% Solve linear transport problem
[Cmin,phi] = LP_transport(C,pn,ph);
temp = Cmin;
indicator = 0;
while indicator == 0
    [C,r_temp,x_temp,P_temp] = cost(phi,h_r,h_x,h_p,model);
    [Cmin,phi] = LP_transport(C,pn,ph);
    if temp - Cmin < 1e-3 && temp >= Cmin
        indicator = 1;
        r_hat = r_temp;
        x_hat = x_temp;
        P_hat = P_temp;
    else
        temp = Cmin;
    end
end


end

