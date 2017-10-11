function [ Tdot ] = tempodes_phi( t, T );
%Climate Model for ODES

%Bring in global vars from main script
global k L ao pcZaverage sigmaB S0 tau asky gamma A phi1 phi2 ...
        phi3 phi4 phi5 phi6;


% Find Closest time in phi
    [dummy, tphi] = min(abs(phi1(:,1) - t));


Tdot = [(1/pcZaverage(1))*(gamma(1)*(1-asky)*(1-ao(1))*S0*phi1(tphi,2)-tau*sigmaB*T(1)^4)+(L(1,2)*k(1,2)/(A(1)*pcZaverage(1)))*(T(2)-T(1));
    
    (1/pcZaverage(2))*(gamma(2)*(1-asky)*(1-ao(2))*S0*phi2(tphi,2)-tau*sigmaB*T(2)^4)+(1/(A(2)*pcZaverage(2)))*(-L(1,2)*k(1,2)*(T(2)-T(1))+L(2,3)*k(2,3)*(T(3)-T(2)));
    
    (1/pcZaverage(3))*(gamma(3)*(1-asky)*(1-ao(3))*S0*phi3(tphi,2)-tau*sigmaB*T(3)^4)+(1/(A(3)*pcZaverage(3)))*(-L(2,3)*k(2,3)*(T(3)-T(2))+L(3,4)*k(3,4)*(T(4)-T(3)));
    
    (1/pcZaverage(4))*(gamma(4)*(1-asky)*(1-ao(4))*S0*phi4(tphi,2)-tau*sigmaB*T(4)^4)+(1/(A(5)*pcZaverage(5)))*(-L(3,4)*k(3,4)*(T(4)-T(3))+L(4,5)*k(4,5)*(T(5)-T(4)));
    
    (1/pcZaverage(5))*(gamma(5)*(1-asky)*(1-ao(5))*S0*phi5(tphi,2)-tau*sigmaB*T(5)^4)+(1/(A(5)*pcZaverage(5)))*(-L(4,5)*k(4,5)*(T(5)-T(4))+L(5,6)*k(5,6)*(T(6)-T(5)));
    
    (1/pcZaverage(6))*(gamma(6)*(1-asky)*(1-ao(6))*S0*phi6(tphi,2)-tau*sigmaB*T(6)^4)-(L(5,6)*k(5,6)/(A(6)*pcZaverage(6)))*(T(6)-T(5))];

end

