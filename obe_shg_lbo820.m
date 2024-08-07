function dE = obe_shg_lbo820(k,t,E)
%define characteristic constants
kappa0 = k(1); %kappa0 = n1/n2;
kappa1 = k(2); %1.0386;%4.1651; kappa1 = 2*pi*n12*E10^2*zI/lambda1

dE = zeros(4,1);    % difference elec. field vector
%E(1) real part of the fundam. field
%E(2) imag. part of the fund. field
%E(3) real part of the sh field
%E(4) imag. part of the sh field

dE(1) = (E(2)*E(3)-E(1)*E(4))-kappa1*(3/4*(E(1)^2+E(2)^2)+3/2*(E(3)^2+E(4)^2))*E(2);
dE(2) = (E(1)*E(3)+E(2)*E(4))+kappa1*(3/4*(E(1)^2+E(2)^2)+3/2*(E(3)^2+E(4)^2))*E(1);
dE(3) = -4*kappa0*E(1)*E(2)-2*kappa1*(3/4*(E(3)^2+E(4)^2)+3/2*(E(1)^2+E(2)^2))*E(4);
dE(4) = 2*kappa0*(E(1)^2-E(2)^2)+2*kappa1*(3/4*(E(3)^2+E(4)^2)+3/2*(E(1)^2+E(2)^2))*E(3);

end