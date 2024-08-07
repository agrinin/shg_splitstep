%splitstep integrates the differential equations for propagation via a
%medium with first, second and third order nonlinearity, i.e. second
%harmonic generation, selfphase and crossphase modulation and simple
%dispersion
%for more details s. "Zhang_Second Harmonic generation from regeneratively
%amplified femto-second laser pulsed in BBO and LBO crystals" eq (4-5)
%Inputs: E1r+iE2i is the fundamental field at the position l. It is a 
%        vector in time domain
%        E2r+iE2i is the second harmonic field at position l.
%        l is the current positon in the crystal
%        dl is the distance step size


%Author: Alexey Grinin
%Date: 19.06.2015


function [zres, Eres, nu] = splitstep(k,tau,Ein,l,dl)
% calculate the fields at the position l+dl using the fields at l. To this
% end integrate the ode's with Runga-Kutta


% propagate through crystal of length l by dl
% 1) solve nonlinear part using l+dl/2 step
% 2) Solve dispersion part using dl and fft/ifft
% 3) solve nonlinear part after dispesion l+dl/2->l+dl
% 4) write down every field and current crystal position

% how many steps through crystal
Npos = ceil(l/dl);
Eres = cell([Npos,1]);
zres = zeros([Npos,1]);
nu = zeros([Npos,1]);
Ntau = length(tau); 
fmax = 1/mean(diff(tau));
f = linspace(-fmax/2,fmax/2,21*Ntau); %21*Ntau
zI = k(1);
lw = k(2);
ld1 = k(3);
ld2 = k(4);
kappa0 = k(5);
kappa1 = k(6);
bar = waitbar(0,'Solving ODEs in progress...');
options = odeset('RelTol',1e-7);

% propage for the first time
Eout = zeros(4,length(Ein));
for j = 1:length(Ein)   
    [z,E] = ode15s(@(t,E) obe_shg_lbo820([kappa0 kappa1],t,E),[0 dl/10],Ein(:,j),options);    
    Eout(1,j) = E(end,1);
    Eout(2,j) = E(end,2);
    Eout(3,j) = E(end,3);
    Eout(4,j) = E(end,4);
end
for i = 1:Npos
    waitbar(i/Npos,bar,'Solving ODEs in progress...');  
    % propage by the first half dl/2 step
    for j = 1:length(Ein)   
        [z,E] = ode15s(@(t,E) obe_shg_lbo820([kappa0 kappa1],t,E),[dl*(i-1) dl*(i-1)+dl/2],Eout(:,j),options);    
        Eout(1,j) = E(end,1);
        Eout(2,j) = E(end,2);
        Eout(3,j) = E(end,3);
        Eout(4,j) = E(end,4);
    end
    % solve the dispersion part with fft
    E1 = padarray(Eout(1,:)+1i.*Eout(2,:),[0 Ntau*10]);   
    E2 =  padarray(Eout(3,:)+1i.*Eout(4,:),[0 Ntau*10]);   
%     fourier transform the solution 
    Ef1 = fftshift(fft(E1)).*exp(2*pi*1i.*((f*zI/lw)+pi*f.^2*zI/(2*ld1))*dl);
    Ef2 = fftshift(fft(E2)).*exp(2*pi*1i.*(pi*f.^2*zI/(2*ld2))*dl);
%     inverse fourier transform
    E1 = ifft(ifftshift(Ef1)); E2 = ifft(ifftshift(Ef2));
    E1 = E1(Ntau*10+1:Ntau*11); E2 = E2(Ntau*10+1:Ntau*11); 
%     create real vector
    Eout = [real(E1); imag(E1); real(E2); imag(E2)];
    % propage by the second half dl/2 step
    for j = 1:length(Ein)   
        options = odeset('RelTol',1e-7);
        [z,E] = ode15s(@(t,E) obe_shg_lbo820(k,t,E),[dl*(i-1)+dl/2 dl*i],Eout(:,j),options);    
        Eout(1,j) = E(end,1);
        Eout(2,j) = E(end,2);
        Eout(3,j) = E(end,3);
        Eout(4,j) = E(end,4);
    end

    zres(i) = z(end);
    Eres{i} = Eout;
    nu(i) = max((Eout(3,:).^2+Eout(4,:).^2))./max((Ein(1,:).^2+Ein(2,:).^2));

    
end
close(bar);
end


       
