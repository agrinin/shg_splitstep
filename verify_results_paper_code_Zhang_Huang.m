%verify the result of the paper in order to chech the correctness of the
%code.
%for more details s. "Zhang_Second Harmonic generation from regeneratively
%amplified femto-second laser pulsed in BBO and LBO crystals" eq (4-5)

%% define fundamental constants 
% BBO crystal constants for lambda1 = 800nm at  T = 298.15K
%natural constants
fpc = fundamentalPhysicalConstantsFromNIST();
c = fpc.speed_of_light_in_vacuum.value;
eps0 = fpc.electric_constant.value;
mu0 = fpc.mag_constant.value;
%% pulse constants
tau0 = sqrt(2)*80e-15;
lambda1 = 410e-9;       %wavelength of the  fundamental beam
omega1 = 2*pi*c/lambda1; %angular frequency of the fundamental wave              
%crystal and SHG process constants
n1 = 1.660;%1.610365;%%;%1.691;%%;          %index of refraction at the fundamental frequency
n2 = 1.660;%1.610365;%1.660;%1.610365;%1.691;%1.660;%1.610365;          %index of refraction at the SH
U0 = 1.5/315e6;
w0x = 10e-6; w0y = 10e-6;
E10 = sqrt(8*sqrt(log(2))*U0/(eps0*c*pi^(3/2)*w0x*w0y*tau0)); %peak electric field
I = 1/2*E10^2*c*eps0;

% I = 7e13;
% E10 = 23.7648*sqrt(I/tau0)*sqrt(tau0);%sqrt(2*I/(c*eps0));%18.9796
% E10 = sqrt(2*I/(c*eps0));
lw = c*tau0/(1.775 - 2.165);% c*tau0/(1.683 -1.737);% c*tau0/(1.631 -1.665);%c*tau0/(1.683 -1.737);c*tau0/(1.775 - 2.165)%;%c*tau0/(1.684-1.742);     %c*tau0/(1.631-1.665);   %walk-off distance
gn1 = 1.03146*10^-25;%BBO 410nm%3.75883*10^-26;%1.0773*10^-25;%;;% 
gn2 = 4.59531*10^-25;%BBO 410nm%7.31414*10^-26;%5.18822*10^-25;%%
ld1 = tau0^2/(4*gn1); %0.086;         %pulse spreading distance for 410nm in BBO
ld2 = tau0^2/(4*gn2); % 0.043;         %pulse spreading distance for 205nm in LBO
deff = 2.64e-13;%0.761e-12;%2e-12;%0.264e-12;%;%7.61E-13; %effective nonlinear coeff. in [m/V]
xi2 = deff;      %second order nonlinear susceptibility
zI = 2*n1*c/(omega1*xi2*E10);%   %interaction length of the SH
n12 = 0.32e-19;     %second order index of refraction of LBO
b = 0.5;              %linear chirp coefficient
%constants for the differential equations
kappa0 = 1;%0.35;
kappa1 = 2*pi*n12*E10^2*zI/lambda1*eps0*c;
k = [zI lw ld1 ld2  kappa0 kappa1];

%% calculate
N = 2^7;
tau = linspace(-30,30,N);
Ein = zeros(4,N);
for i=1:N
Ein(1,i) = real(exp(-(1+1i*b).*tau(i).^2));  %real part of the fundam. field
Ein(2,i) = imag(exp(-(1+1i*b).*tau(i).^2));  %imag. part of the fund. field
Ein(3,i) = 0;                                %real part of the sh field
Ein(4,i) = 0;                                %imag. part of the sh field
end


%% propagate the beam
L  = 1500e-6;
Nl = 100;
l = L/zI;
dl = L/zI*1/Nl;

[zres, Eres, nu] = splitstep(k,tau,Ein,l,dl);
zres = zres*zI;
%% interpolate results
tauq = linspace(min(tau),max(tau),2e3);
[m,n] = size(Eres);

for i = 1:1:m
Eres{i} = [interp1(tau,Eres{i}(1,:),tauq,'spline');...
           interp1(tau,Eres{i}(2,:),tauq,'spline');...
           interp1(tau,Eres{i}(3,:),tauq,'spline');...
           interp1(tau,Eres{i}(4,:),tauq,'spline')];
end


%% plot for diff L
figure; 
subplot(1,2,1);
plot(tau,(Ein(1,:).^2+Ein(2,:).^2));
hold on;
nu1 = zeros([m,n]);



for i = 1:1:m
subplot(1,2,1);
plot(tauq,interp1(tauq,(Eres{i}(1,:).^2+Eres{i}(2,:).^2),tauq,'spline'),'r');  
plot(tauq,interp1(tauq,20*(Eres{i}(3,:).^2+20*Eres{i}(4,:).^2),tauq,'spline'),'m');  
subplot(1,2,2);
plot(tau,angle(Ein(1,:)+1i*Ein(2,:)));
hold on;
% plot(tauq,(angle(Eres{i}(1,:)+1i*Eres{i}(2,:))),'r');  
plot(tauq,(angle(Eres{i}(3,:)+1i*Eres{i}(4,:))),'m');  
nu1(i) = trapz(tauq,(Eres{i}(3,:).^2+Eres{i}(4,:).^2))./trapz(tau,(Ein(1,:).^2+Ein(2,:).^2));
end
xlabel('time in taus');
ylabel('bnorm amplitude');
legend('input pulse','output fundamental','shg pulse');

%%
figure;
plot(zres*1e3,nu1,'*')
%%
ylim([0 1]); xlim([0 1]);
xlabel('Crystal Length in mm');
ylabel('Efficiency');
title('BBO SHG at 820nm and 80fs, Fourier limited')
lgd = legend('5GW/cm^2', '12GW/cm^2', '47GW/cm^2','190GW/cm^2');
lgd.Location = 'best';
set(gca,'FontSize',18);
orient(gcf,'Landscape');
set(gcf,'PaperPositionMode','auto'); 
set(gcf, 'PaperSize', [11.5 6.7]);

%%
print(gcf, '-dpdf', 'BBOb0820nm80fsEff.pdf')

%%
subplot(1,2,1);
xlabel('time in \tau_0');
ylabel('Normalized Intensity [a.u.]');
title('190GW/cm^2,BBO SHG at 820nm and 80fs, FL')
set(gca,'FontSize',18);
orient(gcf,'Landscape');
set(gcf,'PaperPositionMode','auto'); 
set(gcf, 'PaperSize', [11.5 6.7]);
%%
subplot(1,2,2);
xlabel('time in \tau_0');
ylabel('Phase \phi');
title('Phase')
set(gca,'FontSize',18);
orient(gcf,'Landscape');
set(gcf,'PaperPositionMode','auto'); 
set(gcf, 'PaperSize', [11.5 6.7]);

%%
print(gcf, '-dpdf', 'BBOb0820nm80fsAmpPhase5um.pdf')

