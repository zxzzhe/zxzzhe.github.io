clear;%close all;
k=linspace(0, 9/1000,1e3);
n_r=1.5;n_i=0.0;c=1;
omega1=c*linspace(0, 9/1000, 1e5);
[kk,omega]=meshgrid(k,omega1);
result = zeros( length(omega1),length(k) );

epsilon_p = (n_r+1i*n_i)^2;
epsilon_n = (n_r-1i*n_i)^2;

beta = sqrt(kk.^2-omega.^2./c^2);
alpha_p = sqrt( epsilon_p .* omega.^2./c^2 - kk.^2);
alpha_n = sqrt( epsilon_n .* omega.^2./c^2 - kk.^2);

a=150; b=400;

for ii=1:length(omega1)
    for jj=1:length(k)
        F_p=(beta(ii,jj)+1i*alpha_p(ii,jj))/(beta(ii,jj)-1i*alpha_p(ii,jj));
        F_n=(beta(ii,jj)+1i*alpha_n(ii,jj))/(beta(ii,jj)-1i*alpha_n(ii,jj));
        Gamma_p=(F_p*exp(1i*2*alpha_p(ii,jj)*b)-F_p^(-1));
        Gamma_n=(F_n*exp(1i*2*alpha_n(ii,jj)*b)-F_n^(-1));
        
        result(ii,jj)=abs( (Gamma_p*Gamma_n ) / exp(-4*beta(ii,jj)*a) - (exp(1i*2*alpha_n(ii,jj)*b)-1)*(exp(1i*2*alpha_p(ii,jj)*b)-1) );
        result(ii,jj)=log(result(ii,jj));
    end
end

for ii = 1:length(omega1)
    for jj = 1:length(k)
        if omega(ii,jj)>=c*kk(ii,jj) || omega(ii,jj)<=c*kk(ii,jj)/n_r
            omega(ii,jj)=NaN;
            result(ii,jj)=NaN;
        end
    end
end

figure(2)
pcolor(kk,omega,result);
% surf(kk,omega,result);% view(0,90); %equivalent expression
shading interp;
% colorbar; colormap(jet);caxis([0,1]);
xlabel('$k$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('\rm\fontname{Times New Roman} \rm\fontname{Times New Roman}\omega','FontSize',20)
% set(gca,'LooseInset',[0,0,0,0]);% Cancel the white edge of the picture

%% EP point
%ni=0(hermite situation)  k=810------>omega=609 or 631

ny_EP =  882; %1623; % awave vector
figure(ny_EP)
angular_frequency=omega(:,ny_EP);
result1=result(:,ny_EP);
plot(angular_frequency,result1);
xlabel('$\omega$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('\rm\fontname{Times New Roman} \rm\fontname{Times New Roman}\Omega','FontSize',20)
%%
% for n_EP=2*810:2*815
%     figure(n_EP)
%     angular_frequency=omega(:,n_EP);
%     result1=result(:,n_EP);
%     plot(angular_frequency,result1);
% end

%% distribution of Ey
[peaks,locs] = findpeaks(-result1);
num_EP=length(locs); % number of peaks
if num_EP>=3
    nx1_EP = locs(1); % angular frequency
    nx2_EP = locs(2);
elseif num_EP==2
    nx1_EP = locs(1);
end
%%%%%%%% select peak
nx_EP = nx2_EP;

beta_EP=beta(nx_EP,ny_EP);
alpha_p_EP=alpha_p(nx_EP,ny_EP);
alpha_n_EP=alpha_n(nx_EP,ny_EP);


x1 = linspace(-10*a-b, -a-b, 9*a);
x2 = linspace(-a-b, -a, b);
x3 = linspace(-a, a, 2*a);
x4 = linspace(a, a+b, b);
x5 = linspace(a+b, 10*a+b, 9*a);

A  = exp(1i*0.6*pi);
B1 = A*(1 + beta_EP / (1i*alpha_p_EP) ) / (2 * exp(-1i*alpha_p_EP*b) );
B2 = A*(1 - beta_EP / (1i*alpha_p_EP) ) / (2 * exp( 1i*alpha_p_EP*b) );
C1 = (B1+B2 - 1i*alpha_p_EP*(B1-B2)/beta_EP) / (2 * exp( beta_EP*a) );
C2 = (B1+B2 + 1i*alpha_p_EP*(B1-B2)/beta_EP) / (2 * exp(-beta_EP*a) );
D1 = ( C1*exp(-beta_EP*a) + C2*exp(beta_EP*a)...
     - beta_EP/(1i*alpha_n_EP) * ( C1*exp(-beta_EP*a) - C2*exp(beta_EP*a) ) )/2;
D2 = ( C1*exp(-beta_EP*a) + C2*exp(beta_EP*a)...
     + beta_EP/(1i*alpha_n_EP) * ( C1*exp(-beta_EP*a) - C2*exp(beta_EP*a) ) )/2;
F  = D1*exp(1i*alpha_n_EP*b) + D2 * exp(-1i*alpha_n_EP*b);

Ey1 = F*exp(beta_EP.*(x1+a+b));
Ey2 = D1*exp(-1i.*alpha_n_EP.*(x2+a)) + D2*exp(1i.*alpha_n_EP.*(x2+a));
Ey3 = C1*exp(beta_EP.*(x3)) + C2*exp(-beta_EP.*(x3));
Ey4 = B1*exp(-1i.*alpha_p_EP.*(x4-a)) + B2*exp(1i.*alpha_p_EP.*(x4-a));
Ey5 = A*exp(-beta_EP.*(x5-a-b));

x = [x1,x2,x3,x4,x5];
Ey = [Ey1,Ey2,Ey3,Ey4,Ey5];

maxnum = max(abs(Ey));

% 
A = A/maxnum;
B1 = A*(1 + beta_EP / (1i*alpha_p_EP) ) / (2 * exp(-1i*alpha_p_EP*b) );
B2 = A*(1 - beta_EP / (1i*alpha_p_EP) ) / (2 * exp( 1i*alpha_p_EP*b) );
C1 = (B1+B2 - 1i*alpha_p_EP*(B1-B2)/beta_EP) / (2 * exp( beta_EP*a) );
C2 = (B1+B2 + 1i*alpha_p_EP*(B1-B2)/beta_EP) / (2 * exp(-beta_EP*a) );
D1 = ( C1*exp(-beta_EP*a) + C2*exp(beta_EP*a)...
     - beta_EP/(1i*alpha_n_EP) * ( C1*exp(-beta_EP*a) - C2*exp(beta_EP*a) ) )/2;
D2 = ( C1*exp(-beta_EP*a) + C2*exp(beta_EP*a)...
     + beta_EP/(1i*alpha_n_EP) * ( C1*exp(-beta_EP*a) - C2*exp(beta_EP*a) ) )/2;
F  = D1*exp(1i*alpha_n_EP*b) + D2 * exp(-1i*alpha_n_EP*b);
Ey1 = F*exp(beta_EP.*(x1+a+b));
Ey2 = D1*exp(-1i.*alpha_n_EP.*(x2+a)) + D2*exp(1i.*alpha_n_EP.*(x2+a));
Ey3 = C1*exp(beta_EP.*(x3)) + C2*exp(-beta_EP.*(x3));
Ey4 = B1*exp(-1i.*alpha_p_EP.*(x4-a)) + B2*exp(1i.*alpha_p_EP.*(x4-a));
Ey5 = A*exp(-beta_EP.*(x5-a-b));
Ey = [Ey1,Ey2,Ey3,Ey4,Ey5];

figure(5)
subplot(2,2,1)
abs_Ey = abs(Ey);
plot(x,abs_Ey);
xlabel('$x(\mathrm{nm})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$\mathrm{abs}(E_y)$','interpreter','latex','FontName','Times New Roman','FontSize',20)
% ylabel('\rm\fontname{Times New Roman}abs \rm\fontname{Times New Roman}E_y','FontSize',20)

% real part && imag part
figure(5)
subplot(2,2,2)
real_Ey = real(Ey);
plot(x,real_Ey);
xlabel('$x(\mathrm{nm})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$\mathrm{Re}(E_y)$','interpreter','latex','FontName','Times New Roman','FontSize',20)

figure(5)
subplot(2,2,3)
imag_Ey = imag(Ey);
plot(x,imag_Ey);
xlabel('$x(\mathrm{nm})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$\mathrm{Im}(E_y)$','interpreter','latex','FontName','Times New Roman','FontSize',20)

% phase part
figure(5)
subplot(2,2,4)
phase_Ey = angle(Ey)/pi;
plot(x,phase_Ey);
xlabel('$x(\mathrm{nm})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$\mathrm{phase}(\pi)$','interpreter','latex','FontName','Times New Roman','FontSize',20)


%% Poynting vector

Hx1 = -kk(nx_EP,ny_EP)/omega(nx_EP,ny_EP)*(Ey1);
Hx2 = -kk(nx_EP,ny_EP)/omega(nx_EP,ny_EP)*(Ey2);
Hx3 = -kk(nx_EP,ny_EP)/omega(nx_EP,ny_EP)*(Ey3);
Hx4 = -kk(nx_EP,ny_EP)/omega(nx_EP,ny_EP)*(Ey4);
Hx5 = -kk(nx_EP,ny_EP)/omega(nx_EP,ny_EP)*(Ey5);
Hx = -kk(nx_EP,ny_EP)/omega(nx_EP,ny_EP)*Ey;

figure('numbertitle','off','name','Poynting vector'); 

sz = Ey.*conj(Hx);
% abs_Sz  = abs(Sz);
real_sz = real(sz);
imag_sz = imag(sz);

hold on
% plot(x,abs(Sz));
plot(x,real_sz);
plot(x,imag_sz);
% legend('$\mathrm{abs}(S_z)$','$\mathrm{Re}(S_z)$','$\mathrm{Im}(S_z)$','interpreter','latex','FontName','Times New Roman','FontSize',20)
legend('$\mathrm{Re}(S_z)$','$\mathrm{Im}(S_z)$','interpreter','latex','FontName','Times New Roman','FontSize',20)
title('$S_z=E_y\cdot H_x^{*}$','interpreter','latex','FontName','Times New Roman','FontSize',20)
xlabel('$x(\mathrm{nm})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$S_z$','interpreter','latex','FontName','Times New Roman','FontSize',20)
hold off

figure('numbertitle','off','name','Poynting vector'); 

sz = Ey.*(Hx);
% abs_Sz  = abs(Sz);
real_sz = real(sz);
imag_sz = imag(sz);

hold on
% plot(x,abs(Sz));
plot(x,real_sz);
plot(x,imag_sz);
% legend('$\mathrm{abs}(S_z)$','$\mathrm{Re}(S_z)$','$\mathrm{Im}(S_z)$','interpreter','latex','FontName','Times New Roman','FontSize',20)
legend('$\mathrm{Re}(S_z)$','$\mathrm{Im}(S_z)$','interpreter','latex','FontName','Times New Roman','FontSize',20)
title('$S_z=E_y\cdot H_x$','interpreter','latex','FontName','Times New Roman','FontSize',20)
xlabel('$x(\mathrm{nm})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$S_z$','interpreter','latex','FontName','Times New Roman','FontSize',20)
hold off

%% distribution of Hz
Hz1 = beta_EP*F*exp(beta_EP.*(x1+a+b));
Hz2 = -1i*alpha_n_EP*( D1*exp(-1i.*alpha_n_EP.*(x2+a)) - D2*exp(1i.*alpha_n_EP.*(x2+a)) );
Hz3 = beta_EP*( C1*exp(beta_EP.*(x3)) - C2*exp(-beta_EP.*(x3)) );
Hz4 = -1i*alpha_p_EP*( B1*exp(-1i.*alpha_p_EP.*(x4-a)) - B2*exp(1i.*alpha_p_EP.*(x4-a)) );
Hz5 = -beta_EP*A*exp(-beta_EP.*(x5-a-b));
Hz1 = 1i/omega(nx_EP,ny_EP)*Hz1;
Hz2 = 1i/omega(nx_EP,ny_EP)*Hz2;
Hz3 = 1i/omega(nx_EP,ny_EP)*Hz3;
Hz4 = 1i/omega(nx_EP,ny_EP)*Hz4;
Hz5 = 1i/omega(nx_EP,ny_EP)*Hz5;
Hz = [Hz1,Hz2,Hz3,Hz4,Hz5];
real_Hz = real(Hz);
imag_Hz = imag(Hz);
abs_Hz = abs(Hz);
figure(2222)
hold on
% plot(x,real_Hz);
plot(x,imag_Hz);
% plot(x,abs_Hz);
hold off

%% group v
w1 = Ey1.*conj(Ey1) + Hx1.*conj(Hx1);
w2 = epsilon_n*Ey2.*conj(Ey2) + Hx2.*conj(Hx2);
w3 = Ey3.*conj(Ey3) + Hx3.*conj(Hx3);
w4 = epsilon_p*Ey4.*conj(Ey4) + Hx4.*conj(Hx4);
w5 = Ey5.*conj(Ey5) + Hx5.*conj(Hx5);
w = [w1,w2,w3,w4,w5];
W = trapz(x,w);
Sz = trapz(x,sz);
vg = Sz/W;
