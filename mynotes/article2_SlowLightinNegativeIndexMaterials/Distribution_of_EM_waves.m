clear;%close all;
epsilon_p = 4; mu_p = 1;

F = 0.56; omega_p = 2*pi*10e9; omega_0 = 2*pi*4e9; omega = 2*pi*5e9;
epsilon_n = 1 - (omega_p/omega)^2;
mu_n = 1 - F*omega^2/(omega^2-omega_0^2);

mu_r = 1;
mu_0 = 4*pi*1e-7;
epsilon_0 = 8.85e-12;
c = 3e8;

b = 4e-2;
dist_a = linspace(0, 5e-2/2, 1e2+1);
k = linspace(105, 145, 1e2+1);

[aa,kk]=meshgrid(dist_a,k);
result = zeros( length(k),length(dist_a) );

beta = sqrt(kk.^2-omega.^2./c^2);
alpha_p = sqrt( epsilon_p .* mu_p .* omega.^2./c^2 - kk.^2);
alpha_n = sqrt( epsilon_n .* mu_n .* omega.^2./c^2 - kk.^2);

for i=1:length(k)
    for j=1:length(dist_a)
        H=zeros(8,8);
        H(1,1)=1;H(1,2)=-exp(-1i*alpha_n(i,j)*b);H(1,3)=-exp(1i*alpha_n(i,j)*b);
        H(2,2)=1;H(2,3)=1;H(2,4)=-exp(beta(i,j)*aa(i,j));H(2,5)=-exp(-beta(i,j)*aa(i,j));
        H(3,4)=exp(-beta(i,j)*aa(i,j));H(3,5)=exp(beta(i,j)*aa(i,j));H(3,6)=-1;H(3,7)=-1;
        H(4,6)=exp(1i*alpha_p(i,j)*b);H(4,7)=exp(-1i*alpha_p(i,j)*b);H(4,8)=-1;

        H(5,1)=-beta(i,j);H(5,2)=1i*alpha_n(i,j)*exp(-1i*alpha_n(i,j)*b)/mu_n;H(5,3)=-1i*alpha_n(i,j)*exp(1i*alpha_n(i,j)*b)/mu_n;
        H(6,2)=-1i*alpha_n(i,j)/mu_n;H(6,3)=1i*alpha_n(i,j)/mu_n;H(6,4)=-beta(i,j)*exp(beta(i,j)*aa(i,j));H(6,5)=beta(i,j)*exp(-beta(i,j)*aa(i,j));
        H(7,4)=beta(i,j)*exp(-beta(i,j)*aa(i,j));H(7,5)=-beta(i,j)*exp(beta(i,j)*aa(i,j));H(7,6)=1i*alpha_p(i,j)/mu_p;H(7,7)=-1i*alpha_p(i,j)/mu_p;
        H(8,6)=-1i*alpha_p(i,j)*exp(1i*alpha_p(i,j)*b)/mu_p;H(8,7)=1i*alpha_p(i,j)*exp(-1i*alpha_p(i,j)*b)/mu_p;H(8,8)=-beta(i,j);
        
        result(i,j) = log(abs(det(H)));
        result(i,j) = log(result(i,j)+1);
    end
end

figure(2)
pcolor(dist_a*100*2,k/100,result);
% surf(kk,omega,result);% view(0,90); %equivalent expression
shading interp;
% colorbar; colormap(jet);caxis([0,1]);
xlabel('$a\ (\rm{cm})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$k\ (\rm{cm^{-1}})$','interpreter','latex','FontSize',20)
% set(gca,'LooseInset',[0,0,0,0]);% Cancel the white edge of the picture

%% EP point
%%%%% 809
ny = 30 ;
figure(ny)
k_list = kk(:,ny);
result1=result(:,ny);
plot(k_list,result1);
xlabel('$k$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('result','interpreter','latex','FontSize',20)
%%
for n_EP=20:1:30
    figure(n_EP)
    k_list=kk(:,n_EP);
    result1=result(:,n_EP);
    plot(k_list,result1);
end
%% save curve
num_be = 100; % set the begin number of data a
k2a_A = zeros( 1, length(dist_a) );
k2a_B = zeros( 1, length(dist_a) );
k2a_A(1:num_be) = NaN;
k2a_B(1:num_be) = NaN;

for mm = num_be+1:length(dist_a)
    result1=result(:,mm);
    [~,locs] = findpeaks(-result1);
    if any(locs==756)
        locs(locs==756) = [];
%     else
%         print('error')
%         break;
    end
    if isempty(locs)
        k2a_A(mm) = NaN;
        k2a_B(mm) = NaN;
        length_loc = 0;      %%%%%% length_loc, to record last locs's length
    elseif length(locs) == 1
        k2a_A(mm) = k( locs(1) );
        if length_loc == 1 || length_loc == 0
            k2a_B(mm) = NaN;
            length_loc = 1;
        elseif length_loc == 2
            k2a_B(mm) = k( 756 );
            length_loc = 2;
        end
    elseif length(locs) == 2
        k2a_A(mm) = k( locs(1) );
        k2a_B(mm) = k( locs(2) );
        length_loc = 2;
    end
    
end
num_mid = find(k2a_B>0,1) - 1 ;
k2a_B(num_mid) = k2a_A(num_mid);

figure('numbertitle','off','name','save curve');
hold on
plot( dist_a*100*2, k2a_A/100 );
plot( dist_a*100*2, k2a_B/100 );
legend('$\mathrm{curve}\ A$','$\mathrm{curve}\ B$','interpreter','latex','FontName','Times New Roman','FontSize',20,'Location','southwest')
% title('$E_y\ \mathrm{and}\ S_z\ \mathrm{at}\ A$','interpreter','latex','FontName','Times New Roman','FontSize',20)
xlabel('$a\ (\mathrm{cm})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$k\ (\mathrm{cm}^{-1})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
xlim([0,5])
ylim([1.10,1.39])
hold off

%% distribution of Ey at A
[~,locs] = findpeaks(-result1);
num_EP=length(locs); % number of peaks
if num_EP>=3
    nx1_EP = locs(1); % disatance_a
    nx2_EP = locs(2);
elseif num_EP==2
    nx1_EP = locs(1);
end
%%%%%%%% select peak
nx_EP = nx1_EP;

% PIM or NIM
beta_EP=beta(nx_EP, ny);
alpha_p_EP=alpha_p(nx_EP,ny);
alpha_n_EP=alpha_n(nx_EP,ny);
a = aa(nx_EP, ny);

x1 = linspace(-4*a-b, -a-b, 50);
x2 = linspace(-a-b, -a, 200);
x3 = linspace(-a, a, 50);
x4 = linspace(a, a+b, 200);
x5 = linspace(a+b, 4*a+b, 50);

A  = 1;
B1 = A*(1i*alpha_n_EP/mu_n + beta_EP/mu_r ) / (2*1i * alpha_n_EP/mu_n * exp(-1i*alpha_n_EP*b) );
B2 = A*(1i*alpha_n_EP/mu_n - beta_EP/mu_r ) / (2*1i * alpha_n_EP/mu_n * exp(+1i*alpha_n_EP*b) );
C1 = ( ( beta_EP/mu_r - 1i*alpha_n_EP/mu_n )*B1 + ( beta_EP/mu_r + 1i*alpha_n_EP/mu_n )*B2 ) / ( 2*beta_EP/mu_r * exp(+beta_EP * a) );
C2 = ( ( beta_EP/mu_r + 1i*alpha_n_EP/mu_n )*B1 + ( beta_EP/mu_r - 1i*alpha_n_EP/mu_n )*B2 ) / ( 2*beta_EP/mu_r * exp(-beta_EP * a) );
D1 = ( ( 1i*alpha_p_EP/mu_p*exp(-beta_EP*a) - beta_EP/mu_r*exp(-beta_EP*a) )*C1 ...
     + ( 1i*alpha_p_EP/mu_p*exp(+beta_EP*a) + beta_EP/mu_r*exp(+beta_EP*a) )*C2 )...
     / ( 2*1i*alpha_p_EP/mu_p );
D2 = ( ( 1i*alpha_p_EP/mu_p*exp(-beta_EP*a) + beta_EP/mu_r*exp(-beta_EP*a) )*C1 ...
     + ( 1i*alpha_p_EP/mu_p*exp(+beta_EP*a) - beta_EP/mu_r*exp(+beta_EP*a) )*C2 )...
     / ( 2*1i*alpha_p_EP/mu_p );
F  = D1*exp(1i*alpha_p_EP*b) + D2 * exp(-1i*alpha_p_EP*b);

Ey1 = F*exp(beta_EP.*(x1+a+b));
Ey2 = D1*exp(-1i.*alpha_p_EP.*(x2+a)) + D2*exp(1i.*alpha_p_EP.*(x2+a));
Ey3 = C1*exp(beta_EP.*(x3)) + C2*exp(-beta_EP.*(x3));
Ey4 = B1*exp(-1i.*alpha_n_EP.*(x4-a)) + B2*exp(1i.*alpha_n_EP.*(x4-a));
Ey5 = A*exp(-beta_EP.*(x5-a-b));

x = [x1,x2,x3,x4,x5];
Ey = [Ey1,Ey2,Ey3,Ey4,Ey5];

%%%%%% normalize
% maxnum = max(abs(Ey)); 
% A = A/maxnum;
% B1 = A*(1 + beta_EP / (1i*alpha_n_EP) ) / (2 * exp(-1i*alpha_n_EP*b) );
% B2 = A*(1 - beta_EP / (1i*alpha_n_EP) ) / (2 * exp( 1i*alpha_n_EP*b) );
% C1 = (B1+B2 - 1i*alpha_n_EP*(B1-B2)/beta_EP) / (2 * exp( beta_EP*a) );
% C2 = (B1+B2 + 1i*alpha_n_EP*(B1-B2)/beta_EP) / (2 * exp(-beta_EP*a) );
% D1 = ( C1*exp(-beta_EP*a) + C2*exp(beta_EP*a)...
%      - beta_EP/(1i*alpha_p_EP) * ( C1*exp(-beta_EP*a) - C2*exp(beta_EP*a) ) )/2;
% D2 = ( C1*exp(-beta_EP*a) + C2*exp(beta_EP*a)...
%      + beta_EP/(1i*alpha_p_EP) * ( C1*exp(-beta_EP*a) - C2*exp(beta_EP*a) ) )/2;
% F  = D1*exp(1i*alpha_p_EP*b) + D2 * exp(-1i*alpha_p_EP*b);
% 
% Ey1 = F*exp(beta_EP.*(x1+a+b));
% Ey2 = D1*exp(-1i.*alpha_p_EP.*(x2+a)) + D2*exp(1i.*alpha_p_EP.*(x2+a));
% Ey3 = C1*exp(beta_EP.*(x3)) + C2*exp(-beta_EP.*(x3));
% Ey4 = B1*exp(-1i.*alpha_n_EP.*(x4-a)) + B2*exp(1i.*alpha_n_EP.*(x4-a));
% Ey5 = A*exp(-beta_EP.*(x5-a-b));
% Ey = [Ey1,Ey2,Ey3,Ey4,Ey5];

figure(5)
subplot(2,2,1)
abs_Ey = abs(Ey);   abs_Ey_A = abs_Ey;
plot(x*100,abs_Ey);
xlabel('$x(\mathrm{cm})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$\mathrm{abs}(E_y)$','interpreter','latex','FontName','Times New Roman','FontSize',20)
% ylabel('\rm\fontname{Times New Roman}abs \rm\fontname{Times New Roman}E_y','FontSize',20)

% real part && imag part
figure(5)
subplot(2,2,2)
real_Ey = real(Ey);
plot(x*100,real_Ey);
xlabel('$x(\mathrm{cm})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$\mathrm{Re}(E_y)$','interpreter','latex','FontName','Times New Roman','FontSize',20)

figure(5)
subplot(2,2,3)
imag_Ey = imag(Ey);
plot(x*100,imag_Ey);
xlabel('$x(\mathrm{cm})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$\mathrm{Im}(E_y)$','interpreter','latex','FontName','Times New Roman','FontSize',20)

% phase part
figure(5)
subplot(2,2,4)
phase_Ey = angle(Ey)/pi;
plot(x*100,phase_Ey);
xlabel('$x(\mathrm{cm})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$\mathrm{phase}(\pi)$','interpreter','latex','FontName','Times New Roman','FontSize',20)

%% Poynting vector at A

Hx1 = -kk(nx_EP,ny)/omega/mu_0/mu_r .*(Ey1);
Hx2 = -kk(nx_EP,ny)/omega/mu_0/mu_p .*(Ey2);
Hx3 = -kk(nx_EP,ny)/omega/mu_0/mu_r .*(Ey3);
Hx4 = -kk(nx_EP,ny)/omega/mu_0/mu_n .*(Ey4);
Hx5 = -kk(nx_EP,ny)/omega/mu_0/mu_r .*(Ey5);
Hx = [Hx1,Hx2,Hx3,Hx4,Hx5];

figure('numbertitle','off','name','Poynting vector'); 

sz = - Ey.*conj(Hx);
% abs_Sz  = abs(Sz);
real_sz = real(sz);    sz_A = real_sz;
imag_sz = imag(sz);

hold on
% plot(x,abs(Sz));
plot(x*100,real_sz);
plot(x*100,imag_sz);
% legend('$\mathrm{abs}(S_z)$','$\mathrm{Re}(S_z)$','$\mathrm{Im}(S_z)$','interpreter','latex','FontName','Times New Roman','FontSize',20)
legend('$\mathrm{Re}(S_z)$','$\mathrm{Im}(S_z)$','interpreter','latex','FontName','Times New Roman','FontSize',20,'Location','east')
title('$S_z=-E_y\cdot H_x^{*}$','interpreter','latex','FontName','Times New Roman','FontSize',20)
xlabel('$x(\mathrm{cm})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$S_z$','interpreter','latex','FontName','Times New Roman','FontSize',20)
hold off

%% distribution of Ey at B
[peaks,locs] = findpeaks(-result1);
num_EP=length(locs); % number of peaks
if num_EP>=3
    nx1_EP = locs(1); % disatance_a
    nx2_EP = locs(2);
elseif num_EP==2
    nx1_EP = locs(1);
end
%%%%%%%% select peak
nx_EP = nx2_EP;

% PIM or NIM
beta_EP=beta(nx_EP, ny);
alpha_p_EP=alpha_p(nx_EP,ny);
alpha_n_EP=alpha_n(nx_EP,ny);
a = aa(nx_EP, ny);

x1 = linspace(-4*a-b, -a-b, 50);
x2 = linspace(-a-b, -a, 50);
x3 = linspace(-a, a, 50);
x4 = linspace(a, a+b, 50);
x5 = linspace(a+b, 4*a+b, 50);

A  = 1;
B1 = A*(1i*alpha_n_EP/mu_n + beta_EP/mu_r ) / (2*1i * alpha_n_EP/mu_n * exp(-1i*alpha_n_EP*b) );
B2 = A*(1i*alpha_n_EP/mu_n - beta_EP/mu_r ) / (2*1i * alpha_n_EP/mu_n * exp(+1i*alpha_n_EP*b) );
C1 = ( ( beta_EP/mu_r - 1i*alpha_n_EP/mu_n )*B1 + ( beta_EP/mu_r + 1i*alpha_n_EP/mu_n )*B2 ) / ( 2*beta_EP/mu_r * exp(+beta_EP * a) );
C2 = ( ( beta_EP/mu_r + 1i*alpha_n_EP/mu_n )*B1 + ( beta_EP/mu_r - 1i*alpha_n_EP/mu_n )*B2 ) / ( 2*beta_EP/mu_r * exp(-beta_EP * a) );
D1 = ( ( 1i*alpha_p_EP/mu_p*exp(-beta_EP*a) - beta_EP/mu_r*exp(-beta_EP*a) )*C1 ...
     + ( 1i*alpha_p_EP/mu_p*exp(+beta_EP*a) + beta_EP/mu_r*exp(+beta_EP*a) )*C2 )...
     / ( 2*1i*alpha_p_EP/mu_p );
D2 = ( ( 1i*alpha_p_EP/mu_p*exp(-beta_EP*a) + beta_EP/mu_r*exp(-beta_EP*a) )*C1 ...
     + ( 1i*alpha_p_EP/mu_p*exp(+beta_EP*a) - beta_EP/mu_r*exp(+beta_EP*a) )*C2 )...
     / ( 2*1i*alpha_p_EP/mu_p );
F  = D1*exp(1i*alpha_p_EP*b) + D2 * exp(-1i*alpha_p_EP*b);

Ey1 = F*exp(beta_EP.*(x1+a+b));
Ey2 = D1*exp(-1i.*alpha_p_EP.*(x2+a)) + D2*exp(1i.*alpha_p_EP.*(x2+a));
Ey3 = C1*exp(beta_EP.*(x3)) + C2*exp(-beta_EP.*(x3));
Ey4 = B1*exp(-1i.*alpha_n_EP.*(x4-a)) + B2*exp(1i.*alpha_n_EP.*(x4-a));
Ey5 = A*exp(-beta_EP.*(x5-a-b));

x = [x1,x2,x3,x4,x5];
Ey = [Ey1,Ey2,Ey3,Ey4,Ey5];

maxnum = max(abs(Ey));

figure(5)
subplot(2,2,1)
abs_Ey = abs(Ey);   abs_Ey_B = abs_Ey;
plot(x*100,abs_Ey);
xlabel('$x(\mathrm{cm})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$\mathrm{abs}(E_y)$','interpreter','latex','FontName','Times New Roman','FontSize',20)
% ylabel('\rm\fontname{Times New Roman}abs \rm\fontname{Times New Roman}E_y','FontSize',20)

% real part && imag part
figure(5)
subplot(2,2,2)
real_Ey = real(Ey);
plot(x*100,real_Ey);
xlabel('$x(\mathrm{cm})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$\mathrm{Re}(E_y)$','interpreter','latex','FontName','Times New Roman','FontSize',20)

figure(5)
subplot(2,2,3)
imag_Ey = imag(Ey);
plot(x*100,imag_Ey);
xlabel('$x(\mathrm{cm})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$\mathrm{Im}(E_y)$','interpreter','latex','FontName','Times New Roman','FontSize',20)

% phase part
figure(5)
subplot(2,2,4)
phase_Ey = angle(Ey)/pi;
plot(x*100,phase_Ey);
xlabel('$x(\mathrm{cm})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$\mathrm{phase}(\pi)$','interpreter','latex','FontName','Times New Roman','FontSize',20)

%% Poynting vector at B

Hx1 = -kk(nx_EP,ny)/omega/mu_0/mu_r .*(Ey1);
Hx2 = -kk(nx_EP,ny)/omega/mu_0/mu_p .*(Ey2);
Hx3 = -kk(nx_EP,ny)/omega/mu_0/mu_r .*(Ey3);
Hx4 = -kk(nx_EP,ny)/omega/mu_0/mu_n .*(Ey4);
Hx5 = -kk(nx_EP,ny)/omega/mu_0/mu_r .*(Ey5);
Hx = [Hx1,Hx2,Hx3,Hx4,Hx5];

figure('numbertitle','off','name','Poynting vector'); 

sz = - Ey.*conj(Hx);
% abs_Sz  = abs(Sz);
real_sz = real(sz);    sz_B = real_sz;
imag_sz = imag(sz);

hold on
% plot(x,abs(Sz));
plot(x*100,real_sz);
plot(x*100,imag_sz);
% legend('$\mathrm{abs}(S_z)$','$\mathrm{Re}(S_z)$','$\mathrm{Im}(S_z)$','interpreter','latex','FontName','Times New Roman','FontSize',20)
legend('$\mathrm{Re}(S_z)$','$\mathrm{Im}(S_z)$','interpreter','latex','FontName','Times New Roman','FontSize',20,'Location','east')
title('$S_z=E_y\cdot H_x^{*}$','interpreter','latex','FontName','Times New Roman','FontSize',20)
xlabel('$x(\mathrm{cm})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$S_z$','interpreter','latex','FontName','Times New Roman','FontSize',20)
hold off

%% plot Ey and sz
figure('numbertitle','off','name','A'); 
hold on
abs_Ey_A_norm = abs_Ey_A/max(abs_Ey_A);
sz_A_norm = sz_A/max(sz_A);
plot(x*100,abs_Ey_A_norm);
plot(x*100,sz_A_norm);
legend('$\mathrm{abs}(E_y)$','$S_z$','interpreter','latex','FontName','Times New Roman','FontSize',20,'Location','northeast')
title('$E_y\ \mathrm{and}\ S_z\ \mathrm{at}\ A$','interpreter','latex','FontName','Times New Roman','FontSize',20)
xlabel('$x(\mathrm{cm})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
% ylabel('$S_z$','interpreter','latex','FontName','Times New Roman','FontSize',20)
hold off

figure('numbertitle','off','name','B');
hold on
abs_Ey_B_norm = abs_Ey_B/max(abs_Ey_B);
sz_B_norm = sz_B/max(sz_B);
plot(x*100,abs_Ey_B_norm);
plot(x*100,sz_B_norm);
% legend('$\mathrm{abs}(S_z)$','$\mathrm{Re}(S_z)$','$\mathrm{Im}(S_z)$','interpreter','latex','FontName','Times New Roman','FontSize',20)
legend('$\mathrm{abs}(E_y)$','$S_z$','interpreter','latex','FontName','Times New Roman','FontSize',20,'Location','southwest')
title('$E_y\ \mathrm{and}\ S_z\ \mathrm{at}\ B$','interpreter','latex','FontName','Times New Roman','FontSize',20)
xlabel('$x(\mathrm{cm})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
% ylabel('$S_z$','interpreter','latex','FontName','Times New Roman','FontSize',20)
hold off


%% distribution of Hz
Hz1 = beta_EP*F/mu_r*exp(beta_EP.*(x1+a+b));
Hz2 = - 1i*alpha_p_EP*D1/mu_p*exp(-1i.*alpha_p_EP.*(x2+a)) ...
      + 1i*alpha_p_EP*D2/mu_p*exp(+1i.*alpha_p_EP.*(x2+a));
Hz3 = beta_EP*C1/mu_r*exp(beta_EP.*x3) - beta_EP*C2/mu_r*exp(-beta_EP.*x3);
Hz4 = - 1i*alpha_n_EP*B1/mu_n*exp(-1i.*alpha_n_EP.*(x4-a)) ...
      + 1i*alpha_n_EP*B2/mu_n*exp(+1i.*alpha_n_EP.*(x4-a));
Hz5 = -beta_EP*A/mu_r*exp(-beta_EP.*(x5-a-b));
Hz = [Hz1,Hz2,Hz3,Hz4,Hz5];
Hz = Hz/1i*omega/mu_0;
real_Hz = real(Hz);
imag_Hz = imag(Hz);
abs_Hz = abs(Hz);
figure(2222)
hold on
% plot(x,real_Hz);
% plot(x,imag_Hz);
plot(x*100,abs_Hz);
hold off

%% group v

pp_epsilon_p = epsilon_p;
pp_mu_p = mu_p;
pp_epsilon_n = 1 + omega_p^2/omega^2;
pp_mu_n = 1 - F * (omega^4 - 3*omega_0^2*omega^2)/(omega^2 - omega_0^2)^2;

w1 = epsilon_0 * Ey1.*conj(Ey1) + mu_0 * ( Hx1.*conj(Hx1) + Hz1.*conj(Hz1) );
w2 = epsilon_0*pp_epsilon_p * Ey2.*conj(Ey2) + mu_0*pp_mu_p * ( Hx2.*conj(Hx2) + Hz2.*conj(Hz2) );
w3 = epsilon_0 * Ey3.*conj(Ey3) + mu_0 * ( Hx3.*conj(Hx3) + Hz3.*conj(Hz3) );
w4 = epsilon_0*pp_epsilon_n * Ey4.*conj(Ey4) + mu_0*pp_mu_n * ( Hx4.*conj(Hx4) + Hz4.*conj(Hz4) );
w5 = epsilon_0 * Ey5.*conj(Ey5) + mu_0 * ( Hx5.*conj(Hx5) + Hz5.*conj(Hz5) );
w = [w1,w2,w3,w4,w5];
W = 1/4*trapz(x,w);
sz = - Ey.*conj(Hx);
Sz = trapz(x,sz);

vg = Sz/W;
