%%%%% MATLAB2021a
clear; close all;

%%%%%  mu_0 = epsilon_0 = c = 1
mu_0 = 1; epsilon_0 = 1; c = 1;
%%%%% size
d = 1e-3; L = 4*d; % PTBG 4 layers
nb = 1.5; ns = 1.5; A = 0.1; delta = -0.9963; % delta = \pm 0.9963 is EP
n1 = nb + A*(1+1i*delta);
n2 = nb - A*(1-1i*delta);
n3 = nb - A*(1+1i*delta);
n4 = nb + A*(1-1i*delta);
% change sequences of layers
n_temp = n1;
n1 = n2;
n2 = n3;
n3 = n4;
n4 = n_temp;

% space grid
N_x = 1e4+1; % number of sample
x = linspace(-2.5,1.5,N_x); % period for fft
% gaussian
E0 = 1; omega_PBG = pi*c/L/nb;% reduced wavevector at Brillouin zone edge
x0 = -550 * d; w0 = 240*d; % Gaussian wavepackage % x0 = -550d
Ex_G = E0 * exp(-(x-x0).^2./w0^2) .* exp( -1i .* (ns.*omega_PBG)./c .* (x) );
Ex_G(x>x0) = 0;
I_0 = abs(Ex_G).^2;
% frequency grid
N_kn = 1e4+1;
kn_mid = (ns.*omega_PBG)./c;
kn = linspace( (1-5e-1)*kn_mid, (1+5e-1)*kn_mid, N_kn );
% FFT weight of frequency from -(N/2）*(2pi/period)~ +（N/2）*(2pi/period)
Fourier = zeros(1,length(kn));  % weight
for ii = 1:length(x)
    Fourier = Fourier + Ex_G(ii) .* exp(1i .* kn .* x(ii)); % 空间位置为权重，各次谐波进行叠加
end
Fourier = Fourier./N_x; % frequency sprectrum
% select some frequencies to get main part of Gaussian wave
select_list = find(abs(Fourier)>1e-9);
select_Fourier = Fourier(select_list);
select_kn = kn(select_list);
%%%%%% ifft example
N_xx = 2e3+1;
xx = linspace(-2,2,N_xx); % period for caculating distribution
Ex_G_select_ifft = zeros(1,length(xx));
for jj = 1:length(select_kn)
    Ex_G_select_ifft = Ex_G_select_ifft + select_Fourier(jj).*exp(-1i*select_kn(jj).*xx);
end
%%
figure()
plot(x,abs(Ex_G).^2)
figure()
plot(kn,abs(Fourier))
figure()
plot(select_kn,abs(select_Fourier))
figure()
plot(xx,abs(Ex_G_select_ifft).^2)

%% sum of distributions of diffrent frequencies in fiber
%%%% space grid
NN_total = 1e3;% total units L
NN_total_d = 4 * NN_total; %total uniits d
NN_PT = 100; % total PTBG units
N_xxx = 2*NN_total_d +1;
xxx = linspace(-NN_total_d*d/2,NN_total_d*d/2,N_xxx);
Ex_G_select = zeros(size(xxx));
%%%%%  frequency selected
select_omega = c.*select_kn./ns; % angle frequency list
k1 = n1/ns .* select_kn; % layer1 k1 list
k2 = n2/ns .* select_kn; % layer2 k2 list
k3 = n3/ns .* select_kn; % layer3 k3 list
k4 = n4/ns .* select_kn; % layer4 k4 list
ks = select_kn; % surrounding space ks list
% ifft, caculate Ex distribution using selected frequencies
t = 0*d/c;
for ii = 1:length(select_list) % for every selected omega
    [Ex_kn] = sigle_frequency_distribution(ii,d,NN_total,NN_PT,xxx,N_xxx,...
            k1,k2,k3,k4,ks,n1,n2,n3,n4,ns);
    Ex_G_select = Ex_G_select + (select_Fourier(ii)) .* Ex_kn .* exp(1i*select_omega(ii)*t);% k = ns*omega_PBG/c 
end
figure()
plot(xxx(1:end),abs(Ex_G_select).^2/max(abs(Ex_G_select).^2))
%%
t1 = 2200*d/c;
% ifft, caculate Ex distribution using selected frequencies
Ex_G_select_t1 = zeros(size(xxx));
for ii = 1:length(select_list)
    [Ex_kn] = sigle_frequency_distribution(ii,d,NN_total,NN_PT,xxx,N_xxx,...
        k1,k2,k3,k4,ks,n1,n2,n3,n4,ns);
    Ex_G_select_t1 = Ex_G_select_t1 + select_Fourier(ii) .* Ex_kn .* exp(1i*select_omega(ii)*t1);
end
figure()
plot(xxx(1:end),abs(Ex_G_select_t1).^2/max(abs(Ex_G_select).^2));
title('$ t = 2200 \tau_0,\ \Delta = - 0.9963 $','interpreter','latex')
ylim([0,2])

