%%%%% MATLAB2021a
clear; close all;

%%%%%  mu_0 = epsilon_0 = c = 1
mu_0 = 1; epsilon_0 = 1; c = 1;
%%%%% size
d = 1e-3; L = 4*d; % PTBG 4 layers
nb = 1.5; ns = 1.5; A = 0.1; delta =1.2; % index
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

 
w_list = linspace(80*d,600*d,20);
delta_peaks = zeros(size(w_list));
for ww = 1:length(w_list)
    %%%%% gaussian
    E0 = 1; omega_PBG = pi*c/L/nb;% reduced wavevector at Brillouin zone edge
    x0 = -550 * d; w0 = w_list(ww); % Gaussian wavepackage % x0 = -550d
    Ex_G = E0 * exp(-(x-x0).^2./w0^2) .* exp( -1i .* (ns.*omega_PBG)./c .* (x) );
    I_0 = abs(Ex_G).^2;
    % frequency grid
    N_kn = 1e4+1;
    kn_mid = (ns.*omega_PBG)./c;
    kn = linspace( (1-2e-1)*kn_mid, (1+2e-1)*kn_mid, N_kn );
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
    %%%
    % figure()
    % plot(x,abs(Ex_G))
    % figure()
    % plot(kn,abs(Fourier))
    % figure()
    % plot(select_kn,abs(select_Fourier))
    % figure()
    % plot(xx,abs(Ex_G_select_ifft))

    %%sum of distributions of diffrent frequencies in fiber
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
    t = 2400*d/c;
    for ii = 1:length(select_list) % for every selected omega
    %     [Ex_kn] = sigle_frequency_distribution(d,x1,NN_total,NN_PT,N,...
    %         k1(ii),k2(ii),k3(ii),k4(ii),ks(ii),n1,n2,n3,n4,ns);
        %   Ex: | M2 | M4 | M6 | M8 
        [M_be,~,~,~,~] = M1_ReflAndTran(ns, n1 ); 
        [M_nd,~,~,~,~] = M1_ReflAndTran(n4, ns );
        [Ms1,~,~,~,~] = M1_ReflAndTran(ns, ns );
        [Ms2] = M2_propagation(ks(ii), d );

        [M1,~,~,~,~] = M1_ReflAndTran(n4, n1 );
        [M2] = M2_propagation(k1(ii), d );
        [M3,~,~,~,~] = M1_ReflAndTran(n1, n2 );
        [M4] = M2_propagation(k2(ii), d );
        [M5,~,~,~,~] = M1_ReflAndTran(n2, n3 );
        [M6] = M2_propagation(k3(ii), d );
        [M7,~,~,~,~] = M1_ReflAndTran(n3, n4 );
        [M8] = M2_propagation(k4(ii), d );

        MM = M8*M7*M6*M5*M4*M3*M2*M1;
        M_all = M_nd*MM^(NN_PT-1)*M8*M7*M6*M5*M4*M3*M2*M_be;
        r12 = - M_all(2,1)/M_all(2,2);
        r21 = M_all(1,2)/M_all(2,2);
        t12 = M_all(1,1) + r12*M_all(1,2);
        t21 = 1/M_all(2,2);

        Ex_kn = zeros(size(xxx)); % distribution of wave with frequency of omega = c*kn/ns
        lication_G = 0;
        %%%%%%%%%%% method 1 %%%%%%%%%%%%%%%%%
        % before PTBG
        Ex_input = [1;r12];
        xxx0 = xxx( 1: (NN_total/2)*((N_xxx-1)/NN_total) );
        Ex_kn( 1: (NN_total/2)*((N_xxx-1)/NN_total) ) = ...
            Ex_input(1).*exp(-1i .* ks(ii) .* xxx0 ) + Ex_input(2).*exp( 1i .* ks(ii) .* xxx0 );     

        % propagation in PTBG
        delta_xxx = xxx(2) - xxx(1);
        xxx1 = (0:delta_xxx:d-delta_xxx);% all grid points in 1d
        Ex_input = [1;r12];
        for nn = NN_total/2 +1 +lication_G : NN_total/2 + NN_PT % 363~500, 501~600, 601~1000
            if  nn == NN_total/2 +1 % 501
                [unit_Ex,Ex_output] = Unit_Ex_forward(Ex_input,M_be,M2,M3,M4,M5,M6,M7,M8,k1(ii),k2(ii),k3(ii),k4(ii),xxx1);
            elseif nn > NN_total/2 +1 && nn < NN_total/2 +NN_PT +1 % 502~600
                [unit_Ex,Ex_output] = Unit_Ex_forward(Ex_input,M1,M2,M3,M4,M5,M6,M7,M8,k1(ii),k2(ii),k3(ii),k4(ii),xxx1);
            end
            Ex_input = Ex_output;
            Ex_kn( (nn-1)*(N_xxx-1)/NN_total +1 : (nn)*(N_xxx-1)/NN_total ) = unit_Ex;
        end
        % propagation out of PTBG
        Ex_input = M_nd*Ex_input;
        Ex_input(2) = 0; % theoretically [t,0] 
        xxx2 = xxx((NN_total/2+NN_PT)*((N_xxx-1)/NN_total)+1 : end) - xxx((NN_total/2+NN_PT)*((N_xxx-1)/NN_total)+1);
        Ex_kn( (NN_total/2+NN_PT)*((N_xxx-1)/NN_total)+1 : end ) = ...
            Ex_input(1).*exp(-1i .* ks(ii) .* xxx2 ) + Ex_input(2).*exp( 1i .* ks(ii) .* xxx2 );
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Ex_G_select = Ex_G_select + (select_Fourier(ii)) .* Ex_kn .* exp(1i*select_omega(ii)*t);% k = ns*omega_PBG/c 
    end
%     figure()
%     plot(xxx(1:end),abs(Ex_G_select).^2/max(abs(Ex_G_select).^2))
    [~,location_peak] = findpeaks(abs(Ex_G_select).^2);
    location_xxx = xxx(location_peak(end));
    delta_peaks(ww) = location_xxx;
end
%% 
Delta_peaks = delta_peaks - x0 - t*c/nb; % spatial shift of different w0
figure()
plot(w_list./d,Delta_peaks./d)
xlabel('w_0 (d)')
ylabel('\Delta (d)')

