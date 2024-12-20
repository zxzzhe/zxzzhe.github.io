%%%%% MATLAB2021a
clear; close all;

%%%%%  mu_0 = epsilon_0 = c = 1
mu_0 = 1; epsilon_0 = 1; c = 1;
%%%%% size
d = 4e-3; % length of waveguides = gap between waveguides
k_PBG = pi/d; % reduced wavevector at Brillouin zone edge
%%%%% scan field
beta_list = linspace(0.8*k_PBG,1*k_PBG,1e2); %% ii
omega_list = linspace(0.5*pi*c/d,0.85*pi*c/d,1e2); %% jj

nb = 1.5; ns = 1.5; A = 0.1; delta =1.2; % index EP:0.99635
n1 = nb + A*(1-1i*delta);
n2 = nb - A*(1+1i*delta);
n3 = nb - A*(1-1i*delta);
n4 = nb + A*(1+1i*delta);
k1 = n1 .* omega_list./c; % perpendicular polarization []
k2 = n2 .* omega_list./c;
k3 = n3 .* omega_list./c;
k4 = n4 .* omega_list./c;
k_vacuum = omega_list./c;

%%%%%%%%% debug %%%%%%%%%%
% result1 = zeros(1,length(omega_list));
% for jj = 1:length(omega_list)
%     %%%% reflect and transmission %%%%
%     [M1,~,~,~,~] = M1_ReflAndTran(n4, n1 );
%     [M2] = M2_propagation(k1(jj), d/4 );
%     [M3,~,~,~,~] = M1_ReflAndTran(n1, n2 );
%     [M4] = M2_propagation(k2(jj), d/4 );
%     [M5,~,~,~,~] = M1_ReflAndTran(n2, n3 );
%     [M6] = M2_propagation(k3(jj), d/4 );
%     [M7,~,~,~,~] = M1_ReflAndTran(n3, n4 );
%     [M8] = M2_propagation(k4(jj), d/4 );
%     M_exp = [exp(-1i*beta_list(1)*d),0;0,exp(-1i*beta_list(1)*d)];
%     M_all = M8*M7*M6*M5*M4*M3*M2*M1;
%     result1(jj) = log( abs(det(M_all-M_exp)) );
% end
% figure()
% plot(omega_list/(c/d),result1)
%%%%%%%%% debug %%%%%%%%%%%%

%%%% 
result = zeros(length(beta_list),length(omega_list));
for ii = 1:length(beta_list)
    for jj = 1:length(omega_list)
        %%%%reflect and transmission%%%%
        [M_be,~,~,~,~] = M1_ReflAndTran(ns, n2 );
        [M_nd,~,~,~,~] = M1_ReflAndTran(n1, ns );
        [M1,~,~,~,~] = M1_ReflAndTran(n1, n2 );
        [M2] = M2_propagation(k2(jj), d/4 );
        [M3,~,~,~,~] = M1_ReflAndTran(n2, n3 );
        [M4] = M2_propagation(k3(jj), d/4 );
        [M5,~,~,~,~] = M1_ReflAndTran(n3, n4 );
        [M6] = M2_propagation(k4(jj), d/4 );
        [M7,~,~,~,~] = M1_ReflAndTran(n4, n1 );
        [M8] = M2_propagation(k1(jj), d/4 );
%       [M9,~,~,~,~] = M1_ReflAndTran(n1,n2 );
%        M = M8*M7*M6*M5*M4*M3*M2*M1;
        M_exp = [exp(-1i*beta_list(ii)*d),0;0,exp(-1i*beta_list(ii)*d)];
        M_all = M8*M7*M6*M5*M4*M3*M2*M1;
        result(ii,jj) = log( abs(det(M_all-M_exp)) );
    end
end

figure();
[beta,omega] = meshgrid(beta_list,omega_list);
% surf(beta/k_PBG,omega/(c/d),result');
pcolor(beta/k_PBG,omega/pi/(c/d),result');
shading interp;
xlabel('$k\ (k_{PBG})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$\omega\ (\pi c/d)$','interpreter','latex','FontSize',20)

% find delta where is an EP
% figure()
% omega_kPBG_list = result(end,:);
% plot(omega_kPBG_list)

%% Transmission and \delta_t at k_PBG
NN = 100; % number of units
omega_PBG = linspace(0.5*pi*c/d,0.85*pi*c/d,1e3);
k1_PBG = n1 .* omega_PBG./c; % perpendicular polarization
k2_PBG = n2 .* omega_PBG./c;
k3_PBG = n3 .* omega_PBG./c;
k4_PBG = n4 .* omega_PBG./c;

TT = zeros(1,length(omega_PBG));
RR = zeros(1,length(omega_PBG));
Phi = zeros(1,length(omega_PBG));

for ii = 1:length(omega_PBG)
    %%%%reflect and transmission%%%%
    [M_be,~,~,~,~] = M1_ReflAndTran(ns,n2 );
    [M_nd,~,~,~,~] = M1_ReflAndTran(n1,ns );
    [M1,~,~,~,~] = M1_ReflAndTran(n1,n2 );
    [M2] = M2_propagation(k2_PBG(ii),d/4);
    [M3,~,~,~,~] = M1_ReflAndTran(n2,n3 );
    [M4] = M2_propagation(k3_PBG(ii),d/4);
    [M5,~,~,~,~] = M1_ReflAndTran(n3,n4 );
    [M6] = M2_propagation(k4_PBG(ii),d/4);
    [M7,~,~,~,~] = M1_ReflAndTran(n4,n1 );
    [M8] = M2_propagation(k1_PBG(ii),d/4);
    M = M8*M7*M6*M5*M4*M3*M2*M1;
    M_all = M_nd*M^(NN-1)*M8*M7*M6*M5*M4*M3*M2*M_be;
    rr = - M_all(2,1)/M_all(2,2);
    tt = M_all(1,1) + rr*M_all(1,2);
    TT(ii) = abs(tt)^2;
    RR(ii) = abs(rr)^2;
    Phi(ii) = angle(tt);
%     Phi(ii) = angle(tt)/(2*pi);
%     Phi(ii) = acos(real( M(2,2) ))/pi;
end

figure()
subplot(2,2,1)
plot(omega_PBG/(c/d),TT)
title('$Transmission$','interpreter','latex','FontName','Times New Roman','FontSize',20)

% figure()
subplot(2,2,2)
plot(omega_PBG/(c/d),RR)
title('$Reflection$','interpreter','latex','FontName','Times New Roman','FontSize',20)

% figure()
subplot(2,2,3)
plot(omega_PBG/(c/d),Phi)
title('$Phase\ of\ Transmission$','interpreter','latex','FontName','Times New Roman','FontSize',20)

% subplot(2,2,4)
Phi_omega = (Phi(2:end) - Phi(1:end-1))./(omega_PBG(2)-omega_PBG(1));
% figure()
% plot(Phi_omega)
%%%%  aberrant phase changes begin
abr1 = find(Phi_omega > 1); % phase aberrant
for ii = 1:length(abr1)  %  phase change
%     if abr(ii)<length(Phi)
       Phi_omega(abr1(ii)) = (Phi(abr1(ii)+1)-Phi(abr1(ii))-2*pi)./(omega_PBG(2)-omega_PBG(1));
%     end
end
abr2 = find(Phi_omega < -1); % phase aberrant
for ii = 1:length(abr2)  %  phase change
%     if abr(ii)<length(Phi)
       Phi_omega(abr2(ii)) = (Phi(abr2(ii)+1)-Phi(abr2(ii))+2*pi)./(omega_PBG(2)-omega_PBG(1));
%     end
end
%%%%%% aberrant phase changes end

delta_t = - Phi_omega - nb*NN*d/c;
% delta_t = - Phi_omega;
figure()
plot(omega_PBG(1:end-1)/pi/(c/d), delta_t)
title('$\Delta_t$','interpreter','latex','FontName','Times New Roman','FontSize',20)








