%%%%% MATLAB2021a
clear; close all;

%%%%%  mu_0 = epsilon_0 = c = 1
mu_0 = 1; epsilon_0 = 1; c = 1;
%%%%% size
d = 1e-3; L = 4*d; % length of waveguides = gap between waveguides
k_PBG = pi/L; % reduced wavevector at Brillouin zone edge
nb = 1.5; ns = 1.5; A = 0.1;
%%%%% scan field
beta_list = linspace(0.8*k_PBG,1*k_PBG,1e4+1); %% ii
omega_mid = pi*c/(nb*L);
omega_list1 = linspace(omega_mid*(1-1.9e-1),omega_mid*(1-1.9e-2),1e2+1); %% jj
omega_list2 = linspace(omega_mid*(1-1.9e-2),omega_mid*(1+1.9e-2),1e2+1);
omega_list3 = linspace(omega_mid*(1+1.9e-2),omega_mid*(1+1.9e-1),1e2+1);
omega_list = [omega_list1(1:end-1),omega_list2(1:end-1),omega_list3];
 
%% broken phase
delta = 1.2; % index EP:0.99635
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
result1 = zeros(length(beta_list),length(omega_list));
for ii = 1:length(beta_list)
    for jj = 1:length(omega_list)
        %%%%reflect and transmission%%%%
        [M_be,~,~,~,~] = M1_ReflAndTran(ns, n1 );
        [M_nd,~,~,~,~] = M1_ReflAndTran(n4, ns );
        [M1,~,~,~,~] = M1_ReflAndTran(n4, n1 );
        [M2] = M2_propagation(k1(jj), d );
        [M3,~,~,~,~] = M1_ReflAndTran(n1, n2 );
        [M4] = M2_propagation(k2(jj), d );
        [M5,~,~,~,~] = M1_ReflAndTran(n2, n3 );
        [M6] = M2_propagation(k3(jj), d );
        [M7,~,~,~,~] = M1_ReflAndTran(n3, n4 );
        [M8] = M2_propagation(k4(jj), d );
        M_exp = [exp(-1i*beta_list(ii)*L),0;0,exp(-1i*beta_list(ii)*L)];% exp()[1,0,;0,1]
        M_all = M8*M7*M6*M5*M4*M3*M2*M1;
        result1(ii,jj) = log( abs(det(M_all-M_exp)) ); % Bloch mode
    end
end

figure();
[beta,omega] = meshgrid(beta_list,omega_list);
% surf(beta/k_PBG,omega/(c/L),result');
pcolor(beta/k_PBG,omega/(pi*c/L),result1');
shading interp;
xlabel('$k\ (k_{PBG})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$\omega\ (\pi c/L)$','interpreter','latex','FontSize',20)

%% find peaks
% figure()
% beta_kPBG_list = result1(: ,1);
% plot(beta_kPBG_list)
%% exact phase
delta = 0.8; % index EP:0.99635
n1 = nb + A*(1-1i*delta);
n2 = nb - A*(1+1i*delta);
n3 = nb - A*(1-1i*delta);
n4 = nb + A*(1+1i*delta);
k1 = n1 .* omega_list./c; % perpendicular polarization []
k2 = n2 .* omega_list./c;
k3 = n3 .* omega_list./c;
k4 = n4 .* omega_list./c;
k_vacuum = omega_list./c;
 
result2 = zeros(length(beta_list),length(omega_list));
for ii = 1:length(beta_list)
    for jj = 1:length(omega_list)
        %%%%reflect and transmission%%%%
        [M_be,~,~,~,~] = M1_ReflAndTran(ns, n1 );
        [M_nd,~,~,~,~] = M1_ReflAndTran(n4, ns );
        [M1,~,~,~,~] = M1_ReflAndTran(n4, n1 );
        [M2] = M2_propagation(k1(jj), d );
        [M3,~,~,~,~] = M1_ReflAndTran(n1, n2 );
        [M4] = M2_propagation(k2(jj), d );
        [M5,~,~,~,~] = M1_ReflAndTran(n2, n3 );
        [M6] = M2_propagation(k3(jj), d );
        [M7,~,~,~,~] = M1_ReflAndTran(n3, n4 );
        [M8] = M2_propagation(k4(jj), d );
        M_exp = [exp(-1i*beta_list(ii)*L),0;0,exp(-1i*beta_list(ii)*L)];% exp()[1,0,;0,1]
        M_all = M8*M7*M6*M5*M4*M3*M2*M1;
        result2(ii,jj) = log( abs(det(M_all-M_exp)) ); % Bloch mode
    end
end

figure();
[beta,omega] = meshgrid(beta_list,omega_list);
% surf(beta/k_PBG,omega/(c/L),result2');
pcolor(beta/k_PBG,omega/(pi*c/L),result2');
shading interp;
xlabel('$k\ (k_{PBG})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$\omega\ (\pi c/L)$','interpreter','latex','FontSize',20)
%% find peaks
% figure()
% beta_kPBG_list = result2(: ,50);
% [~, locs_temp] = findpeaks(-beta_kPBG_list);
% plot(beta_kPBG_list)
%% dispersion curve
%%%% broken
betas_omega_broken = zeros(1,length(omega_list));
for jj = 1: length(omega_list)
    spectrum_jj = result1(:,jj);
    [~, locs] = findpeaks(-spectrum_jj);
    if length(locs) == 1
        betas_omega_broken(jj) = beta_list(locs(1));
    elseif isempty(locs) == 0
        betas_omega_broken(jj) = nan;
    end 
end
%%%% exact
betas_omega_exact = zeros(1,length(omega_list));
for jj = 1: length(omega_list)
    spectrum_jj = result2(:,jj);
    [~, locs] = findpeaks(-spectrum_jj);
    if length(locs) == 1
        betas_omega_exact(jj) = beta_list(locs(1));
    elseif isempty(locs) == 1
        betas_omega_exact(jj) = nan;
    end 
end
%%%%%%%%%%%%%
% num = 0;
% for jj = 1: length(beta_list)
%     spectrum_jj = result2(jj,:);
%     [~, locs] = findpeaks(-spectrum_jj);
%     if length(locs) == 2
%         betas_omega_exact(jj, 1) = omega_list(locs(1));
%         betas_omega_exact(jj, 2) = omega_list(locs(2));
%     end
%     if length(locs) == 1 && num == 0
%        betas_omega_exact(jj, 1) = omega_list(locs);
%        betas_omega_exact(jj, 2) = omega_list(locs);
%        num =1;
%     elseif length(locs) == 1 && num == 1
%        betas_omega_exact(jj, 1) = nan;
%        betas_omega_exact(jj, 2) = nan;
%     end
% end
%%%%%%%%%%%%%%%
figure()
hold on
scatter(betas_omega_broken./(pi/L),omega_list./(c*pi/L),10,'r');
locs_nan = find(isnan(betas_omega_exact));
betas_omega_exact1 = betas_omega_exact(1:locs_nan(1)-1);
omega_list_exact1 = omega_list(1:locs_nan(1)-1);
betas_omega_exact2 = betas_omega_exact(locs_nan(end)+1 : end);
omega_list_exact2 = omega_list(locs_nan(end)+1 : end);
scatter(betas_omega_exact1./(pi/L),omega_list_exact1./(c*pi/L),10,'b');
scatter(betas_omega_exact2./(pi/L),omega_list_exact2./(c*pi/L),10,'b');
hold off
xlabel('$k\ (k_{PBG})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$\omega\ (\pi c/L)$','interpreter','latex','FontSize',20)
%% group velocity
d_beta_broken = betas_omega_broken(2:end) - betas_omega_broken(1:end-1);
d_beta_exact1 = betas_omega_exact1(2:end) - betas_omega_exact1(1:end-1);
d_beta_exact2 = betas_omega_exact2(2:end) - betas_omega_exact2(1:end-1);
d_omega = omega_list(2) - omega_list(1);
vg_broken = abs( d_omega ./ d_beta_broken );
vg_excat1 = abs( d_omega ./ d_beta_exact1 );
vg_excat2 = abs( d_omega ./ d_beta_exact2 );
figure()
hold on
scatter(vg_broken./c,omega_list(1:end-1)./(c*pi/L),10,'b');
scatter(vg_excat1./c,omega_list_exact1(1:end-1)./(c*pi/L),10,'r');
scatter(vg_excat2./c,omega_list_exact2(1:end-1)./(c*pi/L),10,'r');
hold off
xlabel('$v_g\ (c)$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$\omega\ (\pi c/L)$','interpreter','latex','FontSize',20)
xlim([0.2,1.2])









