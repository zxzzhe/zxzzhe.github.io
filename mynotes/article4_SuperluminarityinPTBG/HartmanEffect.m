%%%%% MATLAB2021a
clear; close all;

%%%%%  mu_0 = epsilon_0 = c = 1
mu_0 = 1; epsilon_0 = 1; c = 1;
%%%%% size
d = 4e-3; % length of waveguides = gap between waveguides
k_PBG = pi/d; % reduced wavevector at Brillouin zone edge
%%%%% scan field
beta_list = linspace(0.5*k_PBG,1*k_PBG,1e2); %% ii
omega_list = linspace(1*c/d,3.5*c/d,1e4); %% jj

nb = 1.5; ns = 1; A = 0.1; delta =1.2; % index: delta=0.8 PT exact phase, 1.2 PT broken phase
n1 = nb + A*(1-1i*delta);
n2 = nb - A*(1+1i*delta);
n3 = nb - A*(1-1i*delta);
n4 = nb + A*(1+1i*delta);

NN_list = [40,80,120]; %number of units
omega_PBG = linspace(1.6*c/d,2.6*c/d,1e3);
k1_PBG = n1 .* omega_PBG./c; % perpendicular polarization
k2_PBG = n2 .* omega_PBG./c;
k3_PBG = n3 .* omega_PBG./c;
k4_PBG = n4 .* omega_PBG./c;

Phi = zeros(3,length(omega_PBG));

for nn = 1:3
    NN = NN_list(nn);
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
        M_exp = [exp(-1i*k_PBG*d),0;0,exp(-1i*k_PBG*d)];
        M_all = M_nd*M^(NN-1)*M8*M7*M6*M5*M4*M3*M2*M_be;
        rr = - M_all(2,1)/M_all(2,2);
        tt = M_all(1,1) + rr*M_all(1,2);
        Phi(nn,ii) = angle(tt);
    end
end


figure()
hold on
plot(omega_PBG/(c/d),Phi(1,:))
plot(omega_PBG/(c/d),Phi(2,:))
plot(omega_PBG/(c/d),Phi(3,:))
title('$Phase\ of\ Transmission$','interpreter','latex','FontName','Times New Roman','FontSize',20)
hold off

figure()
Phi_omega = zeros( 3,length(Phi(1,:))-1 );
Phi_omega(1,:) = (Phi(1,2:end) - Phi(1,1:end-1))./(omega_PBG(2)-omega_PBG(1));
Phi_omega(2,:) = (Phi(2,2:end) - Phi(2,1:end-1))./(omega_PBG(2)-omega_PBG(1));
Phi_omega(3,:) = (Phi(3,2:end) - Phi(3,1:end-1))./(omega_PBG(2)-omega_PBG(1));

for nn = 1:3
    abr = find(Phi_omega(nn,:) >1); % phase aberrant
    for ii = 1:length(abr)  %  phase change
    %     if abr(ii)<length(Phi)
            Phi_omega(nn,abr(ii)) = (Phi(nn,abr(ii)+1)-Phi(nn,abr(ii))-2*pi)./(omega_PBG(2)-omega_PBG(1));
    %     end
    end
end

delta_t = - Phi_omega;
hold on
plot(omega_PBG(1:end-1)/(c/d), delta_t(1,:))
plot(omega_PBG(1:end-1)/(c/d), delta_t(2,:))
plot(omega_PBG(1:end-1)/(c/d), delta_t(3,:))
title('$\Delta_t$','interpreter','latex','FontName','Times New Roman','FontSize',20)
hold off