%%%%% MATLAB2021a
clear; close all;

%%%%%  mu_0 = epsilon_0 = c = 1
mu_0 = 1; epsilon_0 = 1; c = 1;
%%%%% size
d = 4e-3; % length of waveguides = gap between waveguides
k_PBG = pi/d; % reduced wavevector at Brillouin zone edge
%%%%% scan field
beta_list = linspace(0.5*k_PBG,1*k_PBG,1e2); %% ii
omega_list = linspace(1*c/d,3.5*c/d,2e2); %% jj

nb = 1.5; ns = 1; A = 0.1; delta =0.99635; % index
n1 = nb + A*(1-1i*delta);
n2 = nb - A*(1+1i*delta);
n3 = nb - A*(1-1i*delta);
n4 = nb + A*(1+1i*delta);
k1 = n1 .* omega_list./c; % perpendicular polarization []
k2 = n2 .* omega_list./c;
k3 = n3 .* omega_list./c;
k4 = n4 .* omega_list./c;
k_vacuum = omega_list./c;

%%%%%%% forward %%%%%%%%%
NN = 80; % number of units
omega_PBG = linspace(1.6*c/d,2.6*c/d,1e4);
k1_PBG = n1 .* omega_PBG./c; % perpendicular polarization
k2_PBG = n2 .* omega_PBG./c;
k3_PBG = n3 .* omega_PBG./c;
k4_PBG = n4 .* omega_PBG./c;

TT = zeros(1,length(omega_PBG));
RR = zeros(1,length(omega_PBG));
Phi = zeros(1,length(omega_PBG));
for ii = 1:length(omega_PBG)
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
    TT(ii) = abs(tt)^2;
    RR(ii) = abs(rr)^2;
    Phi(ii) = angle(tt);
end
%%%%%%% backward %%%%%%%%%%%%
TT_b = zeros(1,length(omega_PBG));
RR_b = zeros(1,length(omega_PBG));
Phi_b = zeros(1,length(omega_PBG));

for ii = 1:length(omega_PBG)
    [M_be,~,~,~,~] = M1_ReflAndTran(ns,n3 );
    [M_nd,~,~,~,~] = M1_ReflAndTran(n4,ns );
    [M1,~,~,~,~] = M1_ReflAndTran(n4,n3 );
    [M2] = M2_propagation(k3_PBG(ii),d/4);
    [M3,~,~,~,~] = M1_ReflAndTran(n3,n2 );
    [M4] = M2_propagation(k2_PBG(ii),d/4);
    [M5,~,~,~,~] = M1_ReflAndTran(n2,n1 );
    [M6] = M2_propagation(k1_PBG(ii),d/4);
    [M7,~,~,~,~] = M1_ReflAndTran(n1,n4 );
    [M8] = M2_propagation(k4_PBG(ii),d/4);
    M = M8*M7*M6*M5*M4*M3*M2*M1; 
    M_exp = [exp(-1i*k_PBG*d),0;0,exp(-1i*k_PBG*d)];
    M_all = M_nd*M^(NN-1)*M8*M7*M6*M5*M4*M3*M2*M_be;
    rr = - M_all(2,1)/M_all(2,2);
    tt = M_all(1,1) + rr*M_all(1,2);
    TT_b(ii) = abs(tt)^2;
    RR_b(ii) = abs(rr)^2;
    Phi_b(ii) = angle(tt);
end
%%%%%%%% plot %%%%%%%%%%%%%
figure()
hold on
plot(omega_PBG/(c/d),TT)
plot(omega_PBG/(c/d),TT_b)
title('$Transmission$','interpreter','latex','FontName','Times New Roman','FontSize',20)
hold off

figure()
semilogy(omega_PBG/(c/d),RR,omega_PBG/(c/d),RR_b)
title('$Reflection$','interpreter','latex','FontName','Times New Roman','FontSize',20)
