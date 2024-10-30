%%%%% MATLAB2021a
clear; %close all;

%%%%%  mu_0 = epsilon_0 = c = 1
mu_0 = 1; epsilon_0 = 1; c = 1;
%%%%% size
d = 4e-3; % length of waveguides = gap between waveguides
nb = 1.5; ns = 1; A = 0.1; delta =0.8; % index
n1 = nb + A*(1-1i*delta);
n2 = nb - A*(1+1i*delta);
n3 = nb - A*(1-1i*delta);
n4 = nb + A*(1+1i*delta);
% n_temp = n3; %%%% inverse sequence
% n3 = n2;
% n2 = n_temp;
% n_temp = n1;
% n1 = n4;
% n4 = n_temp;

k_PBG = pi/d;  omega_PBG = k_PBG*c/nb;% reduced wavevector at Brillouin zone edge
beta = k_PBG; omega = omega_PBG; 
k1 = n1 .* omega./c; % wavevector
k2 = n2 .* omega./c;
k3 = n3 .* omega./c;
k4 = n4 .* omega./c;
ks = ns .* omega./c;

%%%%reflect and transmission%%%%
[M_be,~,~,~,~] = M1_ReflAndTran(ns, n2 );
[M_nd,~,~,~,~] = M1_ReflAndTran(n1, ns );
[Ms_1,~,~,~,~] = M1_ReflAndTran(ns, ns );
[Ms_2] = M2_propagation(ks, d/4 );

[M1,~,~,~,~] = M1_ReflAndTran(n1, n2 );
[M2] = M2_propagation(k2, d/4 );
[M3,~,~,~,~] = M1_ReflAndTran(n2, n3 );
[M4] = M2_propagation(k3, d/4 );
[M5,~,~,~,~] = M1_ReflAndTran(n3, n4 );
[M6] = M2_propagation(k4, d/4 );
[M7,~,~,~,~] = M1_ReflAndTran(n4, n1 );
[M8] = M2_propagation(k1, d/4 );
[M9,~,~,~,~] = M1_ReflAndTran(n1, n2 );

NNx = 150;% total x in layer digits
NN = 80; % number of layers
d_num = 40; % numbers of one period d
x = linspace(0,NNx*d,d_num*NNx+1);
x1 = x( 1 : d_num/4 );% 4 subunits in one period d
Ey = zeros(size(x));
%%%%% backward and input %%%%%%%
Ey_out = [1,0]';
for nn = NN:-1:1
    if nn < NN
        [unit_Ey,Ey_in] = Unit_Ey_backward(Ey_out,M2,M3,M4,M5,M6,M7,M8,M9,k1,k2,k3,k4,x1);
    elseif nn == NN
        [unit_Ey,Ey_in] = Unit_Ey_backward(Ey_out,M2,M3,M4,M5,M6,M7,M8,M_nd,k1,k2,k3,k4,x1);
    end
    Ey( (nn-1)*d_num+1 : nn*d_num ) = unit_Ey;
    Ey_out = Ey_in;
end
Ey_input_vec = inv(M_be)*Ey_out; %% input vector
x_input = linspace(-10*d,0,1e2);
Ey_input = Ey_input_vec(1)*exp(-1i .* ks .* x_input); %%% only incident part
% Ey_input = Ey_input_vec(1)*exp(-1i .* ks .* x_input)...
%          + Ey_input_vec(2)*exp( 1i .* ks .* x_input); %%% incident and reflection

%%%%%% forward %%%%%%
Ey_in = [1,0]';
for nn = NN+1:NNx
    [unit_Ey,Ey_out] = Unit_Ey_backward(Ey_in,Ms_2,Ms_1,Ms_2,Ms_1,Ms_2,Ms_1,Ms_2,Ms_1,ks,ks,ks,ks,x1);
    Ey( (nn-1)*d_num+1 : nn*d_num ) = unit_Ey;
    Ey_in = Ey_out;
end

figure()
x = [x_input(1:end-1), x];
Ey = [Ey_input(1:end-1), Ey];
abs_Ey = abs(Ey);
plot(x(1:end-1)/d,abs_Ey(1:end-1)./max(abs_Ey));







