clear;%close all;
k=linspace(0, 9/1000,1e3);
n_r=1.5;n_i=0.03;c=1;
omega1=c*linspace(0, 9/1000, 1e5);
[kk,omega]=meshgrid(k,omega1);
result=zeros( length(omega1),length(k) );

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

%%
for n_EP = 880:889
    figure(n_EP)
    angular_frequency = omega(:,n_EP);
    result1 = result(:,n_EP);
    plot(angular_frequency,result1);
end

%% velocity
%ni=0(hermite situation)  k=810------>omega=609 or 631
loc_be = 2;
loc_ed = 881; %811;   %position of EP

% 记录色散曲线和poynting vector得到的群速度
omega2k = zeros(2,loc_ed);
velocity_Sz = zeros(2,loc_ed);
velocity_Sz_conj = zeros(2,loc_ed);

omega2k(:,1:loc_be-1) = NaN;
velocity_Sz(:,1:loc_be-1) = NaN;
velocity_Sz_conj(:,1:loc_be-1) = NaN;

for ny_EP = loc_be : loc_ed
    angular_frequency = omega(:,ny_EP);
    result1 = result(:,ny_EP);
    [peaks,locs] = findpeaks(-result1,'minpeakheight',0);
    num_EP = length(locs); % number of peaks
    
    nx_EPs = [0,0];
    if num_EP == 0
        nx_EPs(1) = find(result1 == min(result1));
        nx_EPs(2) = NaN;
        omega2k(1,ny_EP) = omega(nx_EPs(1),ny_EP);
        omega2k(2,ny_EP) = NaN;
    elseif num_EP == 1
        nx_EPs(1) = locs(1);
        nx_EPs(2) = NaN;
        omega2k(1,ny_EP) = omega(nx_EPs(1),ny_EP);
        omega2k(2,ny_EP) = NaN;
    elseif num_EP >= 2 && ny_EP <= loc_ed
        nx_EPs(1) = locs(1);
        nx_EPs(2) = locs(2);
        omega2k(1,ny_EP) = omega(nx_EPs(1),ny_EP);
        omega2k(2,ny_EP) = omega(nx_EPs(2),ny_EP);
%     elseif num_EP >= 2 && ny_EP >= loc_ed
%         nx_EPs(1) = locs(1);
%         nx_EPs(2) = locs(1);
%         omega2k(1,ny_EP) = omega(nx_EPs(1),ny_EP);
%         omega2k(2,ny_EP) = omega(nx_EPs(2),ny_EP);
    end
    
    % select peak
    for jj = 1:2
        nx_EP = nx_EPs(jj); 
        if isnan(nx_EP)
            velocity_Sz(jj,ny_EP) = NaN;
            velocity_Sz_conj(jj,ny_EP) = NaN;
            break;
        end

        beta_EP = beta(nx_EP,ny_EP);
        alpha_p_EP = alpha_p(nx_EP,ny_EP);
        alpha_n_EP = alpha_n(nx_EP,ny_EP);

        % distribution of Ey
        x1 = linspace(-100*a-b, -a-b, 99*a);
        x2 = linspace(-a-b, -a, b);
        x3 = linspace(-a, a, 2*a);
        x4 = linspace(a, a+b, b);
        x5 = linspace(a+b, 100*a+b, 99*a);

        A  = 1;
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

        % Poynting vector
        Hx1 = -kk(nx_EP,ny_EP)/omega(nx_EP,ny_EP)*(Ey1);
        Hx2 = -kk(nx_EP,ny_EP)/omega(nx_EP,ny_EP)*(Ey2);
        Hx3 = -kk(nx_EP,ny_EP)/omega(nx_EP,ny_EP)*(Ey3);
        Hx4 = -kk(nx_EP,ny_EP)/omega(nx_EP,ny_EP)*(Ey4);
        Hx5 = -kk(nx_EP,ny_EP)/omega(nx_EP,ny_EP)*(Ey5);
        Hx = [Hx1,Hx2,Hx3,Hx4,Hx5];
        sz = - Ey.*Hx;
        Hx_conj = - kk(nx_EP,ny_EP)/omega(nx_EP,ny_EP)*conj(Ey);
        sz_conj = -1/2 * Ey.*Hx_conj;

        % energy density
        w1_conj = Ey1.*conj(Ey1);
%         w2_conj = n_r*Ey2.*conj(Ey2);
        w2_conj = epsilon_n*Ey2.*conj(Ey2);
        w3_conj = Ey3.*conj(Ey3);
%         w4_conj = n_r*Ey4.*conj(Ey4);
        w4_conj = epsilon_p*Ey4.*conj(Ey4);
        w5_conj = Ey5.*conj(Ey5);
        
        w1 = Ey1.*(Ey1);
%         w2 = n_r*Ey2.*(Ey2) ;
        w2 = epsilon_n*Ey2.*(Ey2);
        w3 = Ey3.*(Ey3);
%         w4 = n_r*Ey4.*(Ey4) ;
        w4 = epsilon_p*Ey4.*(Ey4);
        w5 = Ey5.*(Ey5) ;
        
        w_conj = [w1_conj,w2_conj,w3_conj,w4_conj,w5_conj];
        W_conj = 1/2*trapz(x,w_conj);
        w = [w1,w2,w3,w4,w5];
        W = trapz(x,w);
        
        % poynting vector
        Sz = trapz(x,sz);
        vg = (Sz)/(W);
        Sz_conj = trapz(x,sz_conj);
        vg_conj = Sz_conj/W_conj;
        
        % group velocity
        velocity_Sz(jj,ny_EP) = vg;
        velocity_Sz_conj(jj,ny_EP) = vg_conj;
    end
end

%% omega2k
% figure(123)
% subplot(2,2,1)

figure('numbertitle','off','name','omega2k'); 
hold on
plot(k(1:loc_ed),omega2k(1,:));
plot(k(1:loc_ed),omega2k(2,:));
% plot(k(1:4:loc_ed),omega2k(1,1:4:end));
% plot(k(1:4:loc_ed),omega2k(2,1:4:end));
hold off
xlabel('$k$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('\rm\fontname{Times New Roman} \rm\fontname{Times New Roman}\omega','FontSize',20)

%% group velocity of Sz
figure('numbertitle','off','name','group velocity of Sz'); 
% figure(1234)
% subplot(2,2,1)
hold on
abs_velocity_Sz = (velocity_Sz);
delete_num = 20;
plot(k(1:loc_ed-delete_num ),abs_velocity_Sz(1,1:loc_ed-delete_num ));
plot(k(1:loc_ed-delete_num ),abs_velocity_Sz(2,1:loc_ed-delete_num ));

hold off
xlabel('$k$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$v_g$','interpreter','latex','FontName','Times New Roman','FontSize',20)
title('$group \ velocity\ of\ S_z$','interpreter','latex','FontName','Times New Roman','FontSize',20)
text(0.004,1,'$\frac{\int{\vec{E}\times \vec{H} dr}}{\int{\vec{E}\cdot\vec{D} dr}} $','interpreter','latex','fontsize',20,'ho','c')
ylim([0 2])
%% group velocity of Sz(conj)
% figure('numbertitle','off','name','group velocity of Sz(conj)'); 
figure(1234)
subplot(2,2,2)
hold on
abs_velocity_Sz_conj = abs(velocity_Sz_conj);

plot(k(1:loc_ed),abs_velocity_Sz_conj(1,:));
plot(k(1:loc_ed),abs_velocity_Sz_conj(2,:));
hold off
xlabel('$k$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$v_g$','interpreter','latex','FontName','Times New Roman','FontSize',20)
text(0.004,0.5,'$\frac{\int{\vec{E}\times \vec{H} dr}}{\int{\vec{E}\cdot\vec{D}\ dr}} $','interpreter','latex','fontsize',20,'ho','c')
ylim([0 1])
title('$group\ velocity\ of\ S_z(conj)$','interpreter','latex','FontName','Times New Roman','FontSize',20)

%% group velocity from dispersion curves
% figure('numbertitle','off','name','group velocity from dispersion curves');
figure(1234)
subplot(2,2,3)
hold on
diff_omega1 = diff( omega2k(1,1:1:end) )./ diff( k(1:1:loc_ed) );
diff_omega2 = diff( omega2k(2,1:1:end) )./ diff( k(1:1:loc_ed) );
plot(k(2:1:loc_ed),diff_omega1);
plot(k(2:1:loc_ed),diff_omega2);

hold off
xlabel('$k$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$v_g$','interpreter','latex','FontName','Times New Roman','FontSize',20)
title('$group\ velocity\ from\ dispersion\ curves$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylim([0 2])

%% phase velocity from dispersion curves

figure(1234)
subplot(2,2,4)
% figure('numbertitle','off','name','phase velocity from dispersion curves');
hold on

vp  = zeros(2,length(omega2k(1,1:4:end)) );
vp(1,:) = omega2k(1,1:4:end)./k(1:4:loc_ed);
vp(2,:) = omega2k(2,1:4:end)./k(1:4:loc_ed);

plot(k(1:4:loc_ed), vp(1,:));
plot(k(1:4:loc_ed), vp(2,:));
hold off
xlabel('$k$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$v_p$','interpreter','latex','FontName','Times New Roman','FontSize',20)
title('$phase\ velocity\ from\ dispersion\ curves$','interpreter','latex','FontName','Times New Roman','FontSize',20)

%% picture together
figure(6666)
hold on
abs_velocity_Sz = abs(velocity_Sz);
delete_num = 20;
gap = 1;
c1 = plot(k(1:gap:loc_ed-delete_num),abs_velocity_Sz(1,1:gap:loc_ed-delete_num), 'ro', 'linewidth',0.51);
c2 = plot(k(1:gap:loc_ed-delete_num),abs_velocity_Sz(2,1:gap:loc_ed-delete_num), 'ro', 'linewidth',0.51);

diff_omega1 = diff( omega2k(1,1:1:loc_ed) )./ diff( k(1:1:loc_ed) );
diff_omega2 = diff( omega2k(2,1:1:loc_ed) )./ diff( k(1:1:loc_ed) );
% diff_omega1 = [diff_omega1, diff_omega2(loc_ed-1)];
c3 = plot(k(2:1:loc_ed),diff_omega1, 'b', 'linewidth',1.5);
c4 = plot(k(2:1:loc_ed),diff_omega2, 'b', 'linewidth',1.5);

xlabel('$k$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$v_g$','interpreter','latex','FontName','Times New Roman','FontSize',20)
title('$group\ velocity$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylim([0 1.5])
legend([c1,c3],'$v_g^{2\omega}$','$v_g^{D}$','interpreter','latex','FontName','Times New Roman','FontSize',20,'Location','SouthWest')

hold off
