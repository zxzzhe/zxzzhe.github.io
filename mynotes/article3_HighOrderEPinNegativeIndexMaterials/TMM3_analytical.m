%%%%% MATLAB2021a
clear;%close all;
%%%%%  mu_0 = epsilon_0 = c = 1
mu_0 = 1;
epsilon_0 = 1;
c = 1;

%%%%% size
b = 4e-2; % length of WGs123
d = 6e-2; % gap of WGs23
lambda = 6e-2;
lambda_p = 3e-2;
lambda_0 = 7.5e-2;
omega = 2*pi/lambda*c;
omega_p = 2*pi/lambda_p*c;
omega_0 = 2*pi/lambda_0*c;

%%%%% air
mu_r = 1*mu_0; % mu = mu_0*mu_r (in vacuum)
epsilon_r = 1*epsilon_0; % epsilon = epsilon_0*epsilon_r (in vacuum)
%%%%% positive material
epsilon_p = 4.50412933 * epsilon_0; mu_p = 1 * mu_0; % natural unit
%%%%% negative material
F = 0.56; 
epsilon_n = 1 - (omega_p/omega)^2 ;
mu_n = 1 - F*omega^2/(omega^2-omega_0^2);
epsilon_n = epsilon_n * epsilon_0;
mu_n = mu_n * mu_0;

%%%%%  variables
% list_a = linspace(3e-2, 7e-2, 16e2+1).*1e4;
% list_beta = linspace(125, 140, 1e3+1)./1e4; % wavevector_parallel
list_a = linspace(6e-2, 8e-2, 1e3+1);
list_beta = linspace(132.35, 132.65, 1e3+1); % wavevector_parallel
%%%%% generate grid & get z axis(result)
[aa,beta] = meshgrid(list_a,list_beta);

%%%%% wavevector_vertical
% k_m = sqrt( epsilon_r .* mu_r .* omega.^2./c^2 - beta.^2 ); % complex
% k_p = sqrt( epsilon_p .* mu_p .* omega.^2./c^2 - beta.^2 );
% k_n = sqrt( epsilon_n .* mu_n .* omega.^2./c^2 - beta.^2 );
k_m = sqrt( epsilon_r .* mu_r .* omega.^2 - beta.^2 ); % complex
k_p = sqrt( epsilon_p .* mu_p .* omega.^2 - beta.^2 );
k_n = sqrt( epsilon_n .* mu_n .* omega.^2 - beta.^2 );
%%%%% area divided into 7 parts, record k & mu 
k_1 = k_m; mu_1 = mu_r;
k_2 = k_n; mu_2 = mu_n;
k_3 = k_m; mu_3 = mu_r;
k_4 = k_p; mu_4 = mu_p;
k_5 = k_m; mu_5 = mu_r;
k_6 = k_p; mu_6 = mu_p;
k_7 = k_m; mu_7 = mu_r;

result = (exp(b.*k_n.*1i).*((k_m.*mu_n)./(2.*k_n.*mu_r) - 1./2).*(exp(-aa.*k_m.*1i).*((k_n.*mu_r)./(2.*k_m.*mu_n) - 1./2).*(exp(b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*(exp(-d.*k_m.*1i).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(exp(b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2) + exp(-b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1)) + exp(d.*k_m.*1i).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).*(exp(b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1) + exp(-b.*k_p.*1i).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).*(((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1).^2 + 2./((k_p.*mu_r)./(k_m.*mu_p) - 1)))) + exp(-b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1).*(exp(d.*k_m.*1i).*(((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).^2 + 2./((k_m.*mu_p)./(k_p.*mu_r) - 1)).*(exp(b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1) + exp(-b.*k_p.*1i).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).*(((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1).^2 + 2./((k_p.*mu_r)./(k_m.*mu_p) - 1))) + exp(-d.*k_m.*1i).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).*(exp(b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2) + exp(-b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1)))) + exp(aa.*k_m.*1i).*((k_n.*mu_r)./(2.*k_m.*mu_n) - 1./2).*(2./((k_m.*mu_n)./(k_n.*mu_r) - 1) + 1).*(exp(-b.*k_p.*1i).*(((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1).^2 + 2./((k_p.*mu_r)./(k_m.*mu_p) - 1)).*(exp(d.*k_m.*1i).*(((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).^2 + 2./((k_m.*mu_p)./(k_p.*mu_r) - 1)).*(exp(b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1) + exp(-b.*k_p.*1i).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).*(((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1).^2 + 2./((k_p.*mu_r)./(k_m.*mu_p) - 1))) + exp(-d.*k_m.*1i).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).*(exp(b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2) + exp(-b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1))) + exp(b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1).*(exp(-d.*k_m.*1i).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(exp(b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2) + exp(-b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1)) + exp(d.*k_m.*1i).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).*(exp(b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1) + exp(-b.*k_p.*1i).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).*(((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1).^2 + 2./((k_p.*mu_r)./(k_m.*mu_p) - 1)))))) + exp(-b.*k_n.*1i).*((k_m.*mu_n)./(2.*k_n.*mu_r) - 1./2).*(exp(aa.*k_m.*1i).*(((k_n.*mu_r)./(2.*k_m.*mu_n) - 1./2).*(2./((k_m.*mu_n)./(k_n.*mu_r) - 1) + 1).^2 + 2./((k_m.*mu_n)./(k_n.*mu_r) - 1)).*(exp(-b.*k_p.*1i).*(((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1).^2 + 2./((k_p.*mu_r)./(k_m.*mu_p) - 1)).*(exp(d.*k_m.*1i).*(((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).^2 + 2./((k_m.*mu_p)./(k_p.*mu_r) - 1)).*(exp(b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1) + exp(-b.*k_p.*1i).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).*(((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1).^2 + 2./((k_p.*mu_r)./(k_m.*mu_p) - 1))) + exp(-d.*k_m.*1i).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).*(exp(b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2) + exp(-b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1))) + exp(b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1).*(exp(-d.*k_m.*1i).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(exp(b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2) + exp(-b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1)) + exp(d.*k_m.*1i).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).*(exp(b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1) + exp(-b.*k_p.*1i).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).*(((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1).^2 + 2./((k_p.*mu_r)./(k_m.*mu_p) - 1))))) + exp(-aa.*k_m.*1i).*((k_n.*mu_r)./(2.*k_m.*mu_n) - 1./2).*(2./((k_m.*mu_n)./(k_n.*mu_r) - 1) + 1).*(exp(b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*(exp(-d.*k_m.*1i).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(exp(b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2) + exp(-b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1)) + exp(d.*k_m.*1i).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).*(exp(b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1) + exp(-b.*k_p.*1i).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).*(((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1).^2 + 2./((k_p.*mu_r)./(k_m.*mu_p) - 1)))) + exp(-b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1).*(exp(d.*k_m.*1i).*(((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).^2 + 2./((k_m.*mu_p)./(k_p.*mu_r) - 1)).*(exp(b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1) + exp(-b.*k_p.*1i).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).*(((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1).^2 + 2./((k_p.*mu_r)./(k_m.*mu_p) - 1))) + exp(-d.*k_m.*1i).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).*(exp(b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2) + exp(-b.*k_p.*1i).*((k_m.*mu_p)./(2.*k_p.*mu_r) - 1./2).*((k_p.*mu_r)./(2.*k_m.*mu_p) - 1./2).*(2./((k_m.*mu_p)./(k_p.*mu_r) - 1) + 1).*(2./((k_p.*mu_r)./(k_m.*mu_p) - 1) + 1))))).*(2./((k_n.*mu_r)./(k_m.*mu_n) - 1) + 1));
result = log( abs( result ) );

figure(2)
pcolor(list_a*100,list_beta/100,result);
% surf(beta,omega,result);% view(0,90); %equivalent expression
shading interp;
% colorbar; colormap(jet);caxis([0,1]);
xlabel('$a\ (\rm{cm})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$\beta\ (\rm{cm^{-1}})$','interpreter','latex','FontSize',20)
% set(gca,'LooseInset',[0,0,0,0]);% Cancel the white edge of the picture
%% trial
for n_EP=273:1:278 % 71~72
    figure(n_EP)
    k_list=beta(:,n_EP);
    result1=result(:,n_EP);
    plot(k_list,result1);
end
%% data
num_be = 1; % set the begin number of data a
%%%%% record curves
beta2a = zeros( 3, length(list_a) ); % 1--k1, 2--k2, 3--k3
beta2a(:,1:num_be) = NaN;
locs2a = zeros( 3, length(list_a) ); % 1--k1, 2--k2, 3--k3
locs2a(:,1:num_be) = NaN;
%%%%% record velocity
velocity_sz = zeros( 3, length(list_a) );
velocity_sz(:,1:num_be) = NaN;
%%%%% record alpha
alpha2a = zeros( 3, length(list_a) );
alpha2a(:,1:num_be) = NaN;

length_loc = 0;
for mm = num_be+1:length(list_a)
    result1=result(:,mm);
    [peaks,locs] = findpeaks(-result1);
    
    % judgment of order, 1st.:beta2a(2,:) 2ed.:beta2a(3,:) 3rd.:beta2a(1,:)
    if isempty(locs) % length(locs) == 0
        beta2a(1,mm) = NaN;
        beta2a(2,mm) = NaN;
        beta2a(3,mm) = NaN;
        length_loc = 0;      % length_loc, to record last locs's length
    elseif length(locs) == 1
        beta2a(3,mm) = list_beta( locs(1) );
        if length_loc == 1 || length_loc == 0
            beta2a(1,mm) = NaN;
            beta2a(2,mm) = list_beta( locs(1) );
            beta2a(3,mm) = NaN;
            locs2a(1,mm) = NaN;
            locs2a(2,mm) = locs(1);
            locs2a(3,mm) = NaN;
            length_loc = 1;
        end
    elseif length(locs) == 2
        if length_loc == 1
            beta2a(1,mm) = NaN;
            beta2a(2,mm) = list_beta( locs(1) );
            beta2a(3,mm) = list_beta( locs(2) );
            locs2a(1,mm) = NaN;
            locs2a(2,mm) = locs(1);
            locs2a(3,mm) = locs(2);
            length_loc = 2;
        elseif length_loc == 2
            beta2a(1,mm) = NaN;
            beta2a(2,mm) = list_beta( locs(1) );
            beta2a(3,mm) = list_beta( locs(2) );
            locs2a(1,mm) = NaN;
            locs2a(2,mm) = locs(1);
            locs2a(3,mm) = locs(2);
            length_loc = 2;
        end
    elseif length(locs) == 3
        if length_loc == 1
            beta2a(1,mm) = list_beta( locs(1) );
            beta2a(2,mm) = list_beta( locs(2) );
            beta2a(3,mm) = list_beta( locs(3) );
            locs2a(1,mm) = locs(1);
            locs2a(2,mm) = locs(2);
            locs2a(3,mm) = locs(3);
            length_loc = 3;
        elseif length_loc == 2
            beta2a(1,mm) = list_beta( locs(1) );
            beta2a(2,mm) = list_beta( locs(2) );
            beta2a(3,mm) = list_beta( locs(3) );
            locs2a(1,mm) = locs(1);
            locs2a(2,mm) = locs(2);
            locs2a(3,mm) = locs(3);
            length_loc = 3;
        elseif length_loc == 3
            beta2a(1,mm) = list_beta( locs(1) );
            beta2a(2,mm) = list_beta( locs(2) );
            beta2a(3,mm) = list_beta( locs(3) );
            locs2a(1,mm) = locs(1);
            locs2a(2,mm) = locs(2);
            locs2a(3,mm) = locs(3);
            length_loc = 3;
        end
    end
    %%%%% jj represents two part of the curve
    for jj = 1:size(beta2a,1) % first dim of beta2a, second is length(beta2a)
        if isnan( beta2a(jj,mm) )
            velocity_sz(jj,mm) = NaN;
            alpha2a(jj,mm) = NaN;
        else
            %%%%% select wave vector and list_a
            nx_EP = locs2a(jj,mm);
            ny = mm;
            beta_EP=beta(nx_EP, ny);
            a = aa(nx_EP, ny);
            k_1_EP = k_1(nx_EP, ny);
            k_2_EP = k_2(nx_EP, ny);
            k_3_EP = k_3(nx_EP, ny);
            k_4_EP = k_4(nx_EP, ny);
            k_5_EP = k_5(nx_EP, ny);
            k_6_EP = k_6(nx_EP, ny);
            k_7_EP = k_7(nx_EP, ny);
            
            %%%%% Ey Hx and Hz --> sz w and vg
            %%%%% x linspace
            space_num = 20;
            x1 = linspace(0, 3*b, space_num);
            x2 = linspace(3*b, 4*b, space_num);
            x3 = linspace(4*b, 4*b+a, space_num);
            x4 = linspace(4*b+a, 5*b+a, space_num);
            x5 = linspace(5*b+a, 5*b+a+d, space_num);
            x6 = linspace(5*b+a+d, 6*b+a+d, space_num);
            x7 = linspace(6*b+a+d, 9*b+a+d, space_num);
            %%%%% reflected & transmitted
            [M1 ,~,~,~,~] = M1_ReflAndTran(-k_1_EP,k_2_EP,mu_1,mu_2);
            [M2 ] = M2_propagation(k_2_EP,b);
            [M3 ,~,~,~,~] = M1_ReflAndTran(k_2_EP,-k_3_EP,mu_2,mu_3);
            [M4 ] = M2_propagation(-k_3_EP,a);
            [M5 ,~,~,~,~] = M1_ReflAndTran(-k_3_EP,k_4_EP,mu_3,mu_4);
            [M6 ] = M2_propagation(k_4_EP,b);
            [M7 ,~,~,~,~] = M1_ReflAndTran(k_4_EP,-k_5_EP,mu_4,mu_5);
            [M8 ] = M2_propagation(-k_5_EP,d);
            [M9 ,~,~,~,~] = M1_ReflAndTran(-k_5_EP,k_6_EP,mu_5,mu_6);
            [M10] = M2_propagation(k_6_EP,b);
            [M11,~,~,~,~] = M1_ReflAndTran(k_6_EP,-k_7_EP,mu_6,mu_7);
            % M = M11*M10*M9*M8*M7*M6*M5*M4*M3*M2*M1;
            
            A = 1;
            Ey1 = A*exp(1i .* -k_1_EP .* (x1-x1(end)) );
            Ey2_vec = M1 *[0, A].';
            Ey3_vec = M3 *M2 *Ey2_vec;
            Ey4_vec = M5 *M4 *Ey3_vec;
            Ey5_vec = M7 *M6 *Ey4_vec;
            Ey6_vec = M9 *M8 *Ey5_vec;
            Ey7_vec = M11*M10*Ey6_vec;
            %%%%% Ey
            Ey2 = Ey2_vec(1).* exp(-1i.*  k_2_EP .* (x2-x2(1))  ) ...
                + Ey2_vec(2).* exp( 1i.*  k_2_EP .* (x2-x2(1))  );
            Ey3 = Ey3_vec(1).* exp(-1i.* -k_3_EP .* (x3-x3(1))  ) ...
                + Ey3_vec(2).* exp( 1i.* -k_3_EP .* (x3-x3(1))  );
            Ey4 = Ey4_vec(1).* exp(-1i.*  k_4_EP .* (x4-x4(1))  ) ...
                + Ey4_vec(2).* exp( 1i.*  k_4_EP .* (x4-x4(1))  );
            Ey5 = Ey5_vec(1).* exp(-1i.* -k_5_EP .* (x5-x5(1))  ) ...
                + Ey5_vec(2).* exp( 1i.* -k_5_EP .* (x5-x5(1))  );
            Ey6 = Ey6_vec(1).* exp(-1i.*  k_6_EP .* (x6-x6(1))  ) ...
                + Ey6_vec(2).* exp( 1i.*  k_6_EP .* (x6-x6(1))  );
            Ey7 = Ey7_vec(1).* exp(-1i.* -k_7_EP .* (x7-x7(1))  )...
                + Ey7_vec(2).* exp( 1i.* -k_7_EP/1e2 .* (x7-x7(1)) );% Ey7 is not exist in theory
            x = [x1,x2,x3,x4,x5,x6,x7];
            Ey = [Ey1,Ey2,Ey3,Ey4,Ey5,Ey6,Ey7];
            
            %%%%% alpha
            alpha = trapz(x4,Ey4.^2) / ( trapz(x2,Ey2.^2) + trapz(x4,Ey4.^2) + trapz(x6,Ey6.^2) );
            alpha2a(jj,mm) = alpha;
            
            %%%%% Hx
            Hx1 = -beta_EP/omega/mu_1 .*(Ey1);
            Hx2 = -beta_EP/omega/mu_2 .*(Ey2);
            Hx3 = -beta_EP/omega/mu_3 .*(Ey3);
            Hx4 = -beta_EP/omega/mu_4 .*(Ey4);
            Hx5 = -beta_EP/omega/mu_5 .*(Ey5);
            Hx6 = -beta_EP/omega/mu_6 .*(Ey6);
            Hx7 = -beta_EP/omega/mu_7 .*(Ey7);
            Hx = [Hx1,Hx2,Hx3,Hx4,Hx5,Hx6,Hx7];
            
            %%%%% Hz
            Hz1 = k_1_EP/mu_1*Ey1;
            Hz2 = + k_2_EP / mu_2 * Ey2_vec(1).* exp(-1i.*  k_2_EP .* (x2-x2(1))  ) ...
                  - k_2_EP / mu_2 * Ey2_vec(2).* exp( 1i.*  k_2_EP .* (x2-x2(1))  );
            Hz3 = - k_3_EP / mu_3 * Ey3_vec(1).* exp(-1i.* -k_3_EP .* (x3-x3(1))  ) ...
                  + k_3_EP / mu_3 * Ey3_vec(2).* exp( 1i.* -k_3_EP .* (x3-x3(1))  );
            Hz4 = + k_4_EP / mu_4 * Ey4_vec(1).* exp(-1i.*  k_4_EP .* (x4-x4(1))  ) ...
                  - k_4_EP / mu_4 * Ey4_vec(2).* exp( 1i.*  k_4_EP .* (x4-x4(1))  );
            Hz5 = - k_5_EP / mu_5 * Ey5_vec(1).* exp(-1i.* -k_5_EP .* (x5-x5(1))  ) ...
                  + k_5_EP / mu_5 * Ey5_vec(2).* exp( 1i.* -k_5_EP .* (x5-x5(1))  );
            Hz6 = + k_6_EP / mu_6 * Ey6_vec(1).* exp(-1i.*  k_6_EP .* (x6-x6(1))  ) ...
                  - k_6_EP / mu_6 * Ey6_vec(2).* exp( 1i.*  k_6_EP .* (x6-x6(1))  );
            Hz7 = - k_7_EP / mu_7 * Ey7_vec(1).* exp(-1i.* -k_7_EP .* (x7-x7(1))  ) ...
                  + k_7_EP / mu_7 * Ey7_vec(2).* exp( 1i.* -k_7_EP/1e2 .* (x7-x7(1)));
            Hz = [Hz1,Hz2,Hz3,Hz4,Hz5,Hz6,Hz7];
            Hz = Hz/omega;
            
            %%%%% \partial\epsilon / \partial\omega & \partial\mu / \partial\omega
            pp_epsilon_p = epsilon_p;
            pp_mu_p = mu_p;
            pp_epsilon_n = 1 + omega_p^2/omega^2;
            pp_mu_n = 1 - F * (omega^4 - 3*omega_0^2*omega^2)/(omega^2 - omega_0^2)^2;
            pp_epsilon_n = pp_epsilon_n * epsilon_0;
            pp_mu_n = pp_mu_n * mu_0;

            pp_epsilon_1 = epsilon_r;    pp_mu_1 = mu_r;
            pp_epsilon_2 = pp_epsilon_n; pp_mu_2 = pp_mu_n;
            pp_epsilon_3 = epsilon_r;    pp_mu_3 = mu_r;
            pp_epsilon_4 = pp_epsilon_p; pp_mu_4 = pp_mu_p;
            pp_epsilon_5 = epsilon_r;    pp_mu_5 = mu_r;
            pp_epsilon_6 = pp_epsilon_p; pp_mu_6 = pp_mu_p;
            pp_epsilon_7 = epsilon_r;    pp_mu_7 = mu_r;
            
            %%%%% energy density
            w1 = pp_epsilon_1 * Ey1.*conj(Ey1) + pp_mu_1 * ( Hx1.*conj(Hx1) + Hz1.*conj(Hz1) );
            w2 = pp_epsilon_2 * Ey2.*conj(Ey2) + pp_mu_2 * ( Hx2.*conj(Hx2) + Hz2.*conj(Hz2) );
            w3 = pp_epsilon_3 * Ey3.*conj(Ey3) + pp_mu_3 * ( Hx3.*conj(Hx3) + Hz3.*conj(Hz3) );
            w4 = pp_epsilon_4 * Ey4.*conj(Ey4) + pp_mu_4 * ( Hx4.*conj(Hx4) + Hz4.*conj(Hz4) );
            w5 = pp_epsilon_5 * Ey5.*conj(Ey5) + pp_mu_5 * ( Hx5.*conj(Hx5) + Hz5.*conj(Hz5) );
            w6 = pp_epsilon_6 * Ey6.*conj(Ey6) + pp_mu_6 * ( Hx6.*conj(Hx6) + Hz6.*conj(Hz6) );
            w7 = pp_epsilon_7 * Ey7.*conj(Ey7) + pp_mu_7 * ( Hx7.*conj(Hx7) + Hz7.*conj(Hz7) );
            w = [w1,w2,w3,w4,w5,w6,w7];
            W = 1/4*trapz(x,w);
            sz = - 1/2.*Ey.*conj(Hx);
            Sz = trapz(x,sz);

            vg = Sz/W;
            velocity_sz(jj,mm) = vg;
        end
    end
end
%%%%% Data correction, the 1st data in the list is not 1st in fact
num_mid1 = find(beta2a(1,:)>0,1) - 1;
num_mid3 = find(beta2a(3,:)>0,1) - 1;
beta2a(1,num_mid1-3:num_mid1) = beta2a(2,num_mid1-3:num_mid1);
beta2a(3,num_mid3-2:num_mid3) = beta2a(2,num_mid3-2:num_mid3);
alpha2a(1,num_mid1-3:num_mid1) = alpha2a(2,num_mid1-3:num_mid1);
alpha2a(3,num_mid3-2:num_mid3) = alpha2a(2,num_mid3-2:num_mid3);
velocity_sz(1,num_mid1-3:num_mid1) = velocity_sz(2,num_mid1-3:num_mid1);
velocity_sz(3,num_mid3-2:num_mid3) = velocity_sz(2,num_mid3-2:num_mid3);

%% figure
%%%%% curve
figure('numbertitle','off','name','save curve');
hold on
data_gap = 1;
% plot( list_a(1:data_gap:end)*100, beta2a(1,1:data_gap:end), 'o', 'LineWidth', 2, 'MarkerFaceColor', [204/256 0 0]  );
% plot( list_a(1:data_gap:end)*100, beta2a(2,1:data_gap:end), 'o', 'LineWidth', 2, 'MarkerFaceColor', [0 0 153/256]  );
% plot( list_a(1:data_gap:end)*100, beta2a(3,1:data_gap:end), 'o', 'LineWidth', 2, 'MarkerFaceColor', [0 102/256 0]  );
plot( list_a(1:data_gap:end)*100, beta2a(1,1:data_gap:end), '-', 'LineWidth', 2, 'Color', [204/256 0 0]  );
plot( list_a(1:data_gap:end)*100, beta2a(2,1:data_gap:end), '-', 'LineWidth', 2, 'Color', [0 0 153/256]  );
plot( list_a(1:data_gap:end)*100, beta2a(3,1:data_gap:end), '-', 'LineWidth', 2, 'Color', [0 102/256 0]  );
legend('$\mathrm{curve}\ k_1$','$\mathrm{curve}\ k_2$','$\mathrm{curve}\ k_3$','interpreter','latex','FontName','Times New Roman','FontSize',20,'Location','southeast')
% legend('boxoff')
% title('$E_y\ \mathrm{and}\ S_z\ \mathrm{at}\ A$','interpreter','latex','FontName','Times New Roman','FontSize',20)
xlabel('$a\ (\mathrm{cm})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$k\ (\mathrm{cm}^{-1})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
% xlim([0,5])
% ylim([1.23,1.39])
xline( list_a(num_mid3)*100 , '--b', '$a_c$',  'LineWidth', 1,...% order is important, must be "color->text"
    'interpreter','latex','FontName','Times New Roman','FontSize', 20,...
    'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment','bottom',...
    'LabelOrientation', 'horizontal', 'HandleVisibility','off');
hold off
%%
%%%%% ratio alpha
figure('numbertitle','off','name','ratio alpha');

data_gap = 20;
semilogy( list_a(1:data_gap:end)*100, alpha2a(1,1:data_gap:end), 'o', 'LineWidth', 2, 'MarkerEdgeColor', [204/256 0 0]  );
hold on % when using semilogy or loglog, the "hold on" should be used after them
semilogy( list_a(1:data_gap:end)*100, alpha2a(2,1:data_gap:end), 'o', 'LineWidth', 2, 'MarkerEdgeColor', [0 0 153/256]  );
semilogy( list_a(1:data_gap:end)*100, alpha2a(3,1:data_gap:end), 'o', 'LineWidth', 2, 'MarkerEdgeColor', [0 102/256 0]  );

% plot( list_a(1:data_gap:end)*100, alpha2a(1,1:data_gap:end), 'o', 'LineWidth', 2, 'MarkerEdgeColor', [204/256 0 0] );
% plot( list_a(1:data_gap:end)*100, alpha2a(2,1:data_gap:end), 'o', 'LineWidth', 2, 'MarkerEdgeColor', [0 0 153/256]  );
% plot( list_a(1:data_gap:end)*100, alpha2a(3,1:data_gap:end), 'o', 'LineWidth', 2, 'MarkerEdgeColor', [0 102/256 0]  );
legend('$\mathrm{curve}\ k_1$','$\mathrm{curve}\ k_2$','$\mathrm{curve}\ k_3$','interpreter','latex','FontName','Times New Roman','FontSize',20,'Location','east')
% title('$E_y\ \mathrm{and}\ S_z\ \mathrm{at}\ A$','interpreter','latex','FontName','Times New Roman','FontSize',20)
xlabel('$a\ (\mathrm{cm})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$\alpha$','interpreter','latex','FontName','Times New Roman','FontSize',20)
% xlim([0,5])
% ylim([1.23,1.39])
xline( list_a(num_mid3)*100 , '--b', '$a_c$',  'LineWidth', 1,...% order is important, must be "color->text"
    'interpreter','latex','FontName','Times New Roman','FontSize', 20,...
    'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment','bottom',...
    'LabelOrientation', 'horizontal', 'HandleVisibility','off');
hold off
%%
%%%%% velocity
figure('numbertitle','off','name','velocity');
hold on
data_gap = 30;
plot( list_a(1:data_gap:end)*100, velocity_sz(1,1:data_gap:end) , 'o', 'LineWidth', 2, 'MarkerEdgeColor', [204/256 0 0] );
plot( list_a(1:data_gap:end)*100, velocity_sz(2,1:data_gap:end) , 'o', 'LineWidth', 2, 'MarkerEdgeColor', [0 0 153/256] );
plot( list_a(1:data_gap:end)*100, velocity_sz(3,1:data_gap:end) , 'o', 'LineWidth', 2, 'MarkerEdgeColor', [0 102/256 0] );
legend('$\mathrm{velocity}\ k_1$','$\mathrm{velocity}\ k_2$','$\mathrm{velocity}\ k_3$','interpreter','latex','FontName','Times New Roman','FontSize',20,'Location','northwest')
% title('$E_y\ \mathrm{and}\ S_z\ \mathrm{at}\ A$','interpreter','latex','FontName','Times New Roman','FontSize',20)
xlabel('$a\ (\mathrm{cm})$','interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('$v_g\ (c)$','interpreter','latex','FontName','Times New Roman','FontSize',20)
% x_lim = [3,7];
% y_lim = [-0.5e-4,3e-4];
% xlim(x_lim);
% ylim(y_lim);
xline( list_a(num_mid3)*100 , '--b', '$a_c$',  'LineWidth', 1,...% order is important, must be "color->text"
    'interpreter','latex','FontName','Times New Roman','FontSize', 20,...
    'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment','bottom',...
    'LabelOrientation', 'horizontal', 'HandleVisibility','off');
yline(0, ':b', 'LineWidth', 1, 'HandleVisibility','off');
%%%%% notation
% X = [(list_a(num_mid3)*100*1.2 - x_lim(1)) / ( x_lim(2) -  x_lim(1) ),...
%     (list_a(num_mid3)*100 - x_lim(1)) / ( x_lim(2) -  x_lim(1) ) ];
% Y = [ (0.2e-4 - y_lim(1) ) / ( y_lim(2) -  y_lim(1) ),...
%     (0 - y_lim(1) ) / ( y_lim(2) -  y_lim(1) )];
% annotation('textarrow',X,Y,'String', 'EP ')
hold off
%% Coupled Mode Theory
%%%%% EP3
ac = 0.06692;
k0 =  1.324997000000000e+02;
XX = list_a(num_mid3:end);
YY = beta2a(1,num_mid3:end);
YY = (YY-k0).^2;
XX = XX - ac;
myfittype = fittype('B^2*(1-exp(-2*XX/Ld))',...
    'dependent',{'YY'},'independent',{'XX'},...
    'coefficients',{'B','Ld'});
myfit = fit(XX',YY',myfittype);
YY_fit = myfit(XX);
figure;
plot(XX,YY);
hold on
plot(XX,YY_fit);
hold off

%%%%% EP2
BB = 0.1218;
delta = 6e-3;
Ld = 0.0129;
AA = BB * exp(-(list_a - ac)/Ld);
solutions = zeros(3,length(AA));
for ii = 1:length(AA)
    AA_ii = AA(ii);
    if isnan(AA_ii)
        solutions(:,ii) = nan;
        break;
    end
    H1 = [k0+3*delta,-AA_ii,0;
         AA_ii,k0,BB;
         0,BB,k0];
    [~,DD]=eig(H1);
    DD_sort = [DD(1,1);DD(2,2);DD(3,3)];
    DD_sort = sort(DD_sort);
    solutions(:,ii) = DD_sort;
end
figure;
plot(list_a,solutions(1,:),'.');
hold on
plot(list_a,solutions(2,:),'.');
plot(list_a,solutions(3,:),'.');
hold off







