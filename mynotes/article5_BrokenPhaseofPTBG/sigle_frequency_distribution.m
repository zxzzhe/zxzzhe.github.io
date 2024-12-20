function [Ex_kn] = sigle_frequency_distribution(ii,d,NN_total,NN_PT,xxx,N_xxx,k1,k2,k3,k4,ks,n1,n2,n3,n4,ns)
    %SIGLE_FREQUENCY_DISTRIBUTION 
    %   Ex: | M2 | M4 | M6 | M8 
    [M_be,~,~,~,~] = M1_ReflAndTran(ns, n1 ); 
    [M_nd,~,~,~,~] = M1_ReflAndTran(n4, ns );

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
    
end

