function [abs_Ex,Ex_out] = Unit_Ex_forward(Ex_in,M1,M2,M3,M4,M5,M6,M7,M8,k1,k2,k3,k4,x1)
    %  wave propagation Ex: | M2 | M4 | M6 | M8 
    %  "|" is boundary between two layers
    Ex_1_vec = M1*Ex_in;
    Ex_2_vec = M3*M2*Ex_1_vec;
    Ex_3_vec = M5*M4*Ex_2_vec;
    Ex_4_vec = M7*M6*Ex_3_vec;
    Ex_out = M8*Ex_4_vec;   % output
    
    Ex_1 = Ex_1_vec(1)*exp(-1i .* k1 .* x1)...
         + Ex_1_vec(2)*exp( 1i .* k1 .* x1);
    Ex_2 = Ex_2_vec(1)*exp(-1i .* k2 .* x1 )...
         + Ex_2_vec(2)*exp( 1i .* k2 .* x1 );
    Ex_3 = Ex_3_vec(1)*exp(-1i .* k3 .* x1 )...
         + Ex_3_vec(2)*exp( 1i .* k3 .* x1 );
    Ex_4 = Ex_4_vec(1)*exp(-1i .* k4 .* x1 )...
         + Ex_4_vec(2)*exp( 1i .* k4 .* x1 );
    Ex_all = [Ex_1,Ex_2,Ex_3,Ex_4]; % distribution in one unit
    abs_Ex = (Ex_all);
%     abs_Ex = abs(Ex_all); % amplitude
end

