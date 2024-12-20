function [abs_Ex,Ex_output] = Unit_Ex_backward(Ey_input,M1,M2,M3,M4,M5,M6,M7,M8,k1,k2,k3,k4,x1)
    %   | M2 | M4 | M6 | M8 Ey_input
    %   "|" is the boundry between two layers
    %   begin with the output Ey_out, using transform matrix to caculate
    %   former Ey distribution, k_n is propagation constant in part of PTBG, since it
    %   is vertical, k_n equal to wavevector in part of PTBG.
    %   [E1+,E1-]' = M [E0+,E0-]';
    %   for unidirectional wave, E0+ only have reflection E0-
    %   and refraction E1+, which requires E1- = 0. in this case, set
    %   Ey_out = [1,0]',and use TMM to caculate E0+ and E0-
    %%%%reflect and transmission%%%%
    Ex_4_vec = inv(M8)*Ey_input;
    Ex_3_vec = inv(M6)*inv(M7)*Ex_4_vec;
    Ex_2_vec = inv(M4)*inv(M5)*Ex_3_vec;
    Ex_1_vec = inv(M2)*inv(M3)*Ex_2_vec;
    Ex_output = inv(M1)*Ex_1_vec;
    Ex_1 = Ex_1_vec(1)*exp(-1i .* k1 .* x1)...
         + Ex_1_vec(2)*exp( 1i .* k1 .* x1);
    Ex_2 = Ex_2_vec(1)*exp(-1i .* k2 .* x1 )...
         + Ex_2_vec(2)*exp( 1i .* k2 .* x1 );
    Ex_3 = Ex_3_vec(1)*exp(-1i .* k3 .* x1 )...
         + Ex_3_vec(2)*exp( 1i .* k3 .* x1 );
    Ex_4 = Ex_4_vec(1)*exp(-1i .* k4 .* x1 )...
         + Ex_4_vec(2)*exp( 1i .* k4 .* x1 );
    Ex_all = [Ex_1,Ex_2,Ex_3,Ex_4];
    abs_Ex = (Ex_all);
%     abs_Ex = abs(Ex_all);
end

