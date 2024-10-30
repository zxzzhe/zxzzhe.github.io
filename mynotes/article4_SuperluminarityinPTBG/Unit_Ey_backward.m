function [abs_Ey,Ey_in] = Unit_Ey_backward(Ey_out,M2,M3,M4,M5,M6,M7,M8,M9,k1,k2,k3,k4,x1)
    %UNIT_EY_DIST 此处显示有关此函数的摘要
    %   [E1+,E1-]' = M [E0+,E0-]';
    %   for unidirectional wave, E0+ only have reflection E0-
    %   and refraction E1+, which requires E1- = 0. in this case, set
    %   Ey_out = [1,0]',and use TMM to caculate E0+ and E0-
    %%%%reflect and transmission%%%%
    Ey_4_vec = inv(M8)*inv(M9)*Ey_out;
    Ey_3_vec = inv(M6)*inv(M7)*Ey_4_vec;
    Ey_2_vec = inv(M4)*inv(M5)*Ey_3_vec;
    Ey_1_vec = inv(M2)*inv(M3)*Ey_2_vec;
    Ey_in = Ey_1_vec;
    Ey_1 = Ey_1_vec(1)*exp(-1i .* k2 .* x1)...
         + Ey_1_vec(2)*exp( 1i .* k2 .* x1);
    Ey_2 = Ey_2_vec(1)*exp(-1i .* k3 .* x1 )...
         + Ey_2_vec(2)*exp( 1i .* k3 .* x1 );
    Ey_3 = Ey_3_vec(1)*exp(-1i .* k4 .* x1 )...
         + Ey_3_vec(2)*exp( 1i .* k4 .* x1 );
    Ey_4 = Ey_4_vec(1)*exp(-1i .* k1 .* x1 )...
         + Ey_4_vec(2)*exp( 1i .* k1 .* x1 );
    Ey_all = [Ey_1,Ey_2,Ey_3,Ey_4];
    abs_Ey = (Ey_all);
%     abs_Ey = abs(Ey_all);
end

