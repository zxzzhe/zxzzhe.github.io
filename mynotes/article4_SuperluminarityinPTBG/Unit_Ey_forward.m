function [abs_Ey,Ey_out] = Unit_Ey_forward(Ey_in,M2,M3,M4,M5,M6,M7,M8,M9,k1,k2,k3,k4,x1)
    %UNIT_EY_FORWARD 此处显示有关此函数的摘要
    %   此处显示详细说明
    Ey_1_vec = Ey_in;
    Ey_2_vec = M3*M2*Ey_in;
    Ey_3_vec = M5*M4*Ey_2_vec;
    Ey_4_vec = M7*M6*Ey_3_vec;
    Ey_out = M9*M8*Ey_4_vec;
    
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
%     abs_Ey = abs(Ey_all); % amplitude
end

