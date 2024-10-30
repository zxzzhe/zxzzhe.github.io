function [M,rr1,rr2,tt1,tt2] = M1_ReflAndTran(k1,k2,mu1,mu2)
% Reflectance And Transmittance
% k1, input wavevector(vertical)
% k2, output wavevector(vertical)
% mu1,input direction medium
% mu2,output direction medium
tt1 = 2 / ( 1+k2*mu1/k1/mu2 );
rr1 = tt1 - 1;
tt2 = 2 / ( 1+k1*mu2/k2/mu1 );
rr2 = - rr1;
M = [ tt1 - rr2*rr1/tt2 , rr2/tt2 ; -rr1/tt2 , 1/tt2];
end

