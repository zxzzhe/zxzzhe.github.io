function [M,r12,r21,t12,t21] = M1_ReflAndTran(n1,n2)
% Reflectance And Transmittance
% k1, input wavevector(vertical)
% k2, output wavevector(vertical)
% mu1,input direction medium
% mu2,output direction medium
t12 = 2*n1/(n1+n2);
r12 = (n1-n2)/(n1+n2);
t21 = 2*n2/(n1+n2);
r21 = - r12;
M = [ t12 - r21*r12/t21 , r21/t21 ; -r12/t21 , 1/t21];
end


