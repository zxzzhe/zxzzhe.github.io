function [M] = M2_propagation(k,d)
% propagation
% k, Vertical wavevector
% d, propagation distance, between WGs(a) or length of WGs(b)
M = [exp(-1i*k*d) , 0 ; 0 , exp(1i*k*d)];
end

