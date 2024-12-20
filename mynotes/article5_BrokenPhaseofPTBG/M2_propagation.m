function [M] = M2_propagation(nk0,d)
% propagation
% nk0, Vertical wavevector
% d, propagation distance, between WGs(a) or length of WGs(b)
M = [exp(-1i*nk0*d) , 0 ; 0 , exp(1i*nk0*d)];
end

