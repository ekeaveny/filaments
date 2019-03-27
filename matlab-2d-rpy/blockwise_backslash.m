function [x] = blockwise_backslash(A,b,SW_IND)
% BLOCKWISE_BACKSLASH  Performs A\b but expects A in a certain form.
%         
%   x = BLOCKWISE_BACKSLASH(A,b,SW_IND) returns x = A\b where
%   A is expected to be of the form [A1 A2 A3 ...], which is a compact 
%   version of the form A = [A1  O  O ..
%                             O A2  O ..
%                             O  O A3 ..
%                             :  :  :   ]. 
%
%   Each block A1, A2, etc., is of size (3*N_w, 3*N_w), where N_w is found 
%   through the indexing matrix SW_IND.

N_sw_local = size(SW_IND,1);
N_w_local = size(SW_IND,2);
blocksize = N_w_local*3;

x = zeros(N_sw_local*blocksize,1);

for w = 0:N_sw_local-1
    block = (w*blocksize+1:(w+1)*blocksize);
    x(block) = A(:,block)\b(block);
end

end