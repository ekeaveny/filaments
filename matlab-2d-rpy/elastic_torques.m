function [TAUZ] = elastic_torques(TAUZ_IN, TX, TY, KB_array, SW_IND, DL_SW)
% ELASTIC_TORQUES  Adds elastic torques to filaments.
%
%   [TAUZ] = elastic_torques(TAUZ_IN, TX, TY, KB_array, SW_IND, DL_SW)
%   takes in current torques TAUZ_IN, tangent vectors [TX,TY], filament
%   indexing matrix SW_IND and segment separation matrix DL_SW. 
%
%   It also accepts a vector of bending moduli, KB_array, which can be
%   different for different filaments.
%
%   The function returns the elastic torque + any existing torque.
%
%   These expressions correspond to T^E, equation (15), in the paper.

TAUZ = TAUZ_IN;
N_worm = size(SW_IND,2);
N_sw = size(SW_IND,1);

for i = 1:N_sw
    num_filament_types = size(KB_array,2);
    k_index = floor((i-1)/N_sw*num_filament_types) + 1; 
    KB = KB_array(k_index);    
    
    ind_st = SW_IND(i,1);
    
    tx = TX(ind_st+1);
    ty = TY(ind_st+1);
    
    KB_div_DL = KB./DL_SW(i,1);
    
    TAUZ(ind_st) = TAUZ(ind_st) + KB_div_DL*(TX(ind_st)*ty - tx*TY(ind_st));
    
    for j = 2:N_worm-1
        jsw = SW_IND(i,j);
        jswp1 = SW_IND(i,j+1);
        jswm1 = SW_IND(i,j-1);
        
        tx = TX(jswm1);    
        ty = TY(jswm1);
        
        KB_div_DL = KB./DL_SW(i,j-1);
    
        TAUZ(jsw) = TAUZ(jsw) + KB_div_DL*(TX(jsw)*ty - tx*TY(jsw));
        
        tx = TX(jswp1);
        ty = TY(jswp1);
        
        KB_div_DL = KB./DL_SW(i,j);
        
        TAUZ(jsw) = TAUZ(jsw) + KB_div_DL*(TX(jsw)*ty - tx*TY(jsw));
    end
    ind_edm1 = SW_IND(i,N_worm - 1);
    ind_ed = SW_IND(i,N_worm);
    tx = TX(ind_edm1);
    ty = TY(ind_edm1);
    
    KB_div_DL = KB./DL_SW(i,N_worm-1);
    
    TAUZ(ind_ed) = TAUZ(ind_ed) + KB_div_DL*(TX(ind_ed)*ty - tx*TY(ind_ed));
end
