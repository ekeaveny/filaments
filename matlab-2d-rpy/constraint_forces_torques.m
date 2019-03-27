function  [FX, FY, TAUZ] = constraint_forces_torques(FX_IN, FY_IN, ...
                          TAUZ_IN, TX, TY, LAMBDA1, LAMBDA2, SW_IND, DL_SW)
% CONSTRAINT_FORCES_TORQUES  Adds constraint forces and torques to 
%                            filaments.
%
%   [FX, FY, TAUZ] = constraint_forces_torques(FX_IN, FY_IN, TAUZ_IN,...
%                                              TX, TY, LAMBDA1, LAMBDA2,...
%                                              SW_IND, DL_SW)
%   takes in current forces and torques [FX_IN, FY_IN, TAUZ_IN], tangent 
%   vectors [TX,TY], constraints [LAMBDA1, LAMBDA2], filament indexing 
%   matrix SW_IND and segment separation matrix DL_SW. 
%
%   The function returns the constraint forces and torques + any existing 
%   forces and torques.
%
%   These expressions correspond to F^C and T^C, equations (16) and (17)
%   in the paper.
                                                     
                                                     
FX = FX_IN;
FY = FY_IN;
TAUZ = TAUZ_IN;

N_sw = size(SW_IND,1);
N_w = size(SW_IND,2);

N_lam = N_w - 1;

for i=1:N_sw
    jlam_st = (i-1)*N_lam;
    
    jlam = jlam_st + 1; 
    
    jsw = SW_IND(i,1);
    ksw = SW_IND(i,2);
    DL = DL_SW(i,1);
    
    lam1 = LAMBDA1(jlam);
    lam2 = LAMBDA2(jlam);
    
    FX(ksw) = FX(ksw) + lam1;
    FY(ksw) = FY(ksw) + lam2;
    
    FX(jsw) = FX(jsw) - lam1;
    FY(jsw) = FY(jsw) - lam2;
    
    tempx = -0.5*DL*lam1;
    tempy = -0.5*DL*lam2;
    
    TAUZ(ksw) = TAUZ(ksw) + (TX(ksw)*tempy - TY(ksw)*tempx);
    TAUZ(jsw) = TAUZ(jsw) + (TX(jsw)*tempy - TY(jsw)*tempx);
    
    for j=2:N_lam
        jlam = jlam_st + j; 
        
        jsw = SW_IND(i,j);
        ksw = SW_IND(i,j+1);
        DL = DL_SW(i,j);
        
        lam1 = LAMBDA1(jlam);
        lam2 = LAMBDA2(jlam);
        
        FX(ksw) = FX(ksw) + lam1;
        FY(ksw) = FY(ksw) + lam2;
        
        FX(jsw) = FX(jsw) - lam1;
        FY(jsw) = FY(jsw) - lam2;
        
        tempx = -0.5*DL*lam1;
        tempy = -0.5*DL*lam2;
        
        TAUZ(ksw) = TAUZ(ksw) + (TX(ksw)*tempy - TY(ksw)*tempx);
        TAUZ(jsw) = TAUZ(jsw) + (TX(jsw)*tempy - TY(jsw)*tempx);
        
    end
end