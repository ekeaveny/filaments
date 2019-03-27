function [FX, FY] = collision_barrier(X, Y, FX_IN, FY_IN, Lx, Ly, PtoS, map, head, list, RAD)
% COLLISION_BARRIER  Places a short-range steric repulsion force on 
%                    segments closer than a certain cutoff.
%
%   [FX, FY] = collision_barrier(X, Y, FX_IN, FY_IN, Lx, Ly, ...
%                                PtoS, map, head, list, RAD)
%   takes in positions [X,Y], current forces on segments [FX,FY], segment
%   size vector RAD, and returns current forces + steric repulsion forces
%   in [FX,FY]. 
%
%   It is set up for periodic domains of size [Lx,Ly], so for infinite
%   domains, simply give large values of Lx and Ly.
%
%   It is also set up to use a linked list for efficiency, but in the main
%   code, [map, head, list] are hardcoded to one box. For linked list
%   algorithms, see Allen & Tildesley, 'Computer Simulations of Liquids').

ncell = size(head,1);

FX = FX_IN;
FY = FY_IN;

cutoff_radius_multiplier = 1.1;

Lxh = Lx/2.0;
Lyh = Ly/2.0;


for icell = 1:ncell
    i = head(icell);
    while(i > 1e-10)
        
        x_w = X(i);
        y_w = Y(i);
        fxi = FX(i);
        fyi = FY(i);
        isw = PtoS(i);
        
        x_w = x_w - Lx*floor(x_w/Lx);
        y_w = y_w - Ly*floor(y_w/Ly);
        isw = PtoS(i);
        
        j = list(i);
        while(j > 1e-10)
            
            x_o = X(j);
            y_o = Y(j);
            jsw = PtoS(j);
            
            x_o = x_o - Lx*floor(x_o/Lx);
            y_o = y_o - Ly*floor(y_o/Ly);
            
            xij = x_w - x_o;
            yij = y_w - y_o;
            
            xij = xij - Lx*fix(xij/Lxh);
            yij = yij - Ly*fix(yij/Lyh);
            
            rij2 = xij^2 + yij^2;
            rij = rij2^0.5;
            
            rtot = RAD(i) + RAD(j);
            R_ref2 = (cutoff_radius_multiplier*rtot)^2;
            foura2 = (rtot)^2;
            F_ref = 2.5*6.0*pi*rtot;
            if(rij2 < R_ref2 && rij2>1e-10)  
                temp = R_ref2 - foura2;
                temp2 = (R_ref2 - rij2)/temp;
                temp2 = temp2*temp2;
                fxij = F_ref*temp2*temp2*xij/(rij);
                fyij = F_ref*temp2*temp2*yij/(rij);
                
                fxi = fxi + fxij;
                fyi = fyi + fyij;
                
                FX(j) = FX(j) - fxij;
                FY(j) = FY(j) - fyij;
            end
            j = list(j);
        end
        jcello = 4*(icell-1);
        for nghbr=1:4
            jcell = map(jcello + nghbr);
            j = head(jcell);
            while(j > 1e-10)
                
                x_o = X(j);
                y_o = Y(j);
                jsw = PtoS(j);
                
                x_o = x_o - Lx*floor(x_o/Lx);
                y_o = y_o - Ly*floor(y_o/Ly);
                
                xij = x_w - x_o;
                yij = y_w - y_o;
                
                xij = xij - Lx*fix(xij/Lxh);
                yij = yij - Ly*fix(yij/Lyh);
                
                rij2 = xij^2 + yij^2;
                rij = rij2^0.5;
                
                rtot = RAD(i) + RAD(j);
                R_ref2 = (cutoff_radius_multiplier*rtot)^2;
                foura2 = (rtot)^2;
                F_ref = 2.5*6.0*pi*rtot;   
                if(rij2 < R_ref2 && rij2>1e-10)    
                    temp = R_ref2 - foura2;
                    temp2 = (R_ref2 - rij2)/temp;
                    temp2 = temp2*temp2;
                    fxij = F_ref*temp2*temp2*xij/(rij);
                    fyij = F_ref*temp2*temp2*yij/(rij); 
                    
                    if isnan(fyij)
                        keyboard
                    end
                    
                    fxi = fxi + fxij;
                    fyi = fyi + fyij;
                    
                    FX(j) = FX(j) - fxij;
                    FY(j) = FY(j) - fyij;
                end
                j = list(j);
            end
        end
        FX(i) = fxi;
        FY(i) = fyi;
        i = list(i);
    end
end