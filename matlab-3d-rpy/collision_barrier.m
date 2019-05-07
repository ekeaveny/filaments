function [] = collision_barrier(Filaments)
% COLLISION_BARRIER  Places a short-range steric repulsion force on 
%                    segments closer than a certain cutoff.
%
%   collision_barrier(Filaments)  adds the steric repulsion forces to
%   each Filament.F
%
%   It is not set up for periodic domains.

N_sw = length(Filaments);

for i=1:N_sw
    Xi = Filaments(i).X; % Position
    Ri = Filaments(i).R; % Radius
    Fi = Filaments(i).F; % Force
    
    for j=i:N_sw
        Xj = Filaments(j).X;
        Rj = Filaments(j).R;
        Fj = Filaments(j).F;
        
        x = norm(Filaments(i).CentreOfMass - Filaments(j).CentreOfMass);
        d = 0.5*(Filaments(i).N_w * Filaments(i).DL ... 
               + Filaments(j).N_w * Filaments(j).DL);
        
        if x <= d
            for m=1:Filaments(i).N_w
                Xm = Xi(:,m);
                r = Ri(m);
                for n=1:Filaments(j).N_w
                    diff = Xm - Xj(:,n);
                    dist2 = diff(1)*diff(1) + diff(2)*diff(2) ...
                          + diff(3)*diff(3);
                    r2 = (r + Rj(n))^2;
                    chi = 1.21 * r2;
                    if dist2 < chi
                        fac = 15*pi*((chi - dist2)/(chi - r2))^3;
                        Force = fac*diff;
                        Fi(1:3,m) = Fi(1:3,m) + Force;
                        if i~=j % So that we don't cancel out barrier 
                                % forces in the same filament.
                            Fj(1:3,n) = Fj(1:3,n) - Force;
                        end
                    end
                end
            end
        end
        Filaments(j).F = Fj;
    end
    Filaments(i).F = Fi;
end

end

