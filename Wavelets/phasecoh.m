%calculates the phase coherence between 2 matrixes of phases:
%returns a vector of L doubles. one for each freq value

function [coh] = phasecoh(WT1, WT2)
phi1=angle(WT1);
phi2=angle(WT2);
if length(phi2)>length(phi1)
    phi2=phi2(:,1:length(phi1));
elseif length(phi2)<length(phi1)
    phi1=phi1(:,1:length(phi2));
end
dphi = phi2-phi1; 

for k = 1:length(dphi(1,:)) %wrap the angle in +- pi
    for kk = 1:length(dphi(:,1)) 
        if dphi(kk,k) < -pi, dphi(kk,k) = dphi(kk,k) + 2*pi; end;        
        if dphi(kk,k) > pi,  dphi(kk,k) = dphi(kk,k) - 2*pi;  end;
    end
end

cphi = cos(dphi); 
sphi = sin(dphi);
coh = sqrt( (nanmean(cphi,2)).^2 + (nanmean(sphi,2)).^2 ); %mean along rows

end

