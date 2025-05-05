function [dPr_mp,dPhi_mp] = multipath_pseudorange_delay(delay0,delay1,CNR0,CNR1,d_phase,time_spacing,chip_length)

% multipath pseudorange error simulation (unit: chip)
%
% [1] Liu, Liyu, and Moeness G. Amin. "Tracking performance and average 
%     error analysis of GPS discriminators in multipath." 
%     Signal Processing 89.6 (2009): 1224-1239.

% BW = 10*log10(RF_BW);%bandwidth of observation (dB)
% SNR0 = CNR0 - BW;%open-sky snr in dB
% SNR0_ = 10^(SNR0/20); %(dBW)
% SNR1 = CNR1 - BW;%multipath snr in dB
% SNR1_ = 10^(SNR1/20); %(dBW)
alpha = 10^((CNR1-CNR0)/20);
beta = d_phase;%can be extended with complete form
d = time_spacing;
MP_delay = delay1-delay0;
dtao = MP_delay/chip_length;%(chip)

% Code multipath
if alpha<=1
    %by eq (30)
    d1 = ((1+2*alpha*cos(beta)+alpha^2)*d)/(2*(1+alpha*cos(beta)));
    d2 = ((alpha*d*cos(beta)+1)*(2-d)-d*(1-d/2-alpha^2*d/2))/(alpha*d*cos(beta)+2-d);

    if dtao >= 0 && dtao < d1
        dPr_mp = (alpha*(alpha+cos(beta))*dtao)/(1+2*alpha*cos(beta)+alpha^2);
    elseif dtao >= d1 && dtao < d2
        dPr_mp = (((alpha*cos(beta)*(1-dtao)-alpha^2*d/2+1-d/2)^2 + ...
                   2*alpha^2*d*cos(beta)*cos(beta)*(1-d/2) + ...
                   2*alpha^3*d*cos(beta)*(1-dtao))^0.5 - ...
                  (alpha*cos(beta)*(1-dtao)-alpha^2*d/2+1-d/2)) / ...
                 (2*alpha*cos(beta));
    elseif dtao >= d2 && dtao <= (1+d/2)
        dPr_mp = ((alpha*cos(beta)-alpha^2)*dtao+alpha^2*(1+d/2)-alpha*d*cos(beta)-2+d + ...
                  (alpha^2*cos(beta)*cos(beta)*dtao^2+2*alpha*(2-d)*(alpha-cos(beta))*dtao - ...
                   4*alpha^2*cos(beta)*cos(beta)*dtao + ...
                   (2-d)*(2-d+2*alpha*d*cos(beta)) + ...
                   4*alpha^2*(d^2/4-sin(beta)*sin(beta)))^0.5) / ...
                 (2*alpha*cos(beta)-alpha^2);
    elseif dtao > (1+d/2)
        dPr_mp = 0;
    end
    
elseif alpha>1
    %by eq (32)
    d1 = ((1+2*alpha*cos(beta)+alpha^2)*d)/(2*(alpha+cos(beta)));
    d2 = ((alpha+cos(beta))*(2-d)-d*(1-d/2-alpha^2*d/2))/(alpha+cos(beta)+2-d);

    if dtao >= 0 && dtao < d1
        dPr_mp = (alpha*(alpha+cos(beta))*dtao)/(1+2*alpha*cos(beta)+alpha^2);
    elseif dtao >= d1 && dtao < d2
        dPr_mp = (2*alpha*cos(beta)*(1+dtao)+2*alpha^2*(1-d/2)-d-...
                  (4*alpha^2*cos(beta)*cos(beta)*((1-dtao)^2+2*d-d*d)+...
                   (2*alpha^2*(1-d/2)-d)^2+4*alpha*cos(beta)*d*(1-dtao)+...
                   8*alpha^3*cos(beta)*(1-d/2)*(1-dtao))^0.5...
                  )/(4*alpha*cos(beta));
    elseif dtao >= d2 && dtao <= (1+d/2)
        dPr_mp = (2*alpha^2*(1-d/2)+alpha*cos(beta)*(dtao+d)-d/2-1-2*...
                  (4*alpha^4*(1-d/2)^2+alpha^2*cos(beta)*cos(beta)*(2-dtao)^2+...
                   4*alpha^2*(1-d/2)*(dtao-1-d/2)+4*alpha^3*cos(beta)*(1-d/2)*(d-dtao)+...
                   4*alpha^2*(d^2/4-sin(beta)*sin(beta)))^0.5...
                  )/(2*alpha*cos(beta)-1);
    elseif dtao > (1+d/2)
        dPr_mp = 0;
    end
    
end
dPr_mp = dPr_mp * chip_length;
% MP delay threshold
if abs(dPr_mp)>max([delay0,delay1])
    dPr_mp = sign(dPr_mp) * max([delay0,delay1]);
end

% Phase multipath
dPhi_mp = atan(-alpha*sin(beta)/(1+alpha*cos(beta)));
