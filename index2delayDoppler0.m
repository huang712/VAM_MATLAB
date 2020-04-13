function [delay,doppler] = index2delayDoppler(index)
% delay = 1:17;
% Doppler = 1:11
    if(mod(index,17)==0)
          doppler = floor(index/17);
    else
          doppler = floor(index/17)+1;
    end
    delay = index-(doppler-1)*17;
end