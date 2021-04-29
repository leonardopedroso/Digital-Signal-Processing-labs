% obtains the single-sided magnitude and phase spectra of the signal and
% the corresponding normalized frequencies
function [absX,angX,f] = dft_custom(x,N)
    dft = fft(x,N);
    dft = dft(1:N/2+1)/length(x);   % obtains only k between [0,N/2]
    dft(2:end-1) = 2*dft(2:end-1);  % do not double 0 nor N/2
    
    f = 2*pi*(0:length(dft)-1)/N;
    absX = abs(dft);
    angX = angle(dft);
end