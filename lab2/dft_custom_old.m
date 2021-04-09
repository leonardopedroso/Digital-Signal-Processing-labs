function [absX,angX,f] = dft_custom(x,N)
    dft = fft(x,N);
    dft = dft(N/2+1:end)/length(x);
    dft(2:end) = 2*dft(2:end); % do not double 0
    f = 2*pi*(0:N/2-1)'/N;
    absX = abs(dft);
    angX = angle(dft);
end