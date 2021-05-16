% From the index of the peaks gets the amplitudes, phases, and frequencies
% of each peak and computes the reconstructed signal from the peaks

function [x_r,amp_p,phase_p,f_p] = reconstruction(peaks,mag_x_dft,...
                                                    ang_x_dft,N,n)
    amp_p = mag_x_dft(peaks);
    phase_p = ang_x_dft(peaks);
    f_p = (peaks-1)*2*pi/N;     % the frequencies start at 0

    % the reconstructed signals is the sum of the harmonic of each index
    x_r = amp_p(1)*cos(f_p(1)*n+phase_p(1));
    for i=2:length(peaks)
        x_r = x_r + amp_p(i)*cos(f_p(i)*n+phase_p(i));
    end
end

