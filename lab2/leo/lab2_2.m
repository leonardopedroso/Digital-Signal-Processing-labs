% lab2_1.m
% 89683 - José Neves
% 89691 - Leonardo Pedroso

clear
close all

%% R2.a)
% Saves the sound of the Professor reading a sentence as a vector and
% sampling frequency
[sentence, fs] = audioread('How_many_roads.wav');
sentence = sentence(1:end-1,:); % so the duration of the signal is even

%%
% Plays the sound of the sentence
soundsc(sentence,fs)

%% R2.b)
M = 2048;
initial_sample = 48500;

% Gets the indices of the segment
n_segment = initial_sample:initial_sample+M-1;

% Declares the segment as a part of the sentence starting at the index
% 48500 and of length 2048
segment = sentence(initial_sample:initial_sample+M-1,1);

%% R2.c)
N = 2048;

[abs_segment_dft,ang_segment_dft,f] = dft_custom(segment,N);

%%
% Creates the magnitude spectrum
figure('units','normalized','outerposition',[0 0 1 1])
plot(f,abs_segment_dft)
set(gca,'FontSize',24)
xlabel("$2\pi k/N\;[\mathrm{rad}]$",'Interpreter','latex');
ylabel("$ | X(k) |$",'Interpreter','latex');

% Saves it as a PDF
if false
    set(gcf,'Units','Inches');
    pathFigPos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[pathFigPos(3), pathFigPos(4)])
    print(gcf,'figures/mag_segment','-dpdf','-r0')
end

% Creates the phase spectrum
figure('units','normalized','outerposition',[0 0 1 1])
plot(f,ang_segment_dft)
set(gca,'FontSize',24)
xlabel("$2\pi k/N\;[\mathrm{rad}]$",'Interpreter','latex');
ylabel("$\mathrm{arg}(X(k)) \;[rad]$",'Interpreter','latex');

% Saves it as a PDF
if false
    set(gcf,'Units','Inches');
    pathFigPos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[pathFigPos(3), pathFigPos(4)])
    print(gcf,'figures/phase_segment','-dpdf','-r0')
end

%%
% obtains the peaks of the magnitude spectrum
[amp, indexes] = findpeaks(abs_segment_dft);

% chooses the 3 largest
[~,i] = maxk(amp,3);
% obtains the indexes of the 3 largest peaks
coefs = indexes(i);

% reconstructs the signals from these indexes only
[segment_r,amp_p,phase_p,f_p] = reconstruction(coefs,abs_segment_dft,...
                                    ang_segment_dft,N,n_segment);

% chooses the 10 largest
[~,j] = maxk(amp,5);
% obtains the indexes of the 10 largest peaks
coefs = indexes(j);

% reconstructs the signals from these indexes only
[segment_r2,~,~,~] = reconstruction(coefs,abs_segment_dft,...
                                    ang_segment_dft,N,n_segment);

%%
% creates a plot of the segment and its reconstruction
figure('units','normalized','outerposition',[0 0 1 1]);
plot(n_segment,segment,'g')
hold on
plot(n_segment, segment_r,'r')
plot(n_segment, segment_r2,'b')
hold off
set(gca,'FontSize',24)
legend('Segment','Reconstruction with 3','Reconstruction with 5',...
    'Location','best')
xlabel('n')

if false
    set(gcf,'Units','Inches');
    pathFigPos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[pathFigPos(3), pathFigPos(4)])
    print(gcf,'figures/segment_reconstruction','-dpdf','-r0')
end

%% Applies the same method to the complete signal
N = length(sentence);

[mag_song_dft1,ang_song_dft1,~] = dft_custom(sentence(:,1),N);

%%
S = 0.15e-3;
n_song = 1:length(sentence);

% chooses all the peaks above with magnitude larger than S
[~,coefs1] = findpeaks(mag_song_dft1,'MinPeakHeight',S);
song_r1 = reconstruction(coefs1,mag_song_dft1,ang_song_dft1,N,n_song);

%%
soundsc(song_r1,fs)