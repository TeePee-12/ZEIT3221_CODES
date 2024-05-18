%% Generate wave (1) - N=40
Fs = 50000; % [kHz]
Ts = 1/Fs; % [s]
% define time interval
t = 0.:Ts:(39*Ts); % [s]
% create signal
x = 4*sin(35000*pi*t)-2*cos(3750*pi*t);

% calculate FFT of the signal
X  = fft(x);
L = length(X);
% sampling frequency
f = Fs*(0:(L-1))/(L-1);
figure (1)
stem(x) 
title("Lab 1: 1(b) - Stem Plot of Signal (1), N=40")
ylabel("Magnitude")
xlabel("Integer Index [k]")
saveas(figure (1), 'lab_1_1b.png')

figure (2)
tiledlayout(2,1)
nexttile
stem(abs(X))
title("Lab 1: 1(c) - FFT Stem Plot of Signal (1), N=40")
ylabel("Magnitude")
nexttile
stem(angle(X))
ylabel("Phase Shift")
set(gca,'YTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
saveas(figure (2), 'LAB_1_1c.png')

%% Generate wave (1) - N=80
Fs = 50000; % [kHz]
Ts = 1/Fs; % [s]
% define time interval
t = 0.:Ts:(79*Ts); % [s]
% create signal
x = 4*sin(35000*pi*t)-2*cos(3750*pi*t);

% calculate FFT of the signal
X  = fft(x);
L = length(X);
X(abs(X) < 0.0001) = 0;
% sampling frequency
f = Fs*(0:(L-1))/(L-1);
figure (3)
stem(x) 
title("Lab 1: 1(e) - Stem Plot of Signal (1), N=80")
ylabel("Magnitude")
xlabel("Integer Index [k]")
saveas(figure (3), 'LAB_1_1e.png')

figure (4)
tiledlayout(2,1)
nexttile
stem(abs(X))
title("Lab 1: 1(f) - FFT Stem Plot of Signal (1), N=80")
ylabel("Magnitude")
nexttile
stem(angle(X))
ylabel("Phase Shift")
set(gca,'YTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
saveas(figure (4), 'LAB_1_1f.png')

index=(-40:1:39);
figure (5)
tiledlayout(2,1)
nexttile
stem(index,abs(fftshift(X)))
title("Lab 1: 1(h) - Shifted FFT Stem Plot of Signal (1), N=80")
ylabel("Magnitude")
nexttile
stem(index,angle(X))
ylabel("Phase Shift")
set(gca,'YTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
saveas(figure (5), 'LAB_1_1h.png')

%% LAB_1 SECTION 2
TotS = 7787520; %Total Samples 7787520 bits
Fs = 44100; %Sample Rate 44100
Fmax = Fs/2;
Duration = TotS/Fs; %Total song length from file data
Player_Duration = TotS/50000 %Total Song length played on Fs=50 000
FileSizeBytes = TotS*4 %Samples x2 bytes per sample x2 chanels
FileSizeKB    = (TotS*4)/(2^10)
FileSizeMB    = (TotS*4)/(2^20)
tune = audioread('tune.wav');

figure(6)
plot(tune(3100001:3101501,1))
title("Lab 1: 2(c) - 1500 contiguous samples of tune.wav chanel 1")
ylabel("Amplitude");
xlabel("time[s/50000]");
saveas(figure (6), 'LAB_1_2b.png')

DFT_Tune = fft(tune);
a_index = (0:1:(((length(DFT_Tune))/2)-1));
b_index = (-(((length(DFT_Tune))/2)-1):1:0);
index = cat(2,a_index,b_index);
freqs=index*Fs/TotS;

figure (7)
tiledlayout(2,1)
nexttile
ABS_tune = abs(DFT_Tune(:,1));
mag_tune = mag2db(ABS_tune);
plot(freqs, mag_tune);
title("Lab 1: 2(e) - FFT Stem Plot of Tune.wav")
ylabel("Magnitude [dB]")
xlabel("Frequency [Hz]")
nexttile
(plot(freqs,(angle(DFT_Tune(:,1)))));
ylabel("Phase Shift")
xlabel("Frequency [Hz]")
saveas(figure (7), 'LAB_1_2e.png')

%Cutoff Frequency - Digital sample number range
cutoff = round(2000/(Fs/TotS));
mid = round(TotS/2);
Mod_DFT_Tune = DFT_Tune;
Mod_DFT_Tune(cutoff:TotS-cutoff,:) = 0;

Mod_Tune = ifft(Mod_DFT_Tune,'Symmetric');

figure(9)
plot(Mod_Tune(3100001:3101501,1))
title("Lab 1: 2(f) - 1500 contiguous samples of tune.wav chanel 1 after modifying the FFT")
ylabel("Amplitude");
xlabel("time[s/50000]");
saveas(figure (9), 'LAB_1_2f.png')

myplayer = audioplayer(tune, Fs);
modplayer = audioplayer(Mod_Tune, Fs);

%% Section 3 - Own Recording

myRecording = audiorecorder
recordblocking(myRecording, 4)
play(myRecording)
x = getaudiodata(myRecording)
%% 
Fs = 8000;
TotS = 32000;
X=fft(x);
a_index = (0:1:(((length(X))/2)-1));
b_index = (-(((length(X))/2)-1):1:0);
index = cat(2,a_index,b_index);
freqs=index*Fs/TotS;

figure(10)
plot(freqs, mag2db(abs(X)));
title("dB Scale Voice Recording Frequency Spectrum")
xlabel("Frequency [Hz]")
ylabel("Magnitude [dB]")
saveas(figure(10), "Lab_1_3_b.png")

%% Low Freq
recordblocking(myRecording, 4)
play(myRecording)
x_lo = getaudiodata(myRecording)
Fs = 8000;
TotS = 32000;
X_LO=fft(x_lo);
a_index = (0:1:(((length(X_LO))/2)-1));
b_index = (-(((length(X_LO))/2)-1):1:0);
index = cat(2,a_index,b_index);
freqs=index*Fs/TotS;

figure(11)
plot(freqs, abs(X_LO));
title("Low Vocal Tone Frequency Spectrum")
xlabel("Frequency [Hz]")
ylabel("Magnitude")

%% Hi Freq
recordblocking(myRecording, 4)
play(myRecording)
x_hi = getaudiodata(myRecording)
Fs = 8000;
TotS = 32000;
X_HI=fft(x_hi);
a_index = (0:1:(((length(X_HI))/2)-1));
b_index = (-(((length(X_HI))/2)-1):1:0);
index = cat(2,a_index,b_index);
freqs=index*Fs/TotS;

figure(12)
plot(freqs, abs(X_HI));
title("Soprano Vocal Tone Frequency Spectrum")
xlabel("Frequency [Hz]")
ylabel("Magnitude")
%% Hi-Lo PLot
figure(13)
tiledlayout(2,1)
nexttile
plot(freqs, abs(X_LO));
title("Low Vocal Tone Frequency Spectrum")
xlabel("Frequency [Hz]")
ylabel("Magnitude")
hold on
nexttile
plot(freqs, abs(X_HI));
title("High Vocal Tone Frequency Spectrum")
xlabel("Frequency [Hz]")
ylabel("Magnitude")
saveas(figure(13),"Lab_1_3_c.png")
%% Scales
myRecording = audiorecorder(3200,16,1)
recordblocking(myRecording, 10)
play(myRecording)
x = getaudiodata(myRecording)

figure(14)
spectrogram(x)
title("Vocal Scales Spectrogram, Fs=3200Hz")
%xlabel("Frequency [Hz]")
%ylabel("Magnitude")
saveas(figure(14),"Frequency Spectrum.png")
