%% Lab 2 section 1 FIR filtering
load('myFIRFilterCoefficients.mat');

%count number of taps (non-zero values)
taps = 0;
for i = 1:length(h);
    if abs(h(i))>0.000001
        taps = taps+1;
    end
end
taps
figure(1)
stem(h)
title('Lab 2 Q1 A) Plot of FIR Filter Coefficients')
xlim([0 81])
saveas(figure(1), 'LAB_2_Q1_A.svg')

%zero pad the impulse response to 1000
padded = [h,zeros(1,(1000-taps))];
%padded = zeros(1000);
%padded(460:540) = h;
figure(2)
plot(padded)
title('impulse response FIR padded to 1000 samples')
FIRFFT = fft(padded);
% Generate digital frequencies corresponding to the DFT bins
N = length(FIRFFT);
digital_freqs = linspace(0, 2*pi, N);

%plot magnitude and phase of frequency response
figure(3)
tiledlayout(2,1)
nexttile
plot(digital_freqs, mag2db(abs(FIRFFT)))
title('Lab 2 Q1 B Frequency Response of FIR filter')
ylabel('dB')
xlabel('Digital Frequency')
xlim([0 (2*pi)])
nexttile
plot(digital_freqs, unwrap(angle(FIRFFT)))
xlim([0 (2*pi)])
ylabel('Phase Angle')
xlabel('Digital Frequency')
saveas(figure(3), 'LAB_2_Q_1_B.svg')

%% FIR Filter Implementation

%testing the implementation
%g = [1 4 1];
%[inputs, outputs] = myFIRFilter(3, g, 100, 4)

load('myFIRFilterCoefficients.mat');
[inputs, outputs] = myFIRFilter(length(h), h, 100, 1000);
%%
figure(4)
tiledlayout(2,1)
nexttile
plot(inputs)
title('Lab 2 Q2 D Filter Input Signal (Time Domain)')
ylabel('Mag (linear)')
xlabel('Time (n)')
xlim([0 100])
ylim([-1 5])
nexttile
plot(outputs)
title('Lab 2 Q2 D Filter Output Signal (Time Domain)')
ylabel('Mag (linear)')
xlabel('Time (n)')
xlim([0 100])
saveas(figure(4), 'LAB_2_Q_2_D.svg')

%%
fftin = fft(inputs);
fftout = fft(outputs);
N = length(fftin)
digital_freqs = linspace(0, 2*pi, N);

figure(5)
tiledlayout(2,1)
nexttile
plot(digital_freqs, (abs(fftin)))
title('Lab 2 Q2 E Filter Input Signal (Frequency Spectrum)')
ylabel('Mag Linear')
xlabel('Digital Frequency')
nexttile
plot(digital_freqs, (abs(fftout)))
title('Lab 2 Q2 E Filter Output Signal (Frequency Spectrum)')
ylabel('Mag Linear')
xlabel('Digital Frequency')
saveas(figure(5), 'LAB_2_Q_2_E.svg')
%% IIR Filter Analysis
clear all;
clc;
load('myIIRFilterCoefficients.mat')

figure(6)
title('Lab 2 Q3 A')
tiledlayout(2,1)
nexttile
stem(a)
title('IIR Coeficients a')
xlabel('index')
nexttile
stem(b)
title('IIR Coeficients b')
xlabel('index')
saveas(figure(6), 'LAB_2_Q3_A.svg')

figure(7)
freqz(b, a)
saveas(figure(7), 'Lab_2_Q3_B.svg')

%% IIR Filter Implementation
clear all;
clc;
load('myIIRFilterCoefficients.mat');
%Test my implementation
%bk = [0.5 1 0.5];
%ak = [1 1 0.5];
%[inputs, outputs] = myIIRFilter(3, 3, bk, ak, 100, 4)

[inputs, outputs] = myIIRFilter(length(b), length(a), b, a, 100, 1000);

figure(8)
tiledlayout(2,1)
nexttile
plot(inputs)
title('Lab 2 Q3 D Filter Input Signal (Time Domain)')
ylabel('Mag (linear)')
xlabel('Time (n)')
xlim([0 100])
ylim([-3.5 3.5])
nexttile
plot(outputs)
title('Lab 2 Q3 D Filter Output Signal (Time Domain)')
ylabel('Mag (linear)')
xlabel('Time (n)')
xlim([0 100])
saveas(figure(8), 'LAB_2_Q_4_D.svg')

fftin = fft(inputs);
fftout = fft(outputs);
N = length(fftin);
digital_freqs = linspace(0, 2*pi, N);

figure(9)
tiledlayout(2,1)
nexttile
plot(digital_freqs, (abs(fftin)))
title('Lab 2 Q3 E Filter Input Signal (Frequency Spectrum)')
ylabel('Mag Linear')
xlabel('Digital Frequency')
nexttile
plot(digital_freqs, (abs(fftout)))
title('Lab 2 Q3 E Filter Output Signal (Frequency Spectrum)')
ylabel('Mag Linear')
xlabel('Digital Frequency')
saveas(figure(9), 'LAB_2_Q_3_E.svg')