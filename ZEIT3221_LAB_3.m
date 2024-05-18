Ma=25;
tapsa=[];
for n=-(((Ma-1)/2)):(((Ma-1)/2))
    tapsa = [tapsa, (cos((n*pi)/3)*sinc(n/8))];
end
figure(1)
stem(tapsa)
title("FIR Bandpass Filter, 25 Taps")
xlabel("Taps")
xlim([0 Ma])
saveas(figure(1), "FIR FILTER 25 TAPS.svg")

Mb=111;
tapsb=[];
for n=-(((Mb-1)/2)):(((Mb-1)/2))
    tapsb = [tapsb, (cos((n*pi)/3)*sinc(n/8))];
end
figure(2)
stem(tapsb)
title("FIR Bandpass Filter, 111 Taps")
xlabel("Taps")
xlim([0 Mb])
saveas(figure(2), "FIR FILTER 111 TAPS.svg")
%%
Fs = 160; %Hz
padded_tapsa = [tapsa,zeros(1,(1000-Ma))];
padded_tapsb = [tapsb,zeros(1,(1000-Mb))];

fft_tapsa = (fft(padded_tapsa));
freqsa = linspace(0, Fs, length(fft_tapsa));
fft_tapsb = (fft(padded_tapsb));
freqsb = linspace(0, Fs, length(fft_tapsb));

figure(3)
tiledlayout(2,1)
nexttile
plot(freqsa, abs(fft_tapsa))
title("Phase Resposne, 25 Tap Filter")
xlabel("Frequency [Hz]")
ylabel("Magnitude (linear)")

nexttile
plot(freqsa, unwrap(angle(fft_tapsa)))
xlabel("Frequency [Hz]")
ylabel("Phase")
saveas(figure(3),"25 Taps Freq Response.svg")

figure(4)
tiledlayout(2,1)
nexttile
plot(freqsb, abs(fft_tapsb))
title("Phase Resposne, 111 Tap Filter")
xlabel("Frequency [Hz]")
ylabel("Magnitude (linear)")
nexttile
plot(freqsb, unwrap(angle(fft_tapsb)))
xlabel("Frequency [Hz]")
ylabel("Phase")
saveas(figure(4),"111 Taps Freq Response.svg")
%%
window = blackman(length(tapsb));
windowtaps = zeros(1,length(tapsb));
for n = 1 : length(tapsb)
    windowtaps(n) = window(n)*tapsb(n);
end
windowtapspadded = [windowtaps,zeros(1,(1000-length(windowtaps)))];
windowedfft = fft(windowtapspadded);
figure(5)
tiledlayout(2,1)
nexttile
plot(freqsb, abs(windowedfft))
title("Frequency Resposne, 111 Tap Filter, Blackman Window Application")
xlabel("Frequency [Hz]")
ylabel("Magnitude (linear)")
nexttile
plot(freqsb, unwrap(angle(windowedfft)))
xlabel("Frequency [Hz]")
ylabel("Phase")
saveas(figure(5),"111 Taps Windowed Freq Response.svg")

%% Minimax method filter design
Fs = 160;%Hz
M = 251;
f = [0 0.23 0.25 0.5 0.52 1];
a = [0 0 4 4 0 0];
w = [1,6,1];
[b,err] = firpm(M-1,f,a,w);

figure(6)
stem(b)
xlim([0 M])
title("Filter taps as determined by the MATLAB firpm method")
saveas(figure(6), "FIRPM Filter Taps.svg")

b = [b,zeros(1,(1000-length(b)))];
B = fft(b);
freqsb = linspace(0, Fs, length(B));

figure(7)
tiledlayout(2,1)
nexttile
plot(f*80,a)
hold on;
plot(freqsb, abs(B))
legend('Ideal','firpm Design')
title("Frequency Resposne, 109 Tap Filter, FIRPM Method Application")
xlabel("Frequency [Hz]")
ylabel("Magnitude (linear)")
xlim([0 80])
nexttile
plot(freqsb, unwrap(angle(B)))
xlabel("Frequency [Hz]")
ylabel("Phase")
xlim([0 80])
saveas(figure(7),"Ideal vs firpm FREQ RESP.svg")
err

%% Minimax filter run 2
Fs = 160;%Hz
M = 675;
f = [0 0.295 0.305 0.695 0.705 1.695 1.705 2.095 2.105 pi]/pi;
a = [0 0     2      2     0     0     1     1     0     0];
[b,err] = firpm(M-1,f,a);

figure(8)
stem(b)
xlim([0 M])
title("Filter taps as determined by the MATLAB firpm method - 675 Taps")
saveas(figure(8), "VERSION 2 FIRPM Filter Taps.svg")

b = [b,zeros(1,(1000-length(b)))];
B = fft(b);
freqsb = linspace(0, Fs, length(B));

figure(9)
tiledlayout(2,1)
nexttile
plot(f*80,a)
hold on;
plot(freqsb, abs(B))
legend('Ideal','firpm Design')
title("Frequency Resposne, 675 Tap Filter, FIRPM Method Application")
xlabel("Frequency [Hz]")
ylabel("Magnitude (linear)")
xlim([0 80])
nexttile
plot(freqsb, unwrap(angle(B)))
xlabel("Frequency [Hz]")
ylabel("Phase")
xlim([0 80])
saveas(figure(9),"VERSION 2 Ideal vs firpm FREQ RESP.svg")
err
%% continuous time filter
Fs = 160;%Hz
wp = 30*2*pi;
ws = 40*2*pi;
rp = 3;
rs = 15;
N = cheb1ord(wp,ws,rp,rs,'s');
[NUM,DEN] = cheby1(N,rp,wp,'s');
systemc = tf(NUM,DEN)
figure(10)
bode(systemc)
title("4th-Order CT Chebyshev Filter Prototype")

%%IIR filter Bilinear Transform
figure(11)
[NUMd, DENd] = bilinear(NUM,DEN,Fs);
freqz(NUMd,DENd)
title('Biliear Transform of the Prototype Chebyshev Filter, Freq Response')
saveas(figure(11), "Bilinear Transform Q4B.svg")

figure(12)
systemd = filt(NUMd,DENd,1/Fs);
bode(systemd)
%% Section 5

[NUMd, DENd] = cheby1(4,3,30/80);
figure(13)
freqz(NUMd,DENd)
systemd = filt(NUMd,DENd,1/160)
%% Section 6 Cascaded IIR Filter
Fs = 160; %Hz

[NUMd, DENd] = cheby1(10,1,8/80);
figure(14)
freqz(NUMd,DENd)
title({"IIR 10^{th} Order Chebyshev LPF Frequency Reposne"; "1dB Passband Ripple, 0.1\pi Cutoff, 0dB Passband Gain"})

NUMd = NUMd*db2mag(20);
figure(15)
freqz(NUMd,DENd)
title({"IIR 10^{th} Order Chebyshev LPF Frequency Reposne"; "1dB Passband Ripple, 0.1\pi Cutoff, 20dB Passband Gain"})

systemd = filt(NUMd,DENd,1/Fs);
figure(16)
pzmap(systemd)
title({"IIR 10^{th} Order Chebyshev 20dB Gain Pole Zero Plot"})

[zeros,poles,b0]=tf2zpk(NUMd,DENd)
format long 
[sos,b0]=tf2sos(NUMd,DENd)

%% section 7 cacade IIR implementation
clear all;
clear;
clc;
tic
Fs = 160; %Hz
L = 5;
lengthInput = 1000;
[NUMd, DENd] = cheby1(10,1,8/80);
NUMd = NUMd*db2mag(20);
[sos,b0]=tf2sos(NUMd,DENd);
[inputs, outputs] = myCascadedIIRFilter(L, sos, b0, Fs, lengthInput);

toc
figure(17)
plot (inputs)
hold on;
plot (outputs)
legend("input", "output")
title("Cascaded IIR Filter Implementation.")
ylabel("Signal Amplitude")
xlabel("Time [n]")
%saveas(figure(17), "IIR Filter performance.svg")

figure(18)
IN = fft(inputs);
OUT = fft(outputs);
digital_freqs = linspace(0, 2*pi, length(IN));
freqs = (digital_freqs*Fs)/2*pi;
plot((digital_freqs/(2*pi))*160, mag2db(abs(IN)))
hold on
plot((digital_freqs/(2*pi))*160,mag2db(abs(OUT)))
xlim([0 80])
ylabel("Magnitude [dB]")
xlabel("Frequency [Hz]")
title("Cascaded IIR Filter Implementation - Frequency Spectrum.")
%saveas(figure(18), "IIR FILTER FREQUENCY SPECTRUM.svg")
