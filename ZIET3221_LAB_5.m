
xrealization = 2*randn(200000, 1) + 3;
figure(1)

subplot(2, 2, 1);
stem(xcorr(xrealization,'biased'))
title("Standard Biased Autocorrelation Estimate")
xlabel('Lag');
ylabel('Autocorrelation');
subplot(2, 2, 2);
stem(xcorr(xrealization,'unbiased'))
title("Standard Unbiased Autocorrelation Estimate")
xlabel('Lag');
ylabel('Autocorrelation');
subplot(2, 2, 3);
stem(xcov(xrealization,'biased'))
title("Standard Biased Autocovariance Estimate")
xlabel('Lag');
ylabel('Autocovariance');
subplot(2, 2, 4);
stem(xcov(xrealization,'unbiased'))
title("Standard Unbiased Autocovariance Estimate")
xlabel('Lag');
ylabel('Autocovariance');


%% part 2

white_noise = sqrt(2) .* randn(100000, 1);
RxxEstimate = xcorr(white_noise, 'biased');

figure(2);
plot(RxxEstimate);
title('2a) Autocorrelation of Gaussian White Noise');
xlabel('Lag');
ylabel('Autocorrelation');

Sc = abs(fft(RxxEstimate));
freqs = linspace(-0.5, 0.5, length(Sc));

% Plot
figure(3);
plot(freqs, fftshift(Sc));
title('2b) Correlogram');
xlabel('Digital Frequency');
ylabel('Magnitude');

RxxEstimate(abs(RxxEstimate) < 0.05) = 0;
Sc = abs(fft(RxxEstimate));
figure(4);
plot(freqs, fftshift(Sc));
title('2b) Correlogram (Spurious small estimates removed)');
xlabel('Digital Frequency');
ylabel('Magnitude');

%% 2C
h = [3 1 3];
xrealisation = sqrt(2) .* randn(100000, 1);
yrealisation = filter(h, 1, xrealisation);


Ryx_biased = xcorr(yrealisation,xrealisation,'biased',10)
Ryx_unbiased = xcorr(yrealisation,xrealisation,'unbiased',10)


% Plotting the biased and unbiased estimates
figure;
subplot(2,1,1);
stem(-10:10, Ryx_biased);
title('Biased Cross-Correlation Estimate');
xlabel('Lag n'); ylabel('ˆRyx[n]');

subplot(2,1,2);
stem(-10:10, Ryx_unbiased);
title('Unbiased Cross-Correlation Estimate');
xlabel('Lag n'); ylabel('ˆRyx[n]')

%% 2d 

Syx = abs(fft(Ryx_unbiased));

figure;
subplot(2,1,1);
plot(linspace(0, pi, length(Syx)), Syx);
title('Cross-Spectral Density Estimate');
xlabel('Digital Frequency (radians)'); ylabel('ˆS_y_x(e^jω)');
xlim([0 pi/2])

Syx(Syx < max(Syx)*0.05) = 0;
subplot(2,1,2);
plot(linspace(0, pi, length(Syx)), Syx);
title('Cross-Spectral Density Estimate (Spurious values removed)');
xlabel('Digital Frequency (radians)'); ylabel('ˆS_y _x(e^jω)');
xlim([0 pi/2])

%% 2e
%PSD of y[n]
%yrealisation(abs(yrealisation) < 0.05) = 0;
figure;
subplot(2,1,1)
plot(xcorr(yrealisation,'biased'))
title("Biased Autocorrelation Estimate of Filter Output y[n]")
xlabel('Lag');
ylabel('Autocorrelation');
Pyy = abs(fft(yrealisation)).^2 / length(yrealisation);

subplot(2,1,2)
plot(linspace(0, pi, length(Pyy)), 10*log10(Pyy));
title('Power Spectral Density Estimate Filter Output y[n]');
xlabel('Digital Frequency (radians)');
ylabel('Power/Frequency (dB/Hz)');
xlim([0 pi])

%% 3a

xr = sqrt(16)*randn(1, 16384);
% Calculate the periodogram
Spk = fft(xr) .* conj(fft(xr)) / 16384;

% Plot the periodogram
figure;
plot(linspace(0, pi, 16384/2), Spk(1:16384/2));
title('Periodogram of Gaussian White Noise');
xlabel('Digital Frequency (radians)');
ylabel('Magnitude');
xlim([0 pi])

% Calculate the average value of the periodogram
avg_periodogram = mean(Spk(1:16384/2));

disp(['Average Periodogram Value: ', num2str(avg_periodogram)]);
disp(['Analytical PSD: ', num2str(sigma2)]);

%% 3b
[Spk_builtin, w] = periodogram(xr);

figure;
plot(w, Spk_builtin);
title('Periodogram using MATLAB function');
xlabel('Digital Frequency (radians)');
ylabel('Magnitude');
xlim([0 pi])

avg_periodogram_builtin = mean(Spk_builtin)

%% 3c
periodogram(xr, rectwin(length(xr)), length(xr), 16000);
title('Periodogram in dB against Real Frequency');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');

%% 3d
[Sw, w] = pwelch(xr);
figure;
plot(w, Sw);
title('Spectral Estimate using Welch’s method');
xlabel('Digital Frequency (radians)');
ylabel('Magnitude');
xlim([0 pi])

avg_pwelch = mean(Sw)
