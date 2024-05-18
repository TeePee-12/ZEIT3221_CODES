function [inputs,outputs] = myFIRFilter(numTaps,filterCoeffs,Fs,lengthInput)

% this function implements a finite impulse response discrete-time filter,  
% calling on another function, delayLineAndOutputFIR, to be written by you

% numTaps is number of filter taps (also equal to length of delay line)
% filterCoeffs is a vector of filter coefficients
% Fs is sampling frequency in Hz
% lengthInput is length of the input signal

inputs = zeros(lengthInput,1); % keep a record of filter inputs 
outputs = zeros(lengthInput,1); % keep a record of filter outputs, for length of input signal

x = zeros(numTaps,1); % initialize delay line

Ts = 1/Fs; % [s]

sigFs = 50000; % [kHz]
sigTs = 1/Fs; % [s]
% define time interval
t = 0.:sigTs:(lengthInput*sigTs); % [s]
% create signal
sig = sin(50*pi*t)+cos(500*pi*t)+cos(5000*pi*t)+cos(35000*pi*t);

for n=1:lengthInput 
       inputs(n) = sig(n);
       %inputs(n+1) = 3*n+1; % simple input for testing
       tic; % start a timer for subsequent iterations
       [outputs(n),x] = delayLineAndOutputFIR(numTaps,filterCoeffs,x,inputs(n)); % update delay line and filter output
       while toc<Ts
       % generate approximately real-time execution
       end
end



