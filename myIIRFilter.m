
function [inputs,outputs]=myIIRFilter(M,N,b,a,Fs,lengthInput)

% this function implements an infinite impulse response discrete-time filter, 
% calling on another function delayLinesAndOutputIIR, to be written by you

% lengthInput is the length of the input signal
% Fs is sampling frequency in Hz
% M is the number of feedforward coefficients 
% N is the number of feedback coefficients
% b is the vector of feedforward coefficients
% a is the vector of feedback coefficients, including 1 as the first entry

inputs = zeros(lengthInput,1); % keep a record of filter inputs (so we can look at data after execution)
outputs = zeros(lengthInput,1); % keep a record of filter outputs (so we can look at data after execution)

x = zeros(M,1); % initialize feedforward delay line
y = zeros(N,1); % initialize feedback delay line

Ts=1/Fs; % sampling period

sigFs = 50000; % [kHz]
sigTs = 1/Fs; % [s]
% define time interval
t = 0.:sigTs:(lengthInput*sigTs); % [s]
% create signal
sig = sin(5*pi*t)+cos(10*pi*t)+cos(50*pi*t)+cos(500*pi*t);

for n=0:lengthInput-1
       inputs(n+1) = sig(n+1);
       tic; % start a timer for subsequent iterations
       [outputs(n+1),x,y] = delayLinesAndOutputIIR(M,N,b,a,x,y,inputs(n+1)); % update delay line and filter output
       while toc<Ts
       % generate approximately real-time execution
       end   
end



