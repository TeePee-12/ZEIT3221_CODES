
function [inputs,outputs]=myCascadedIIRFilter(L,sos,b0,Fs,lengthInput)

% this function implements an infinite impulse response discrete-time filter, 
% via a cascade of second-order sections, implemented by another function
% mySOS, to be written by you

% L is the number of second-order sections in the filter
% sos is the matrix of filter coefficients for each second-order section
% b0 is the z^0 coefficient in the numerator of the transfer function
% Fs is sampling frequency in Hz
% lengthInput is the length of the input signal

inputs = zeros(lengthInput,1); % keep a record of filter inputs (so we can look at data after execution)
outputs = zeros(lengthInput,1); % keep a record of filter outputs (so we can look at data after execution)

x = zeros(L,3); % initialize feedforward delay lines for each of L second-order sections
y = zeros(L,3); % initialize feedback delay lines for each of L second-order sections

Ts=1/Fs; % sampling period

%Create our input signal for testing
sigTs = 1/Fs; % [s]
% define time interval
t = 0.:sigTs:(lengthInput*sigTs); % [s]
% create signal
sig = sin(2*pi*t)+(0.25*cos(50*pi*t));

for n=0:lengthInput-1
       inputs(n+1) = sig(n+1);
        
       tic; % start a timer for subsequent iterations
       
       SOSinput=inputs(n+1);
       for r=1:L
            [SOSoutput,x,y] = mySOS(SOSinput,x,y,sos,L,r); % evaluate output of each second-order section
            SOSinput=SOSoutput;
       end
       
       outputs(n+1) = b0*SOSoutput; % multiple by b0 for final output
             
       %outputs(n+1) % print filter outputs (you may wish to comment this out)
       
       while toc<Ts
       % generate approximately real-time execution
       end   
end



