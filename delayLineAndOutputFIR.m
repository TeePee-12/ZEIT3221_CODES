function [output, x] = delayLineAndOutputFIR(numTaps,filterCoeffs,x,input)

%Initialize outputs to zero
output = 0;

%Set delay line index zero to new input, increment others by 1 index
x = [input; x(1:end-1)];

%Calculate filter outputs
for i=1:numTaps
    output = (filterCoeffs(i)*x(i))+output;
end



