function [output, x ,y ] = delayLinesAndOutputIIR(M,N,b,a,x,y,input)

%Set feedforward delay line index zero to new input, increment others by 1 index
%Set feedback delay line index zero to 0, increment others up
x = [input; x(1:end-1)];
y = [0; y(1:end-1)];

%First summation block of the IIR formula
bx = 0;
for i=1:M
    bx = bx + (b(i)*x(i));
end  

%Second summation block of the IIR formula
ay = 0;
for j=2:N
    ay = ay + (a(j)*y(j));
end 
%Difference of summations - CCDE
output = bx - ay;
y(1) = output;


end