function [SOSoutput,xnew,ynew]=mySOS(input,x,y,sos,L,r)

    % r is the number of the second-order section in the cascade IIR filter
    % sos is matrix of filter coefficients for all second-order sections
    % x is matrix of samples in feedforward delay lines
    % y is matrix of samples in feedback delay lines
    % input is the input to the r^(th) second-order section

    % assign old delay line values to new delay line values
    % shift delay line for r^(th) row of x
    % new input value for r^(th) row of x
    % shift delay line for r^(th) row of y
    xnew = x;
    xr = xnew(r,:);
    xshift = zeros(1,length(xr));
    for i = length(xr):-1:2
        xshift(i) = xr(i-1);
    end
    xnew(r,:) = xshift;
    xnew(r,1) = input;
   
    ynew = y;
    yr = ynew(r,:);
    yshift = zeros(1,length(yr));
    for i = length(yr):-1:2
        yshift(i) = yr(i-1);
    end
    ynew(r,:) = yshift;
    
    % calculate second-order section output 
    SOSoutput = 0;
    for i = 1:3
        SOSoutput = SOSoutput + sos(r,i) * xnew(r,i);
        SOSoutput = SOSoutput - sos(r,i+3) * ynew(r,i);
    end
    
    % new output value for r^(th) row of y
    ynew(r,1) = SOSoutput;
end
