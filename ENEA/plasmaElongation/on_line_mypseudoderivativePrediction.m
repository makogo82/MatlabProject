function out = on_line_mypseudoderivativePrediction(indata,c,d,Ts)
%x  the domain vector
%y  the image vector
%c  the length of the two subwindows where the mean is evaluated
%d  is total number of data considered to evaluate the pseudo-derivative, d>=2c

persistent mybuffer2 j

% (------d------)
% (-c-)-----(-c-)
% [---]-----[---]


if(isempty(mybuffer2) )
    mybuffer2 = zeros(d,1);
    j = 1;
end

if(j < d+1)
    mybuffer2(j) = indata;
    out = 0.0;
else
    %Buffer update
    for k=1:d-1
        mybuffer2(k) = mybuffer2(k+1);
    end
    mybuffer2(d) = indata;
    
    %ora calcolo le medie
    temp1 = 0;
    for k=1:c
        temp1 = temp1 + mybuffer2(k);
    end
    temp1 = temp1/c;
    
    temp2 = 0;
    for k=d-c+1:d
        temp2 = temp2 + mybuffer2(k);
    end
    temp2 = temp2/c;
    
    out = (temp2-temp1) / (Ts*(d-c));
end
j = min(j +1,d+1);