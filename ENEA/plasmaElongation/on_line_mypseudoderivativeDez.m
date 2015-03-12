function out = on_line_mypseudoderivativeDez(indata,c,d,Ts)
%x  the domain vector
%y  the image vector
%c  the length of the two subwindows where the mean is evaluated
%d  is total number of data considered to evaluate the pseudo-derivative, d>=2c

persistent mybuffer j

% (------d------)
% (-c-)-----(-c-)
% [---]-----[---]

if(isempty(mybuffer) )
    mybuffer = zeros(d,1);
    j = 1;
end

if(j < d+1)
    mybuffer(j) = indata;
    out = 0.0;
else
    %Buffer update
    for k=1:d-1
        mybuffer(k) = mybuffer(k+1);
    end
    mybuffer(d) = indata;
    
    %ora calcolo le medie
    temp1 = 0;
    for k=1:c
        temp1 = temp1 + mybuffer(k);
    end
    temp1 = temp1/c;
    
    temp2 = 0;
    for k=d-c+1:d
        temp2 = temp2 + mybuffer(k);
    end
    temp2 = temp2/c;
    
    out = (temp2-temp1) / (Ts*(d-c));
end
j = min(j +1,d+1);