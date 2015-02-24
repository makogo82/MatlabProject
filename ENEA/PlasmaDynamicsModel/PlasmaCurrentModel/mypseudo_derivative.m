function out  = mypseudo_derivative(x,y,c,d)

%x  the domain vector
%y  the image vector
%c  the length of the two subwindows where the mean is evaluated
%d  is total number of data considered to evaluate the pseudo-derivative, d>=2c

% (------d------)
% (-c-)-----(-c-)
% [---]-----[---] 


if(c>= 1 && d>=2*c && d>=2 && length(x)==length(y))
   
    
    T = mean(x(2:end)-x(1:end-1));
    temp_y = zeros(length(y),1);
    for j=d:length(y),
        temp_y(j) = (mean(y(j-c+1:j)) - mean(y(j-d+1:j-d+c)))/(T*(d-c));
    end
    for j=1:d-1,
        temp_y(j) = temp_y(d);
    end
    out = temp_y;
else
    disp('ERROR: c<1 or d < max(2*c,2) or  length(x) != length(y)')
    out = y*0;
end


%% Axample: how to use it
% 
% c = 10; %try also c=1
% d = 20; %d = 2 and see the difference
% 
% tempx = [0:0.001:1];
% tempy = 3*tempx + 2*sin(4*tempx) + 0.1*sin(2*pi*50*tempx) + 0.05*sin(2*pi*150*tempx);
% 
% tempy_derivative = mypseudo_derivative(tempx,tempy,c,d);
% 
% figure(1)
% subplot(2,1,1)
% plot(tempx,tempy_derivative)
% subplot(2,1,2)
% plot(tempx,tempy)