close all; 

yv = lsim(tf([1],[0.03 1]),vloopN/((2*pi*0.86)),newTime,0.2)
plot(newTime,yv,newTime,EllN,'r');
