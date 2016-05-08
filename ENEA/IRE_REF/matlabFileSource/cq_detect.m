function y = cq_detect(activate, ipl, dipl, minIPL, dipl_threshold, counterThreshold)
%% CurrentQuench and RE palteau detection
% INPUT:
%          activate = boolean var for eneabled the alghoritm
%                      0 : active
%                      1 : non active
%           ipl = absolute value of the  plasma current
%           dipl, minIPL, dipl_threshold, counterThreshold : min
%
% OUTPUT:  y in R
%         -1 : detection non active
%          0 : detection active
%          1 : CQ detect
%          2 : RE palteau detect
%%

persistent  startDetection;

if(activate)
    if(isempty(startDetection))
        startDetection = 0;
    end
    
    %STEP1: FIND CQ
    if( ipl >= minIPL)
        
        if( dipl <= dipl_threshold  && mod(startDetection,2)  == 0)
            startDetection = startDetection+1;
        end
        
        if( dipl >=  dipl_threshold  && mod(startDetection,2) == 1)
            startDetection = startDetection+1;
        end
        
    else
        
        if(startDetection > 0)
            startDetection =-1;
        end
        
    end
    
    y = startDetection;
    
else
    y = -2;
end
