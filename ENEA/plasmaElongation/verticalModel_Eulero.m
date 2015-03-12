%% DOCUMENT TITLE
% INTRODUCTORY TEXT
%%

function [x1New , x2New] = verticalModel_Eulero(time, Ts, x1 , x2 , IPLmis,Hmis,Fmis, c_2, c_3 , c_4 , c_5  )
    f  = verticalModel(time,x1,IPLmis,Hmis,Fmis, c_2 , c_3 , c_4 , c_5);
    x2New = x2 + Ts*f;
    x1New = x1 + Ts*x2; 
end

