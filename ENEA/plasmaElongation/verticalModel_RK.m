%% DOCUMENT TITLE
% INTRODUCTORY TEXT
%%

function [x1New , x2New ] = verticalModel_RK(time, Ts, x1 , x2 , IPLmis,Hmis,Fmis, c_2, c_3 , c_4 , c_5  )
 
    f  = verticalModel(time,x1,IPLmis,Hmis,Fmis, c_2 , c_3 , c_4 , c_5);
    k1 = Ts * f;
    m1 = Ts * x2;
    f = verticalModel(time+Ts/2,x1+m1/2,IPLmis,Hmis,Fmis, c_2 , c_3 , c_4 , c_5);
    k2 = Ts * f;
    m2=  Ts * (x2+k2/2);
    f = verticalModel(time+Ts/2,x1+m2/2,IPLmis,Hmis,Fmis, c_2 , c_3 , c_4 , c_5);
    k3 = Ts * f;
    m3=  Ts * (x2+k3/2);
    f = verticalModel(time+Ts,x1+m3,IPLmis,Hmis,Fmis, c_2 , c_3 , c_4 , c_5);
    k4 = Ts * f;
    m4 = Ts * (x2+k4); 

    x2New = x2 + ( k1 + 2*k2 + 2*k3 + k4)/6;
    x1New = x1 + ( m1 + 2*m2 + 2*m3 + m4)/6;

    %x2New = x2 + Ts*f;
    %x1New = x1 + Ts*x2; 
end

