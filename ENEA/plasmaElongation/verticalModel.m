%% DOCUMENT TITLE
% INTRODUCTORY TEXT
%%

function [t ,t1 ,t2 ] = verticalModel(time , x1 , IPLmis,Hmis,Fmis, c_2, c_3 , c_4 , c_5  )
 

t1 =   c_2*IPLmis*Hmis*(1 /(c_3-1.5*(x1))    - 1/(-c_3-1.5*(x1)) );

t2 =  -c_4*Fmis*IPLmis*(2*( 1/(c_5-1.5*(x1)) + 1/(-c_5-1.5*(x1))))^2;

t = t1 +t2;

%t = - time * sin(time) + 2* cos(time);

end

