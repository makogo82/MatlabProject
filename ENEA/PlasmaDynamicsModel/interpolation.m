function output = interpolation(temp,t_i,t_f,Ts)

j1 = find(temp.x >= t_i,1);
j2 = find(temp.x >= t_f,1);   

new_time = [temp.x(j1):Ts:temp.x(j2)];

output.time = new_time;
output.data = interp1(temp.x,temp.y,new_time);


end

