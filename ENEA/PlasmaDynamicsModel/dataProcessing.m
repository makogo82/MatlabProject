function Data = dataProcessing(t_i,t_f,Ts,shot, n , wn) 
%  clc;clear all;
%  shot = 38570
%  t_i = -4;
%  t_f = 2;
%  Ts = 0.0005;
%  n=2;
%  wn=0.1;


disp(strcat('ELABORAZIONE SHOT: ', num2str(shot)));
%%
tempus = load([ 'shot_' num2str(shot) '.mat']);
data = tempus.data;
time = data.itime;
temp = data.ip;
%%
[blow,alow] = butter(n,wn,'low');
[bhigh,ahigh] = butter(n,wn,'high');

%new vector time
j1 = find(time >= t_i,1);
j2 = find(time >= t_f,1);       
new_time = [time(j1):Ts:time(j2)];
Data.time = new_time;
%% clc
Data.IPLmis = interp1(temp.x,temp.y,new_time);
Data.IPLmisLowFilter = interp1(temp.x,temp.rft_lf,new_time);
Data.IPLmisHighFilter = interp1(temp.x,temp.rft_hf,new_time);
%%
temp = data.if;
Data.IFmis = interp1(temp.x,temp.y,new_time);
Data.IFmisLowFilter =  interp1(time,temp.rft_lf,new_time);
Data.IFmisHighFilter = interp1(time,temp.rft_hf,new_time);

%%
temp = data.iv;
%%
Data.IVmis = interp1(temp.x,temp.y,new_time);
Data.IVmisLowFilter = interp1(time, temp.rft_lf, new_time);
Data.IVmisHighFilter = interp1(time, temp.rft_hf, new_time);

temp = data.it;
Data.ITmis = interp1(temp.x,temp.y,new_time);
Data.ITmisLowFilter =  interp1(time, temp.rft_lf, new_time);
Data.ITmisHighFilter =  interp1(time, temp.rft_hf, new_time);

temp = data.ih;
Data.IHmis = interp1(temp.x,temp.y,new_time);
Data.IHmisLowFilter =  interp1(time,temp.rft_lf,new_time);
Data.IHmisHighFilter =  interp1(time,temp.rft_hf,new_time);
 
data.vbclean.y
data.vsclean.y

for k = 1:1:size(data.vb_rft_out_lf,2)
    
    Data.vblow(k,:) = interp1(data.vb(k).x,  data.vb_rft_out_lf(k).y, new_time);
    Data.vbhigh(k,:) = interp1(data.vb(k).x,  data.vb_rft_out_hf(k).y, new_time);
    
    Data.vslow(k,:) = interp1(data.vs(k).x,  data.vs_rft_out_lf(k).y, new_time);
    Data.vshigh(k,:) = interp1(data.vb(k).x,  data.vb_rft_out_hf(k).y, new_time);
    
end

Data.VBclean = interp1( data.vb(1).x ,  data.vbclean(1).y, new_time);
Data.VSclean = interp1( data.vb(1).x ,  data.vsclean(1).y, new_time);

% N_txN_u
% manca il campo magnetico toroidale
Data.inputLow = [Data.IFmisLowFilter' Data.IVmisLowFilter' Data.ITmisLowFilter' Data.IHmisLowFilter'];
Data.inputHigh = [Data.IFmisHighFilter' Data.IVmisHighFilter' Data.ITmisHighFilter' Data.IHmisHighFilter'];
% N_txN_u
% TODO inserire Vloop come una possibile uscita
Data.outputLow = [Data.vblow' Data.vslow'];
Data.outputHigh = [Data.vbhigh' Data.vshigh'];

%iddata low
TsIddata = (Data.time(2)-Data.time(1));

Data.iddDataLow = iddata(Data.outputLow, Data.inputLow, TsIddata, 'Name', 'Input Low Frequency Magnetic Model');
Data.iddDataLowSim = iddata([], Data.inputLow , TsIddata, 'Name', 'Output Low Frequency Magnetic Model');

%iddata high
Data.iddDataHigh = iddata(Data.outputHigh, Data.inputHigh, TsIddata, 'Name', 'Input Low Frequency Magnetic Model');
Data.iddDataHighSim = iddata([], Data.inputHigh , TsIddata, 'Name', 'Output Low Frequency Magnetic Model');

