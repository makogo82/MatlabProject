function Data = CleanFTUdataLOCAL_v4(t_i,t_f,Ts,shot,radiusNum,radiusIntPosition) 


globalVar; 

disp(strcat('ELABORAZIONE SHOT: ', num2str(shot)));

tempus = load([ '..\shot_' num2str(shot) '.mat']);

DATA = tempus.Data;

[valjj,jj] = max(abs(DATA.IPLprep.y));
CurrentSignIp = sign(DATA.IPLprep.y(jj));
Data.CurrentSignIp =CurrentSignIp;


%Ip
if(isstruct(DATA.IPLmis)==0)
    DATA.IPLmis = DATA.oldALFmis;
end
temp = DATA.IPLmis;

j1 = find(temp.x >= t_i,1);
j2 = find(temp.x >= t_f,1);       
new_time = [temp.x(j1):Ts:temp.x(j2)];

Data.time = new_time;
Data.IPLmis = interp1(temp.x,temp.y,new_time);


%Ip preprogrammed
temp = DATA.IPLprep;
Data.Ipprep = interp1(temp.x,temp.y,new_time);

%T coil
%We first leave the offset of specific signal on the measured values
if(isstruct(DATA.ALTmis)==0)
    DATA.ALTmis = DATA.oldALTmis;
end

temp = DATA.ALTmis;
tj1 = find(temp.x>=-1.4,1);
tj2 = find(temp.x>=-1.1,1);
temp.y = temp.y - mean(temp.y(tj1:tj2));
SignTemp = sign(temp.y(find(temp.x>=0,1)));
if(SignTemp ~= CurrentSignIp)
  temp.y = -(temp.y);
end  
Data.Tmis = interp1(temp.x,temp.y,new_time);

% coil T required current 
if(isstruct(DATA.ALTcalc)==0)
    DATA.ALTcalc = DATA.oldALTcalc;
end
temp = DATA.ALTcalc;
SignTemp = sign(temp.y(find(temp.x>=0,1)));
if(SignTemp ~= CurrentSignIp)
  temp.y = -(temp.y);
end    
Data.Tcalc = interp1(temp.x,temp.y,new_time);


temp = DATA.ALTprep;
SignTemp = sign(temp.y(find(temp.x>=0,1)));
if(SignTemp ~= CurrentSignIp)
  temp.y = -(temp.y);
end    
Data.Tprep = interp1(temp.x,temp.y,new_time);

%COIL V
if(isstruct(DATA.ALVmis)==0)
    DATA.ALVmis = DATA.oldALVmis;
end
temp = DATA.ALVmis;
tj1 = find(temp.x>=-1.4,1);
tj2 = find(temp.x>=-1.1,1);
temp.y = temp.y - mean(temp.y(tj1:tj2));
temp.y = abs(temp.y)*(-CurrentSignIp);  
Data.Vmis = interp1(temp.x,temp.y,new_time);

if(isstruct(DATA.ALVcalc)==0)
    DATA.ALVcalc = DATA.oldALVcalc;
end
temp = DATA.ALVcalc;
tj1 = find(temp.x>=-1.4,1);
tj2 = find(temp.x>=-1.1,1);
temp.y = temp.y - mean(temp.y(tj1:tj2));
temp.y = abs(temp.y)*(-CurrentSignIp);  
Data.Vcalc = interp1(temp.x,temp.y,new_time);

temp = DATA.ALVprep;
temp.y = abs(temp.y)*(-CurrentSignIp);  
Data.Vprep = interp1(temp.x,temp.y,new_time);


if(isstruct(DATA.ALFmis)==0)
    DATA.ALFmis = DATA.oldALFmis;
end
temp = DATA.ALFmis;
tj1 = find(temp.x>=-1.4,1);
tj2 = find(temp.x>=-1.1,1);
temp.y = temp.y - mean(temp.y(tj1:tj2));
Data.Fmis = interp1(temp.x,temp.y,new_time);

if(isstruct(DATA.ALFcalc)==0)
    DATA.ALFcalc = DATA.oldALFcalc;
end
temp = DATA.ALFcalc;
Data.Fcalc = interp1(temp.x,temp.y,new_time);

temp = DATA.ALFprep;
Data.Fprep = interp1(temp.x,temp.y,new_time);

% COIL H
if(isstruct(DATA.ALHmis)==0)
    DATA.ALHmis = DATA.oldALHmis;
end
temp = DATA.ALHmis;
tj1 = find(temp.x>=-1.4,1);
tj2 = find(temp.x>=-1.1,1);
temp.y = temp.y - mean(temp.y(tj1:tj2));
Data.Hmis = interp1(temp.x,temp.y,new_time);

if(isstruct(DATA.ALHcalc)==0)
    DATA.ALHcalc = DATA.oldALHcalc;
end
temp = DATA.ALHcalc;
Data.Hcalc = interp1(temp.x,temp.y,new_time);

temp = DATA.ALHprep;
Data.Hprep = interp1(temp.x,temp.y,new_time);


%temp = densita;
%density = interp1(temp.x,temp.y,new_time);
if(isstruct(DATA.DEP)==0)
    DATA.DEP = DATA.oldDEP;
end
temp = DATA.DEP;
Data.dep = interp1(temp.x,temp.y,new_time);

if(isstruct(DATA.DEZ)==0)
    DATA.DEZ = DATA.oldDEZ;
end
temp = DATA.DEZ;
Data.dez = interp1(temp.x,temp.y,new_time);


temp=0;
if(isstruct(DATA.RS2))
    temp = DATA.RS2;
end
if(isstruct(temp)==0)
    temp = DATA.oldRS2;
    if(isstruct(temp)==0)
        clear temp;
        temp.x = new_time;
        temp.y = Data.IPLmis*0;
    end
end
Data.rs2 = interp1(temp.x,max(0,min(1.5,temp.y)),new_time);


temp=0;
if(isstruct(DATA.RS1))
    temp = DATA.RS1;
end
if(isstruct(temp)==0)
    temp = DATA.oldRS1;
    if(isstruct(temp)==0)
        clear temp;
        temp.x = new_time;
        temp.y = Data.IPLmis*0;
    end
end
Data.rs1 = interp1(temp.x,max(-0,min(1.2,temp.y)),new_time);

temp=0;
if(isstruct(DATA.ZS2))
    temp = DATA.ZS2;
end
if(isstruct(temp)==0)
    temp = DATA.oldZS2;
    if(isstruct(temp)==0)
        clear temp;
        temp.x = new_time;
        temp.y = Data.IPLmis*0;
    end
end
Data.zs2 = interp1(temp.x,max(-0.5,min(0.5,temp.y)),new_time);

temp=0;
if(isstruct(DATA.ZS1))
    temp = DATA.ZS1;
end
if(isstruct(temp)==0)
    temp = DATA.oldZS1;
    if(isstruct(temp)==0)
        clear temp;
        temp.x = new_time;
        temp.y = Data.IPLmis*0;
    end
end
Data.zs1 = interp1(temp.x,max(-0.5,min(0.5,temp.y)),new_time);

temp = DATA.R_ext_prep;
Data.RS2prep = interp1(temp.x,max(-2,min(2,temp.y)),new_time);

temp = DATA.R_int_prep;
Data.RS1prep = interp1(temp.x,max(-2,min(2,temp.y)),new_time);

temp = DATA.Z_up_prep;
Data.ZS2prep = interp1(temp.x,max(-2,min(2,temp.y)),new_time);

temp = DATA.Z_down_prep;
Data.ZS1prep = interp1(temp.x,max(-2,min(2,temp.y)),new_time);


temp = DATA.marteferunaway;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.marteferunaway = interp1(temp.x,temp.y,new_time);


temp = DATA.runawayplatea;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.runawayplatea = interp1(temp.x,temp.y,new_time);

temp = DATA.AIFref;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.aifref = interp1(temp.x,temp.y,new_time);


temp = DATA.AdeltaIF;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.adeltaif = interp1(temp.x,temp.y,new_time);


temp = DATA.runipprep;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.runipprep = interp1(temp.x,temp.y,new_time);

temp = DATA.MHDFSTCH04;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.mhdfstc04 = temp.y;
Data.mhdfstc04_time = temp.x;


temp = DATA.TZXAMPN;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.tzxampn = interp1(temp.x,temp.y,new_time);

temp = DATA.Temperatura;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.temperatura = interp1(temp.x,temp.y,new_time);

temp = DATA.neutrnfc144;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.neutrnfc144 = interp1(temp.x,temp.y,new_time);


temp = DATA.disruption;
if(isstruct(temp)==0) 
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.disruption = interp1(temp.x,max(-2,min(2,temp.y)),new_time);


temp = DATA.oldVloop;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.vloop = interp1(temp.x,temp.y,new_time);


if(exist('DATA.martefe_xi'))
    temp = DATA.martefe_xi;
    if(isstruct(temp)==0)
        clear temp;
        temp.x = new_time;
        temp.y = Data.IPLmis*0;
    end
else
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.xi = interp1(temp.x,temp.y,new_time);


if(exist('DATA.martefe_xif'))
    temp = DATA.martefe_xif;
    if(isstruct(temp)==0)
        clear temp;
        temp.x = new_time;
        temp.y = Data.IPLmis*0;
    end
else
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.xif = interp1(temp.x,temp.y,new_time);


if(exist('DATA.martefe_alt_new'))
    temp = DATA.martefe_alt_new;
    if(isstruct(temp)==0)
        clear temp;
        temp.x = new_time;
        temp.y = Data.IPLmis*0;
    end
else
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.alt_new = interp1(temp.x,temp.y,new_time);


if(exist('DATA.martefe_neufc02'))
    temp = DATA.martefe_neufc02;
    if(isstruct(temp)==0)
        clear temp;
        temp.x = new_time;
        temp.y = Data.IPLmis*0;
    end
else
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.neufc02 = interp1(temp.x,temp.y,new_time);


if(exist('DATA.martefe_rerif'))
    temp = DATA.martefe_rerif;
    if(isstruct(temp)==0)
        clear temp;
        temp.x = new_time;
        temp.y = Data.IPLmis*0;
    end
else
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.rerif = interp1(temp.x,temp.y,new_time);

if(exist('DATA.martefe_rerif1'))
    temp = DATA.martefe_rerif1;
    if(isstruct(temp)==0)
        clear temp;
        temp.x = new_time;
        temp.y = Data.IPLmis*0;
    end
else
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.rerif1 = interp1(temp.x,temp.y,new_time);


if(exist('DATA.martefe_rerif2'))
    temp = DATA.martefe_rerif2;
    if(isstruct(temp)==0)
        clear temp;
        temp.x = new_time;
        temp.y = Data.IPLmis*0;
    end
else
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.rerif2 = interp1(temp.x,temp.y,new_time);


temp = DATA.FBKP_T;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.FBKP_T = interp1(temp.x,temp.y,new_time);

temp = DATA.FBKP_F;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.FBKP_F = interp1(temp.x,temp.y,new_time);

temp = DATA.FBKP_V;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.FBKP_V = interp1(temp.x,temp.y,new_time);

temp = DATA.FBKP_H;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.FBKP_H = interp1(temp.x,temp.y,new_time);

temp = DATA.FBKI_T;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.FBKI_T = interp1(temp.x,temp.y,new_time);

temp = DATA.FBKI_F;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.FBKI_F = interp1(temp.x,temp.y,new_time);

temp = DATA.FBKI_V;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.FBKI_V = interp1(temp.x,temp.y,new_time);

temp = DATA.FBKI_H;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.FBKI_H = interp1(temp.x,temp.y,new_time);

temp = DATA.FBKD_T;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.FBKD_T = interp1(temp.x,temp.y,new_time);

temp = DATA.FBKD_F;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.FBKD_F = interp1(temp.x,temp.y,new_time);

temp = DATA.FBKD_V;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.FBKD_V = interp1(temp.x,temp.y,new_time);

temp = DATA.FBKD_H;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.FBKD_H = interp1(temp.x,temp.y,new_time);

c_dez = 2;
d_dez = 9;
c_ddez = 2;
d_ddez = 10;

dez_d_error = mypseudo_derivative(Data.time,Data.dez,c_dez,d_dez);
dez_dd_error = mypseudo_derivative(Data.time,dez_d_error,c_ddez,d_ddez);
Data.dez_d_error = dez_d_error;
Data.dez_dd_error = dez_dd_error;
Data.mytime = Data.time ;
Data.Itu = Data.Tcalc;
Data.Ivu = Data.Vcalc;
Data.Ifu = Data.Fcalc;
Data.Ihu = Data.Hcalc;
Data.Itm = Data.Tmis;
Data.Ivm = Data.Vmis;
Data.Ifm = Data.Fmis;
Data.Ihm = Data.Hmis;
Data.Ip = Data.IPLmis;
Data.Tprep = Data.Tprep;
Data.Vprep = Data.Vprep;
Data.Fprep = Data.Fprep;
Data.Hprep = Data.Hprep;
Data.Ipprep = Data.Ipprep;
Data.myRS2prep = Data.RS2prep;
Data.myRS1prep = Data.RS1prep;
Data.myZS2prep = Data.ZS2prep;
Data.myZS1prep = Data.ZS1prep;
Data.vloop = Data.vloop;
Data.neutrnfc144=Data.neutrnfc144;

Data.inputIddData=[Data.IPLmis' Data.Hmis' Data.Vmis'];
Data.outputIddData=[Data.dez'];
zIddData = iddata(Data.outputIddData, Data.inputIddData , (Data.mytime(2)-Data.mytime(1)), 'Name', strcat(['Identificazione HCoil t_{plateau}=[',num2str(t_i),'-',num2str(t_f),']. Utilizzo segnale FC SHOT: ', num2str(shot)]));
zIddDataSim = iddata([], Data.inputIddData , (Data.mytime(2)-Data.mytime(1)), 'Name', 'H coil');
Data.zIddDataSim=zIddDataSim;
Data.zIddData=zIddData;


