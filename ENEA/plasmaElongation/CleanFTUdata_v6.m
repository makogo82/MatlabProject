function Data = CleanFTUdata_v6(t_i, t_f, Ts, shot, local) 

%% t_i : initial time
%  t_f : final time
% Ts : sampling time
% shot : shot number
% local : one if the local files have to be loaded

if (local == 0)
    DATA = GetShotData_v1(shot);
else
    tempus = load(['shot_' num2str(shot) '.mat']);
    DATA = tempus.Data;
end

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
if(isempty(j2))
    disp('End data time too large for the considered shot.')
end
new_time = [temp.x(j1):Ts:temp.x(j2)];

Data.time = new_time;
Data.IPLmis = interp1(temp.x,temp.y,new_time);
Data.names.IPLmis = 'Plasma current';

%Ip preprogrammed
temp = DATA.IPLprep;
Data.Ipprep = interp1(temp.x,temp.y,new_time);
Data.names.Ipprep = 'Desired Plasma current';

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
Data.names.Tmis = 'T coil current measured';

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
Data.names.Tcalc = 'T coil current desired';

temp = DATA.ALTprep;
SignTemp = sign(temp.y(find(temp.x>=0,1)));
if(SignTemp ~= CurrentSignIp)
  temp.y = -(temp.y);
end    
Data.Tprep = interp1(temp.x,temp.y,new_time);
Data.names.Tprep = 'T coil current preprogrammed';

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
Data.names.Vmis = 'V coil current measured';

if(isstruct(DATA.ALVcalc)==0)
    DATA.ALVcalc = DATA.oldALVcalc;
end
temp = DATA.ALVcalc;
tj1 = find(temp.x>=-1.4,1);
tj2 = find(temp.x>=-1.1,1);
temp.y = temp.y - mean(temp.y(tj1:tj2));
temp.y = abs(temp.y)*(-CurrentSignIp);  
Data.Vcalc = interp1(temp.x,temp.y,new_time);
Data.names.Vcalc = 'V coil current desired';

temp = DATA.ALVprep;
temp.y = abs(temp.y)*(-CurrentSignIp);  
Data.Vprep = interp1(temp.x,temp.y,new_time);
Data.names.Vprep = 'V coil current preprogrammed';

if(isstruct(DATA.ALFmis)==0)
    DATA.ALFmis = DATA.oldALFmis;
end
temp = DATA.ALFmis;
tj1 = find(temp.x>=-1.4,1);
tj2 = find(temp.x>=-1.1,1);
temp.y = temp.y - mean(temp.y(tj1:tj2));
Data.Fmis = interp1(temp.x,temp.y,new_time);
Data.names.Fmis = 'F coil current measured';

if(isstruct(DATA.ALFcalc)==0)
    DATA.ALFcalc = DATA.oldALFcalc;
end
temp = DATA.ALFcalc;
Data.Fcalc = interp1(temp.x,temp.y,new_time);
Data.names.Fcalc = 'F coil current desired';

temp = DATA.ALFprep;
Data.Fprep = interp1(temp.x,temp.y,new_time);
Data.names.Fprep = 'F coil current preprogrammed';


% COIL H
if(isstruct(DATA.ALHmis)==0)
    DATA.ALHmis = DATA.oldALHmis;
end
temp = DATA.ALHmis;
tj1 = find(temp.x>=-1.4,1);
tj2 = find(temp.x>=-1.1,1);
temp.y = temp.y - mean(temp.y(tj1:tj2));
Data.Hmis = interp1(temp.x,temp.y,new_time);
Data.names.Hmis = 'H coil current measured';

if(isstruct(DATA.ALHcalc)==0)
    DATA.ALHcalc = DATA.oldALHcalc;
end
temp = DATA.ALHcalc;
Data.Hcalc = interp1(temp.x,temp.y,new_time);
Data.names.Hcalc = 'H coil current desired';

temp = DATA.ALHprep;
Data.Hprep = interp1(temp.x,temp.y,new_time);
Data.names.Hprep = 'H coil current preprogrammed';

temp = DATA.esidens;
if(isstruct(temp)==0)  
    if(isstruct(DATA.densita)==0)
        temp.x = DATA.IPLmis.x;
        temp.y = 0*DATA.IPLmis.y;
    end
end
Data.densita = interp1(temp.x,temp.y,new_time);
Data.names.densita = 'Central plasma density';

if(isstruct(DATA.DEP)==0)
    DATA.DEP = DATA.oldDEP;
end
temp = DATA.DEP;
Data.dep = interp1(temp.x,temp.y,new_time);
Data.names.dep = 'DEP, flux horizontal position error';

if(isstruct(DATA.DEZ)==0)
    DATA.DEZ = DATA.oldDEZ;
end
temp = DATA.DEZ;
Data.dez = interp1(temp.x,temp.y,new_time);
Data.names.dez = 'DEZ, flux vertical position error';

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
Data.names.rs2 = 'External horizontal plasma radius';

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
Data.names.rs1 = 'Internal horizontal plasma radius';

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
Data.names.zs2 = 'Lower vertical plasma radius';

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
Data.names.zs1 = 'Upper vertical plasma radius';

temp = DATA.R_ext_prep;
Data.RS2prep = interp1(temp.x,max(-2,min(2,temp.y)),new_time);
Data.names.RS2prep = 'Desired external horizontal plasma radius';

temp = DATA.R_int_prep;
Data.RS1prep = interp1(temp.x,max(-2,min(2,temp.y)),new_time);
Data.names.RS1prep = 'Desired internal horizontal plasma radius';

temp = DATA.Z_up_prep;
Data.ZS2prep = interp1(temp.x,max(-2,min(2,temp.y)),new_time);
Data.names.ZS2prep = 'Desired lower vertical plasma radius';

temp = DATA.Z_down_prep;
Data.ZS1prep = interp1(temp.x,max(-2,min(2,temp.y)),new_time);
Data.names.ZS1prep = 'Desired upper vertical plasma radius';

temp = DATA.marteferunaway;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.marteferunaway = interp1(temp.x,temp.y,new_time);
Data.names.marteferunaway = 'Runaway on';

temp = DATA.runawayplatea;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.runawayplatea = interp1(temp.x,temp.y,new_time);
Data.names.runawayplatea = 'Runaway plateau detected';

temp = DATA.AIFref;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.aifref = interp1(temp.x,temp.y,new_time);
Data.names.aifref = 'Desired F current during runaway control';

temp = DATA.AdeltaIF;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.adeltaif = interp1(temp.x,temp.y,new_time);
Data.names.adeltaif = 'Variation F current during runaway control';


temp = DATA.runipprep;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.runipprep = interp1(temp.x,temp.y,new_time);
Data.names.runipprep = 'Desired plasma current during runaway control';

temp = DATA.MHDFSTCH04;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.mhdfstc04 = temp.y;
Data.mhdfstc04_time = temp.x;
Data.names.mhdfstc04 = 'MHD signal from mirnov coil N.4';

temp = DATA.TZXAMPN;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.tzxampn = interp1(temp.x,temp.y,new_time);
Data.names.tzxampn = 'Estimated MHD amplitude';

temp = DATA.Temperatura;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
if(length(temp.x)<=1)
       temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.temperatura = interp1(temp.x,temp.y,new_time);
Data.names.temperatura = 'Plasma temperature';

temp = DATA.neutrnfc144;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.neutrnfc144 = interp1(temp.x,temp.y,new_time);
Data.names.neutrnfc144 = 'Fission Chamber signal, Hard-X >= 6MeV, neutrn.fc144_.statn2, counts each 1 ms';

temp = DATA.disruption;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.disruption = interp1(temp.x,max(-2,min(2,temp.y)),new_time);
Data.names.disruption = 'Disruption detected';

%temp = DATA.Vloop;
%temp = DATA.e_Vloop;
%temp = DATA.vl_Vloop;
temp = DATA.oldVloop;
if(isstruct(temp)==0)
    temp = oldVloop
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
else
    if(max(abs(temp.y))<0.1)
       disp('Vloop odd...');
    end        
end
Data.vloop = interp1(temp.x,temp.y,new_time);
Data.names.vloop = 'Loop voltage';

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
Data.names.xi = 'martefe xi variable';

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
Data.names.xif = 'martefe xif variable';


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
Data.names.alt_new = 'alt_new variable';

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
Data.names.neufc02 = 'Real-time FC signal Martefe.neufc02, >= 6MeV, counts each 0.5ms';

temp = DATA.rerif;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.rerif = interp1(temp.x,temp.y,new_time);
Data.names.rerif = 'Desired RS2';

temp = DATA.rerif1;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.rerif1 = interp1(temp.x,temp.y,new_time);
Data.names.rerif1 = 'Desired1 RS2 when RE control active ';

temp = DATA.rerif2;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.rerif2 = interp1(temp.x,temp.y,new_time);
Data.names.rerif2 = 'Desired2 RS2 when RE control active ';

temp = DATA.zurif;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.zurif = interp1(temp.x,temp.y,new_time);
Data.names.zurif = 'Desired zs1';


temp = DATA.zlrif;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.zlrif = interp1(temp.x,temp.y,new_time);
Data.names.zlrif = 'Desired zs2';


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


temp = DATA.avalvecmd;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.avalvecmd = interp1(temp.x,temp.y,new_time);
Data.names.avalvecmd = 'Command signal of deuterium valve';

temp = DATA.FBINA_M_RIS;
if(isstruct(temp)==0)
        temp = DATA.old_FBINA_M_RIS;
        if(isstruct(temp)==0)
            temp = DATA.old_FBINA_M_RIS;
            clear temp;
            temp.x = new_time;
            temp.y = Data.IPLmis*0;
        end
end
Data.fbina_m_ris = interp1(temp.x,temp.y,new_time);
Data.names.fbina_m_ris = 'Hard-X signal over 200KeV';

temp = DATA.NEU213;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.neu213 = interp1(temp.x,temp.y,new_time);
Data.names.neu213 = 'Scitillator sensible to neutron and gamma rays';

temp = DATA.BF3;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.bf3 = interp1(temp.x,temp.y,new_time);
Data.names.bf3 = 'Neutron detector, not sensible to gamma rays';

temp = DATA.elong;
if(isstruct(temp)==0)
        temp = DATA.old_elong;
        if(isstruct(temp)==0)
            clear temp;
            temp.x = new_time;
            temp.y = Data.IPLmis*0;
        end
end
Data.elong = interp1(temp.x,temp.y,new_time);
Data.names.elong = 'Plasma elongation';

temp = DATA.esoftx;
if(isstruct(temp)==0)
    clear temp;
    temp.x = new_time;
    temp.y = Data.IPLmis*0;
end
Data.esoftx = interp1(temp.x,temp.y,new_time);
Data.names.esoftx = 'Soft-X sensor  up to 10keV';


