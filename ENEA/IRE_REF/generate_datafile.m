function generate_datafile_1(shots)

delete data_actual.dat;

signals = {
    'MARTEFE.POTPER' ...
    'MARTEFE.DCN_F' ...
    'MARTEFE.IPL' ...
    'MARTEFE.VM' ...
    'MARTEFE.VL' ...
    'MARTEFE.VS01' ...
    'MARTEFE.VS02' ...
    'MARTEFE.VS03' ...
    'MARTEFE.VS04' ...
    'MARTEFE.VS05' ...
    'MARTEFE.VS06' ...
    'MARTEFE.VS07' ...
    'MARTEFE.VS08' ...
    'MARTEFE.VS09' ...
    'MARTEFE.VS10' ...
    'MARTEFE.VS11' ...
    'MARTEFE.VS12' ...
    'MARTEFE.VS13' ...
    'MARTEFE.VS14' ...
    'MARTEFE.VS15' ...
    'MARTEFE.VS16' ...
    'MARTEFE.VB01' ...
    'MARTEFE.VB02' ...
    'MARTEFE.VB03' ...
    'MARTEFE.VB04' ...
    'MARTEFE.VB05' ...
    'MARTEFE.VB06' ...
    'MARTEFE.VB07' ...
    'MARTEFE.VB08' ...
    'MARTEFE.VB09' ...
    'MARTEFE.VB10' ...
    'MARTEFE.VB11' ...
    'MARTEFE.VB12' ...
    'MARTEFE.VB13' ...
    'MARTEFE.VB14' ...
    'MARTEFE.VB15' ...
    'MARTEFE.VB16' ...
    '0_MACC.FVCOR-ALF.A' ...
    '0_MACC.FVCOR-ALH.A' ...
    '0_MACC.FVCOR-ALT.A' ...
    '0_MACC.FVCOR-ALV.A' ...
    '0_MACC.FBCOR-IP.A' ...
    '0_MACC.FBTEN-VP.A' ...
    '0_MACC.FBADM-BP.A' ...
    'MARTEFE.FBINA_M_RIS' ...
    'MARTEFE.ZFM=I' ...
    'MARTEFE.VALVECMD' ...
    'MARTEFE.FVCOR-ALV' ...
    'MARTEFE.FBIP_M_I' ...
    'MARTEFE.FBINA_M_RIS' ...
    '%E.VLOOP10' ... %martefe.vloop
    'martefe.neufc02' ... %martefe.neufc02
    'MARTEFE.FBMEA_ALX' ...
    'MARTEFE.FBMEA_ALF' ...
    'MARTEFE.FBMEA_ALV' ...
    'MARTEFE.AVALVECMD' ...
    'MARTEFE.AVALVEMEAS' ...
    'MARTEFE.AVALVEIMP' ...
    'MARTEFE.IMPURITY_REF' ...
    'MARTEFE.NOISEPOLM5' ...
    'MARTEFE.AVALPOLM5' ...
    'MARTEFE.APLASMA_DENS' ...
    'MARTEFE.RIRIF' ...
    'MARTEFE.RERIF' ...
    'MARTEFE.ZLRIF' ...
    'MARTEFE.ZURIF' ...
    'MARTEFE.ZTM=I'
    };

varnames = {
    'POTPER' ...
    'DCN' ...
    'IPLA' ...
    'VM' ...
    'VL' ...
    'VS01' ...
    'VS02' ...
    'VS03' ...
    'VS04' ...
    'VS05' ...
    'VS06' ...
    'VS07' ...
    'VS08' ...
    'VS09' ...
    'VS10' ...
    'VS11' ...
    'VS12' ...
    'VS13' ...
    'VS14' ...
    'VS15' ...
    'VS16' ...
    'VB01' ...
    'VB02' ...
    'VB03' ...
    'VB04' ...
    'VB05' ...
    'VB06' ...
    'VB07' ...
    'VB08' ...
    'VB09' ...
    'VB10' ...
    'VB11' ...
    'VB12' ...
    'VB13' ...
    'VB14' ...
    'VB15' ...
    'VB16' ...
    'MFVCOR_ALF' ...
    'MFVCOR_ALH' ...
    'MFVCOR_ALT' ...
    'MFVCOR_ALV' ...
    'MFBCOR_IP' ...
    'MFBTEN_VP' ...
    'MFBADM_BP' ...
    'FBINA_M_RIS' ...
    'FVCOR_ALH' ...
    'FVCOR_ALT' ...
    'FVCOR_ALV' ...
    'FBCOR_IP' ...
    'FBINA_M_RIS' ...
    'VLOOP' ...
    'NEUFC02' ...
    'FBMEA_ALX' ...
    'FBMEA_ALF' ...
    'FBMEA_ALV' ...
    'AVALVECMD' ...
    'AVALVEMEAS' ...
    'AVALVEIMP' ...
    'AVALPOLM3' ...
    'AVALPOLM2F' ...
    'AVALPOLM5' ...
    'APLASMA_DENS' ...
    'RIRIF' ...
    'RERIF' ...
    'ZLRIF' ...
    'ZURIF' ...
    'ZTM_I'
    };




for j = 1:length(shots)
    %fid = fopen(['data_' num2str(shots(j)) '.dat'], 'w');
    fid = fopen(['data_actual.dat'], 'w');
    numofwaves(1)=length(signals);
    fwrite(fid, numofwaves, 'int32');
    numofsamples(1)= 30000;
    disp(['processing shot ' num2str(shots(j))]);
   
    for i = 1:length(signals)
        disp(['i:' num2str(i) ' SHOT: ' num2str(shots(j)) '--> CHANNEL: '  signals{i} ''])
        eval([varnames{i} ' = ftudata(shots(j), signals{i});']);
        
        if(eval(['isstruct( ' varnames{i} ')==0']))
            disp(['Channel ' signals{i} ' not found use zero vector!'])
            numofsamples(1)= 30000;
            time=[-9.9995:.0005:5.];
            fwrite(fid,numofsamples,'int32');
            fwrite(fid,time,'float');
            fwrite(fid,zeros(size(time)),'float');
            
        else
            
            numofsamples(1)= eval(['length(' varnames{i} '.x)']);
            fwrite(fid,numofsamples,'int32');
            fwrite(fid,eval([varnames{i} '.x']),'float');
            fwrite(fid,eval([varnames{i} '.y']),'float');
        end
    end
    
    fclose(fid);
    disp(['create Zero File for shot ' num2str(shots(j))]);
    searchStart=shots(j)-1
    searchEnd=searchStart-100
    current_discharge=ftudata(shots(j),'000001');
    toroidalField=current_discharge.yl(1:3)
    toroidalField(1)='Z'; %force to look zerop
    
    if(toroidalField=='ZAC')
        toroidalRamp=current_discharge.yl(5:8)
    end
    
    for k=searchStart:-1:searchEnd
        clear pickupVB
        clear pickupVS
        
        discard=0;
        discharge=ftudata(k,'000001');
        if(isfield(discharge,'yl')==1)
            if(toroidalField==discharge.yl(1:3))
                %if ZAC..then we have to chk toroidalRamp
                disp([num2str(k) ' is a possible Zero shot for shot ' num2str(shots(j))]);
                discharge.yl
                zero=k;
                
                for l=1:16
                    if(l<10)
                        chbmarte=sprintf('MARTEFE.VB0%dF1',l);
                        chsmarte=sprintf('MARTEFE.VS0%dF1',l);
                    else
                        chbmarte=sprintf('MARTEFE.VB%dF1',l);
                        chsmarte=sprintf('MARTEFE.VS%dF1',l);
                    end
                    pickupVB(l)=ftudata(k,chbmarte);
                    pickupVS(l)=ftudata(k,chsmarte);
                    if(isstruct(pickupVB(l))==0 || isstruct(pickupVS(l))==0)
                        disp(['discard shot ' num2str(k) ' no  MARTEFE data  for' chsmarte ' or '  chbmarte]);
                        discard=1;
                        break;
                    end
                end
                if(discard==0)
                    %fid = fopen(['zero_' num2str(shots(j)) '.dat'], 'w');
                    fid = fopen(['zero_actual.dat'], 'w');
                    numofwaves(1)=32
                    fwrite(fid, numofwaves, 'int32');
                    numofsamples(1)= length(pickupVS(1).x);
                    fwrite(fid, numofsamples, 'int32');
                    
                    fwrite(fid,pickupVS(1).x,'float');
                    
                    for l=1:16
                        fwrite(fid,pickupVS(l).y,'float');
                    end
                    
                    for l=1:16
                        fwrite(fid,pickupVB(l).y,'float');
                    end
                    
                    fclose(fid);
                    
                    disp([num2str(k) ' selected as Zero shot for shot ' num2str(shots(j))]);
                    
                    break;
                end
            end
            
        end
        
    end
   if(k==searchEnd)
        disp('ATTENTION: Do not find a suitable zero shot change UseZeroShot to 0 in MARTe cfg File!!!!');
    end
    
end






