function Spri = Dreicer( ne,Zeff,Te,vloop)
    % Drecier formula 
    globalVar;    
    global c e me lnD eps0 R_torus;
    e=abs(e);
    
    %lnD = log(4.9*1e7*( (Te.^(3/2)) ./ ((ne*1e-20).^(1/2)) ));    
    %The value for lnLambda I use is 18, because they told me that 17 or 18 is more or less typical for FTU. 
    lnD = 18;    
    J=Te*1e3*e;                    %we must convert keV in J
    Ell=vloop;%/(2*pi*R_torus);
    k=0.21+0.11*Zeff;
    Spri.k=k;
    vTe=sqrt(2*J/me);
    vee=(ne*e^4.*lnD)./(4*pi*eps0^2*me^2*(vTe.^3));
    Spri.vee=vee;
    ED=(ne*e^3.*lnD)./(4*pi*eps0^2*J);
    eD=Ell./ED;    
    Spri.eD=eD;
    
    
    DreicerExpFirst = -(1/(4*eD))'-sqrt((1+Zeff)./eD); 
    Spri.DreicerExpFirst=exp(DreicerExpFirst);
    DreicerExpSecond = -(J/(me*c^2)).*((1./(8*(eD.^2)))+2/3*sqrt((1+Zeff)./(eD.^3)));
    Spri.DreicerExpSecond=exp(DreicerExpSecond);
    DreicerExp = exp(DreicerExpFirst+DreicerExpSecond);
    Spri.DreicerExp= exp(DreicerExp);
    eDZeff = (eD.^((-3/16)*(1+Zeff)));
    Spri.eDZeff = eDZeff;    
    Spri.Spri=k.*ne.*vee.*eDZeff.*DreicerExp;
    
end

