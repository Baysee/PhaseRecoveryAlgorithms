function [xt,Si]=phaseRecovLoop(nIncsInterp,windowCentersInterp,lent,winInterp,winLen,t,dt,xt0,tspgm,fspgm,spgm,SUT,analysisWin,Fs,maxIteration)

% % % Plotting parameters
h1=figure;
h1.Position=[55 117 990 861];%[-1528 112 821 865];
xlimsZoom=t(round(lent/2))+winLen*dt*[-2 2];
diCLabel='Normalized inconsistency (dB)';
diRLabel='Frog Error (dB)';
SUTnorm=SUT/max(real(SUT));

%%% Initialization parameters
S0=abs(sqrt(spgm));
xt=xt0;

%%% Over correction parameters
% FROM FGLA
nominalAlpha=0.8;%0.8;%.8;
alpha=nominalAlpha;
tn=S0; tnm1=S0;

% FROM FORG OVERCorrection
beta=0

i=1;
while i<maxIteration+1
   

    %     Si=get_stft_fullSigLen(nIncsInterp,winIndsInterp,windowCentersInterp,lent,winInterp,dt,winLenInterp,xt);         % Get spgm of present signal guess xt
    Si=get_stft_fullSigLen(nIncsInterp,windowCentersInterp,lent,winInterp,dt,xt);         % Get spgm of present signal guess xt
%     Si=Si+alpha*(Si-tn)+beta*(spgm-abs(Si).^2)./abs(Si).^2;
  
if i>2
Si=Si+alpha*(Si-tn)+beta*(spgm-abs(Si).^2)./abs(Si).^2;
end
    Sip1=S0.*exp(1j*angle(Si));%.*Si./abs(Si);%.*exp(1j*angle(Si));%.*Si./abs(Si);                           % Enforce magnitude along with calculated phase from Si
    Sip1=S0.*Si./abs(Si);%.*exp(1j*angle(Si));%.*Si./abs(Si);                           % Enforce magnitude along with calculated phase from Si
    Sip1(isnan(Sip1))=0;
%     Sip1=Sip1+alpha*(Sip1-tn);
% Sip1=Sip1+alpha*(Sip1-tn)+beta*(spgm-abs(Sip1).^2)./abs(Sip1).^2;
% Sip1=Sip1+alpha*(Si-tn)+beta*Sip1.*(spgm-abs(Si).^2)./abs(Si).^2;

    xt=get_istft_fullSigLen(lent,windowCentersInterp,analysisWin,Fs,nIncsInterp,Sip1);


    diC(i)=sqrt(    sum(sum( abs(  abs(S0) - abs(Si)  ).^2))  /    sum(sum(  abs(S0).^2 )) );%)/norm(abs(Si).^2)

    %     xc=max(abs(xcorr(abs(SUT),abs(xt),50)))/(sqrt(sum(abs(SUT).^2)*sum(abs(xt).^2)));
    % %     diR(i)=(1-xc);%)/norm(abs(Si).^2)

    trebinoErrorspgm=spgm/max(max(spgm)); trebinoErrorSip1=abs(Si).^2/max(max(abs(Si).^2));
    diR(i)=sqrt(  1/numel(S0) * sum(sum(  (trebinoErrorspgm - trebinoErrorSip1).^2 ))  );

    %)/norm(abs(Si).^2)

    % plot(t,real(xt)); drawnow();



    if mod(i,5)==1
        diCPlot=10*log10(diC);
        diCLabel='Normalized inconsistency (dB)';
        diRPlot=10*log10(diR);
        diRLabel='Frog Error (dB)';
        updatePlot(h1,xlimsZoom,t,tspgm,fspgm,xt0,real(SUTnorm*max(real(xt))),xt,spgm,abs(Si).^2,diCPlot,diCLabel,diRPlot,diRLabel)
    end



%     if 10*log10(diC(i))<-8;
    if (diC(i))<(1-nominalAlpha) && alpha~=0;
        convDiff=10*log10((diC(i-1)-diC(i))/diC(i))
        alpha=1-diC(i);
        
        if convDiff<-20
            alpha=0.999
        end
% 
%         if abs(convDiff)>nominalAlpha
%             alpha=nominalAlpha
%         end
% 
%         if abs(convDiff)<nominalAlpha
%             alpha=1-abs(convDiff)
%         end
    end


            tn=Si;



    i=i+1;
end

diCPlot=10*log10(diC);
diRPlot=10*log10(diR);

updatePlot(h1,xlimsZoom,t,tspgm,fspgm,xt0,real(SUTnorm*max(real(xt))),xt,spgm,abs(Si).^2,diCPlot,diCLabel,diRPlot,diRLabel)

