function [xt,Si]=phaseRecovLoopTLS(nIncs,windowCenters,lent,winLen,t,f,dt,xt0,tspgm,fspgm,spgm,SUT,Fs,maxIteration,phit,phiw)

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
alpha1=0;%.7;%8;%0.8;%.8;
alpha2=0.4;
alpha3=0.8;
% FROM FORG OVERCorrection
beta1=1.1;
beta2=1.3;
beta3=1.5;

plotFilt=0;
fMaxTLS=fspgm(end)-fspgm(1);
filtm=50;
xt=filterSig(filtm,fMaxTLS,xt0,plotFilt,f);
% % 
% alpha1=0; alpha2=alpha1;alpha3=alpha1;
% beta1=1; beta2=beta1;beta3=beta1;

alpha=alpha1;
tn=S0; tnm1=S0;



beta=beta1;

i=1;
while i<maxIteration+1
   

    %     Si=get_stft_fullSigLen(nIncsInterp,winIndsInterp,windowCentersInterp,lent,winInterp,dt,winLenInterp,xt);         % Get spgm of present signal guess xt
    Si=stft_TLS(xt,winLen,phit,phiw);%get_stft_fullSigLen(nIncsInterp,windowCentersInterp,lent,winInterp,dt,xt);         % Get spgm of present signal guess xt
%     Si=Si+alpha*(Si-tn)+beta*(spgm-abs(Si).^2)./abs(Si).^2;
  
['current parameters: alpha: ' num2str(alpha)  ' beta: ' num2str(beta)]
% if i>2
Si=Si+alpha*(Si-tn);%+beta*(spgm-abs(Si).^2)./abs(Si).^2;
% Si=Si.*(S0./Si).^beta; Sip1(isnan(Sip1))=0;
% end
%     Sip1=S0.*exp(1j*angle(Si));%.*Si./abs(Si);%.*exp(1j*angle(Si));%.*Si./abs(Si);                           % Enforce magnitude along with calculated phase from Si
   Sip1=Si.*abs(S0./Si).^beta; Sip1(isnan(Sip1))=0;

% Sip1=S0.*Si./abs(Si); Sip1(isnan(Sip1))=0;                          % Enforce magnitude along with calculated phase from Si
   
%    Sip1=S0.*exp(1j*angle(Si));%.*Si./abs(Si); 
%     Sip1=Sip1+alpha*(Sip1-tn);
% Sip1=Sip1+alpha*(Sip1-tn)+beta*(spgm-abs(Sip1).^2)./abs(Sip1).^2;
% Sip1=Sip1+alpha*(Si-tn)+beta*Sip1.*(spgm-abs(Si).^2)./abs(Si).^2;

    xt=istft_TLS(Sip1,winLen,phit,phiw);%get_istft_fullSigLen(lent,windowCenters,analysisWin,Fs,nIncs,Sip1);
xt=filterSig(filtm,fMaxTLS,xt,plotFilt,f);


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



    if 10*log10(diC(i))<-5 && i>1;
%     if (diC(i))<(1-nominalAlpha) && alpha~=0;
  
    beta=beta2;
   alpha=alpha2;
        convDiff=10*log10(abs(diC(i-1)-diC(i))/diC(i))
%         alpha=1-diC(i);
        
        if convDiff<-20 && 10*log10(diC(i))<-7;
            beta=beta3;
            alpha=alpha3;
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

