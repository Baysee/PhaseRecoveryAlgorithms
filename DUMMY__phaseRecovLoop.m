function [xt,Si,diCPlot,diRPlot]=DUMMY__phaseRecovLoop(nIncsInterp,windowCentersInterp,lent,winInterp,winLen,t,dt,...
    xt0,tspgm,fspgm,spgm,SUT,analysisWin,Fs,maxIteration,metrics,plotIter,filtm,fMaxStft,plotFilt,f)
%% metrics is a 2-element cell, corresponding to the left and right axis content of the convergence curve

% if plotIter;end
    % % % % % % Make video
    makeVid=0;
    if makeVid
        v = VideoWriter('../GLA2');
        v.FrameRate=2; open(v);
    end
    
    % % % Plotting parameters
    h1=figure;
    h1.Position=[55 117 990 861];%[-1528 112 821 865];
    xlimsZoom=t(round(lent/2))+winLen*dt*[-5 5];
    diCLabel='Normalized inconsistency (dB)';
    diRLabel='Frog Error (dB)';
    SUTnorm=SUT/max(real(SUT));


%%% Initialization parameters
S0=abs(sqrt(spgm));
xt=xt0;
xt=filterSig(filtm,fMaxStft,xt0,plotFilt,f);

%%% Over correction parameters
% FROM FGLA
alpha1=0.1;%.7;%8;%0.8;%.8;
alpha2=0.2;
alpha3=0.4;
% FROM FORG OVERCorrection
beta1=1.15;
beta2=1.2;
beta3=1.4;

%
%
%
% alpha1=0; alpha2=alpha1;alpha3=alpha1;
% beta1=1; beta2=beta1;beta3=beta1;

alpha=alpha1;
tn=S0; tnm1=S0;

Si=S0;

beta=beta1;

i=1;
while i<maxIteration+1
    
    tn=Si;
    
    alpha
    beta
    Si=get_stft_fullSigLen(nIncsInterp,windowCentersInterp,lent,winInterp,dt,xt);         % Get spgm of present signal guess xt
    
    Si=Si+alpha*(Si-tn);
    
    %     Sip1=S0.*exp(1j*angle(Si));%                   % Enforce magnitude along with calculated phase from Si
   SiDiv=Si; 
   SiDiv(abs(SiDiv)<max(max(abs(SiDiv)))*1e-9)=1;
    Sip1=Si.*abs(S0./SiDiv).^beta; Sip1(isnan(Sip1))=0;
%     Sip1=Si.*abs(S0./Si).^beta; Sip1(isnan(Sip1))=0;
    

    xt=get_istft_fullSigLen(lent,windowCentersInterp,analysisWin,Fs,nIncsInterp,Sip1);
    xt=filterSig(filtm,fMaxStft,xt,plotFilt,f);
    
    trebinoErrorspgm=spgm/max(max(spgm)); trebinoErrorSip1=abs(Si).^2/max(max(abs(Si).^2));
    diR(i)=sqrt(  1/numel(S0) * sum(sum(  (trebinoErrorspgm - trebinoErrorSip1).^2 ))  );
    diRPlot=10*log10(diR);
    diRLabel='Frog Error (dB)';
    
    
    
%     if plotIter
        if strcmp(metrics{1},'Inconsistency')
            diC(i)=sqrt(    sum(sum( abs(  abs(S0) - abs(Si)  ).^2))  /    sum(sum(  abs(S0).^2 )) );%)/norm(abs(Si).^2)
            diCPlot=10*log10(diC);
            diCLabel='Normalized inconsistency (dB)';
        else if strcmp(metrics{1},'cc')
                %         xc=max(abs(xcorr(abs(SUT).^2   ,abs(xt).^2,50)))/(sqrt(sum(abs(SUT).^4)*sum(abs(xt).^4)))
                xc4=corrcoef(abs(SUT).^2,abs(xt).^2);
                xc=xc4(1,2);
                
                diC(i)=(1-xc);%)/norm(abs(Si).^2)
                diCPlot=10*log10(diC);
                
                diCLabel='10Log(1-cc)';
                
            else if strcmp(metrics{1},'mse')
                    %         diC(i)=mse(abs(SUT).^2,abs(xt).^2,'normalization','standard');
                    %         diCPlot=10*log10(diC); diCLabel='10Log(MSE)';
                    
                    diC(i)=sqrt(mean((abs(SUT).^2/max(abs(SUT).^2)-abs(xt).^2/max(abs(xt).^2)).^2));
                    diCPlot=10*log10(diC); diCLabel='10Log(MSE)';
                end
            end
        end
         if plotIter
        if mod(i,5)==1
            updatePlot(h1,xlimsZoom,t,tspgm,fspgm,xt0,real(SUTnorm*max(real(xt))),xt,spgm,abs(Si).^2,diCPlot,diCLabel,diRPlot,diRLabel)
            if makeVid
                writeVideo(v,getframe(gcf));% M(i)=getframe;
            end
        end
        
    end
    
    
    
    if diRPlot(i)<-10 && i>1
        %     if (diC(i))<(1-nominalAlpha) && alpha~=0;
        
        beta=beta2;
        alpha=alpha2;
        convDiff=10*log10(abs(diR(i-1)-diR(i))/diR(i))
        %         alpha=1-diC(i);
        
        if convDiff<-15 && diRPlot(i)<-7;
            beta=beta3;
            alpha=alpha3;
        end
     
    end
    
    
    
    
    
    i=i+1;
end


 if strcmp(metrics{1},'Inconsistency')
            diC(i)=sqrt(    sum(sum( abs(  abs(S0) - abs(Si)  ).^2))  /    sum(sum(  abs(S0).^2 )) );%)/norm(abs(Si).^2)
            diCPlot=10*log10(diC);
            diCLabel='Normalized inconsistency (dB)';
        else if strcmp(metrics{1},'cc')
                %         xc=max(abs(xcorr(abs(SUT).^2   ,abs(xt).^2,50)))/(sqrt(sum(abs(SUT).^4)*sum(abs(xt).^4)))
                xc4=corrcoef(abs(SUT).^2,abs(xt).^2);
                xc=xc4(1,2);
                
                diC(i)=(1-xc);%)/norm(abs(Si).^2)
                diCPlot=10*log10(diC);
                
                diCLabel='10Log(1-cc)';
                
            else if strcmp(metrics{1},'mse')
                    %         diC(i)=mse(abs(SUT).^2,abs(xt).^2,'normalization','standard');
                    %         diCPlot=10*log10(diC); diCLabel='10Log(MSE)';
                    
                    diC(i)=sqrt(mean((abs(SUT).^2/max(abs(SUT).^2)-abs(xt).^2/max(abs(xt).^2)).^2));
                    diCPlot=10*log10(diC); diCLabel='10Log(MSE)';
                end
            end
        end
% if plotIter
    diCPlot=10*log10(diC);
    diRPlot=10*log10(diR);
    
    updatePlot(h1,xlimsZoom,t,tspgm,fspgm,xt0,real(SUTnorm*max(real(xt))),xt,spgm,abs(Si).^2,diCPlot,diCLabel,diRPlot,diRLabel)
    if makeVid
        writeVideo(v,getframe(gcf));% M(i)=getframe;
        close(v)
    end
% end