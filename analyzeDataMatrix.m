% datDir='C:\Users\Lord Photon\Documents\MATLAB\time-frequency analysis\outputFigs\'
% datMatName='dataMatrix21h58m.mat';

datDir='/Users/ben/Documents/MATLAB/timeFrequencyAnalysis/phaseRecoveryAlgorithm_outputFigs/outputFigs/'
datMatName='dataMatrixAll01h50m.mat';

load([datDir datMatName])


% % define loop parameters
clips_l=[20,40,60];
nptPerWin_l=[16,32,64];
bwIRF_l=[30:5:70]*1e9;
interpAmounts=[1,2,4,8,16];
% (icl, inptWin, ibwIRF, iInter);
b21=-12758; b22=5082;


% Check only unInterpolated data
alldiRs=alldiRs(:,:,:,1);
allMeanb2s1=allMeanb2s1(:,:,:,1);
allMeanb2s2=allMeanb2s2(:,:,:,1);
allStdb2s1=allStdb2s1(:,:,:,1);
allStdb2s2=allStdb2s2(:,:,:,1);

% Normalize data
b21Error=10*log10(abs(allMeanb2s1-b21)/abs(b21));
b22Error=10*log10(abs(allMeanb2s2-b22)/abs(b22));
CV1=10*log10(abs(allStdb2s1/b21));
CV2=10*log10(abs(allStdb2s2/b22));


varNames={'alldiRs','b21Error','b22Error','CV1','CV2'};
poses=zeros(5,numel(CV1)); % Here I will store the linear index of each metric in the order of best performance

for iv=1:numel(varNames);
    
    varNow=eval(varNames{iv}); % Select metric to analyze
    [minVals, poses(iv,:)] = sort(varNow(:));
    [i,j,k,l] = ind2sub(size(varNow),poses(iv,:));

    figure;
    subplot(3,1,1)
    plot(minVals)
        title(varNames{iv})

    subplot(3,1,2)
    plot(i); ylabel('clipVal')
    yyaxis right; plot(j); ylabel('nptPerWin');
    subplot(3,1,3);
    plot(k); ylabel('bw IRF');
    yyaxis right
    plot(l); ylabel('interpAmount')
    
    
end

noneFound=0
iter=1;
nParas=4;%numel(poses(:,1));
% nParas=numel(poses(:,1));
while noneFound==0;
    
%     currentMatUnique=unique(poses(:,iter));
%     nReps=numel(poses(:,iter))-numel(currentMatUnique);
    matNow=poses(2:5,1:iter);
%     matNow=poses(:,1:iter);
    A=matNow; A=A(:);
    [x, b, uidx] = unique(A);
% counts = accumarray(uidx, 1);
N=numel(x);

   count = zeros(N,1);
   for k = 1:N
      count(k) = sum(A==x(k));
   end
   
[MaxCount,Ind]=max(count);


    if MaxCount==nParas
        bestIndexOverall=x(Ind);
        noneFound=1;
        
        [iOpt,jOpt,kOpt,lOpt] = ind2sub(size(varNow),bestIndexOverall);
    end
    iter=iter+1;
end
    
    
%# finds the max of A and its position, when A is viewed as a 1D array
% [max_valb22, posMeanB22] = sort(b22Error(:)); 
% [max_valb21, posMeanB21] = sort(b21Error(:)); 
% [max_valCV2, posCV2] = sort(CV2(:)); 
% [max_valCV1, posCV1] = sort(CV1(:)); 
% [i,j,k,l] = ind2sub(size(CV1),posCV1);
% 













%#transform the index in the 1D view to 4 indices, given the size of A

% 
% 
% dim1=2;dim2=2;
% datToPlot1=squeeze(alldiRs(dim1,dim2,:,:));plotTitle1='Frog error';
% datToPlot2=squeeze(b21Error(dim1,dim2,:,:));plotTitle2='b2 error';
% datToPlot3=squeeze((CV(dim1,dim2,:,:)));plotTitle3='std(b2)';
% datToPlot4=squeeze(b22Error(dim1,dim2,:,:));plotTitle2='b2 error';
% datToPlot5=squeeze((CV(dim1,dim2,:,:)));plotTitle3='std(b2)';
% 
% 
% 
% figure;
% subplot(3,2,1:2);
% imagesc(interpAmounts,bwIRF_l,datToPlot1);colorbar(); title(plotTitle1)
% subplot(3,2,3);
% imagesc(interpAmounts,bwIRF_l,datToPlot2);colorbar(); title(plotTitle2)
% subplot(3,2,5);
% imagesc(interpAmounts,bwIRF_l,datToPlot3);colorbar(); title(plotTitle3)
% 
% subplot(3,2,4);
% imagesc(interpAmounts,bwIRF_l,datToPlot4);colorbar(); title(plotTitle2)
% subplot(3,2,6);
% imagesc(interpAmounts,bwIRF_l,datToPlot5);colorbar(); title(plotTitle3)