%% 60Hz removal
function Filter=Filter60Hz(RawFile, OutFile, Filter)
    %Want to get rid of powerline noise in data. Can assume that powerline
    %noise is reasonably close to 60Hz.
    
    %Procedure:
    %1)Estimate exact frequency
    %cut the recording in blocks of 6x6 periods (with slightly varying
    %periods). Average over columns and maximize L1 norm of resulting
    %vector to determine current frequency. Use a prior for the
    %maximixation and also smoothen the signal to avoid spurious jumps.
    %2)Estimate phase and average
    %get new phase via frequency estimate
    %refine by correlating with past traces
    %median over past traces (5, 108 periods, need to remember all for (3)) 
    %as correction use phase and frequency to
    %determine indices for corrective array. apply correction to raw trace,
    %need to fit amplitude for individual channels.
    %use a slow (AR filtered) avg. to bias phase estimate, in case of low
    %amplitude of oscillation
    %3)(optional) correct single channel noise by averaging past and future
    %of a single channel and correlating with its trace. use only
    %a single block for that average (low amplitude expected), i.e. 36 osc.
    
    Nch=length(Filter.ChMask);
    Nchx=sum(Filter.ChMask,2);
    ChMask=Filter.ChMask;
    ChMap=Filter.ChMap;
    iGroup=0;
    
    bitVolt=h5read(RawFile,['/recordings/' num2str(iGroup) '/application_data/channel_bit_volts']);
    sRate=h5read(RawFile,['/recordings/' num2str(iGroup) '/application_data/channel_sample_rates']);
    time_stamps=h5read(RawFile,['/recordings/' num2str(iGroup) '/application_data/timestamps']);
    %DataSize=size(h5read(RawFile,'/recordings/' num2str(iGroup) '/data'));%probably need to find a better solution here for large files
    a=h5info(RawFile);
    DataSize=a.Groups.Groups(iGroup+1).Datasets.Dataspace.Size;
    if Filter.tEnd==-1
        Filter.nEnd=double(floor(DataSize(1,2)));
    else
        Filter.nEnd=double(floor(Filter.tEnd*60*sRate(1,1)));
    end
    if Filter.tStart==0
        Filter.nStart=1;
    else
        Filter.nStart=double(floor(Filter.tStart*60*sRate(1,1))+1);
    end
    
    version_name=h5readatt(RawFile,'/','kwik_version');
    %rec_method=h5readatt(RawFile,'/recordings/0','name');
    start_time=h5readatt(RawFile,'/recordings/0','start_time');
    start_sample=h5readatt(RawFile,'/recordings/0','start_sample');
    sample_rate=h5readatt(RawFile,'/recordings/0','sample_rate');
    bit_depth=h5readatt(RawFile,'/recordings/0','bit_depth');
    
    
    h5create(OutFile,'/recordings/0/data',[Nchx Inf],'ChunkSize',[1 2048],'Datatype','int16','FillValue',int16(0))
    h5writeatt(OutFile,'/','kwik_version',version_name);
    %h5writeatt(OutFile,'/recordings/0','name',rec_method{:});
    h5writeatt(OutFile,'/recordings/0','start_time',start_time);
    h5writeatt(OutFile,'/recordings/0','start_sample',start_sample);
    h5writeatt(OutFile,'/recordings/0','sample_rate',sample_rate);
    h5writeatt(OutFile,'/recordings/0','bit_depth',bit_depth);
    h5create(OutFile,'/recordings/0/application_data/channel_bit_volts',Nchx);
    h5create(OutFile,'/recordings/0/application_data/channel_sample_rates',Nchx);
    h5create(OutFile,'/recordings/0/application_data/timestamps',[Nch Inf],'ChunkSize',[Nch 16]);
    h5writeatt(OutFile,'/recordings/0/application_data','is_multiSampleRate_data',0);
    
    h5write(OutFile,'/recordings/0/application_data/channel_bit_volts',bitVolt(ChMap(ChMask),:));
    h5write(OutFile,'/recordings/0/application_data/channel_sample_rates',sRate(ChMap(ChMask),1));
    h5write(OutFile,'/recordings/0/application_data/timestamps',time_stamps(1:Nch,:),[1 1],[Nch size(time_stamps,2)]);
    
    bitVolt=median(bitVolt(ChMap(ChMask),1));
    sRate=median(sRate(ChMap(ChMask),1));
    
    %Parameters
    Filter.LenRec=Filter.nEnd-Filter.nStart+1;
    Filter.bitVolt=bitVolt;
    Filter.sRate=sRate;
    Filter.Nch=Nchx;
    
    npx=12;
    nx=double(round(sRate*npx/60));
    nsx=double(round(sRate/60));
    nsxh=floor(nsx/2);
    nT=2000;%bins for time warping (1 period)
    nTh=nT/2;
    npy=2;
    npyh=npy/2;
    %nBlk=nx*npy;
    nBlkh=nx*npyh;
    %nBlkq=nBlkh/2;
    %nBlkd=2*nx*npy;
    npBlk=npx*npy;
    npBlkh=npx*npyh;
    
    %sx=double(sRate/30);
    Nfft=2^14;
    fftBins=1024;
    
    deltaNh=double(round(sRate/5000));
    deltaN=2*deltaNh;%want this to be even
    fPrior=hanning(deltaN+1);
    fPrior=fPrior*npyh/sum(fPrior);
    fP=hanning(5);
    fP=fP/sum(fP);
    
    fPriorLong=1*fPrior;
    fTau=20;
    fA=exp(-1/fTau);
    fScale=0.1*(1-fA);
    fAR=zeros(2*deltaN+1,1);    
    
    nhAll=2;
    nBlkAll=2*nhAll+1;
    nhSingle=2;
    nBlkSingle=2*nhSingle+1;
    %nBlkRaw=nBlkAll+nhSingle;
    nBlkX=nBlkh+deltaNh*npyh;
    nDelay=nhAll+nhSingle+1;%time lag after which two step processing finishes and data can be written   
    
    dtPhase=double(deltaNh*npyh*nT*60/sRate(1,1));
    %Variables
    
    xCut=zeros(4*nBlkX,Nchx);
    fPow=zeros(15,2*deltaN+1);
    
    %for phase alignment use sine wave
    xSine=reshape(sin(2*pi*(1:(nT+1)*nT)/nT),nT+1,nT);
    xSine=xSine(1:nT,:);
    


    xFilt=zeros(4*nBlkX,nhSingle+1,Nchx);
    
    xOsc=zeros(nT,nBlkAll);%modulo
    sOsc=zeros(nT,nBlkSingle,Nchx);%modulo
    
    skipBlk=zeros(nDelay+1,1);%rolling
    sepFrame=zeros(nDelay+1,1);%rolling
    
    Ftmp=complex(zeros(fftBins,nhSingle+1,Nchx),0);%modulo
    Fctmp=complex(zeros(fftBins,nhSingle+1,Nchx),0);%modulo
    
    cAll=zeros(nBlkX,nhSingle+1,Nchx);%global correction,modulo %don't remember as long
    %cSingle=zeros(nBlkX,nhSingle+1,Nch);%single channel
    %correction,modulo%no need to remember
    wPhase=zeros(nDelay,1);%modulo
    fXa=zeros(nDelay,1);%modulo
    %phasePrior=zeros(1,nT);
    stepTimes=zeros(1);
    stepThr=30/bitVolt(1,1);
    
    %warmup
    for k=10:-1:1%need to fix
        xRaw=double(h5read(RawFile,['/recordings/' num2str(iGroup) '/data'],[1,Filter.nStart+k*nBlkX-nsx],[Nch,2*nBlkX+2*nsx+1])');
        xRaw=xRaw(:,ChMap(ChMask));
        xCum=cumsum(xRaw,1);
        xB=(xCum(2*nsx+2:end,:)-xCum(2:end-2*nsx,:))/(4*nsx)...
            +(xCum(nsx+nsxh+2:end-nsx+nsxh,:)-xCum(2+nsxh:end-2*nsx+nsxh,:))/(2*nsx);
        xCut(1:2:end-1,:)=(xCum(nsx+3:end-nsx+1,:)-xCum(nsx:end-nsx-2,:))/3-xB;
        xCut(2:2:end,:)=(xCum(nsx+3:end-nsx+1,:)-xCum(nsx+1:end-nsx-1,:))/2-xB;
        for j=0:2*deltaN
            xBlk=reshape(xCut(1:end-j*npy,:),2*nx+deltaN-j,npy,Nchx);
            fPow(end,j+1)=median(mean(abs(mean(xBlk,2)),1),3);
        end
        fMean=median(fPow,1)'+fScale*fAR;
        fAR=fA*fAR+fMean-fScale*fAR;%in case of small osc.
        fCM=conv(fMean,fPrior,'same');%insensitive to individual channel phase
        [~,fXarg]=max(fCM(deltaNh+1:end-deltaNh));%need to constrain for prior
        %create an updated prior
        fConv=conv(fMean,fP,'same');
        fConv=fConv(fXarg:fXarg+deltaN)-min(fConv(fXarg:fXarg+deltaN));
        fP0=(fConv+fConv(end:-1:1))*0.5/sum(fConv);
        fPriorLong=fA*fPriorLong+fScale*fP0;
        fPrior=0.1*fPriorLong+0.9*fP0;%need a separate, long-term component as a bias term
        fPow=circshift(fPow,-1,1);
    end

    %conversion to warped time (nT frames per osc.), always use 1 block)
    Ind0=cell(1,2*deltaNh+1);
    fracInd0=cell(1,2*deltaNh+1);
    for i=1:2*deltaNh+1
        Ind0{i}=floor(((1:(npBlk)*nT)'-1)*(2*nx+i-deltaNh-1)/(npx*nT))+1;%need to adjust this conversion (npx*nT),(npx*npy+1)*nT
        fracInd0{i}=1-mod(((1:(npBlk)*nT)'-1)*(2*nx+i-deltaNh-1)/(npx*nT),1);
    end
    %conversion to real time (convert 1 block)
    Ind1=cell(1,2*deltaNh+1);
    Ind2=cell(1,2*deltaNh+1);
    Ind3=cell(1,2*deltaNh+1);
    fracInd1=cell(1,2*deltaNh+1);
    for i=1:2*deltaNh+1
        %Ind3{i}=floor(((1:(2*nx+i-deltaNh-1))-1)*2*nx/(2*nx+i-deltaNh-1))+1;
        Ind3{i}=floor(((1:4*nBlkX)'-1)*(npx*nT)/(2*nx+i-deltaNh-1))+1;
        Ind1{i}=mod(Ind3{i}-1,nT)+1;
        Ind2{i}=mod(Ind3{i},nT)+1;
        fracInd1{i}=1-mod(((1:4*nBlkX)'-1)*(npx*nT)/(2*nx+i-deltaNh-1),1);
    end
    
    numRec=floor((Filter.LenRec-nsx)/nBlkX)-1;%number of iterations
    
    xRaw=zeros(nBlkX+2*nsx,nDelay+1,Nchx);%rolling
    xCum=zeros(nBlkX+1+2*nsx,1,Nchx);
    xSmth=zeros(2*nBlkX,nDelay+1,Nchx);%double resolution,rolling
    
    %first read in necessary data (nBlkAll+nBlkSingle/2 blocks)
    xRaw0=double(h5read(RawFile,['/recordings/' num2str(iGroup) '/data'],[1,Filter.nStart],[Nch,nBlkX+nsx])');
    xRaw(1+nsx:end,end,:)=xRaw0(:,ChMap(ChMask));
    xRaw(nsx,end,:)=xRaw(nsx+1,end,:);%boundary...
    xRaw(1:nsx-1,end,:)=xRaw(nsx+1:2*nsx-1,end,:);
    xCum(2:end,1,:)=cumsum(xRaw(:,end,:),1);
    xB=(xCum(2*nsx+2:end,1,:)-xCum(2:end-2*nsx,1,:))/(4*nsx)...
        +(xCum(nsx+nsxh+2:end-nsx+nsxh,1,:)-xCum(2+nsxh:end-2*nsx+nsxh,1,:))/(2*nsx);
    xSmth(1:2:end-1,end,:)=(xCum(nsx+3:end-nsx+1,1,:)-xCum(nsx:end-nsx-2,1,:))/3-xB;
    xSmth(2:2:end,end,:)=(xCum(nsx+3:end-nsx+1,1,:)-xCum(nsx+1:end-nsx-1,1,:))/2-xB;
    xRaw=circshift(xRaw,-1,2);
    xSmth=circshift(xSmth,-1,2);
    for k=1:nDelay%need to fix
        %modulo indices
        mG=mod(k-1,nBlkAll)+1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        m1G=mod(k-1-nhAll,nBlkAll)+1;
        mAll=mod(k-1,nDelay)+1;
        m1All=mod(k-1-nhAll,nDelay)+1;
        m2All=mod(k,nDelay)+1;
        mSingle=mod(k-1,nhSingle+1)+1;
        m1Single=mod(k,nhSingle+1)+1;
        
        kAll=k-nhAll;%index for which global filtered avg is becoming available
        kSingle=k-nhAll-nhSingle;%index for which single filtered avg is becoming available
        xRaw0=double(h5read(RawFile,['/recordings/' num2str(iGroup) '/data'],[1,k*nBlkX-nsx+Filter.nStart],[Nch,nBlkX+2*nsx])');
        xRaw(:,end,:)=xRaw0(:,ChMap(ChMask));
        xCum(2:end,1,:)=cumsum(xRaw(:,end,:),1);
        xB=(xCum(2*nsx+2:end,1,:)-xCum(2:end-2*nsx,1,:))/(4*nsx)...
            +(xCum(nsx+nsxh+2:end-nsx+nsxh,1,:)-xCum(2+nsxh:end-2*nsx+nsxh,1,:))/(2*nsx);
        xSmth(1:2:end-1,end,:)=(xCum(nsx+3:end-nsx+1,1,:)-xCum(nsx:end-nsx-2,1,:))/3-xB;
        xSmth(2:2:end,end,:)=(xCum(nsx+3:end-nsx+1,1,:)-xCum(nsx+1:end-nsx-1,1,:))/2-xB;
        xCut=reshape(xSmth(:,end-1:end,:),4*nBlkX,Nchx);
        for j=0:2*deltaN
            xBlk=reshape(xCut(1:end-j*npy,:),2*nx+deltaN-j,npy,Nchx);
            fPow(end,j+1)=median(mean(abs(mean(xBlk,2)),1),3);
        end
        fMean=mean(fPow,1)'+fScale*fAR;
        fAR=fA*fAR+fMean-fScale*fAR;%in case of small osc.
        fCM=conv(fMean,fPrior,'same');%insensitive to individual channel phase
        [~,fXarg]=max(fCM(deltaNh+1:end-deltaNh));%need to constrain for prior
        fXa(mAll)=fXarg;
        %create an updated prior
        fConv=conv(fMean,fP,'same');
        fConv=fConv(fXarg:fXarg+deltaN)-min(fConv(fXarg:fXarg+deltaN));
        fP0=(fConv+fConv(end:-1:1))*0.5/sum(fConv);
        fPriorLong=fA*fPriorLong+fScale*fP0;
        fPrior=0.1*fPriorLong+0.9*fP0;%need a separate, long-term component as a bias term
        
        
        %interpolate at 60 Hz grid
        xWarp=reshape(repmat(fracInd0{fXarg},1,Nchx).*xCut(Ind0{fXarg},:)...
            +repmat(1-fracInd0{fXarg},1,Nchx).*xCut(Ind0{fXarg}+1,:),nT,npBlk,Nchx);
        %median after normalization of amplitude
        xStd=std(xCut,0,1);
        xWmed=median(squeeze(median(xWarp,2))./repmat(xStd,nT,1),2);
        
        %determine phase
        %phasePrior=0.2*circshift(phasePrior,deltaNh*npy,2)+xWmed'*xSine;
        [~,hPhase]=max(xWmed'*xSine);
        wPhase(mAll,1)=hPhase-1;%phase for the current avg.
        
        %%global avg.
        %align
        xOsc(:,mG)=0.5*circshift(xWmed,wPhase(mAll,1),1)-0.5*circshift(xWmed,wPhase(mAll,1)+nTh,1);%aligned with sine wave
        %fix stuff with no data
        if k==1
            wPhase(2:end,1)=wPhase(1,1);
            xOsc(:,2:end)=repmat(xOsc(:,1),1,nBlkAll-1);
        end
        %median oscillation (should have a lag of (nBlkAll+1)/2 blocks here
        if kAll>0%k==3
            %avg. and add past phase
            xOmed=circshift(mean(xOsc,2)+xOsc(:,m1G),-wPhase(m1All,1),1)/2;
        
            xBase=(fracInd1{fXa(m1All)}.*xOmed(Ind1{fXa(m1All)},1)+...
            (1-fracInd1{fXa(m1All)}).*xOmed(Ind2{fXa(m1All)},1))*ones(1,Nchx);%need to fix...
            xCutD=reshape(xSmth(:,end-nhAll-1:end-nhAll,:),4*nBlkX,Nchx);
            
            
            
            xProj=repmat((sum(xCutD(1:4*nBlkX,:).*xBase,1)./(sum(xBase.^2,1)+1)),4*nBlkX,1).*xBase;
            %smoothing to compensate for phase shift -- later
            %estimate global 60Hz
            if kAll==1
                xProjLast=xProj(1,:);
                cAll(:,mSingle,:)=(xProj(2*nBlkX:2:end-2,:)+2*xProj(2*nBlkX+1:2:end-1,:)+xProj(2*nBlkX+2:2:end,:))/8;
            end
            %cAll(:,mSingle,:)=squeeze(cAll(:,mSingle,:))+(xProj(2*nBlkX:2:end-2,:)+2*xProj(2*nBlkX+1:2:end-1,:)+xProj(2*nBlkX+2:2:end,:))/8;
            cAll(:,mSingle,:)=squeeze(cAll(:,mSingle,:))+([xProjLast;xProj(2:2:2*nBlkX-2,:)]+2*xProj(1:2:2*nBlkX-1,:)+xProj(2:2:2*nBlkX,:))/8;
            %remainder
            xFilt(:,mSingle,:)=xCutD-xProj;%Do remember!
            
            %%single channel avg.
            %convert to warped time
            wFilt=reshape(repmat(fracInd0{fXa(m1All)},1,Nchx).*squeeze(xFilt(Ind0{fXa(m1All)},mSingle,:))...
                    +repmat(1-fracInd0{fXa(m1All)},1,Nchx).*squeeze(xFilt(Ind0{fXa(m1All)}+1,mSingle,:)),nT,npBlk,Nchx);
            %correct for phase and avg.
            sOh=circshift(median(wFilt,2),wPhase(m1All,1),1);
            sOsc(:,mSingle,:)=sOh-repmat(mean(sOh,1),nT,1,1);
            if kAll==1
                sOsc=repmat(sOsc(:,mSingle,:),1,nhSingle+1,1);
            end
            F=fft(complex(xCutD(1:Nfft,:)));
            Fc=fft(complex(xProj(1:Nfft,:)));
            Ftmp(:,mSingle,:)=F(1:fftBins,:);
            Fctmp(:,mSingle,:)=Fc(1:fftBins,:);
            %projection on single traces
            if k>nDelay-1
                %add past phase
                sOmed=circshift(squeeze(median(sOsc,2)),-wPhase(m2All,1),1);
                sBase=((fracInd1{fXa(m2All)}*ones(1,Nchx)).*sOmed(Ind1{fXa(m2All)},:)+...
                (1-fracInd1{fXa(m2All)}*ones(1,Nchx)).*sOmed(Ind2{fXa(m2All)},:));
                sProj=repmat(sum(squeeze(xFilt(:,m1Single,:)).*sBase,1)./(sum(sBase.^2,1)+1),4*nBlkX,1).*sBase;
                if kSingle==1
                    sProjLast=sProj(1,:);
                    cSingle=(sProj(2*nBlkX:2:end-2,:)+2*sProj(2*nBlkX+1:2:end-1,:)+sProj(2*nBlkX+2:2:end,:))/8;
                end
                cSingle=cSingle+([sProjLast;sProj(2:2:2*nBlkX-2,:)]+2*sProj(1:2:2*nBlkX-1,:)+sProj(2:2:2*nBlkX,:))/8;
                %cSingle=cSingle+(sProj(2*nBlkX:2:end-2,:)+2*sProj(2*nBlkX+1:2:end-1,:)+sProj(2*nBlkX+2:2:end,:))/8;
                h5write(OutFile,'/recordings/0/data',int16(round((squeeze(xRaw(nsx+1:nBlkX+nsx,end-nDelay,:))...
                    -squeeze(cAll(:,m1Single,:))-cSingle)')),[1,(k-nDelay)*nBlkX+1],[Nchx,nBlkX]);
                %cSingle=([sProjLast;sProj(2:2:2*nBlkX-2,:)]+2*sProj(1:2:2*nBlkX-1,:)+sProj(2:2:2*nBlkX,:))/8;
                cSingle=(sProj(2*nBlkX:2:end-2,:)+2*sProj(2*nBlkX+1:2:end-1,:)+sProj(2*nBlkX+2:2:end,:))/8; 
                sProjLast=sProj(2*nBlkX,:);
            end
            %cAll(:,m1Single,:)=([xProjLast;xProj(2:2:2*nBlkX-2,:)]+2*xProj(1:2:2*nBlkX-1,:)+xProj(2:2:2*nBlkX,:))/8;
            cAll(:,m1Single,:)=(xProj(2*nBlkX:2:end-2,:)+2*xProj(2*nBlkX+1:2:end-1,:)+xProj(2*nBlkX+2:2:end,:))/8;
            xProjLast=xProj(2*nBlkX,:);%%%later
        end
        fPow=circshift(fPow,-1,1);
        xRaw=circshift(xRaw,-1,2);
        xSmth=circshift(xSmth,-1,2);
    end

    xOscSgle=zeros(nT,Nchx);
    xOscAll=zeros(nT,Nchx);
    xOscOld=zeros(nT,1);
    %xProjAvg=zeros(nBlkh,Nch);
    %xRawAvg=zeros(nBlkh,Nch);
    xPhase=zeros(numRec,1);%ok
    xRawV=zeros(numRec,Nchx);
    xCV=zeros(numRec,Nchx);
    xSV=zeros(numRec,Nchx);
    xNV=zeros(numRec,Nchx);
    %xShiftRes=0;
    %dxPhase=0;
    QS=complex(zeros(fftBins,numRec),0);
    QSc=complex(zeros(fftBins,1),0);
    QSs=complex(zeros(fftBins,1),0);
    QSr=complex(zeros(fftBins,1),0);%raw
    CSx=complex(zeros(fftBins,Nchx),0);
    nSkip=0;
    
    for k=nDelay+1:numRec
        %modulo indices
        mG=mod(k-1,nBlkAll)+1;
        m1G=mod(k-1-nhAll,nBlkAll)+1;
        %mxG=mod(k-2,nBlkAll)+1;
        mAll=mod(k-1,nDelay)+1;
        mlAll=mod(k-2,nDelay)+1;
        %mllAll=mod(k-3,nDelay)+1;
        m1All=mod(k-1-nhAll,nDelay)+1;
        m2All=mod(k,nDelay)+1;
        mSingle=mod(k-1,nhSingle+1)+1;
        m1Single=mod(k,nhSingle+1)+1;
        
        %kAll=k-nhAll;%index for which global filtered avg is becoming available
        %kSingle=k-nhAll-nhSingle;%index for which single filtered avg is becoming available
        xRaw0=double(h5read(RawFile,['/recordings/' num2str(iGroup) '/data'],[1,k*nBlkX-nsx+Filter.nStart],[Nch,nBlkX+2*nsx])');
        %test whether there are any jumps in voltage, corresponding to
        %recalibration events of the amplifiers
        skipBlk(end)=false;
        avgRaw0=mean(xRaw0,2);
        [stepMx,stepArg]=max(abs(avgRaw0(4:end)-avgRaw0(1:end-3)));
        %see whether surround is comparable
        if abs(stepArg+1.1-nBlkX/2-nsx)<=nBlkX/2
            wideAvg=abs(mean(avgRaw0(stepArg-28:stepArg+1))-mean(avgRaw0(stepArg+2:stepArg+31)));
            %similar magnituse as step?
            if wideAvg<1.5*stepMx
                if wideAvg>0.5*stepMx
                    if stepMx>stepThr
                        skipBlk(end)=true;
                        sepFrame(end)=stepArg+1-nsx;%between 0,nBlkX-1
                        stepTimes(end+1)=k*nBlkX+sepFrame(end);
                    end
                end
            end
        end
        
        %continue processing
        xRaw(:,end,:)=xRaw0(:,ChMap(ChMask));
        xCum(2:end,1,:)=cumsum(xRaw(:,end,:),1);
        xB=(xCum(2*nsx+2:end,1,:)-xCum(2:end-2*nsx,1,:))/(4*nsx)...
            +(xCum(nsx+nsxh+2:end-nsx+nsxh,1,:)-xCum(2+nsxh:end-2*nsx+nsxh,1,:))/(2*nsx);
        xSmth(1:2:end-1,end,:)=(xCum(nsx+3:end-nsx+1,1,:)-xCum(nsx:end-nsx-2,1,:))/3-xB;
        xSmth(2:2:end,end,:)=(xCum(nsx+3:end-nsx+1,1,:)-xCum(nsx+1:end-nsx-1,1,:))/2-xB;
        xCut=reshape(xSmth(:,end-1:end,:),4*nBlkX,Nchx);
        for j=0:2*deltaN
            xBlk=reshape(xCut(1:end-j*npy,:),2*nx+deltaN-j,npy,Nchx);
            fPow(end,j+1)=median(mean(abs(mean(xBlk,2)),1),3);
        end
        fMean=median(fPow,1)'+fScale*fAR;
        fAR=fA*fAR+fMean-fScale*fAR;%in case of small osc.
        fCM=conv(fMean,fPrior,'same');%insensitive to individual channel phase
        [~,fXarg]=max(fCM(deltaNh+1:end-deltaNh));%need to constrain for prior
        fXa(mAll)=fXarg;
        %create an updated prior
        fConv=conv(fMean,fP,'same');
        fConv=fConv(fXarg:fXarg+deltaN)-min(fConv(fXarg:fXarg+deltaN));
        fP0=(fConv+fConv(end:-1:1))*0.5/sum(fConv);
        fPriorLong=fA*fPriorLong+fScale*fP0;
        fPrior=0.1*fPriorLong+0.9*fP0;%need a separate, long-term component as a bias term
        
        
        %interpolate at 60 Hz grid
        xWarp=reshape(repmat(fracInd0{fXarg},1,Nchx).*xCut(Ind0{fXarg},:)...
            +repmat(1-fracInd0{fXarg},1,Nchx).*xCut(Ind0{fXarg}+1,:),nT,npBlk,Nchx);
        
        %median after normalization of amplitude
        xStd=std(xCut,0,1);
        if skipBlk(end-1)
            %use second half (first half affected by jump)
            xWmed=median(squeeze(median(xWarp(:,npBlkh+1:end,:),2))./repmat(xStd,nT,1),2);
        elseif skipBlk(end)
            %use first half (second half affected by jump)
            xWmed=median(squeeze(median(xWarp(:,1:npBlkh,:),2))./repmat(xStd,nT,1),2);
        else
            xWmed=median(squeeze(median(xWarp,2))./repmat(xStd,nT,1),2);
        end
        %xWBias=circshift(median(xOsc,2),-round((mean(fXa)+deltaNh-1)*npyh*nTh/nsx),1)/2;
        
                %determine phase
        %dPhase=round((mod(wPhase(m1All)-wPhase(m2All)+nT/2,nT)-nT/2)/4+(mod(hPhase-1-wPhase(m1All)+nT/2,nT)-nT/2)/2);
        %phasePrior=0.5*circshift(phasePrior,dPhase,2)+xWmed'*xSine;
        [~,hPhase]=max(xWmed'*xSine);
        %avoid jumps to nowhere? -- better:see which surrounding fits
        %better
%         if abs(mod(hPhase-1-wPhase(mlAll,1)-dtPhase+nTh,nT)-nTh)>2*deltaN
%             %make sure last phase wasn't an outlier
%             if abs(mod(hPhase-1-wPhase(mod(k-4,nDelay)+1,1)-dtPhase+nTh,nT)-nTh)>3*deltaN
%                 if nSkip<2
%                     nSkip=nSkip+1;
%                     hPhase=round(wPhase(mlAll,1)+dtPhase+1);%same phase
%                 else
%                     nSkip=0;
%                 end
%             end
%         else
%             nSkip=0;
%         end
        xOscOld=0.96*xOscOld+median(xOsc,2);
        if abs(mod(hPhase-1-wPhase(mlAll,1)-dtPhase+nTh,nT)-nTh)>2*deltaN
            oPow=zeros(4*deltaN+1,1);
            for j=1:2*deltaN+1
                oPow(j,1)=sum(circshift(xWmed,round(wPhase(mlAll,1)+dtPhase+1)-2*deltaN-2+j,1).*xOscOld,1);
            end
            [oMx1,dxShiftInd1]=max(oPow);
            oPow=zeros(4*deltaN+1,1);
            for j=1:2*deltaN+1
                oPow(j,1)=sum(circshift(xWmed,hPhase-2*deltaN-2+j,1).*xOscOld,1);
            end
            [oMx2,dxShiftInd]=max(oPow);
            if oMx1>oMx2
                dxShiftInd=dxShiftInd1;
                hPhase=round(wPhase(mlAll,1)+dtPhase+1);%same phase
            end
        else
            oPow=zeros(4*deltaN+1,1);
            for j=1:2*deltaN+1
                oPow(j,1)=sum(circshift(xWmed,hPhase-2*deltaN-2+j,1).*xOscOld,1);
            end
            [~,dxShiftInd]=max(oPow);
        end
        wPhase(mAll,1)=hPhase+dxShiftInd-2*deltaN-2;%phase for the current avg.
        if abs(dxShiftInd-2*deltaN-1)>2
            oPhase=round((dxShiftInd-2*deltaN-1)/(nBlkAll-1));
            xOsc=circshift(xOsc,-oPhase,1);
            xOscOld=circshift(xOscOld,-oPhase,1);
            sOsc=circshift(sOsc,-oPhase,1);
            wPhase=wPhase-oPhase;
        end
        %for plotting
        xOscAll=xOscAll+circshift(squeeze(median(xWarp,2)),wPhase(mAll,1),1);
        xPhase(k,1)=mod(hPhase+dxShiftInd-deltaNh-2-k*dtPhase,nT);%not integer
        %%global avg.
        %align
        xOsc(:,mG)=0.5*circshift(xWmed,wPhase(mAll,1),1)-0.5*circshift(xWmed,wPhase(mAll,1)+nTh,1);%aligned with sine wave

        %fix stuff with no data
        %median oscillation (should have a lag of (nBlkAll+1)/2 blocks here
        %avg. and add past phase
        xOmed=circshift(median(xOsc,2)+xOsc(:,m1G),-wPhase(m1All,1),1)/2;
        
        xBase=(fracInd1{fXa(m1All)}.*xOmed(Ind1{fXa(m1All)},1)+...
            (1-fracInd1{fXa(m1All)}).*xOmed(Ind2{fXa(m1All)},1))*ones(1,Nchx);%need to fix...
        xCutD=reshape(xSmth(:,end-nhAll-1:end-nhAll,:),4*nBlkX,Nchx);
        %add first half projection
        if skipBlk(end-nhAll-1)
            %first half affected, use second half
            xProj=repmat((sum(xCutD(2*nBlkX+1:4*nBlkX,:).*xBase(2*nBlkX+1:4*nBlkX,:),1)./...
                (sum(xBase(2*nBlkX+1:4*nBlkX,:).^2,1)+1)),4*nBlkX,1).*xBase;
            %can maybe avoid this distinction
            if sepFrame(end-nhAll-1)>0
                %second part new
                cAll(sepFrame(end-nhAll-1)+1:end,mSingle,:)=(xProj(2*sepFrame(end-nhAll-1):2:2*nBlkX-2,:)...
                    +2*xProj(2*sepFrame(end-nhAll-1)+1:2:2*nBlkX-1,:)+xProj(2*sepFrame(end-nhAll-1)+2:2:2*nBlkX,:))/4;
            else
                cAll(1:end,mSingle,:)=([xProj(2,:);xProj(2:2:2*nBlkX-2,:)]+2*xProj(1:2:2*nBlkX-1,:)+xProj(2:2:2*nBlkX,:))/4;
            end
        elseif skipBlk(end-nhAll)
            %second half affected, use first half
            xProj=repmat((sum(xCutD(1:2*nBlkX,:).*xBase(1:2*nBlkX,:),1)./...
                (sum(xBase(1:2*nBlkX,:).^2,1)+1)),4*nBlkX,1).*xBase;
            %only second part affected, do later
            cAll(:,mSingle,:)=squeeze(cAll(:,mSingle,:))+([xProjLast;xProj(2:2:2*nBlkX-2,:)]+2*xProj(1:2:2*nBlkX-1,:)+xProj(2:2:2*nBlkX,:))/8;
        else
            xProj=repmat((sum(xCutD(1:4*nBlkX,:).*xBase,1)./(sum(xBase.^2,1)+1)),4*nBlkX,1).*xBase;
            %smoothing to compensate for phase shift -- later
            %estimate global 60Hz
            cAll(:,mSingle,:)=squeeze(cAll(:,mSingle,:))+([xProjLast;xProj(2:2:2*nBlkX-2,:)]+2*xProj(1:2:2*nBlkX-1,:)+xProj(2:2:2*nBlkX,:))/8; 
        end
        %remainder
        xFilt(:,mSingle,:)=xCutD-xProj;%Do remember!
        %%single channel avg.
        %convert to warped time
        wFilt=reshape(repmat(fracInd0{fXa(m1All)},1,Nchx).*squeeze(xFilt(Ind0{fXa(m1All)},mSingle,:))...
            +repmat(1-fracInd0{fXa(m1All)},1,Nchx).*squeeze(xFilt(Ind0{fXa(m1All)}+1,mSingle,:)),nT,npBlk,Nchx);
        %correct for phase and avg.
        if skipBlk(end-nhAll-1)
            %first half affected, use second half
            sOh=circshift(median(wFilt(:,1:npBlkh,:),2),wPhase(m1All,1),1); 
        elseif skipBlk(end-nhAll)
            %second half affected, use first half
            sOh=circshift(median(wFilt(:,npBlkh+1:end,:),2),wPhase(m1All,1),1); 
        else
            sOh=circshift(median(wFilt,2),wPhase(m1All,1),1); 
        end
        sOsc(:,mSingle,:)=sOh-repmat(mean(sOh,1),nT,1,1);
        %for plotting
        xOscSgle=xOscSgle+circshift(squeeze(sOsc(:,mSingle,:)),wPhase(m1All,1),1);
        %projection on single traces
        %add past phase
        sOmed=circshift(squeeze(median(sOsc,2)),-wPhase(m2All,1),1);
        sBase=((fracInd1{fXa(m2All)}*ones(1,Nchx)).*sOmed(Ind1{fXa(m2All)},:)+...
            (1-fracInd1{fXa(m2All)}*ones(1,Nchx)).*sOmed(Ind2{fXa(m2All)},:));
        %add first half projection
        if skipBlk(end-nDelay)
            %first half affected, use cut first half, double
            sProj=repmat((sum(squeeze(xFilt(2*nBlkX+1:4*nBlkX,m1Single,:)).*sBase(2*nBlkX+1:4*nBlkX,:),1)./...
                (sum(sBase(2*nBlkX+1:4*nBlkX,:).^2,1)+1)),4*nBlkX,1).*sBase;
            if sepFrame(end-nDelay)>0
                %initial part untouched
                cSingle(sepFrame(end-nDelay)+1:end,:)=(sProj(2*sepFrame(end-nDelay):2:2*nBlkX-2,:)...
                    +2*sProj(2*sepFrame(end-nDelay)+1:2:2*nBlkX-1,:)+sProj(2*sepFrame(end-nDelay)+2:2:2*nBlkX,:))/4;
            else
                cSingle=([sProj(2,:);sProj(2:2:2*nBlkX-2,:)]+2*sProj(1:2:2*nBlkX-1,:)+sProj(2:2:2*nBlkX,:))/4;
            end
        elseif skipBlk(end-nhAll-nhSingle)
            %second half affected, use first half
            sProj=repmat((sum(squeeze(xFilt(1:2*nBlkX,m1Single,:)).*sBase(1:2*nBlkX,:),1)./...
                (sum(sBase(1:2*nBlkX,:).^2,1)+1)),4*nBlkX,1).*sBase;
            %only second part affected, do later
            cSingle=cSingle+([sProjLast;sProj(2:2:2*nBlkX-2,:)]+2*sProj(1:2:2*nBlkX-1,:)+sProj(2:2:2*nBlkX,:))/8;
        else
            sProj=repmat(sum(squeeze(xFilt(:,m1Single,:)).*sBase,1)./(sum(sBase.^2,1)+1),4*nBlkX,1).*sBase;
            cSingle=cSingle+([sProjLast;sProj(2:2:2*nBlkX-2,:)]+2*sProj(1:2:2*nBlkX-1,:)+sProj(2:2:2*nBlkX,:))/8;
        end
        h5write(OutFile,'/recordings/0/data',int16(round((squeeze(xRaw(nsx+1:nBlkX+nsx,end-nDelay,:))...
            -squeeze(cAll(:,m1Single,:))-cSingle)')),[1,(k-nDelay)*nBlkX+1],[Nchx,nBlkX]);
        %for plotting
        xRawV(k,:)=std(xCutD,0,1);
        xCV(k,:)=std(xProj,0,1);
        xSV(k,:)=std(sProj,0,1);
        xNV(k,:)=std(squeeze(xFilt(:,m1Single,:))-sProj,0,1);
        F=fft(complex(xCutD(1:Nfft,:)));
        Fc=fft(complex(xProj(1:Nfft,:)));
        Fs=fft(complex(sProj(1:Nfft,:)));
        Ftmp(:,mSingle,:)=F(1:fftBins,:);
        Fctmp(:,mSingle,:)=Fc(1:fftBins,:);
        QSr=QSr+mean(squeeze(Ftmp(:,m1Single,:).*conj(Ftmp(:,m1Single,:))),2);
        QSc=QSc+mean(squeeze(Fctmp(:,m1Single,:).*conj(Fctmp(:,m1Single,:))),2);
        QSs=QSs+mean(Fs(1:fftBins,:).*conj(Fs(1:fftBins,:)),2);
        QSx=(squeeze(Ftmp(:,m1Single,:)-Fctmp(:,m1Single,:))-Fs(1:fftBins,:)).*...
            conj((squeeze(Ftmp(:,m1Single,:)-Fctmp(:,m1Single,:))-Fs(1:fftBins,:)));
        CSx=CSx+QSx;
        QS(:,k)=mean(QSx,2);
        
        %remember stuff for next iteration
        %put second half projection
        if skipBlk(end-nhAll-nhSingle)
            %second half affected, cut second half, double
            if sepFrame(end-nhAll-nhSingle)>0
                cSingle(1:sepFrame(end-nhAll-nhSingle),:)=(sProj(2*nBlkX:2:2*(nBlkX+sepFrame(end-nhAll-nhSingle))-2,:)...
                    +2*sProj(2*nBlkX+1:2:2*(nBlkX+sepFrame(end-nhAll-nhSingle))-1,:)+sProj(2*nBlkX+2:2:2*(nBlkX+sepFrame(end-nhAll-nhSingle)),:))/4;
            end
        else
            cSingle=(sProj(2*nBlkX:2:end-2,:)+2*sProj(2*nBlkX+1:2:end-1,:)+sProj(2*nBlkX+2:2:end,:))/8; 
        end
        %put second half projection
        if skipBlk(end-nhAll)
            %second half affected, cut second half, double
            if sepFrame(end-nhAll)>0
                cAll(1:sepFrame(end-nhAll),mSingle,:)=(xProj(2*nBlkX:2:2*(nBlkX+sepFrame(end-nhAll))-2,:)...
                    +2*xProj(2*nBlkX+1:2:2*(nBlkX+sepFrame(end-nhAll))-1,:)+xProj(2*nBlkX+2:2:2*(nBlkX+sepFrame(end-nhAll)),:))/4;
            end
        else
            cAll(:,m1Single,:)=(xProj(2*nBlkX:2:end-2,:)+2*xProj(2*nBlkX+1:2:end-1,:)+xProj(2*nBlkX+2:2:end,:))/8;
        end
        sProjLast=sProj(2*nBlkX,:);
        xProjLast=xProj(2*nBlkX,:);
        fPow=circshift(fPow,-1,1);
        xRaw=circshift(xRaw,-1,2);
        xSmth=circshift(xSmth,-1,2);
        skipBlk=circshift(skipBlk,-1,1);
        sepFrame=circshift(sepFrame,-1,1);
    end
    %treat ending
    xTmp=double(h5read(RawFile,['/recordings/' num2str(iGroup) '/data'],[1,(numRec+1)*nBlkX+Filter.nStart],[Nch,Filter.LenRec-(numRec+1)*nBlkX])');
    for k=1:nhAll
        m1All=mod(k-1+numRec-nhAll,nDelay)+1;
        m2All=mod(k+numRec,nDelay)+1;
        mSingle=mod(k-1+numRec,nhSingle+1)+1;
        m1Single=mod(k+numRec,nhSingle+1)+1;
        xOmed=circshift(median(xOsc,2),-wPhase(m1All,1),1);
        
        xBase=(fracInd1{fXa(m1All)}.*xOmed(Ind1{fXa(m1All)},1)+...
            (1-fracInd1{fXa(m1All)}).*xOmed(Ind2{fXa(m1All)},1))*ones(1,Nchx);%need to fix...
        xCutD=reshape(xSmth(:,end-nhAll-1-k:end-nhAll-k,:),4*nBlkX,Nchx);
        if skipBlk(end-nhAll-1)
            %first half affected, use second half
            xProj=repmat((sum(xCutD(2*nBlkX+1:4*nBlkX,:).*xBase(2*nBlkX+1:4*nBlkX,:),1)./...
                (sum(xBase(2*nBlkX+1:4*nBlkX,:).^2,1)+1)),4*nBlkX,1).*xBase;
            %can maybe avoid this distinction
            if sepFrame(end-nhAll-1)>0
                %second part new
                cAll(sepFrame(end-nhAll-1)+1:end,mSingle,:)=(xProj(2*sepFrame(end-nhAll-1):2:2*nBlkX-2,:)...
                    +2*xProj(2*sepFrame(end-nhAll-1)+1:2:2*nBlkX-1,:)+xProj(2*sepFrame(end-nhAll-1)+2:2:2*nBlkX,:))/4;
            else
                cAll(1:end,mSingle,:)=([xProj(2,:);xProj(2:2:2*nBlkX-2,:)]+2*xProj(1:2:2*nBlkX-1,:)+xProj(2:2:2*nBlkX,:))/4;
            end
        elseif skipBlk(end-nhAll)
            %second half affected, use first half
            xProj=repmat((sum(xCutD(1:2*nBlkX,:).*xBase(1:2*nBlkX,:),1)./...
                (sum(xBase(1:2*nBlkX,:).^2,1)+1)),4*nBlkX,1).*xBase;
            %only second part affected, do later
            cAll(:,mSingle,:)=squeeze(cAll(:,mSingle,:))+([xProjLast;xProj(2:2:2*nBlkX-2,:)]+2*xProj(1:2:2*nBlkX-1,:)+xProj(2:2:2*nBlkX,:))/8;
        else
            xProj=repmat((sum(xCutD(1:4*nBlkX,:).*xBase,1)./(sum(xBase.^2,1)+1)),4*nBlkX,1).*xBase;
            %smoothing to compensate for phase shift -- later
            %estimate global 60Hz
            cAll(:,mSingle,:)=squeeze(cAll(:,mSingle,:))+([xProjLast;xProj(2:2:2*nBlkX-2,:)]+2*xProj(1:2:2*nBlkX-1,:)+xProj(2:2:2*nBlkX,:))/8; 
        end
        %xProj=repmat((sum(xCutD(1:4*nBlkX,:).*xBase,1)/(sum(xBase.^2,1)+1)),4*nBlkX,1).*xBase;
        %smoothing to compensate for phase shift -- later
        %estimate global 60Hz
        %cAll(:,mSingle,:)=squeeze(cAll(:,mSingle,:))+(xProj(2*nBlkX:2:end-2,:)+2*xProj(2*nBlkX+1:2:end-1,:)+xProj(2*nBlkX+2:2:end,:))/8;
        %remainder
        xFilt(:,mSingle,:)=xCutD-xProj;%Do remember!
        
        %%single channel avg.
        %convert to warped time
        wFilt=reshape(repmat(fracInd0{fXa(m1All)},1,Nchx).*squeeze(xFilt(Ind0{fXa(m1All)},mSingle,:))...
            +repmat(1-fracInd0{fXa(m1All)},1,Nchx).*squeeze(xFilt(Ind0{fXa(m1All)}+1,mSingle,:)),nT,npBlk,Nchx);
        %correct for phase and avg.
        if skipBlk(end-nhAll-1)
            %first half affected, use second half
            sOh=circshift(median(wFilt(:,1:npBlkh,:),2),wPhase(m1All,1),1); 
        elseif skipBlk(end-nhAll)
            %second half affected, use first half
            sOh=circshift(median(wFilt(:,npBlkh+1:end,:),2),wPhase(m1All,1),1); 
        else
            sOh=circshift(median(wFilt,2),wPhase(m1All,1),1); 
        end
        %sOh=circshift(median(wFilt,2),wPhase(m1All,1),1);
        sOsc(:,mSingle,:)=sOh-repmat(mean(sOh,1),nT,1,1);
        %for plotting
        xOscSgle=xOscSgle+squeeze(sOsc(:,mSingle,:));
        %projection on single traces
        %add past phase
        sOmed=circshift(squeeze(median(sOsc,2)),-wPhase(m2All,1),1);
        sBase=((fracInd1{fXa(m2All)}*ones(1,Nchx)).*sOmed(Ind1{fXa(m2All)},:)+...
        (1-fracInd1{fXa(m2All)}*ones(1,Nchx)).*sOmed(Ind2{fXa(m2All)},:));
        %add first half projection
        if skipBlk(end-nDelay)
            %first half affected, use cut first half, double
            sProj=repmat((sum(squeeze(xFilt(2*nBlkX+1:4*nBlkX,m1Single,:)).*sBase(2*nBlkX+1:4*nBlkX,:),1)./...
                (sum(sBase(2*nBlkX+1:4*nBlkX,:).^2,1)+1)),4*nBlkX,1).*sBase;
            if sepFrame(end-nDelay)>0
                %initial part untouched
                cSingle(sepFrame(end-nDelay)+1:end,:)=(sProj(2*sepFrame(end-nDelay):2:2*nBlkX-2,:)...
                    +2*sProj(2*sepFrame(end-nDelay)+1:2:2*nBlkX-1,:)+sProj(2*sepFrame(end-nDelay)+2:2:2*nBlkX,:))/4;
            else
                cSingle=([sProj(2,:);sProj(2:2:2*nBlkX-2,:)]+2*sProj(1:2:2*nBlkX-1,:)+sProj(2:2:2*nBlkX,:))/4;
            end
        elseif skipBlk(end-nhAll-nhSingle)
            %second half affected, use first half
            sProj=repmat((sum(squeeze(xFilt(1:2*nBlkX,m1Single,:)).*sBase(1:2*nBlkX,:),1)./...
                (sum(sBase(1:2*nBlkX,:).^2,1)+1)),4*nBlkX,1).*sBase;
            %only second part affected, do later
            cSingle=cSingle+([sProjLast;sProj(2:2:2*nBlkX-2,:)]+2*sProj(1:2:2*nBlkX-1,:)+sProj(2:2:2*nBlkX,:))/8;
        else
            sProj=repmat(sum(squeeze(xFilt(:,m1Single,:)).*sBase,1)./(sum(sBase.^2,1)+1),4*nBlkX,1).*sBase;
            cSingle=cSingle+([sProjLast;sProj(2:2:2*nBlkX-2,:)]+2*sProj(1:2:2*nBlkX-1,:)+sProj(2:2:2*nBlkX,:))/8;
        end
        %sProj=repmat(sum(squeeze(xFilt(:,m1Single,:)).*sBase,1)./(sum(sBase.^2,1)+1),4*nBlkX,1).*sBase;
        %cSingle=cSingle+(sProj(2*nBlkX:2:end-2,:)+2*sProj(2*nBlkX+1:2:end-1,:)+sProj(2*nBlkX+2:2:end,:))/8;
        h5write(OutFile,'/recordings/0/data',int16(round((squeeze(xRaw(nsx+1:nBlkX+nsx,end-nDelay+k,:))...
        -squeeze(cAll(:,m1Single,:))-cSingle)')),[1,(k+numRec-nDelay)*nBlkX+1],[Nchx,nBlkX]);
        %put second half projection
        if skipBlk(end-nhAll-nhSingle)
            %second half affected, cut second half, double
            if sepFrame(end-nhAll-nhSingle)>0
                cSingle(1:sepFrame(end-nhAll-nhSingle),:)=(sProj(2*nBlkX:2:2*(nBlkX+sepFrame(end-nhAll-nhSingle))-2,:)...
                    +2*sProj(2*nBlkX+1:2:2*(nBlkX+sepFrame(end-nhAll-nhSingle))-1,:)+sProj(2*nBlkX+2:2:2*(nBlkX+sepFrame(end-nhAll-nhSingle)),:))/4;
            end
        else
            cSingle=(sProj(2*nBlkX:2:end-2,:)+2*sProj(2*nBlkX+1:2:end-1,:)+sProj(2*nBlkX+2:2:end,:))/8; 
        end
        %put second half projection
        if skipBlk(end-nhAll) && ~(k==nhAll)
            %second half affected, cut second half, double
            if sepFrame(end-nhAll)>0
                cAll(1:sepFrame(end-nhAll),mSingle,:)=(xProj(2*nBlkX:2:2*(nBlkX+sepFrame(end-nhAll))-2,:)...
                    +2*xProj(2*nBlkX+1:2:2*(nBlkX+sepFrame(end-nhAll))-1,:)+xProj(2*nBlkX+2:2:2*(nBlkX+sepFrame(end-nhAll)),:))/4;
            end
        else
            cAll(:,m1Single,:)=(xProj(2*nBlkX:2:end-2,:)+2*xProj(2*nBlkX+1:2:end-1,:)+xProj(2*nBlkX+2:2:end,:))/8;
        end
        sProjLast=sProj(2*nBlkX,:);
        xProjLast=xProj(2*nBlkX,:);
        skipBlk=circshift(skipBlk,-1,1);
        sepFrame=circshift(sepFrame,-1,1);
        %cSingle=([sProjLast;sProj(2:2:2*nBlkX-2,:)]+2*sProj(1:2:2*nBlkX-1,:)+sProj(2:2:2*nBlkX,:))/8;
        %sProjLast=sProj(1,:);
        %cAll(:,m1Single,:)=([xProjLast;xProj(2:2:2*nBlkX-2,:)]+2*xProj(1:2:2*nBlkX-1,:)+xProj(2:2:2*nBlkX,:))/8;
        %xProjLast=xProj(2*nBlkX,:);%%%later
    end
    cAll(:,m1Single,:)=2*cAll(:,m1Single,:);
    for k=nhAll:(nhAll+nhSingle)
        m2All=mod(k+numRec,nDelay)+1;
        m1Single=mod(k+numRec,nhSingle+1)+1;
        %add past phase
        sOmed=circshift(squeeze(median(sOsc,2)),-wPhase(m2All,1),1);
        sBase=((fracInd1{fXa(m2All)}*ones(1,Nchx)).*sOmed(Ind1{fXa(m2All)},:)+...
        (1-fracInd1{fXa(m2All)}*ones(1,Nchx)).*sOmed(Ind2{fXa(m2All)},:));
        %add first half projection
        if skipBlk(end-nDelay)
            %first half affected, use cut first half, double
            sProj=repmat((sum(squeeze(xFilt(2*nBlkX+1:4*nBlkX,m1Single,:)).*sBase(2*nBlkX+1:4*nBlkX,:),1)./...
                (sum(sBase(2*nBlkX+1:4*nBlkX,:).^2,1)+1)),4*nBlkX,1).*sBase;
            if sepFrame(end-nDelay)>0
                %initial part untouched
                cSingle(sepFrame(end-nDelay)+1:end,:)=(sProj(2*sepFrame(end-nDelay):2:2*nBlkX-2,:)...
                    +2*sProj(2*sepFrame(end-nDelay)+1:2:2*nBlkX-1,:)+sProj(2*sepFrame(end-nDelay)+2:2:2*nBlkX,:))/4;
            else
                cSingle=([sProj(2,:);sProj(2:2:2*nBlkX-2,:)]+2*sProj(1:2:2*nBlkX-1,:)+sProj(2:2:2*nBlkX,:))/4;
            end
        elseif skipBlk(end-nhAll-nhSingle)
            %second half affected, use first half
            sProj=repmat((sum(squeeze(xFilt(1:2*nBlkX,m1Single,:)).*sBase(1:2*nBlkX,:),1)./...
                (sum(sBase(1:2*nBlkX,:).^2,1)+1)),4*nBlkX,1).*sBase;
            %only second part affected, do later
            cSingle=cSingle+([sProjLast;sProj(2:2:2*nBlkX-2,:)]+2*sProj(1:2:2*nBlkX-1,:)+sProj(2:2:2*nBlkX,:))/8;
        else
            sProj=repmat(sum(squeeze(xFilt(:,m1Single,:)).*sBase,1)./(sum(sBase.^2,1)+1),4*nBlkX,1).*sBase;
            cSingle=cSingle+([sProjLast;sProj(2:2:2*nBlkX-2,:)]+2*sProj(1:2:2*nBlkX-1,:)+sProj(2:2:2*nBlkX,:))/8;
        end
        %sProj=repmat(sum(squeeze(xFilt(:,m1Single,:)).*sBase,1)./(sum(sBase.^2,1)+1),4*nBlkX,1).*sBase;
        %cSingle=cSingle+(sProj(2*nBlkX:2:end-2,:)+2*sProj(2*nBlkX+1:2:end-1,:)+sProj(2*nBlkX+2:2:end,:))/8;
        h5write(OutFile,'/recordings/0/data',int16(round((squeeze(xRaw(nsx+1:nBlkX+nsx,end-nDelay+k,:))...
        -squeeze(cAll(:,m1Single,:))-cSingle)')),[1,(k+numRec-nDelay)*nBlkX+1],[Nchx,nBlkX]);
        if skipBlk(end-nhAll-nhSingle) && ~(k==nhAll+nhSingle)
            %second half affected, cut second half, double
            if sepFrame(end-nhAll-nhSingle)>0
                cSingle(1:sepFrame(end-nhAll-nhSingle),:)=(sProj(2*nBlkX:2:2*(nBlkX+sepFrame(end-nhAll-nhSingle))-2,:)...
                    +2*sProj(2*nBlkX+1:2:2*(nBlkX+sepFrame(end-nhAll-nhSingle))-1,:)+sProj(2*nBlkX+2:2:2*(nBlkX+sepFrame(end-nhAll-nhSingle)),:))/4;
            end
        else
            cSingle=(sProj(2*nBlkX:2:end-2,:)+2*sProj(2*nBlkX+1:2:end-1,:)+sProj(2*nBlkX+2:2:end,:))/8; 
        end
        sProjLast=sProj(2*nBlkX,:);
        skipBlk=circshift(skipBlk,-1,1);
        sepFrame=circshift(sepFrame,-1,1);
        %cSingle=([sProjLast;sProj(2:2:2*nBlkX-2,:)]+2*sProj(1:2:2*nBlkX-1,:)+sProj(2:2:2*nBlkX,:))/8;
        %sProjLast=sProj(1,:);
    end
    cSingle=2*cSingle;
    m1Single=mod(nDelay+numRec,nhSingle+1)+1;
    %add past phase
    h5write(OutFile,'/recordings/0/data',int16(round((squeeze(xRaw(nsx+1:nBlkX+nsx,end,:))...
    -squeeze(cAll(:,m1Single,:))-cSingle)')),[1,numRec*nBlkX+1],[Nchx,nBlkX]);
    newPhase=(deltaN-fXarg)*npyh;
    newC=circshift(squeeze(cAll(:,m1Single,:))+cSingle,newPhase,1);
    if Filter.LenRec-(numRec+1)*nBlkX>nBlkX
        h5write(OutFile,'/recordings/0/data',int16(round((xTmp(1:nBlkX,ChMap(ChMask))-newC)'))...
            ,[1,(numRec+1)*nBlkX+1],[Nchx,nBlkX]);
        newPhase=2*(deltaN-fXarg)*npyh;
        newC=circshift(squeeze(cAll(:,m1Single,:))+cSingle,newPhase,1);
        h5write(OutFile,'/recordings/0/data',int16(round((xTmp(nBlkX+1:end,ChMap(ChMask))-newC(1:Filter.LenRec-(numRec+2)*nBlkX,:))'))...
            ,[1,(numRec+2)*nBlkX+1],[Nchx,Filter.LenRec-(numRec+2)*nBlkX]);
    else
        h5write(OutFile,'/recordings/0/data',int16(round((xTmp(:,ChMap(ChMask))-newC(1:Filter.LenRec-(numRec+1)*nBlkX,:))'))...
        ,[1,(numRec+1)*nBlkX+1],[Nchx,Filter.LenRec-(numRec+1)*nBlkX]);
    end
    

    %save parameters -- see which ones necessary and add to struct (Filter)
    Filter.pwl.Nch=Nch;
    Filter.pwl.Nchx=Nchx;
    Filter.pwl.ChMask=ChMask;
    Filter.pwl.ChMap=ChMap;
    Filter.pwl.npx=npx;
    Filter.pwl.npy=npy;
    Filter.pwl.nT=nT;
    Filter.pwl.Nfft=Nfft;
    Filter.pwl.fftBins=fftBins;
    Filter.pwl.fTau=fTau;
    Filter.pwl.nhAll=nhAll;
    Filter.pwl.nhSingle=nhSingle;
    Filter.pwl.deltaNh=deltaNh;
    Filter.pwl.sRate=sRate(1,1);
    Filter.pwl.bitVolt=bitVolt(1,1);
    
    Filter.pwl.plot.xOscAll=xOscAll'/(numRec-nDelay)*bitVolt(1,1);%mV
    Filter.pwl.plot.xOscSgle=xOscSgle'/(numRec-nDelay)*bitVolt(1,1);%mV
    Filter.pwl.plot.timeMinutes=(nDelay+1:numRec)'*nBlkX/sRate/60;
    Filter.pwl.plot.rawStdev=mean(xRawV(nDelay+1:numRec,:),2)*bitVolt(1,1);
    Filter.pwl.plot.cAllStdev=mean(xCV(nDelay+1:numRec,:),2)*bitVolt(1,1);
    Filter.pwl.plot.cSingleStdev=mean(xSV(nDelay+1:numRec,:),2)*bitVolt(1,1);
    Filter.pwl.plot.newStdev=mean(xNV(nDelay+1:numRec,:),2)*bitVolt(1,1);
    Filter.pwl.plot.Phase=xPhase(nDelay+1:numRec,1)*100/(6*nT);%ms
    Filter.pwl.plot.rawPower=abs(QSr/(numRec-nDelay)*2000*bitVolt(1,1)^2/(sRate(1,1)*Nfft));
    Filter.pwl.plot.cAllPower=abs(QSc/(numRec-nDelay)*2000*bitVolt(1,1)^2/(sRate(1,1)*Nfft));
    Filter.pwl.plot.cSinglePower=abs(QSs/(numRec-nDelay)*2000*bitVolt(1,1)^2/(sRate(1,1)*Nfft));
    Filter.pwl.plot.newPower=mean(abs(QS(:,nDelay+1:end)),2)*2000*bitVolt(1,1)^2/(sRate(1,1)*Nfft);
    Filter.pwl.plot.powerScale=(1:fftBins)*sRate(1,1)/(Nfft/2*1000);
    Filter.pwl.plot.stepTimes=stepTimes(2:end);
    
    if Filter.pwl.plotResults
        %plotting
        stepsPerMin=60*sRate(1,1)/nBlkX;
        
        fig1=figure('Position',[0 0 1600 1200]);
        ax1 = axes('OuterPosition',[0   3/4 1/3 1/4]);
        ax2 = axes('OuterPosition',[0   1/2 1/3 1/4]);
        ax3 = axes('OuterPosition',[0   1/4 1/3 1/4]);
        ax4 = axes('OuterPosition',[1/3 3/4 1/3 1/4]);
        ax5 = axes('OuterPosition',[1/3 2/4 1/3 1/4]);
        ax6 = axes('OuterPosition',[1/3 1/4 2/3 1/4]);
        ax4a = axes('OuterPosition',[2/3 3/4 1/3 1/4]);
        ax5a = axes('OuterPosition',[2/3 2/4 1/3 1/4]);
        axQSa= axes('OuterPosition',[0   0   1/4 1/4]);
        axQSb= axes('OuterPosition',[1/4 0   1/4 1/4]);
        axQSc= axes('OuterPosition',[2/4 0   1/4 1/4]);
        axQSd= axes('OuterPosition',[3/4 0   1/4 1/4]);
        
        plot(ax1,(1:nT)*100/(6*nT),mean(xOscSgle,2)/(numRec-nDelay)*bitVolt(1,1),'c-')
        hold(ax1,'on')
        plot(ax1,(1:nT)*100/(6*nT),mean(xOscSgle,2)*10/(numRec-nDelay)*bitVolt(1,1),'c:')
        plot(ax1,(1:nT)*100/(6*nT),mean(xOscAll,2)/(numRec-nDelay)*bitVolt(1,1),'b-')
        xlabel(ax1,'time/ms')
        ylabel(ax1,'µV')
        title(ax1,'average')
        xlim(ax1,[0 100/6])
        
        imagesc(ax2,[0 100/6],[1 Nchx],Filter.pwl.plot.xOscAll)
        xlabel(ax2,'time/ms')
        ylabel(ax2,'channel id')
        title(ax2,'avg. correction')
        hC=colorbar(ax2);
        hC.Label.String = 'µV';
        
        
        imagesc(ax3,[0 100/6],[1 Nchx],Filter.pwl.plot.xOscSgle)
        xlabel(ax3,'time/ms')
        ylabel(ax3,'channel id')
        title(ax3,'single correction')
        hC=colorbar(ax3);
        hC.Label.String = 'µV';
        
        xtEV=1-(xNV./xRawV).^2;
        xtEVm=median(median(xtEV,2),1);
        imagesc(ax4,[nDelay+1 nDelay+500]/stepsPerMin,[1 Nchx],...
            xtEV(nDelay+1:nDelay+500,:)','AlphaData',xtEV(nDelay+1:nDelay+500,:)'>0,[0 min(1,2*xtEVm)])
        xlabel(ax4,'time/min')
        ylabel(ax4,'channel id')
        %set(ax4,'Xtick',[1*stepsPerMin 2*stepsPerMin 3*stepsPerMin]+nDelay,'XTicklabel',{'1 min' '2 min' '3 min'});
        title(ax4,'fraction of explained variance')
        hC=colorbar('peer',ax4);
        hC.Label.String = 'explained variance';
        set(hC,'Ytick',0:1,'YTicklabel',{'0' '1'});
        
        imagesc(ax4a,[nDelay+1 numRec]/stepsPerMin,[1 Nchx],...
            xtEV(nDelay+1:end,:)','AlphaData',xtEV(nDelay+1:end,:)'>0,[0 min(1,2*xtEVm)])
        xlabel(ax4a,'time/min')
        ylabel(ax4a,'channel id')
        %set(ax4a,'Xtick',[10*stepsPerMin 20*stepsPerMin 30*stepsPerMin]+nDelay,'XTicklabel',{'10 min' '20 min' '30 min'});
        %xlim(ax4a,[0 5000])
        title(ax4a,'fraction of explained variance')
        hC=colorbar(ax4a);
        hC.Label.String = 'explained variance';
        set(hC,'Ytick',0:1,'YTicklabel',{'0' '1'});
        
        plot(ax5,(nDelay+1:nDelay+500)'*nBlkX/sRate/60,mean(xRawV(nDelay+1:nDelay+500,:),2)*bitVolt(1,1),'k-')%seconds
        hold(ax5,'on')
        plot(ax5,(nDelay+1:nDelay+500)'*nBlkX/sRate/60,mean(xCV(nDelay+1:nDelay+500,:),2)*bitVolt(1,1),'b-')%seconds
        plot(ax5,(nDelay+1:nDelay+500)'*nBlkX/sRate/60,mean(xSV(nDelay+1:nDelay+500,:),2)*bitVolt(1,1),'c-')%seconds
        plot(ax5,(nDelay+1:nDelay+500)'*nBlkX/sRate/60,mean(xNV(nDelay+1:nDelay+500,:),2)*bitVolt(1,1),'r-')%seconds
        xlabel(ax5,'time/min')
        %set(ax5,'Xtick',[1 2 3],'XTicklabel',{'1 min' '2 min' '3 min'});
        ylabel(ax5,'standard deviation/µV')
        xlim(ax5,[nDelay+1 nDelay+500]/stepsPerMin)
        
        plot(ax5a,Filter.pwl.plot.timeMinutes,Filter.pwl.plot.rawStdev,'k-')%seconds
        hold(ax5a,'on')
        plot(ax5a,Filter.pwl.plot.timeMinutes,Filter.pwl.plot.cAllStdev,'b-')%seconds
        plot(ax5a,Filter.pwl.plot.timeMinutes,Filter.pwl.plot.cSingleStdev,'c-')%seconds
        plot(ax5a,Filter.pwl.plot.timeMinutes,Filter.pwl.plot.newStdev,'r-')%seconds
        xlabel(ax5a,'time/min')
        %set(ax5a,'Xtick',[10 20 30],'XTicklabel',{'10 min' '20 min' '30 min'});
        ylabel(ax5a,'standard deviation/µV')
        xlim(ax5a,[nDelay+1 numRec]/stepsPerMin)
        
        plot(ax6,[1 1]*Filter.pwl.plot.stepTimes(1)/(60*sRate(1,1)),[17 18],'m-')
        hold(ax6,'on')
        for i=2:length(Filter.pwl.plot.stepTimes)
            plot(ax6,[1 1]*Filter.pwl.plot.stepTimes(i)/(60*sRate(1,1)),[17 18],'m-')
        end
        plot(ax6,Filter.pwl.plot.timeMinutes,Filter.pwl.plot.Phase,'k.','MarkerSize',1)%ok
        xlabel(ax6,'time/min')
        ylabel(ax6,'phase/ms')
        %set(ax6,'Xtick',[10 20 30],'XTicklabel',{'10 min' '20 min' '30 min'});
        ylim(ax6,[0 18])
        xlim(ax6,[nDelay+1 numRec]/stepsPerMin)
        
        semilogy(axQSa,Filter.pwl.plot.powerScale,Filter.pwl.plot.rawPower,'k-')
        hold(axQSa,'on')
        semilogy(axQSa,Filter.pwl.plot.powerScale,Filter.pwl.plot.cAllPower,'b-')
        semilogy(axQSa,Filter.pwl.plot.powerScale,Filter.pwl.plot.cSinglePower,'c-')
        semilogy(axQSa,Filter.pwl.plot.powerScale,Filter.pwl.plot.newPower,'r-')
        axQSa.YScale='log';%bug in MATLAB
        xlim(axQSa,[0 2])
        ylim(axQSa,[1 1e8])
        xlabel(axQSa,'frequency/kHz')
        ylabel(axQSa,'power/µV^2/kHz')
        title(axQSa,'power spectrum')
        
        imagesc(axQSb,[0.001 fftBins/1000]*sRate(1,1)/(Nfft/2),[1 Nchx],log10(CSx'*2000*bitVolt(1,1)^2/(sRate(1,1)*Nfft*(numRec-nDelay))))
        hC=colorbar(axQSb,'Ticks',[2 4 6],...
            'TickLabels',{'10','100','1000'});
        hC.Label.String = 'voltage rms/µV/kHz';
        xlabel(axQSb,'frequency/kHz')
        ylabel(axQSb,'channel id')
        title(axQSb,'power avg. (corrected)')
        
        imagesc(axQSc,[nDelay+1 nDelay+500]*nBlkX/sRate/60,[0.001 fftBins/1000]*sRate(1,1)/(Nfft/2),log10(QS(:,nDelay+1:nDelay+500)...
            *2000*bitVolt(1,1)^2/(sRate(1,1)*Nfft)))
        hC=colorbar(axQSc,'Ticks',[2 4 6],...
            'TickLabels',{'10','100','1000'});
        hC.Label.String = 'voltage rms/µV/kHz';
        xlabel(axQSc,'time/min')
        ylabel(axQSc,'frequency/kHz')
        title(axQSc,'power (corrected)')
        
        imagesc(axQSd,[nDelay+1 numRec]*nBlkX/sRate/60,[0.001 fftBins/1000]*sRate(1,1)/(Nfft/2),log10(QS(:,nDelay+1:end)...
            *2000*bitVolt(1,1)^2/(sRate(1,1)*Nfft)))
        hC=colorbar(axQSd,'Ticks',[2 4 6],...
            'TickLabels',{'10','100','1000'});
        hC.Label.String = 'voltage rms/µV/kHz';
        xlabel(axQSd,'time/min')
        ylabel(axQSd,'frequency/kHz')
        title(axQSd,'power (corrected)')
        
        saveas(fig1,[Filter.pwl.plotFolder filesep Filter.pwl.plotName])
    close(fig1)
    end
    if ~Filter.pwl.returnPlotData
        Filter.pwl.plot=[];
    end    
end
  