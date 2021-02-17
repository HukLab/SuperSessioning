classdef filterBase
    properties
        nChBoard
        nBoards
        nCycle
        dCycle
        nStepSmth
        StepSmth
        StepWeights
        nOsc
        nAvgh
        nAvgF
        PhaseTau
        nPhaseTau
        nSmth
        SmthKer
        nOverlap
        nBlk
        countBlk
        recParameter
        templateSine
        plots
%         Nfft
%         dtCS
%         f0
%         f1
%         fb0
%         fb1
%         CS
%         PS
%         QS
%         PHch
%         PH
%         QH
%         PHdc
%         phiN
    end
    methods
        function obj = filterBase(recParameter,nBlk)
            obj.recParameter=recParameter;
            if obj.recParameter.Nch==64
                obj.nChBoard=32;
                obj.nBoards=2;
            else
                obj.nChBoard=obj.recParameter.Nch;
                obj.nBoards=1;
            end
            obj.nCycle=500;%assume even
            obj.dCycle=obj.nCycle/10;
            %class properties
            obj.nStepSmth=60;
            obj.StepSmth=hanning(2*obj.nStepSmth+1)/(obj.nStepSmth+1);
            obj.StepWeights=[1:obj.nCycle obj.nCycle-1:-1:1]'/obj.nCycle;%relative weights for smoothed version
            obj.nOsc=12;%200ms,6k frames,even!
            obj.nAvgh=150;%to determine kernel
            obj.nAvgF=30*obj.nCycle;%to determine match with raw trace,even
            obj.plots.PHaseTau=300;%timescale for determining powerline phase (#osc.)
            obj.nPhaseTau=2*obj.plots.PHaseTau;
            obj.nSmth=2;%even!
            obj.SmthKer=hamming(obj.nSmth*obj.nCycle-1);
            obj.SmthKer=[obj.SmthKer(1:obj.nSmth/2*obj.nCycle,1); ...
                obj.SmthKer(obj.nSmth/2*obj.nCycle:end,1)]/(sum(obj.SmthKer)+1);
            obj.nOverlap=300*obj.nCycle;%add nPhaseTau*obj.nCycle for phase estimate (CA only)
            obj.nBlk=nBlk;%3e6? (100 sec of recording)
            %scale sine wave with logistic function (want to estimate phase!)
            obj.templateSine=reshape(2./(1+exp(-3*sin((0:obj.nOsc*obj.nCycle*(obj.nCycle+1)-1)*...
                2*pi/obj.nCycle)))-1,(obj.nCycle+1),obj.nOsc*obj.nCycle);
            obj.templateSine=obj.templateSine(1:end-1,:)./sum(obj.templateSine(1:end-1,:).^2,2);
            if obj.recParameter.LenRec<=obj.nBlk
                obj.nBlk=obj.recParameter.LenRec;
                obj.countBlk=1;
            else
                obj.countBlk=ceil((obj.recParameter.LenRec-3*obj.nOverlap)/...
                    (obj.nBlk-3*obj.nOverlap));
            end
            %plotting stuff
            obj.plots.dtCS=2^14;
            obj.plots.Nfft=floor((obj.nBlk-2*obj.nOverlap)./obj.plots.dtCS);
            %frequency band for visualization
            obj.plots.f0=500;%Hz
            obj.plots.f1=2000;%Hz
            obj.plots.fb0=floor(obj.plots.dtCS*obj.plots.f0/obj.recParameter.sRate);%bin index in spectrum
            obj.plots.fb1=floor(obj.plots.dtCS*obj.plots.f1/obj.recParameter.sRate);
            obj.plots.CS={};
            obj.plots.PS=complex(zeros(obj.plots.dtCS,obj.nBoards*2),0);%power spectrum (ignore zero frequency component)
            obj.plots.QS=complex(zeros(obj.plots.dtCS,obj.nBoards*2),0);
            obj.plots.PH=zeros(200,obj.nBoards*2);
            obj.plots.QH=zeros(200,obj.nBoards*2);
            obj.plots.PHch={};
            for iii=1:obj.nBoards*2
                obj.plots.CS{iii}=complex(zeros(obj.nChBoard,obj.nChBoard),0);%coherence spectrum
                obj.plots.PHch{iii}=((iii-1)*obj.nChBoard:iii*obj.nChBoard-1)'*200;
            end
            obj.plots.PHdc=zeros(200*obj.nChBoard*obj.nBoards*2,1);
            obj.plots.PHiN=2000*obj.recParameter.bitVolt^2/(sum(hanning(obj.plots.dtCS).^2)*obj.recParameter.sRate*(obj.plots.fb1-obj.plots.fb0));
        end
        function obj=filterAll(obj, RawFile, OutFile)
            obj.plots.cPhaseAll=zeros(floor(obj.recParameter.LenRec/obj.nCycle)+2,2);
            obj.plots.medPwl=zeros(obj.nCycle,floor(obj.recParameter.LenRec/(60*obj.nCycle))+1,2);
            obj.plots.iqrAvg=zeros(floor(obj.recParameter.LenRec/obj.dCycle)+1,2);
            obj.plots.projWeightCA=zeros(floor(obj.recParameter.LenRec/obj.nCycle),2);
            for Board=1:obj.nBoards
                %initial segment
                x=double(h5read(RawFile,obj.recParameter.HdfRawDataPath,...
                    [1+(Board-1)*obj.nChBoard,1],[obj.nChBoard,obj.nBlk])');
                %reduce voltage steps in individual traces
                z=obj.reduceRecalibrationSteps(x);
                %use smoothed (hamming kernel, 2 osc.) version subtracted common avg. signal
                zSmth0=conv2(z,obj.SmthKer,'valid');
                zSmth=[repmat(zSmth0(1,:),obj.nSmth/2*obj.nCycle,1); ...
                    zSmth0; repmat(zSmth0(end,:),obj.nSmth/2*obj.nCycle-1,1)];
                filteredZ=z-zSmth;
                obj.plots=obj.PSD(obj.plots,filteredZ,Board);
                %next, compute common average
                CommonAverage=reshape(mean(median(reshape(z,[],obj.nChBoard/8,8),3),2),[],1); %better:group into 4x8 channels, then median, mean
                %need to subtract a smoothed version
                CASmth0=conv(CommonAverage,obj.SmthKer,'valid');
                CASmth=[CASmth0(1,1)*ones(obj.nSmth/2*obj.nCycle,1); ...
                    CASmth0; CASmth0(end,1)*ones(obj.nSmth/2*obj.nCycle-1,1)];
                filteredCA=CommonAverage-CASmth;
                %now, before I do common average referencing, I should try to estimate 60 Hz
                %phase and noise.
                continuousPhase=obj.estimatePhase(filteredCA);
                %plotting
                obj.plots.medPwl0=obj.medianPwl(filteredCA,continuousPhase);
                obj.plots.medPwl(:,1:size(obj.plots.medPwl0,2),Board)=obj.plots.medPwl0;
%                 lenCA=floor(size(filteredCA,1)/obj.nCycle)*obj.nCycle;
%                 cPhaseTemp=interp1(0.5+(0:floor(lenCA/obj.nCycle)+1)*obj.nCycle,continuousPhase,...
%                 1:lenCA)-obj.nCycle;
%                 chWarped=reshape(interp1(cPhaseTemp(1,1:end)', ...
%                     filteredCA(1:floor(lenCA/obj.nCycle)*obj.nCycle,1)',...
%                     1:floor(lenCA/obj.nCycle)*obj.nCycle),obj.nCycle,[]);
%                 imagesc(min(max(reshape(chWarped,obj.nCycle,[])*obj.recParameter.bitVolt,-100),100))
                %
                obj.plots.cPhaseAll(1:floor(obj.nBlk/obj.nCycle)+2,Board)=continuousPhase;%for plotting
                [nfilteredZ, filterWeights]=obj.CAR(filteredCA,filteredZ);
                obj.plots.projWeightCA(1:floor(obj.nBlk/obj.nCycle),Board)=filterWeights;
                nfilteredPwl=obj.subtractResidualPowerline(nfilteredZ,continuousPhase);
                obj.plots.iqrAvg(1:floor(obj.nBlk/obj.dCycle),Board)=...
                    sqrt(mean(reshape(iqr(nfilteredPwl(1:floor(obj.nBlk/obj.dCycle)*obj.dCycle,:),2).^2,...
                    obj.dCycle,[]),1));
                obj.plots=obj.PSD(obj.plots,nfilteredPwl,Board+obj.nBoards);
                %plotting
                %imagesc(min(max(nfilteredZ'*obj.recParameter.bitVolt,-100),100))%not needed.
                %if this is the only segment, save data, else only save the
                %first segment and continue processing
                if obj.countBlk==1
                    h5write(OutFile,'/data',int16(round(nfilteredPwl')),...
                        [1+(Board-1)*obj.nChBoard,1],[obj.nChBoard,obj.nBlk]);
                else
                    h5write(OutFile,'/data',int16(round(...
                        nfilteredPwl(1:end-2*obj.nOverlap,:)')),...
                        [1+(Board-1)*obj.nChBoard,1],...
                        [obj.nChBoard,obj.nBlk-2*obj.nOverlap]);
                    %further iterations
                    for k=2:obj.countBlk-1
                        x=double(h5read(RawFile,obj.recParameter.HdfRawDataPath,...
                            [1+(Board-1)*obj.nChBoard,1+(k-1)*(obj.nBlk-3*obj.nOverlap)],...
                            [obj.nChBoard,obj.nBlk])');
                        %reduce voltage steps in individual traces
                        z=obj.reduceRecalibrationSteps(x);
                        %use smoothed (hamming kernel, 2 osc.) version subtracted common avg. signal
                        zSmth0=conv2(z,obj.SmthKer,'valid');
                        zSmth=[repmat(zSmth0(1,:),obj.nSmth/2*obj.nCycle,1); ...
                            zSmth0; repmat(zSmth0(end,:),obj.nSmth/2*obj.nCycle-1,1)];
                        filteredZ=z-zSmth;
                        obj.plots=obj.PSD(obj.plots,filteredZ,Board);
                        %next, compute common average
                        CommonAverage=reshape(mean(median(reshape(z,[],obj.nChBoard/8,8),3),2),[],1); %better:group into 4x8 channels, then median, mean
                        %need to subtract a smoothed version
                        CASmth0=conv(CommonAverage,obj.SmthKer,'valid');
                        CASmth=[CASmth0(1,1)*ones(obj.nSmth/2*obj.nCycle,1); ...
                            CASmth0; CASmth0(end,1)*ones(obj.nSmth/2*obj.nCycle-1,1)];
                        filteredCA=CommonAverage-CASmth;
                        %now, before I do common average referencing, I should try to estimate 60 Hz
                        %phase and noise.
                        continuousPhase=obj.estimatePhase(filteredCA);
                        k0=(k-1)*(obj.nBlk-3*obj.nOverlap)/obj.nCycle;
                        k1=k0+obj.nOverlap/obj.nCycle;%
                        k2=k1+obj.nOverlap/obj.nCycle;
                        k3=floor(k1/60);
                        k4=floor(obj.nOverlap/(60*obj.nCycle));
                        obj.plots.medPwl0=obj.medianPwl(filteredCA,continuousPhase);
                        obj.plots.medPwl(:,k3+1:k3-k4+size(obj.plots.medPwl0,2),Board)=obj.plots.medPwl0(:,k4+1:end);
                        cPhaseOffset=round((obj.plots.cPhaseAll(k1+1,Board)-...
                            continuousPhase(obj.nOverlap/obj.nCycle+1))/obj.nCycle)*obj.nCycle;
                        obj.plots.cPhaseAll(k1+1:k2,Board)=...
                            obj.plots.cPhaseAll(k1+1:k2,Board).*linspace(1,0,obj.nOverlap/obj.nCycle)'+...
                            linspace(0,1,obj.nOverlap/obj.nCycle)'.*(cPhaseOffset+...
                            continuousPhase(obj.nOverlap/obj.nCycle+1:2*obj.nOverlap/obj.nCycle))';
                        obj.plots.cPhaseAll(k2+1:k0+obj.nBlk/obj.nCycle+2,Board)=...
                            cPhaseOffset+continuousPhase(2*obj.nOverlap/obj.nCycle+1:end);
                        [nfilteredZ, filterWeights]=obj.CAR(filteredCA,filteredZ);
                        obj.plots.projWeightCA(k1+1:k2,Board)=...
                            obj.plots.projWeightCA(k1+1:k2,Board).*linspace(1,0,obj.nOverlap/obj.nCycle)'+...
                            linspace(0,1,obj.nOverlap/obj.nCycle)'.*...
                            filterWeights(obj.nOverlap/obj.nCycle+1:2*obj.nOverlap/obj.nCycle,1);
                        obj.plots.projWeightCA(k2+1:k0+obj.nBlk/obj.nCycle,Board)=...
                            filterWeights(2*obj.nOverlap/obj.nCycle+1:end);
                        nfilteredPwl=obj.subtractResidualPowerline(nfilteredZ,...
                            obj.plots.cPhaseAll(k0+1:k0+obj.nBlk/obj.nCycle+2,Board));
                        obj.plots.iqrAvg(k1*10+1:k1*10+(obj.nBlk-3*obj.nOverlap)/obj.dCycle,Board)=...
                            sqrt(mean(reshape(iqr(nfilteredPwl(obj.nOverlap+1:end-...
                            2*obj.nOverlap,:),2).^2,obj.dCycle,[]),1));
                        obj.plots=obj.PSD(obj.plots,nfilteredPwl,Board+obj.nBoards);
                        h5write(OutFile,'/data',int16(round(...
                        nfilteredPwl(obj.nOverlap+1:end-2*obj.nOverlap,:)')),...
                        [1+(Board-1)*obj.nChBoard,k1*obj.nCycle+1],...
                        [obj.nChBoard,obj.nBlk-3*obj.nOverlap]);
                    end
                    %last block (reduced length)
                    lastN=obj.recParameter.LenRec-(obj.countBlk-1)*(obj.nBlk-3*obj.nOverlap);
                    x=double(h5read(RawFile,obj.recParameter.HdfRawDataPath,...
                        [1+(Board-1)*obj.nChBoard,1+(obj.countBlk-1)*(obj.nBlk-3*obj.nOverlap)],...
                            [obj.nChBoard,lastN])');
                    %reduce voltage steps in individual traces
                    z=obj.reduceRecalibrationSteps(x);
                    %use smoothed (hamming kernel, 2 osc.) version subtracted common avg. signal
                    zSmth0=conv2(z,obj.SmthKer,'valid');
                    zSmth=[repmat(zSmth0(1,:),obj.nSmth/2*obj.nCycle,1); ...
                        zSmth0; repmat(zSmth0(end,:),obj.nSmth/2*obj.nCycle-1,1)];
                    filteredZ=z-zSmth;
                    %next, compute common average
                    CommonAverage=reshape(mean(median(reshape(z,[],obj.nChBoard/8,8),3),2),[],1); %better:group into 4x8 channels, then median, mean
                    %need to subtract a smoothed version
                    CASmth0=conv(CommonAverage,obj.SmthKer,'valid');
                    CASmth=[CASmth0(1,1)*ones(obj.nSmth/2*obj.nCycle,1); ...
                        CASmth0; CASmth0(end,1)*ones(obj.nSmth/2*obj.nCycle-1,1)];
                    filteredCA=CommonAverage-CASmth;
                    %now, before I do common average referencing, I should try to estimate 60 Hz
                    %phase and noise.
                    continuousPhase=obj.estimatePhase(filteredCA);
                    k0=(obj.countBlk-1)*(obj.nBlk-3*obj.nOverlap)/obj.nCycle;
                    k1=k0+obj.nOverlap/obj.nCycle;%
                    k2=k1+obj.nOverlap/obj.nCycle;
                    k3=floor(k1/60);
                    k4=floor(obj.nOverlap/(60*obj.nCycle));
                    obj.plots.medPwl0=obj.medianPwl(filteredCA,continuousPhase);
                    obj.plots.medPwl(:,k3+1:k3-k4+size(obj.plots.medPwl0,2),Board)=obj.plots.medPwl0(:,k4+1:end);
                    cPhaseOffset=round((obj.plots.cPhaseAll(k1+1,Board)-...
                        continuousPhase(obj.nOverlap/obj.nCycle+1))/obj.nCycle)*obj.nCycle;
                    obj.plots.cPhaseAll(k1+1:k2,Board)=...
                        obj.plots.cPhaseAll(k1+1:k2,Board).*linspace(1,0,obj.nOverlap/obj.nCycle)'+...
                        linspace(0,1,obj.nOverlap/obj.nCycle)'.*(cPhaseOffset+...
                        continuousPhase(obj.nOverlap/obj.nCycle+1:2*obj.nOverlap/obj.nCycle))';
                    obj.plots.cPhaseAll(k2+1:end,Board)=...
                        cPhaseOffset+continuousPhase(2*obj.nOverlap/obj.nCycle+1:end);
                    [nfilteredZ, filterWeights]=obj.CAR(filteredCA,filteredZ);
                    obj.plots.projWeightCA(k1+1:k2,Board)=...
                        obj.plots.projWeightCA(k1+1:k2,Board).*linspace(1,0,obj.nOverlap/obj.nCycle)'+...
                        linspace(0,1,obj.nOverlap/obj.nCycle)'.*...
                        filterWeights(obj.nOverlap/obj.nCycle+1:2*obj.nOverlap/obj.nCycle,1);
                    obj.plots.projWeightCA(k2+1:end,Board)=...
                        filterWeights(2*obj.nOverlap/obj.nCycle+1:end,1);
                    obj.plots.projWeightCA(end,Board)=obj.plots.projWeightCA(end-1,Board);
                    nfilteredPwl=obj.subtractResidualPowerline(nfilteredZ,...
                        obj.plots.cPhaseAll(k0+1:end,Board));
                    obj.plots.iqrAvg(10*k1+1:end-1,Board)=...
                        sqrt(mean(reshape(iqr(nfilteredPwl(obj.nOverlap+...
                        1:floor(lastN/obj.dCycle)*obj.dCycle,:),2).^2,obj.dCycle,[]),1));
                    obj.plots.iqrAvg(end,Board)=obj.plots.iqrAvg(end-1,Board);
                    h5write(OutFile,'/data',int16(round(...
                        nfilteredPwl(obj.nOverlap+1:end,:)')),...
                        [1+(Board-1)*obj.nChBoard,k1*obj.nCycle+1],...
                        [obj.nChBoard,lastN-obj.nOverlap]);
                end
            end
        end
        function z=reduceRecalibrationSteps(obj,x)
            y=cumsum([zeros(1,obj.nChBoard); x],1);
            y0=y(obj.nCycle+1:end,:)-y(1:end-obj.nCycle,:);%start at (obj.nCycle)/2, bin edges
            y1=median(y0,2);
            y2=y0-y1;
            %baseline
            ybl=(y2(1:end-obj.nCycle-1,:)+y2(obj.nCycle+2:end,:));%start at obj.nCycle,
            yDiff=abs(y2(1:end-obj.nCycle-1,:)-(y2(obj.nCycle+2:end,:)));%start at obj.nCycle
            
            %need to know how the difference is after 500 frames, use as weight?
            %--better:find isolated maxima, with a distance of at least 500 frames
            yDiffS=mean(yDiff,2);
            StepV=ones(size(yDiffS),'logical');
            for i=1:obj.nCycle
                StepV(i+1:end,:)=StepV(i+1:end,:) &(yDiffS(i+1:end,:)>yDiffS(1:end-i,:));
                StepV(1:end-i,:)=StepV(1:end-i,:) &(yDiffS(1:end-i,:)>=yDiffS(i+1:end,:));
            end
            %make sure that peaks are in region of interest
            StepV(1:obj.nCycle/2,:,:)=false;
            StepV(end-(obj.nCycle/2-1):end,:,:)=false;
            %determine which peaks can be considered outliers
            yThr=4*prctile(yDiffS(StepV,1),25);
            %all traces on a common voltage level.
            StepV=find(StepV.*(yDiffS-yThr)>0);
            for qq=1:length(StepV)
                yPk=StepV(qq);
                %linear extrapolation
                StepTriang=[ybl(yPk-obj.nCycle+(-obj.nStepSmth+1:0),:,:);
                    y2(yPk-obj.nCycle+(1:obj.nCycle-1),:,:)+2*y2(yPk,:,:)-y2(yPk-1:-1:yPk-obj.nCycle+1,:,:); ...
                    ybl(yPk,:); ...
                    2*y2(yPk+obj.nCycle,:,:)+y2(yPk+obj.nCycle+(1:obj.nCycle-1),:,:)-...
                    y2(yPk+obj.nCycle+(obj.nCycle-1:-1:1),:,:); ...
                    ybl(yPk+obj.nCycle+(0:obj.nStepSmth-1),:,:)];
                StepTriangSmth=convn(StepTriang,obj.StepSmth,'full');
                
                ybl(yPk-obj.nCycle+1:yPk+obj.nCycle-1,:)=obj.StepWeights.*StepTriangSmth(2*obj.nStepSmth+1:end-2*obj.nStepSmth,:)+...
                    (1-obj.StepWeights).*StepTriang(obj.nStepSmth+1:end-obj.nStepSmth,:);
            end
            z=x-[repmat(ybl(1,:),obj.nCycle,1); ...
                ybl; repmat(ybl(end,:),obj.nCycle,1)]/1000;%same length output
        end
        function continuousPhase=estimatePhase(obj,filteredCA)
            lenCA=floor(size(filteredCA,1)/obj.nCycle)*obj.nCycle;
            lenCAx=floor(size(filteredCA,1)/obj.nCycle-obj.nOsc+1)*obj.nCycle;
            %lenCA=floor(obj.nBlk/obj.nCycle)*obj.nCycle;%make sure data only includes full cycles
            %lenCAx=floor(obj.nBlk/obj.nCycle-obj.nOsc+1)*obj.nCycle;%length for overlapping blocks with length nOsc
            oscCycles=reshape(repmat(filteredCA(1:lenCA,1),1,obj.nOsc+1),[],1);
            oscCycles=reshape(oscCycles(1:obj.nOsc*lenCA+obj.nOsc*obj.nCycle,1),lenCA+obj.nCycle,obj.nOsc);
            oscCycles=reshape(permute(reshape(oscCycles(1:lenCA-(obj.nOsc-1)*obj.nCycle,:),...
                obj.nCycle,[],obj.nOsc),[1 3 2]),obj.nCycle*obj.nOsc,[]);%(CCxN)yields blocks of nOsc cycles, shifted by one cycle
            oscMean=mean(oscCycles,1);%(1xN)
            oscStd=std(oscCycles,1);
            %threshold for low fluctuations
            b=prctile(oscStd,10:70);
            oscThr=prctile(b(1:end-1)+diff(b).*(90:-1:31),15);
            oscGood=find(oscStd<oscThr);
            %nGood=length(oscGood);
            oscCyclesGood=oscCycles(:,oscGood)-oscMean(1,oscGood);
            
            %roughly estimate phase using sine wave
            %get some number for explained Std, plus phase
            [explainedStd,Phase0]=max(obj.templateSine*oscCyclesGood,[],1);
            explainedStd=explainedStd./oscStd(oscGood);
            %map to unit circle? -- no. locally determine phase drift (allow upto 0.1%
            %variation)
            
            oscPhase=complex(zeros(1,size(oscMean,2)));
            oscPhase(1,oscGood)=explainedStd.*exp(1i*2*pi*Phase0/obj.nCycle);
            
            %fit model with variable slope, average 900 cycles (about 10 sec)
            %do a hanning window, for weighting?--no. rather two timescales (do 450 cycles
            %as well)
            %better forward-backward AR[1] with tau=300 cycles
            oscSlope=-0.001:0.0002:0.001;
            nSlope=length(oscSlope);
            qSlope=zeros(nSlope,lenCAx/obj.nCycle);
            for jj=1:nSlope
                q=oscPhase.*exp(1i*2*pi*oscSlope(jj)*(1:(lenCAx/obj.nCycle)));
                qSlope(jj,:)=abs(filtfilt(1-exp(-1/300),[1 -exp(-1/300)],q))+1e-12*(jj==6);
            end
            %best phase drift
            [~,indSlope0]=max(qSlope,[],1);
            %smoothen a bit, want median across 10 reasonably independent cycles
            %increase resolution and determine absolute phase
            
            indSlope0=conv([indSlope0(100:-1:1) indSlope0 indSlope0(end:-1:end-99)],...
                hanning(201)/101,'valid')';
            indSlope0=[indSlope0(1)*ones(obj.nOsc/2,1); indSlope0; indSlope0(end)*ones(obj.nOsc/2,1)];
            %(extrapolating here)

            %second round of fitting, wider kernel, correct for phase drift; aim to
            %explain about 3/4 of phase drift in first round
            oscPhaseFit0=cumsum(0.00015*(indSlope0-6),1)';
            oscPhaseFit0=oscPhaseFit0-oscPhaseFit0(obj.nOsc/2);
            
            oscPhaseC=oscPhase.*exp(1i*2*pi*oscPhaseFit0(obj.nOsc/2+1:end-obj.nOsc/2));
            %do reflective boundaries (300 bins)
            oscPhaseC=[oscPhaseC(601:-1:2) oscPhaseC oscPhaseC(end-1:-1:end-600)];
            
            oscSlope=-0.0005:0.00005:0.0005;
            nSlope=length(oscSlope);
            qSlope=zeros(nSlope,lenCAx/obj.nCycle+1200);
            qPhase=zeros(nSlope,lenCAx/obj.nCycle+1200);
            for jj=1:nSlope
                q=oscPhaseC.*exp(1i*2*pi*oscSlope(jj)*([601:-1:2 ...
                    1:lenCAx/obj.nCycle lenCAx/obj.nCycle-1:-1:lenCAx/obj.nCycle-600]));
                qSlope(jj,:)=abs(filtfilt(1-exp(-1/300),[1 -exp(-1/300)],q))+1e-12*(jj==6);
                qPhase(jj,:)=angle(filtfilt(1-exp(-1/300),[1 -exp(-1/300)],q))/(2*pi);
            end
            %make one osc longer than needed.
            qPhase=qPhase(:,601-obj.nOsc/2:end-600+obj.nOsc/2)-oscSlope'*(-obj.nOsc/2+1:lenCAx/obj.nCycle+obj.nOsc/2);
            %best phase
            [~,indSlope1]=max(qSlope(:,601-obj.nOsc/2:end-600+obj.nOsc/2),[],1);
            
            qMask=accumarray([reshape(indSlope1,[],1) (1:lenCA/obj.nCycle+1)'],...
                ones(lenCA/obj.nCycle+1,1),[nSlope lenCA/obj.nCycle+1]);
            
            deltaPhase=mod(diff((sum(qMask.*qPhase,1)-oscPhaseFit0)*obj.nCycle,1,2)+...
                obj.nCycle/2,obj.nCycle)-obj.nCycle/2;
            %distribute outliers (smoothing jumps)
            deltaPhase=conv(deltaPhase,hanning(201)'/101,'full');
            deltaPhase(1,101:200)=deltaPhase(1,101:200)+deltaPhase(1,100:-1:1);
            deltaPhase(1,end-199:end-100)=deltaPhase(1,end-199:end-100)+deltaPhase(1,end:-1:end-99);
            deltaPhase=deltaPhase(1,101:end-100);
            Phase0=mod((qPhase(indSlope1(1),1)-oscPhaseFit0(1))*obj.nCycle,obj.nCycle)+0.5;
            %reverse phase mapping --> don't want individual frames, need
            %to do interpolation anyway.
            continuousPhase=(0:lenCA/obj.nCycle+1)*obj.nCycle+Phase0+...
                cumsum([0 deltaPhase deltaPhase(end)],2);
        end
        function [nfilteredZ, filterWeights]=CAR(obj,filteredCA,filteredZ)
            %%Common average referencing (function that does not clip)
            %Aim: remember phase, then do a common average referencing, then get rid
            %of of 60 Hz component. Do a bit of smoothing, and use a projection for CAR.

            %want a normalized common average referencing kernel (use about 12 osc) to
            %control gain
            lenCA=floor(size(filteredCA,1)/obj.nCycle)*obj.nCycle;
            %cfilteredCA=cumsum(filteredCA,1);
            cfilteredCA2=cumsum([0; filteredCA.^2],1);
            %normalization for weights
            cfilteredZ2=cumsum([zeros(1,obj.nChBoard); reshape(sum(reshape(filteredZ(1:lenCA,:),...
                obj.nCycle,[]).^2,1),lenCA/obj.nCycle,[])],1);
            %sliding dot product
            cfilteredZ=cumsum([zeros(1,obj.nChBoard); filteredCA.*filteredZ],1);
            %gain
            filteredGain=(cfilteredZ(obj.nOsc*obj.nCycle+1:end,:)-...
                cfilteredZ(1:end-obj.nOsc*obj.nCycle,:))./...
                (cfilteredCA2(obj.nOsc*obj.nCycle+1:end,1)-...
                cfilteredCA2(1:end-obj.nOsc*obj.nCycle,1));
            %weights
            filterWeights=mean((cfilteredZ(obj.nOsc*obj.nCycle+1:obj.nCycle:end,:)-...
                cfilteredZ(1:obj.nCycle:end-obj.nOsc*obj.nCycle,:))./...
                sqrt(cfilteredCA2(obj.nOsc*obj.nCycle+1:obj.nCycle:end,1)-...
                cfilteredCA2(1:obj.nCycle:end-obj.nOsc*obj.nCycle,1))./...
                sqrt(cfilteredZ2(obj.nOsc+1:end,:)-...
                cfilteredZ2(1:end-obj.nOsc,:)),2);
            filterWeights=[repmat(filterWeights(1,1),obj.nOsc/2-1,1); filterWeights; ...
                repmat(filterWeights(end,1),obj.nOsc/2,1)];
            %average and normalize, subtract from filtered raw traces
            nfilteredZ=filteredZ-filteredCA.*...
                [repmat(filteredGain(1,:),obj.nOsc/2*obj.nCycle-1,1); filteredGain; ...
                repmat(filteredGain(end,:),obj.nOsc/2*obj.nCycle,1)];
        end
%         function nfilteredZ=CAR(obj,filteredCA,z)
%             %%Common average referencing (function that does clip)
%             %use smoothed (hamming kernel, 2 osc.) version subtracted common avg. signal
%             zSmth=conv2(z,obj.SmthKer,'valid');
%             filteredZ=z-zSmth;
%             %want a normalized common average referencing kernel (use about 12 osc) to
%             %control gain
%             %cfilteredCA=cumsum(filteredCA,1);
%             cfilteredCA2=cumsum(filteredCA.^2,1);
%             %sliding dot product
%             cfilteredZ=cumsum(filteredCA.*filteredZ,1);
%             %average and normalize, subtract from filtered raw traces
%             nfilteredZ=filteredZ(obj.nOsc/2*obj.nCycle+1:end-obj.nOsc/2*obj.nCycle,:)-...
%                 filteredCA(obj.nOsc/2*obj.nCycle+1:end-obj.nOsc/2*obj.nCycle,1).*...
%                 (cfilteredZ(obj.nOsc*obj.nCycle+1:end,:)-cfilteredZ(1:end-obj.nOsc*obj.nCycle,:))./...
%                 (cfilteredCA2(obj.nOsc*obj.nCycle+1:end,1)-cfilteredCA2(1:end-obj.nOsc*obj.nCycle,1));
%         end
        function nfilteredPwl=subtractResidualPowerline(obj,nfilteredZ,continuousPhase)
            %%Residual powerline for each channel (no clipping)
            %need to do channelwise (interpolation), can do parallel?!
            lenCA=floor(size(nfilteredZ,1)/obj.nCycle)*obj.nCycle;
            nfilteredPwl=1*nfilteredZ;
            cPhaseTemp=interp1(0.5+(0:floor(lenCA/obj.nCycle)+1)*obj.nCycle,continuousPhase,...
                1:size(nfilteredZ,1));
            i0=floor(min(cPhaseTemp-1)/obj.nCycle+1)*obj.nCycle;
            i1=floor(max(cPhaseTemp)/obj.nCycle-1)*obj.nCycle;
            nMed=floor((i1-i0)/obj.nCycle/6)*6;
            nOffset=(i1-i0)/obj.nCycle-nMed;
            nO0=floor(nOffset/2);
            nO1=floor((nOffset+1)/2);
            for jj=1:obj.nChBoard
                %forward mapping; reshape to obj.nCycle,[];
                chWarped0=reshape(interp1(cPhaseTemp(1,1:end)'-obj.nCycle, ...
                    nfilteredZ(:,jj)',...
                    i0+1:i1),obj.nCycle,[]);
                %sliding median + average (about 300 cycles?)
                chWarped=cumsum([zeros(obj.nCycle,1) chWarped0(:,1:nO0) ...
                    reshape(repmat(median(reshape(chWarped0(:,nO0+1:end-nO1),...
                    obj.nCycle,6,[]),2),1,6,1),obj.nCycle,[])  chWarped0(:,end-nO1+1:end)],2);
                chAvg0=chWarped(:,2*obj.nAvgh+2:end)-chWarped(:,1:end-2*obj.nAvgh-1);
                chAvg=reshape([repmat(chAvg0(:,1),1,obj.nAvgh+1) chAvg0 repmat(chAvg0(:,end),1,obj.nAvgh+1)],[],1);
                %backward mapping
                chAverage=interp1(i0-obj.nCycle+1:i1+obj.nCycle, chAvg', cPhaseTemp');
                %normalization
                chAvgc2=cumsum(chAverage.^2,1);
                %projection and subtraction
                cPwlZ=cumsum(chAverage.*nfilteredZ(:,jj),1);
                filteredGain=(cPwlZ(obj.nAvgF+1:end,1)-cPwlZ(1:end-obj.nAvgF,1))./...
                    (chAvgc2(obj.nAvgF+1:end,1)-chAvgc2(1:end-obj.nAvgF,1));
                nfilteredPwl(:,jj)=nfilteredZ(:,jj)-chAverage(:,1).*...
                    [repmat(filteredGain(1,:),obj.nAvgF/2,1); filteredGain; ...
                    repmat(filteredGain(end,:),obj.nAvgF/2,1)];
            end
        end
        function medPwl=medianPwl(obj,filteredCA,continuousPhase)
            %%want a one second median for powerline noise
            lenCA=floor(size(filteredCA,1)/obj.nCycle)*obj.nCycle;
            cPhaseTemp=interp1(0.5+(0:floor(lenCA/obj.nCycle)+1)*obj.nCycle,continuousPhase,...
                1:lenCA);
            i0=floor(min(cPhaseTemp-1)/obj.nCycle+1)*obj.nCycle;
            i1=floor(max(cPhaseTemp)/obj.nCycle-1)*obj.nCycle;
            nMed=floor((i1-i0)/obj.nCycle/60)*60;
            %nOffset=(i1-i0)/obj.nCycle-nMed;
            chWarped0=reshape(interp1(cPhaseTemp(1,1:end)'-obj.nCycle, ...
                filteredCA(1:floor(lenCA/obj.nCycle)*obj.nCycle,1)',...
                i0+1:i1),obj.nCycle,[]);
            %sliding median + average (about 300 cycles?)
            %if nOffset>5
            %    obj.plots.medPwl=[reshape(median(reshape(chWarped0(:,1:nMed),...
            %        obj.nCycle,60,[]),2),obj.nCycle,[])  median(chWarped0(:,nMed+1:end),2)];
            %else
            medPwl=reshape(median(reshape(chWarped0(:,1:nMed),...
                obj.nCycle,60,[]),2),obj.nCycle,[]);
            %end
        end
        %%plotting
        function Xplots=PSD(obj,Xplots,Xnew,iii)
            for j=1:Xplots.Nfft
                %fft
                F=fft(complex(diag(hanning(Xplots.dtCS))*Xnew((j-1)*Xplots.dtCS+1:j*Xplots.dtCS,:)));
                %cross spectrum in frequency interval (need to normalize)
                Xplots.CS{iii}=Xplots.CS{iii}+F(Xplots.fb0+1:Xplots.fb1,:)'*F(Xplots.fb0+1:Xplots.fb1,:);
                Xplots.QS(:,iii)=Xplots.QS(:,iii)+mean(F.*conj(F),2);
                %PowerHisto
                cs=F(Xplots.fb0+1:Xplots.fb1,:)'*F(Xplots.fb0+1:Xplots.fb1,:);
                phi=reshape(cs,1,obj.nChBoard^2);
                phi=reshape(phi(1:end-1),obj.nChBoard-1,obj.nChBoard+1);
                phi=mean(mean(abs(phi(5:end-3,:))));
                pInd=floor(max(min((log10(phi*Xplots.PHiN))*50,200),1));
                Xplots.QH(pInd,iii)=Xplots.QH(pInd,iii)+1;
                phi=abs(diag(cs));
                pInd=floor(max(min((log10(mean(phi)*Xplots.PHiN))*50,200),1));
                Xplots.PH(pInd,iii)=Xplots.PH(pInd,iii)+1;
                pInd=floor(max(min((log10(phi*Xplots.PHiN))*50,200),1))+Xplots.PHch{iii};
                Xplots.PHdc(pInd,1)=Xplots.PHdc(pInd,1)+1;
                %power spectrum
                P=fft(complex(mean(Xnew((j-1)*Xplots.dtCS+1:j*Xplots.dtCS,:),2).*hanning(Xplots.dtCS)));
                Xplots.PS(:,iii)=Xplots.PS(:,iii)+P.*conj(P);
            end
        end
        function plotResults(obj,Filter)
            %want to plot phase, amplitude of commonAverage, amplitude of
            %single electrode powerline (and phase-locked temporal average)
            %auto and cross-spectrum, covariance matrix.
            figVT=figure('Position',[0 0 1600 1200]);
            %covariance matrices
            axCS1 = axes('OuterPosition',[0   3/4 1/4 1/4]);
            axCS2 = axes('OuterPosition',[0   1/2 1/4 1/4]);
            axCS3 = axes('OuterPosition',[0   1/4 1/4 1/4]);
            %axCS4 = axes('OuterPosition',[0   0 1/4 1/4]);
            %phase
            axPhase = axes('OuterPosition',[1/4 3/4 1/4 1/4]);
            %power spectrum
            axPS = axes('OuterPosition',[1/4 2/4 0.35 1/4]);
            axQS = axes('OuterPosition',[1/4 1/4 0.35 1/4]);
            %phase aligned
            axPhaseAligned = axes('OuterPosition',[0.5 3/4 0.5 1/4]);
            %noise power histos on single electrodes
            axPSHsingle = axes('OuterPosition',[0.6 2/4 0.4 1/4]);
            % avged histos
            axPSH = axes('OuterPosition',[0.6 1/4 0.4 1/4]);
            %
            axWt= axes('OuterPosition',[0   0   0.45 1/4]);
            axIqr= axes('OuterPosition',[0.5 0   0.45 1/4]);
            axWtH= axes('OuterPosition',[0.4   0   0.1 1/4]);
            axIqrH= axes('OuterPosition',[0.9 0   0.1 1/4]);
            %Voltage traces:
            %xlabel(axRawV,'time/ms')
            %ylabel(axRawV,'voltage/µV')
            %xlabel(axFilt1V,'time/ms')
            %ylabel(axFilt1V,'voltage/µV')
            
            %ylabel(axPSt,'0.5-2 kHz power (global)/(µV^2/kHz)')
            %xlabel(axPSt,'time/s')
            %ylabel(axPStd,'0.5-2 kHz power (channel)/(µV^2/kHz)')
            %xlabel(axPStd,'time/s')
            %Coherence
            hold(axPS,'on');
            hold(axQS,'on');
            xlabel(axPS,'frequency/Hz')
            ylabel(axPS,'power (global)/(µV^2/kHz)')
            xlabel(axQS,'frequency/Hz')
            ylabel(axQS,'power (channel)/(µV^2/kHz)')

            
            Snorm=2000*obj.recParameter.bitVolt^2/(sum(hanning(obj.plots.dtCS).^2)*...
                obj.recParameter.sRate*(obj.countBlk-1)*obj.plots.Nfft);
            %need 4 times
            imagesc(axCS1,log10(min(max(abs([obj.plots.CS{1}; obj.plots.CS{2}])*Snorm/(obj.plots.fb1-obj.plots.fb0),1),1e4)))
            %norm=mplt.colors.LogNorm(5e1,1e5))
            hC=colorbar(axCS1);
            set(hC,'Ytick',0:3,'YTicklabel',{'1' '1e1' '1e2' '1e3'});
            hC.Label.String = 'power*kHz/µV^2';
            ylim(axCS1,[0.5 2*obj.nChBoard+0.5])
            xlim(axCS1,[0.5 obj.nChBoard+0.5])
            xticks(axCS1,8:8:obj.nChBoard)
            xticklabels(axCS1,(8:8:obj.nChBoard))
            yticks(axCS1,8:8:2*obj.nChBoard)
            yticklabels(axCS1,(8:8:2*obj.nChBoard))
            imagesc(axCS2,log10(min(max(abs([obj.plots.CS{3}; obj.plots.CS{4}])*Snorm/(obj.plots.fb1-obj.plots.fb0),1),1e4)))
            %norm=mplt.colors.LogNorm(5e1,1e5))
            hC=colorbar(axCS2);
            set(hC,'Ytick',0:3,'YTicklabel',{'1' '1e1' '1e2' '1e3'});
            hC.Label.String = 'power*kHz/µV^2';
            ylim(axCS2,[0.5 2*obj.nChBoard+0.5])
            xlim(axCS2,[0.5 obj.nChBoard+0.5])
            xticks(axCS2,8:8:obj.nChBoard)
            xticklabels(axCS2,(8:8:obj.nChBoard))
            yticks(axCS2,8:8:2*obj.nChBoard)
            yticklabels(axCS2,(8:8:2*obj.nChBoard))
            %correction 
            imagesc(axCS3,log10(min(max(abs([obj.plots.CS{3}; obj.plots.CS{4}]-...
                [obj.plots.CS{1}; obj.plots.CS{2}])*Snorm/(obj.plots.fb1-obj.plots.fb0),1),1e4)))
            %norm=mplt.colors.LogNorm(5e1,1e5))
            hC=colorbar(axCS3);
            set(hC,'Ytick',0:3,'YTicklabel',{'1' '1e1' '1e2' '1e3'});
            hC.Label.String = 'power*kHz/µV^2';
            ylim(axCS3,[0.5 2*obj.nChBoard+0.5])
            xlim(axCS3,[0.5 obj.nChBoard+0.5])
            xticks(axCS3,8:8:obj.nChBoard)
            xticklabels(axCS3,(8:8:obj.nChBoard))
            yticks(axCS3,8:8:2*obj.nChBoard)
            yticklabels(axCS3,(8:8:2*obj.nChBoard))
            title(axCS1,'raw')
            title(axCS2,'corrected')
            title(axCS3,'correction')
            %
            h=zeros(4,1);
            h(1)=semilogy(axPS,(1:obj.plots.dtCS/2+1)*obj.recParameter.sRate/obj.plots.dtCS,...
                abs(obj.plots.PS(1:obj.plots.dtCS/2+1,1))*Snorm,'k-');
            h(2)=semilogy(axPS,(1:obj.plots.dtCS/2+1)*obj.recParameter.sRate/obj.plots.dtCS,...
                mean(abs(obj.plots.PS(1:obj.plots.dtCS/2+1,1:obj.nBoards)),2)*Snorm,'b-');
            h(3)=semilogy(axPS,(1:obj.plots.dtCS/2+1)*obj.recParameter.sRate/obj.plots.dtCS,...
                abs(obj.plots.PS(1:obj.plots.dtCS/2+1,obj.nBoards+1))*Snorm,'r-');
            h(4)=semilogy(axPS,(1:obj.plots.dtCS/2+1)*obj.recParameter.sRate/obj.plots.dtCS,...
                mean(abs(obj.plots.PS(1:obj.plots.dtCS/2+1,obj.nBoards+1:2*obj.nBoards)),2)*Snorm,'m-');
            legend(h,'initial, board 1','initial, other', 'filtered, board 1', 'filtered, other')
            axPS.YScale='log';
            xlim(axPS,[100 2000])
            ylim(axPS,[1e-1 1e5])
            semilogy(axQS,(1:obj.plots.dtCS/2+1)*obj.recParameter.sRate/obj.plots.dtCS,...
                abs(obj.plots.QS(1:obj.plots.dtCS/2+1,1))*Snorm,'k-')
            semilogy(axQS,(1:obj.plots.dtCS/2+1)*obj.recParameter.sRate/obj.plots.dtCS,...
                mean(abs(obj.plots.QS(1:obj.plots.dtCS/2+1,2:obj.nBoards)),2)*Snorm,'b-')
            semilogy(axQS,(1:obj.plots.dtCS/2+1)*obj.recParameter.sRate/obj.plots.dtCS,...
                abs(obj.plots.QS(1:obj.plots.dtCS/2+1,obj.nBoards+1))*Snorm,'r-')
            semilogy(axQS,(1:obj.plots.dtCS/2+1)*obj.recParameter.sRate/obj.plots.dtCS,...
                mean(abs(obj.plots.QS(1:obj.plots.dtCS/2+1,obj.nBoards+2:2*obj.nBoards)),2)*Snorm,'m-')
            axQS.YScale='log';%bug in MATLAB?
            xlim(axQS,[100 2000])
            ylim(axQS,[1e-1 1e5])
            %normalize to density
            %PH=PH*diag(1./sum(PH{iii},1));
            %plot
            loglog(axPSH,logspace(0,4,200),obj.plots.QH(:,1)/...
                sum(obj.plots.QH(:,1)),'k-')
            hold(axPSH,'on');
            loglog(axPSH,logspace(0,4,200),sum(obj.plots.QH(:,2:obj.nBoards),2)/...
                sum(obj.plots.QH(:,2:obj.nBoards),'all'),'b-')
            loglog(axPSH,logspace(0,4,200),obj.plots.QH(:,obj.nBoards+1)/...
                sum(obj.plots.QH(:,obj.nBoards+1)),'r-')
            loglog(axPSH,logspace(0,4,200),sum(obj.plots.QH(:,obj.nBoards+2:2*obj.nBoards),2)/...
                sum(obj.plots.QH(:,obj.nBoards+2:2*obj.nBoards),'all'),'m-')
            loglog(axPSH,logspace(0,4,200),obj.plots.PH(:,1)/...
                sum(obj.plots.PH(:,1)),'k:')
            loglog(axPSH,logspace(0,4,200),sum(obj.plots.PH(:,2:obj.nBoards),2)/...
                sum(obj.plots.PH(:,2:obj.nBoards),'all'),'b:')
            loglog(axPSH,logspace(0,4,200),obj.plots.PH(:,obj.nBoards+1)/...
                sum(obj.plots.PH(:,obj.nBoards+1)),'r:')
            loglog(axPSH,logspace(0,4,200),sum(obj.plots.PH(:,obj.nBoards+2:2*obj.nBoards),2)/...
                sum(obj.plots.PH(:,obj.nBoards+2:2*obj.nBoards),'all'),'m:')
            axPSH.YScale='log';
            axPSH.XScale='log';
            xlim(axPSH,[1 1e4])
            ylim(axPSH,[1e-4 1])
            xlabel(axPSH,'0.5-2 kHz power/(µV^2/kHz)')
            ylabel(axPSH,'density')
            PHdx=reshape(obj.plots.PHdc,200,2*obj.nChBoard*obj.nBoards);
            PHdx=PHdx(:,2:2:end);
            PHdx=PHdx*diag(1./sum(PHdx,1));
            imagesc(axPSHsingle,[0 4],[0 obj.nChBoard*obj.nBoards],log10(max(PHdx',1e-4)))
            set(axPSHsingle,'Xtick',0:4,'XTicklabel',{'1' '1e1' '1e2' '1e3' '1e4'});
            hC=colorbar(axPSHsingle);
            set(hC,'Ytick',-4:-1,'YTicklabel',{'1e-4' '1e-3' '1e-2' '1e-1'});
            hC.Label.String = 'density';
            %saveas(figVT,[Filter.car.PlotBase filesep Filter.car.plotFolder filesep Filter.car.plotName])
            %close(figVT)
            plot(axPhase,(1:floor(obj.recParameter.LenRec/obj.nCycle)+1)/(60*obj.recParameter.sRate/obj.nCycle),...
                mod(obj.plots.cPhaseAll(1:end-1,1),obj.nCycle)*1000/obj.recParameter.sRate,'k.','MarkerSize',1)
            hold(axPhase,'on')
            plot(axPhase,(1:floor(obj.recParameter.LenRec/obj.nCycle)+1)/(60*obj.recParameter.sRate/obj.nCycle),...
                mod(obj.plots.cPhaseAll(1:end-1,2),obj.nCycle)*1000/obj.recParameter.sRate,'b.','MarkerSize',0.5)
            xlabel(axPhase,'time [min]')
            ylabel(axPhase,'phase [ms]')
            ylim(axPhase,[0 obj.nCycle*1000/obj.recParameter.sRate])
            imagesc(axPhaseAligned,[squeeze(obj.plots.medPwl(:,:,1)); nan*zeros(50,size(obj.plots.medPwl,2));
            squeeze(obj.plots.medPwl(:,:,2))])
            xticks(axPhaseAligned,(0:3600*obj.nCycle:obj.recParameter.LenRec)/(60*obj.nCycle))
            xticklabels(axPhaseAligned,(0:3600*obj.nCycle:obj.recParameter.LenRec)/(60*obj.recParameter.sRate))
            yticks(axPhaseAligned,[0:5*obj.recParameter.sRate/1000:obj.nCycle ...
                obj.nCycle+50:5*obj.recParameter.sRate/1000:2*obj.nCycle+50])
            yticklabels(axPhaseAligned,[0:5:obj.nCycle*1000/obj.recParameter.sRate ...
                0:5:obj.nCycle*1000/obj.recParameter.sRate])
            hC=colorbar(axPhaseAligned);
            hC.Label.String = 'phase aligned voltage [µV]';
            ylabel(axPhaseAligned,'phase [ms]')
            xlabel(axPhaseAligned,'time [min]')
            plot(axWt,(1:size(obj.plots.projWeightCA,1))/(60*obj.recParameter.sRate/obj.nCycle),...
                obj.plots.projWeightCA(:,1),'k.','MarkerSize',1)
            hold(axWt,'on')
            plot(axWt,(1:size(obj.plots.projWeightCA,1))/(60*obj.recParameter.sRate/obj.nCycle),...
                obj.plots.projWeightCA(:,2),'b.','MarkerSize',1)
            ylim(axWt,[0 1])
            xlabel(axWt,'time [min]')
            ylabel(axWt,'correlation coefficient')
            hBins=linspace(0,1,201);
            semilogx(axWtH,histcounts(obj.plots.projWeightCA(:,1),hBins(1:2:end)),...
                hBins(2:2:end),'k-')
            hold(axWtH,'on')
            semilogx(axWtH,histcounts(obj.plots.projWeightCA(:,2),hBins(1:2:end)),...
                hBins(2:2:end),'b-')
            yticks(axWtH,[])
            xticks(axWtH,[])
            %
            semilogy(axIqr,(1:size(obj.plots.iqrAvg,1))/(600*obj.recParameter.sRate/obj.nCycle),...
                obj.plots.iqrAvg(:,1)*obj.recParameter.bitVolt,'k.','MarkerSize',1)
            hold(axIqr,'on')
            semilogy(axIqr,(1:size(obj.plots.iqrAvg,1))/(600*obj.recParameter.sRate/obj.nCycle),...
                obj.plots.iqrAvg(:,2)*obj.recParameter.bitVolt,'b.','MarkerSize',1)
            ylim(axIqr,[10 1e4])
            xlabel(axIqr,'time [min]')
            ylabel(axIqr,'interquartile range [µV]')
            hBins=logspace(1,4,201);
            loglog(axIqrH,histcounts(obj.plots.iqrAvg(:,1)*obj.recParameter.bitVolt,hBins(1:2:end)),...
                hBins(2:2:end),'k-')
            hold(axIqrH,'on')
            loglog(axIqrH,histcounts(obj.plots.iqrAvg(:,2)*obj.recParameter.bitVolt,hBins(1:2:end)),...
                hBins(2:2:end),'b-')
            yticks(axIqrH,[])
            xticks(axIqrH,[])
            %save figure
            saveas(figVT,[Filter.PlotBase filesep Filter.plotFolder filesep Filter.plotName])
            close(figVT)
        end
    end
    %%Synchronized stuff?, optional
    %want to get rid of stuff with zero common avg., but high variance across the
    %array (in a sense that more than half of the electrodes reflect outliers
    %or so
    %find(sum(abs(nfilteredPwl),2)-4*mean(sum(abs(nfilteredPwl),2))-10*std(sum(abs(nfilteredPwl),2))>0)
end