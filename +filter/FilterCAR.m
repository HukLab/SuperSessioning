function [Filter,recParameter]=FilterCAR(RawFile,OutFile,Filter,recParameter)

%%Common noise reduction

%This script aims at reducing high-frequency noise found on all channels.
%It can be applied to raw data (without 300Hz high-pass filtering).

%Noise reduction is done subtractively, in two steps:

%First to reduce relatively broad-band noise, it does
% a leaky integration (timescale 15 ms, see ARtau)
% of the smoothed (dt, 3 frames) derivative of all signals. To avoid outliers, these
% derivatives are transformed by a saturating function (1-exp(-x),
% with a maximum amplitude of vSat=50 muV/dt).
%These differentiated, clipped and integrated signals are then averaged
% and the correlation coefficient with each of these signals and their average
% (minus its own contribution to the average) is determined.
%That corrected averaged signal times its correlation coefficient
% (across the whole recording) is subtracted from the raw signal.

%The second step takes the corrected signals from the first step
% and does the same procedure, with a shorter timescale (5 ms)
% and lower saturation amplitude of the differentiated signal (10 muV/dt).
%The hope here is that it could partially correct for errors from the first step
% that were due to signals with a large amplitude but on few channels.

%Additional notes:
%I believe that the smoothing step (dt, 3 frames) is required
% to avoid inducing high frequency noise.
%Differentiating, clipping and leaky integration are done
% to circumvent high-pass filtering at this stage.
%Therefore, this script will barely affect low frequency noise components.
%Using a causal filter, would be better to run it in forward and backward direction (TODO).
%TODO: Boundary conditions (few frames, so ignore for now)

%TODO: Mask for use directly with raw data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nch=sum(recParameter.ChMask,2);

%%Parameters (to be set in this script):
%maximum number of frames in one batch
nBatchMax=1e5;

%frequency band for visualization
f0=500;%Hz
f1=2000;%Hz

%number of channels
%Nch=64
%time step for computing voltage derivative (in frames)
%same for both iterations
dt=3;%odd
dtshift=floor((dt+1)/2);
dtCS=2^14;

%clipping of voltage derivative
%units: muV/dt or rather 1./bitVolt/dt (bitVolt is about 0.2)
%first iteration (allow for large steps in voltage (movement artefacts?))
vSaturation1=100;
%second iteration (no large steps in voltage on individual channels, accounting for spikes)
vSaturation2=100;

%timescales of the AR filters
%note that this is the maximum time window for the history of the averaged fluctuations.
%units: ms
%first iteration (depending on high pass filter setting)
Tau1=15;
%second iteration (faster fluctuations, timescale including spike repolarization phase)
Tau2=15;

%bias median estimate of difference towards 0 (percent channels)
BiasPerc=0.05;%a bit more than one in 32 (to ignore spikes)
BiasPercs=[0.5-BiasPerc 0.5+BiasPerc];


%whether signal power histograms shall be generated
PowerHisto=true;

PlotChannel=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialization (do not change without debugging)
if recParameter.readfromKWIK && ~isfield(Filter,'pwl')
    iGroup=recParameter.iGroup;
    bitVolt=h5read(RawFile,['/recordings/' num2str(iGroup) '/application_data/channel_bit_volts']);
    sRate=h5read(RawFile,['/recordings/' num2str(iGroup) '/application_data/channel_sample_rates']);
    %time_stamps=h5read(RawFile,['/recordings/' num2str(iGroup) '/application_data/timestamps']);
    %DataSize=size(h5read(RawFile,'/recordings/' num2str(iGroup) '/data'));%probably need to find a better solution here for large files
    %version_name=h5readatt(RawFile,'/','kwik_version');
    %rec_method=h5readatt(RawFile,'/recordings/0','name');
    %start_time=h5readatt(RawFile,'/recordings/0','start_time');
    %start_sample=h5readatt(RawFile,'/recordings/0','start_sample');
    %sample_rate=h5readatt(RawFile,'/recordings/0','sample_rate');
    %bit_depth=h5readatt(RawFile,'/recordings/0','bit_depth');
    
    recParameter.bitVolt=median(double(bitVolt(ChMap(ChMask),1)));
    recParameter.sRate=median(double(sRate(ChMap(ChMask),1)));
    recParameter.bitVoltCh=double(bitVolt(ChMap(ChMask),1));
    recParameter.sRateCh=double(sRate(ChMap(ChMask),1));
    HdfDataPath=['/recordings/' num2str(iGroup) '/data'];
    a=h5info(RawFile,HdfDataPath);
    DataSize=a.Dataspace.Size;
    if recParameter.tEnd==-1
        recParameter.nEnd=double(floor(DataSize(1,2)));
    else
        recParameter.nEnd=double(floor(recParameter.tEnd*60*recParameter.sRate));
    end
    if recParameter.tStart==0
        recParameter.nStart=1;
    else
        recParameter.nStart=double(floor(recParameter.tStart*60*recParameter.sRate)+1);
    end
    recParameter.LenRec=double(floor(DataSize(1,2)));
    recParameter.Nch=sum(recParameter.ChMask,2);
else
    HdfDataPath='/data';
end
if ~isfield(Filter,'pwl')
    assert(recParameter.ChMap==1:length(recParameter.ChMap),'remapping not implemented at this stage')
    assert(all(recParameter.ChMask),'masking of channels not implemented at this stage')
end
h5create(OutFile,'/data',[Nch Inf],'ChunkSize',[1 2048],'Datatype','int16','FillValue',int16(0))

vSat=vSaturation1./recParameter.bitVoltCh';
vSat2=vSaturation2./recParameter.bitVoltCh';
%frequency band for visualization
fb0=floor(dtCS*f0/recParameter.sRate);%bin index in spectrum
fb1=floor(dtCS*f1/recParameter.sRate);

%assert recParameter.LenRec>10e3
nBatch=floor(nBatchMax/dtCS)*dtCS;
%use reflective boundaries (warmup time)?
%better: estimate forward+backward and correct half for each. (too much computing time?)
%create a batch length compatible with the length of the time series
while mod((recParameter.LenRec-dt-2),nBatch)>3*dtCS
    nBatch=nBatch-dtCS;
end
%assert nBatch>1001
Nfft=nBatch/dtCS;

%nBatch=10000*2%set to multiple of dt (saves a lot of modulo operations
%assert (nBatch%2)==0
nRuns=floor((recParameter.LenRec-dt-2)/nBatch);


%want to take a spatial average of differential changes
%and integrate them in a leaky way (3ms or 100 frames)

%do for each channel separately first (3 frames diff, 10muV)
%compute corelation coefficients with individual traces (subtract (1/Nch of contrib. of channel itself)

ARtau=exp(-1000/recParameter.sRate/Tau1);
ARtau2=exp(-1000/recParameter.sRate/Tau2);
%ARslow=exp(-0.1/recParameter.sRate);%want to prevent accumulation of numerical errors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%first run: estimate variances

%initial condition
A=zeros(1,Nch);
V=zeros(1,Nch);
%compute variances (in filtered component) for weights or scaling
%iterate
for i=1:nRuns
    x=double(h5read(RawFile,HdfDataPath,[1,(i-1)*nBatch+1],[Nch,nBatch+dt])');
    h=(1-exp(-abs(x(1+dt:end,:)-x(1:end-dt,:))*diag(1./vSat))).*sign(x(1+dt:end,:)-x(1:end-dt,:));
    %treat as driving noise of AR process (can assume zero mean)
    VAi=zeros(1,Nch);
    for t=1:nBatch
        A=A*ARtau;
        A=A+h(t,:);
        %need variances
        VAi=VAi+A.^2;
    end
    V=V+VAi;
end

%factor for normalization
%another fix: don't allow any signal to have
% less than 1/2 of the standard deviation of the median across all channels
Sdev=max(sqrt(V),0.5*median(sqrt(V)));
Sdev=Sdev/mean(Sdev);

%initial condition
A=zeros(1,Nch);
Am=0;
VA=zeros(1,Nch);
VAm=0;
Pm=zeros(1,Nch);
Cm=zeros(1,nRuns);
Cstd=zeros(1,nRuns);


%iterate
for i=1:nRuns
    x=double(h5read(RawFile,HdfDataPath,[1,(i-1)*nBatch+1],[Nch,nBatch+dt])');
    h=(1-exp(-abs(x(1+dt:end,:)-x(1:end-dt,:))*diag(1./vSat))).*sign(x(1+dt:end,:)-x(1:end-dt,:));
    prc=quantile(h*diag(1./Sdev),BiasPercs,2);
    xM=min(abs(prc),[],2)*0.5.*sum(sign(prc),2)*Nch;
    %treat as driving noise of AR process (can assume zero mean)
    VAi=zeros(1,Nch);
    VAmi=0;
    Pmi=zeros(1,Nch);
    for t=1:nBatch
        A=A*ARtau;
        Am=Am*ARtau;
        A=A+h(t,:);
        Am=Am+xM(t,1);
        %correlation coefficient:
        %need variances and avg. product
        VAi=VAi+A.^2;
        VAmi=VAmi+Am^2;
        Pmi=Pmi+A.*(Am-A./Sdev);%subtracting each channel's own contribution to mean
    end
    %want to see the temporal development
    CC=Pmi./sqrt(VAmi*VAi);
    Cm(1,i)=mean(CC);
    Cstd(1,i)=std(CC);
    Pm=Pm+Pmi;
    VAm=VAm+VAmi;
    VA=VA+VAi;
end

%compute correlation coefficients
CCx=Pm./sqrt(VAm*VA);
%Sx=sqrt(VA);
%Sxm=sqrt(VAm);
PCx=CCx.*sqrt(VA/VAm).*vSat;%to convert from avg. to single var

%%%%%%%%%%%%%%%%%%%%%
%second run: estimate for shorter timescale from residuals
%initial condition
VB=zeros(1,Nch);
VBm=0;
Am=zeros(1,2);
A=zeros(2,Nch);
Bm=0;
B=zeros(1,Nch);
Za=zeros(nBatch+1,Nch);
Qm=zeros(1,Nch);
Dm=zeros(1,nRuns);
Dstd=zeros(1,nRuns);
%ACM=zeros(Nch,Nch);
%BCM=zeros(Nch,Nch);

%iterate
%initialize with first 2*dt frames
x=double(h5read(RawFile,HdfDataPath,[1,1],[Nch,dt+1])');
h=(1-exp(-abs(x(dt+1:end,:)-x(1:end-dt,:))*diag(1./vSat))).*sign(x(dt+1:end,:)-x(1:end-dt,:));
prc=quantile(h*diag(1./Sdev),BiasPercs,2);
xM=min(abs(prc),[],2)*0.5.*sum(sign(prc),2)*Nch;
Am(1,1)=Am(1,2)*ARtau;
Am(1,1)=Am(1,1)+xM(1,1);
A(1,:)=A(2,:)*ARtau;
A(1,:)=A(1,:)+h(1,:);
%correction
Za(end,:)=(Am(1,1)-A(1,:)./Sdev).*PCx;%need to be shifted in time below!
for i=1:nRuns
    x=double(h5read(RawFile,HdfDataPath,[1,(i-1)*nBatch+1],[Nch,nBatch+dt+1])');
    %clipped derivative (time shifted)
    h=(1-exp(-abs(x(1+dt:end,:)-x(1:end-dt,:))*diag(1./vSat))).*sign(x(1+dt:end,:)-x(1:end-dt,:));
    prc=quantile(h*diag(1./Sdev),BiasPercs,2);
    xM=min(abs(prc),[],2)*0.5.*sum(sign(prc),2)*Nch;
    %treat as driving noise of AR process (can assume zero mean)
    VBi=zeros(1,Nch);
    VBmi=0;
    %ACMi=zeros(Nch,Nch);
    %BCMi=zeros(Nch,Nch);
    Qmi=zeros(1,Nch);
    Za(1,:)=Za(end,:);
    for t=1:nBatch
        t2=mod(t,2)+1;
        t1=3-t2;
        Am(1,t2)=Am(1,t1)*ARtau;
        Am(1,t2)=Am(1,t2)+xM(t,1);
        A(t2,:)=A(t1,:)*ARtau;
        A(t2,:)=A(t2,:)+h(t,:);
        %ACMi=ACMi+A(t2,:)'*A(t2,:);
        %correction
        Za(t+1,:)=(Am(1,t2)-A(t2,:)./Sdev).*PCx;%estimates integrated version of x
        %Y will be used to correct x in a time-shifted manner (always tshift bins before)
        %need to scale with vSat (done above)
    end
    q=(x(dt+2:end,:)-x(2:end-dt,:))-(Za(2:end,:)-Za(1:end-1,:));
    hq=(1-exp(-abs(q)./vSat2)).*sign(q);
    hprc=quantile(hq*diag(1./Sdev),BiasPercs,2);
    qM=min(abs(hprc),[],2)*0.5.*sum(sign(hprc),2)*Nch;
    for t=1:nBatch
        B=B*ARtau2;
        Bm=Bm*ARtau2;
        B=B+hq(t,:);
        Bm=Bm+qM(t,1);
        %correlation coefficient:
        %need variances and avg. product
        VBi=VBi+B.^2;
        VBmi=VBmi+Bm^2;
        Qmi=Qmi+B.*(Bm-B./Sdev);
        %BCMi=BCMi+B'*B;
    end
    %want to see the temporal development
    CC=Qmi./sqrt(VBmi*VBi);
    Dm(1,i)=mean(CC);
    Dstd(1,i)=std(CC);
    Qm=Qm+Qmi;
    VBm=VBm+VBmi;
    VB=VB+VBi;
    %ACM=ACM+ACMi;
    %BCM=BCM+BCMi;
end
%see compute correlation coefficients
CCy=Qm./sqrt(VBm*VB);
%Sy=sqrt(VB);
%Sym=sqrt(VBm);
PCy=CCy.*sqrt(VB/VBm).*vSat2;%to convert from avg. to single var

%%%%%%%%%%%%%%%%%%%%%%%%%

if Filter.car.plotResults
    figVT=figure('Position',[0 0 1600 1200]);
    axCSr = axes('OuterPosition',[0   3/4 1/4 1/4]);
    axCS = axes('OuterPosition',[0   1/2 1/4 1/4]);
    axCSc = axes('OuterPosition',[0   1/4 1/4 1/4]);
    axRawV = axes('OuterPosition',[1/4 3/4 0.375 1/4]);
    axPS = axes('OuterPosition',[1/4 2/4 0.35 1/4]);
    axQS = axes('OuterPosition',[1/4 1/4 0.35 1/4]);
    axFilt1V = axes('OuterPosition',[0.625 3/4 0.375 1/4]);
    axPSHsingle = axes('OuterPosition',[0.6 2/4 0.4 1/4]);
    axPSH = axes('OuterPosition',[0.6 1/4 0.4 1/4]);
    axPSt= axes('OuterPosition',[0   0   1/2 1/4]);
    axPStd= axes('OuterPosition',[2/4 0   1/2 1/4]);
    %Voltage traces:
    hold(axRawV,'on');
    hold(axFilt1V,'on');
    hold(axPSH,'on');
    hold(axPSt,'on');
    hold(axPStd,'on');
    xlabel(axRawV,'time/ms')
    ylabel(axRawV,'voltage/µV')
    xlabel(axFilt1V,'time/ms')
    ylabel(axFilt1V,'voltage/µV')
    xlabel(axPSH,'0.5-2 kHz power/(µV^2/kHz)')
    ylabel(axPSH,'density')
    ylabel(axPSt,'0.5-2 kHz power (global)/(µV^2/kHz)')
    xlabel(axPSt,'time/s')
    ylabel(axPStd,'0.5-2 kHz power (channel)/(µV^2/kHz)')
    xlabel(axPStd,'time/s')
    %Coherence
    hold(axCSr,'on');
    hold(axCS,'on');
    hold(axCSc,'on');
    hold(axQS,'on');
    hold(axPS,'on');
    xlabel(axPS,'frequency/Hz')
    ylabel(axPS,'power (global)/(µV^2/kHz)')
    xlabel(axQS,'frequency/Hz')
    ylabel(axQS,'power (channel)/(µV^2/kHz)')
    title(axCSr,'raw')
    title(axCS,'corrected')
    title(axCSc,'correction')
end
%FFCM=zeros(Nch,Nch);
%FCM=zeros(Nch,Nch);
NPow=zeros(1,nRuns*Nfft);
%initial condition
%VCM=zeros((Nch,Nch))
if Filter.car.plotResults
    CS=complex(zeros(Nch,Nch),0);%coherence spectrum
    CSr=complex(zeros(Nch,Nch),0);%raw
    CSc=complex(zeros(Nch,Nch),0);%correction
    PS=complex(zeros(dtCS,1),0);%power spectrum (ignore zero frequency component)
    PSa=complex(zeros(dtCS,1),0);
    PSb=complex(zeros(dtCS,1),0);
    PSr=complex(zeros(dtCS,1),0);%raw
    QS=complex(zeros(dtCS,1),0);
    QSc=complex(zeros(dtCS,1),0);
    QSr=complex(zeros(dtCS,1),0);%raw
    TS=zeros(1,2000);%temporal development, nondiagonal elements
    TSc=zeros(1,2000);
    TSr=zeros(1,2000);
    TQ=zeros(1,2000);%temporal development, diagonal elements
    TQc=zeros(1,2000);
    TQr=zeros(1,2000);
    PH=zeros(200,4);
    PHdc=zeros(200*Nch,1);
    PHch=(0:(Nch-1))'*200;
    phiN=2000*recParameter.bitVolt^2/(sum(hanning(dtCS).^2)*recParameter.sRate*(fb1-fb0));
end
Am=zeros(1,2);
A=zeros(2,Nch);
Bm=0;
B=zeros(1,Nch);
Za=zeros(nBatch+1,Nch);
Zb=zeros(nBatch,Nch);

%initialize with first 3*dt frames
x=double(h5read(RawFile,HdfDataPath,[1,1],[Nch,dt+1])');
h=(1-exp(-abs(x(dt+1,:)-x(1,:))*diag(1./vSat))).*sign(x(dt+1,:)-x(1,:));
prc=quantile(h*diag(1./Sdev),BiasPercs);
Am(1,1)=min(abs(prc))*0.5*sum(sign(prc))*Nch;
%Am(0)=sum(h*1./Sdev)
A(1,:)=A(1,:)+h;%may need to reshape
%correction
Za(end,:)=(Am(1,1)-A(1,:)./Sdev).*PCx;%need to save later
for i=1:nRuns
    x=double(h5read(RawFile,HdfDataPath,[1,(i-1)*nBatch+1],[Nch,nBatch+dt+1])');
    %clipped derivative (time shifted)
    h=(1-exp(-abs(x(dt+1:end,:)-x(1:end-dt,:))*diag(1./vSat))).*sign(x(dt+1:end,:)-x(1:end-dt,:));
    prc=quantile(h*diag(1./Sdev),BiasPercs,2);
    xM=min(abs(prc),[],2)*0.5.*sum(sign(prc),2)*Nch;
    %treat as driving noise of AR process (can assume zero mean)
    Za(1,:)=Za(end,:);
    for t=1:nBatch
        t2=mod(t,2)+1;
        t1=3-t2;
        Am(1,t2)=Am(1,t1)*ARtau;
        Am(1,t2)=Am(1,t2)+xM(t,1);
        A(t2,:)=A(t1,:)*ARtau+h(t,:);
        %correction
        Za(t+1,:)=(Am(1,t2)-A(t2,:)./Sdev).*PCx;
        %Y will be used to correct x in a time-shifted manner (always tshift bins before)
        %need to scale with vSat (done above)
        %clipped derivative (time shifted by dt) --> (dt:-dt)
    end
    q=(x(dt+2:end,:)-x(2:end-dt,:))-(Za(2:end,:)-Za(1:end-1,:));
    hq=(1-exp(-abs(q)./vSat2)).*sign(q);
    hprc=quantile(hq*diag(1./Sdev),BiasPercs,2);
    qM=min(abs(hprc),[],2)*0.5.*sum(sign(hprc),2)*Nch;
    for t=1:nBatch
        B=B*ARtau2;
        Bm=Bm*ARtau2;
        B=B+hq(t,:);
        Bm=Bm+qM(t,1);
        Zb(t,:)=(Bm-B./Sdev).*PCy;
    end
    %Shared noise:
    Xcorrection=(Za(2:end,:)+Zb)/dt;
    Xnew=x(dtshift+2:dtshift+1+nBatch,:)-Xcorrection;
    %VCM+=sum(Xnew(:,:,None)*Xnew(:,None,:),axis=0)%useless
    
    if Filter.car.plotResults
        %compute periodograms
        for j=1:Nfft
            %fft
            F=fft(complex(diag(hanning(dtCS))*Xnew((j-1)*dtCS+1:j*dtCS,:)));
            %cross spectrum in frequency interval (need to normalize)
            CS=CS+F(fb0+1:fb1,:)'*F(fb0+1:fb1,:);
            QS=QS+mean(F.*conj(F),2);
            %correction term
            %fft
            Fc=fft(complex(diag(hanning(dtCS))*Xcorrection((j-1)*dtCS+1:j*dtCS,:)));
            %absolute cross spectrum in frequency interval (need to normalize)
            CSc=CSc+Fc(fb0+1:fb1,:)'*Fc(fb0+1:fb1,:);
            QSc=QSc+mean(Fc.*conj(Fc),2);
            %Raw signal
            CSr=CSr+(Fc(fb0+1:fb1,:)'+F(fb0+1:fb1,:)')*(Fc(fb0+1:fb1,:)+F(fb0+1:fb1,:));
            QSr=QSr+mean((F+Fc).*conj(F+Fc),2);
            %temporal properties
            if ((i-1)*Nfft+j)<=2000
                %roughly 1 min of data
                %average over non-diagonal elements
                cs=F(fb0+1:fb1,:)'*(F(fb0+1:fb1,:));
                phi=reshape(cs,1,Nch^2);
                phi=reshape(phi(1:end-1),Nch-1,Nch+1);
                TS((i-1)*Nfft+j)=mean(mean(abs(phi(5:end-3,:))));
                TQ((i-1)*Nfft+j)=mean(abs(diag(cs)));
                cs=Fc(fb0+1:fb1,:)'*(Fc(fb0+1:fb1,:));
                phi=reshape(cs,1,Nch^2);
                phi=reshape(phi(1:end-1),Nch-1,Nch+1);
                TSc((i-1)*Nfft+j)=mean(mean(abs(phi(5:end-3,:))));
                TQc((i-1)*Nfft+j)=mean(abs(diag(cs)));
                cs=(F(fb0+1:fb1,:)'+Fc(fb0+1:fb1,:)')*(F(fb0+1:fb1,:)+Fc(fb0+1:fb1,:));
                phi=reshape(cs,1,Nch^2);
                phi=reshape(phi(1:end-1),Nch-1,Nch+1);
                TSr((i-1)*Nfft+j)=mean(mean(abs(phi(5:end-3,:))));
                TQr((i-1)*Nfft+j)=mean(abs(diag(cs)));
            end
            if PowerHisto
                cs=F(fb0+1:fb1,:)'*F(fb0+1:fb1,:);
                phi=reshape(cs,1,Nch^2);
                phi=reshape(phi(1:end-1),Nch-1,Nch+1);
                phi=mean(mean(abs(phi(5:end-3,:))));
                pInd=floor(max(min((log10(phi*phiN))*50,200),1));
                PH(pInd,2)=PH(pInd,2)+1;
                phi=abs(diag(cs));
                pInd=floor(max(min((log10(mean(phi)*phiN))*50,200),1));
                PH(pInd,4)=PH(pInd,4)+1;
                pInd=floor(max(min((log10(phi*phiN))*50,200),1))+PHch;
                PHdc(pInd,1)=PHdc(pInd,1)+1;
                cs=(F(fb0+1:fb1,:)'+Fc(fb0+1:fb1,:)')*(F(fb0+1:fb1,:)+Fc(fb0+1:fb1,:));
                phi=reshape(cs,1,Nch^2);
                phi=reshape(phi(1:end-1),Nch-1,Nch+1);
                phi=mean(mean(abs(phi(5:end-3,:))));
                pInd=floor(max(min((log10(phi*phiN))*50,200),1));
                PH(pInd,1)=PH(pInd,1)+1;
                phi=mean(abs(diag(cs)));
                pInd=floor(max(min((log10(phi*phiN))*50,200),1));
                PH(pInd,3)=PH(pInd,3)+1;
            end
            cs=Fc(fb0+1:fb1,:)'*Fc(fb0+1:fb1,:);
            phi=reshape(cs,1,Nch^2);
            phi=reshape(phi(1:end-1),Nch-1,Nch+1);
            NPow((i-1)*Nfft+j)=mean(mean(abs(phi(5:end-3,:))));
            %power spectrum
            %correction terms
            Pa=fft(complex(mean(Za(2+(j-1)*dtCS:1+j*dtCS,:)/dt,2).*hanning(dtCS)));
            PSa=PSa+Pa.*conj(Pa);
            Pb=fft(complex(mean(Zb(1+(j-1)*dtCS:j*dtCS,:)/dt,2).*hanning(dtCS)));
            PSb=PSb+Pb.*conj(Pb);
            %raw
            Pr=fft(complex(mean(x(dtshift+2+(j-1)*dtCS:dtshift+1+j*dtCS,:),2).*hanning(dtCS)));
            PSr=PSr+Pr.*conj(Pr);
            %new
            P=fft(complex(mean(Xnew((j-1)*dtCS+1:j*dtCS,:),2).*hanning(dtCS)));
            PS=PS+P.*conj(P);
        end
        %FCM=FCM+Za(2:nBatch+1,:)'*Za(2:nBatch+1,:);
        %FFCM=FFCM+Zb(1:nBatch,:)'*Zb(1:nBatch,:);
    end
    
    h5write(OutFile,'/data',int16(round(Xnew)'),[1,(i-1)*nBatch+dtshift+2],[Nch,nBatch]);
    
    
    if (i==3) && Filter.car.plotResults
        %voltage traces
        %Raw
        plot(axRawV,(1:1000)*1000/recParameter.sRate,x(dtshift+2:dtshift+1001,PlotChannel)*recParameter.bitVolt+30,'k-')
        %axRawV.plot(x(dtshift+1:dtshift+1001,PlotChannel)-Za(1:1001,PlotChannel)*1./dt,'b-')
        plot(axRawV,(1:1000)*1000/recParameter.sRate,Xnew(1:1000,PlotChannel)*recParameter.bitVolt-30,'r-')
        plot(axFilt1V,(1:1000)*1000/recParameter.sRate,mean(x(dtshift+2:dtshift+1001,:),2)*recParameter.bitVolt,'g')
        plot(axFilt1V,(1:1000)*1000/recParameter.sRate,Za(2:1001,PlotChannel)*recParameter.bitVolt/(dt),'b')
        plot(axFilt1V,(1:1000)*1000/recParameter.sRate,Zb(1:1000,PlotChannel)*recParameter.bitVolt/(dt),'c')
    end
end

if i*nBatch+dtshift+1<recParameter.LenRec
    x=double(h5read(RawFile,HdfDataPath,[1,i*nBatch+dtshift+2],[Nch,recParameter.LenRec-i*nBatch-dtshift-1])');
    h5write(OutFile,'/data',int16(round(x)'),[1,i*nBatch+dtshift+2],[Nch,recParameter.LenRec-i*nBatch-dtshift-1]);
end

if Filter.car.plotResults
    %degrees of freedom
    %%%%%nu=sum(hanning(dtCS).^2).^2/sum(hanning(dtCS).^4)*2*(fb1-fb0)*nRuns*Nfft;
    %normalization factor
    %recParameter.sRate/2000/dtCS: kHz per bin, need factor of 2 for negative symmetrical part
    Snorm=2000*recParameter.bitVolt^2/(sum(hanning(dtCS).^2)*recParameter.sRate*nRuns*Nfft);
    SnormT=2000*recParameter.bitVolt^2/(sum(hanning(dtCS).^2)*recParameter.sRate*(fb1-fb0));
    %Perc=(1./scipy.stats.chi2.isf(array((0.005,0.995)),nu))*nu
    %Bonferroni correction
    %%%%%q=0.01/(Nch^2+Nch);
    %%%%%PercB=(1./scipy.stats.chi2.isf(array((q,1.-q)),nu))*nu
    %%%%%figCS.text(0.75,0.6,'confidence interval',ha='left',va='center',fontsize=10)
    %%%%%figCS.text(0.75,0.55,'(%1.7f, %1.7f)'%(PercB(0),PercB(1)),ha='left',va='center',fontsize=10)
    imagesc(axCS,log10(min(max(abs(CS)*Snorm/(fb1-fb0),1),1e4)))
    %norm=mplt.colors.LogNorm(5e1,1e5))
    hC=colorbar(axCS);
    set(hC,'Ytick',0:3,'YTicklabel',{'1' '1e1' '1e2' '1e3'});
    hC.Label.String = 'power*kHz/µV^2';
    ylim(axCS,[0.5 Nch+0.5])
    xlim(axCS,[0.5 Nch+0.5])
    imagesc(axCSr,log10(min(max(abs(CSr)*Snorm/(fb1-fb0),1),1e4)))
    %norm=mplt.colors.LogNorm(5e1,1e5))
    hC=colorbar(axCSr);
    set(hC,'Ytick',0:4,'YTicklabel',{'1' '1e1' '1e2' '1e3' '1e4'});
    hC.Label.String = 'power*kHz/µV^2';
    ylim(axCSr,[0.5 Nch+0.5])
    xlim(axCSr,[0.5 Nch+0.5])
    imagesc(axCSc,log10(min(max(abs(CSc)*Snorm/(fb1-fb0),0.1),1e2)))
    %norm=mplt.colors.LogNorm(1e1,1e4))
    hC=colorbar(axCSc);
    set(hC,'Ytick',-1:2,'YTicklabel',{'0.1' '1' '10' '100'});
    hC.Label.String = 'power*kHz/µV^2';
    ylim(axCSc,[0.5 Nch+0.5])
    xlim(axCSc,[0.5 Nch+0.5])
    semilogy(axPS,(1:dtCS/2+1)*recParameter.sRate/dtCS,abs(PS(1:dtCS/2+1))*Snorm,'r-')
    semilogy(axPS,(1:dtCS/2+1)*recParameter.sRate/dtCS,abs(PSa(1:dtCS/2+1))*Snorm,'b-')
    semilogy(axPS,(1:dtCS/2+1)*recParameter.sRate/dtCS,abs(PSb(1:dtCS/2+1))*Snorm,'c-')
    semilogy(axPS,(1:dtCS/2+1)*recParameter.sRate/dtCS,abs(PSr(1:dtCS/2+1))*Snorm,'k-')
    axPS.YScale='log';
    xlim(axPS,[100 4000])
    ylim(axPS,[1e-1 1e4])
    semilogy(axQS,(1:dtCS/2+1)*recParameter.sRate/dtCS,abs(QS(1:dtCS/2+1))*Snorm,'r-')
    semilogy(axQS,(1:dtCS/2+1)*recParameter.sRate/dtCS,abs(QSc(1:dtCS/2+1))*Snorm,'b-')
    semilogy(axQS,(1:dtCS/2+1)*recParameter.sRate/dtCS,abs(QSr(1:dtCS/2+1))*Snorm,'k-')
    axQS.YScale='log';%bug in MATLAB
    xlim(axQS,[100 4000])
    ylim(axQS,[1 1e5])
    semilogy(axPSt,(1:2000)*dtCS/recParameter.sRate,TSr*SnormT,'k-')
    semilogy(axPSt,(1:2000)*dtCS/recParameter.sRate,TS*SnormT,'r-')
    axPSt.YScale='log';
    ylim(axPSt,[2 2e3])
    semilogy(axPStd,(1:2000)*dtCS/recParameter.sRate,TQr*SnormT,'k-')
    semilogy(axPStd,(1:2000)*dtCS/recParameter.sRate,TQ*SnormT,'r-')
    semilogy(axPStd,(1:2000)*dtCS/recParameter.sRate,TQc*SnormT,'b-')
    axPStd.YScale='log';
    ylim(axPStd,[1 5e3])
    if PowerHisto
        %normalize to density
        PH=PH*diag(1./sum(PH,1));
        %plot
        loglog(axPSH,logspace(0,4,200),PH(:,1),'k-')
        loglog(axPSH,logspace(0,4,200),PH(:,2),'r-')
        loglog(axPSH,logspace(0,4,200),PH(:,3),'k:')
        loglog(axPSH,logspace(0,4,200),PH(:,4),'r:')
        axPSH.YScale='log';
        axPSH.XScale='log';
        xlim(axPSH,[1 1e4])
        ylim(axPSH,[1e-4 1])
        PHdc=reshape(PHdc,200,Nch);
        PHdc=PHdc*diag(1./sum(PHdc,1));
        imagesc(axPSHsingle,[0 4],[0 Nch],log10(max(PHdc',1e-4)))
        set(axPSHsingle,'Xtick',0:4,'XTicklabel',{'1' '1e1' '1e2' '1e3' '1e4'});
        hC=colorbar(axPSHsingle);
        set(hC,'Ytick',-4:-1,'YTicklabel',{'1e-4' '1e-3' '1e-2' '1e-1'});
        hC.Label.String = 'density';
    end
    saveas(figVT,[Filter.car.PlotBase filesep Filter.car.plotFolder filesep Filter.car.plotName])
    close(figVT)
    Filter.car.powerHist=PH;
    Filter.car.powerHistSingle=PHdc;
    Filter.car.NoisePower=NPow';
end
end
            
