function blindTemplateMatchingGPU(Origin,Target,bTM,recParameter)
%threshold based spike detection and characterization
%read file in batches (single channel)
g={};
g.bitVolt=recParameter.bitVolt;
g.LenRec=recParameter.LenRec;
g.Nch=recParameter.Nch;


%h=load(MatFileCAR);
%[~,Vmx]=max(conv2(Filter.car.powerHistSingle,hanning(7),'same'),[],1);

%g.varCh=10.^(Vmx/100);
g.varCh=bTM.Noise.sigmaADC;
g.lowVar=bTM.Noise.lowVar;
g.lowVarDT=bTM.Noise.lowVarDT;
g.nBatch=900000;%increase if GPU memory allows 18
g.Sampling=30;%kHz
nRuns=ceil(g.LenRec/g.nBatch);



%do everything in variance-normalized units, but convert back the
%on-channel data to µV
g.LMC=zeros(floor(g.LenRec/5),g.Nch);
g.LME=cell(g.Nch,1);
GPUfilt=blindTM.FilterGPUnl(g.nBatch,g.varCh,g.LenRec,g.Sampling,bTM.Noise.ACF,g.lowVar,g.lowVarDT);
g.Ntail=GPUfilt.Ntail;
g.Nwidth=GPUfilt.Nwidth;
g.Namp=GPUfilt.Namp;
g.Nphase=GPUfilt.Nphase;
g.Ttail=GPUfilt.Ttail;
g.Twidth=GPUfilt.Twidth;
g.Xamp=GPUfilt.Xamp;
g.templateBank=GPUfilt.templateBank;
g.BiasTail=GPUfilt.BiasTail;
g.templateScale=sqrt(pi./g.Twidth);
SHist=zeros(g.Namp*g.Nwidth*g.Ntail,g.Nch);
SHist0=zeros(g.Namp*g.Nwidth*g.Ntail,g.Nch);
g.dt0=GPUfilt.dt0;
%g.bV=zeros(ceil(g.LenRec/g.dt0),g.Nch);
for ii=1:g.Nch
    disp(ii)
    %GPUfilt = GPUfilt.scaleWithACF(ii);
    %GPUfilt.varCh=g.varCh(ii);
    SpkOld=zeros(GPUfilt.dt2h,3,'double','gpuArray');
    g.LME{ii}=[];
    for i=1:nRuns-1
        [LocMaxContinous,LocMaxEvents,LocMaxHist,LocMaxHist0,SpkOld]=GPUfilt.FiltSection(Origin,ii,i,SpkOld);
        g.LMC((i-1)*g.nBatch/5+1:i*g.nBatch/5,ii)=LocMaxContinous;
        g.LME{ii}=[g.LME{ii};LocMaxEvents];
        SHist(:,ii)=SHist(:,ii)+LocMaxHist;
        SHist0(:,ii)=SHist0(:,ii)+LocMaxHist0;
    end
    i=nRuns;
    [LocMaxContinous,LocMaxEvents,LocMaxHist,LocMaxHist0,~]=GPUfilt.FiltSection(Origin,ii,i,SpkOld);
    g.LMC(end-length(LocMaxContinous)+1:end,ii)=LocMaxContinous;
    g.LME{ii}=[g.LME{ii};LocMaxEvents];
    SHist(:,ii)=SHist(:,ii)+LocMaxHist;
    SHist0(:,ii)=SHist0(:,ii)+LocMaxHist0;
end

g.SpkHist=reshape(SHist,g.Ntail,g.Nwidth,g.Namp,g.Nch);
g.SpkHist0=reshape(SHist0,g.Ntail,g.Nwidth,g.Namp,g.Nch);

if bTM.plotResults
    ChMask=bTM.ChMap;
    xMap=[3 3 3 3 4 4 4 4 6 6 6 6 5 5 5 5 4 4 4 4 3 3 3 3 1 1 1 1 2 2 2 2];
    xMap=[xMap+4 xMap];
    yMap=[11 13 15 17 16 18 20 22 22 20 18 16 17 15 13 11 12 10 8 6 7 5 3 1 1 3 5 7 6 8 10 12];
    yMap=[yMap yMap];

    xPos=(0:9)*0.1+0.005;
    yPos=(21:-1:0)*0.042+0.01;
    fig1=figure('Position',[0 0 1600 1000]);
    %width
    for i=1:g.Nch
        ax1 = axes('Position',[xPos(xMap(ChMask(i))) yPos(yMap(ChMask(i))) 0.05 0.082]);
        imagesc(ax1,log10(squeeze(sum(g.SpkHist0(2:end-1,:,:,i),1)))');
        ax1.CLim=[0 5];
        xticks(ax1,[])
        yticks(ax1,[])
        ax1.YDir='normal';
    end
    %tail
    for i=1:g.Nch
        ax1 = axes('Position',[xPos(xMap(ChMask(i)))+0.05 yPos(yMap(ChMask(i))) 0.04 0.082]);
        imagesc(ax1,log10(squeeze(sum(g.SpkHist0(:,2:end-1,:,i),2)))');
        colormap(ax1,'hot')
        ax1.CLim=[0 5];
        xticks(ax1,[])
        yticks(ax1,[])
        ax1.YDir='normal';
    end
    ax1 = axes('OuterPosition',[xPos(8) yPos(5)+0.05 0.28 0.21]);
    imagesc(ax1,log10(squeeze(sum(sum(g.SpkHist0(2:end-1,:,:,:),4),1)))');
    ax1.CLim=[0 5];
    xticks(ax1,2.5:8:26.5)
    xticklabels(ax1,round(g.Twidth(3:8:27)*10)/10)
    xticks(ax1,interp1(g.Twidth,1:g.Nwidth,[0.1 0.2 0.3 0.5]))
    xticklabels(ax1,[0.1 0.2 0.3 0.5])
    xlabel(ax1,'width/ms')
    yticks(ax1,interp1([g.Xamp 10000],1:length(g.Xamp)+1,[2:10 15:5:50 60:10:100 125:25:200])-0.5)
    xTl={'2', '', '', '5', '', '', '', '', '10', '', '', '25', '', '', '', '', '50',...
        '', '', '', '', '100', '', '', '', '200'};
    yticklabels(ax1,xTl)
    ylabel(ax1,'amplitude/SD')
    ax1.YDir='normal';
    hC=colorbar(ax1,'Ticks',[1 3 5],...
        'TickLabels',{'10','1e3','1e5'});
    hC.Label.String = 'count';
    %tail
    ax1 = axes('OuterPosition',[xPos(1)-0.005 yPos(22)-0.01 0.28 0.21]);
    imagesc(ax1,log10(squeeze(sum(sum(g.SpkHist0(:,2:end-1,:,:),4),2)))');
    colormap(ax1,'hot')
    ax1.CLim=[0 5];
    xticks(ax1,1:5:16)
    xticklabels(ax1,round(g.Ttail(1:5:16)*100)/100)
    xlabel(ax1,'frac. of tail')
    %yticks(ax1,interp1(g.Xamp(13:end),1:length(g.Xamp),[3 10 25 80 300])-0.5)
    %yticklabels(ax1,[3 10 25 80 300])
    yticks(ax1,interp1([g.Xamp 10000],1:length(g.Xamp)+1,[2:10 15:5:50 60:10:100 125:25:200])-0.5)
    xTl={'2', '', '', '5', '', '', '', '', '10', '', '', '25', '', '', '', '', '50',...
        '', '', '', '', '100', '', '', '', '200'};
    yticklabels(ax1,xTl)
    ylabel(ax1,'amplitude/SD')
    ax1.YDir='normal';
    hC=colorbar(ax1,'Ticks',[1 3 5],...
        'TickLabels',{'10','1e3','1e5'});
    hC.Label.String = 'count';
    saveas(fig1,[bTM.PlotBase filesep bTM.plotFolder filesep bTM.plotFile])
    close(fig1)
end
%save
save(Target,'g','-v7.3')
end

