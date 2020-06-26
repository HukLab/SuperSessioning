function blindTemplateMatchingGPU(Origin,Target,bTM)
%threshold based spike detection and characterization
%read file in batches (single channel)
g={};
g.bitVolt=h5read(Origin,'/recordings/0/application_data/channel_bit_volts');
a=h5info(Origin);
DataSize=a.Groups.Groups.Datasets.Dataspace.Size;
g.LenRec=floor(DataSize(1,2));
g.Nch=floor(DataSize(1,1));


%h=load(MatFileCAR);
%[~,Vmx]=max(conv2(Filter.car.powerHistSingle,hanning(7),'same'),[],1);

%g.varCh=10.^(Vmx/100);
g.varCh=bTM.Noise;
g.nBatch=3600000;
g.Sampling=30;%kHz
nRuns=ceil(g.LenRec/g.nBatch);



%do everything in variance-normalized units, but convert back the
%on-channel data to µV
g.LMC=zeros(floor(g.LenRec/5),g.Nch);
g.LME=cell(g.Nch,1);
GPUfilt=FilterGPUdd(g.nBatch,g.varCh(1),g.LenRec,g.Sampling);
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
g.bV=zeros(ceil(g.LenRec/g.dt0),g.Nch);
for ii=1:g.Nch
    disp(ii)
    GPUfilt.varCh=g.varCh(ii);
    g.LME{ii}=[];
    for i=1:nRuns-1
        [LocMaxContinous,LocMaxEvents,LocMaxHist,LocMaxHist0,xbV]=GPUfilt.FiltSection(Origin,ii,i);
        g.LMC((i-1)*g.nBatch/5+1:i*g.nBatch/5,ii)=LocMaxContinous;
        g.LME{ii}=[g.LME{ii};LocMaxEvents];
        g.bV((i-1)*g.nBatch/g.dt0+1:i*g.nBatch/g.dt0,ii)=xbV;
        SHist(:,ii)=SHist(:,ii)+LocMaxHist;
        SHist0(:,ii)=SHist0(:,ii)+LocMaxHist0;
    end
    i=nRuns;
    [LocMaxContinous,LocMaxEvents,LocMaxHist,LocMaxHist0,xbV]=GPUfilt.FiltSection(Origin,ii,i);
    g.LMC(end-length(LocMaxContinous)+1:end,ii)=LocMaxContinous;
    g.LME{ii}=[g.LME{ii};LocMaxEvents];
    g.bV((i-1)*g.nBatch/g.dt0+1:end,ii)=xbV(1:end-1);
    SHist(:,ii)=SHist(:,ii)+LocMaxHist;
    SHist0(:,ii)=SHist0(:,ii)+LocMaxHist0;
end

g.SpkHist=reshape(SHist,g.Ntail,g.Nwidth,g.Namp,g.Nch);
g.SpkHist0=reshape(SHist0,g.Ntail,g.Nwidth,g.Namp,g.Nch);

if bTM.plotResults
    ChMask=1:g.Nch;
    xMap=[4 4 4 4 3 3 3 3 1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 6 6 6 6 5 5 5 5];
    xMap=[xMap+4 xMap];
    yMap=[12 10 8 6 7 5 3 1 1 3 5 7 6 8 10 12 11 13 15 17 16 18 20 22 22 20 18 16 17 15 13 11];
    yMap=[yMap yMap];

    xPos=(0:9)*0.1+0.005;
    yPos=(21:-1:0)*0.042+0.01;
    fig1=figure('Position',[0 0 1600 1000]);
    %width
    for i=1:g.Nch
        ax1 = axes('Position',[xPos(xMap(ChMask(i))) yPos(yMap(ChMask(i))) 0.05 0.082]);
        imagesc(ax1,log10(squeeze(sum(g.SpkHist(2:end-1,:,:,i),1)))');
        ax1.CLim=[0 5];
        xticks(ax1,[])
        yticks(ax1,[])
    end
    %tail
    for i=1:g.Nch
        ax1 = axes('Position',[xPos(xMap(ChMask(i)))+0.05 yPos(yMap(ChMask(i))) 0.04 0.082]);
        imagesc(ax1,log10(squeeze(sum(g.SpkHist(:,2:end-1,:,i),2)))');
        colormap(ax1,'hot')
        ax1.CLim=[0 5];
        xticks(ax1,[])
        yticks(ax1,[])
    end
    ax1 = axes('OuterPosition',[xPos(8) yPos(5)+0.05 0.28 0.21]);
    imagesc(ax1,log10(squeeze(sum(sum(g.SpkHist(2:end-1,:,:,:),4),1)))');
    ax1.CLim=[0 5];
    xticks(ax1,2.5:8:26.5)
    xticklabels(ax1,round(g.Twidth(3:8:27)*10)/10)
    xlabel(ax1,'width/ms')
    yticks(ax1,[0.5 32.5 64.5])
    yticklabels(ax1,g.Xamp(1:32:65))
    ylabel(ax1,'power/stdev')
    hC=colorbar(ax1,'Ticks',[1 3 5],...
        'TickLabels',{'10','1e3','1e5'});
    hC.Label.String = 'count';
    %tail
    ax1 = axes('OuterPosition',[xPos(1)-0.005 yPos(22)-0.01 0.28 0.21]);
    imagesc(ax1,log10(squeeze(sum(sum(g.SpkHist(:,2:end-1,:,:),4),2)))');
    colormap(ax1,'hot')
    ax1.CLim=[0 5];
    xticks(ax1,1:5:16)
    xticklabels(ax1,round(g.Ttail(1:5:16)*100)/100)
    xlabel(ax1,'frac. of tail')
    yticks(ax1,[0.5 32.5 64.5])
    yticklabels(ax1,g.Xamp(1:32:65))
    ylabel(ax1,'power/stdev')
    hC=colorbar(ax1,'Ticks',[1 3 5],...
        'TickLabels',{'10','1e3','1e5'});
    hC.Label.String = 'count';
    saveas(fig1,[bTM.plotFolder filesep bTM.plotFile])
    close(fig1)
end
%save
save(Target,'g','-v7.3')
end

