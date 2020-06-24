function clusterSession(Origin, Target, Clust)

disp(File)
%Parameter
%smoothing (half kernel)
nsih=10;% ISI, samples
nsth=2;% tail
nswh=2;% width
nsah=2;% amplitude

nsi=2*nsih+1;
nst=2*nsth+1;
nsw=2*nswh+1;
nsa=2*nsah+1;

hwi=hanning(nsi);
hwt=hanning(nst);
hww=hanning(nsw);
hwa=hanning(nsa);

%regularization
Bsc=0.000000001;%baseline spike rate (false positives!)

Jthreshold=0.3*log(2);

if ~isfolder(plotBase)
    mkdir(plotBase)
end

%% load spike data from template matching
p=load(Origin);
g=p.g;
spikeTimes=g.LME;
spikeMask=convn([g.bV;g.bV(end,:)],[1;1;1]/3,'same');%floor(t/g.dt0)+1,g.Nch

Ckernel=permute(hwt,[1 2 3]).*permute(hww,[2 1 3]).*permute(hwa,[2 3 1]);
Cks=sum(Ckernel,'all');


%% prepare output
m=struct();
m.RecId=File;
m.RawFile=[File '_CAR.kwd'];
m.KSunits=(1:g.Nch)';
m.Sampling=g.Sampling*1000;
m.LenRec=g.LenRec;
m.Bsc=Bsc;
m.nsTWA=[nst nsw nsa];
m.nsTWAh=[nsth nswh nsah];
m.Jthreshold=Jthreshold;
m.Ckernel=Ckernel;
m.Cks=Cks;
m.HzNorm=m.Sampling/m.LenRec;
m.xPwr=g.Xamp';
m.rPwr=g.Xamp(2)/g.Xamp(1);
m.Nch=g.Nch;
m.Namp=g.Namp;
m.Nwidth=g.Nwidth;
m.Ntail=g.Ntail;
m.noiseThr=25;
m.noisedt=g.dt0;

m.aPwr=reshape(repmat(permute(conv(m.xPwr(3:end,1),hwa,'full'),[2 3 1]),m.Ntail+2*nsth,m.Nwidth+2*nswh,1),1,[]);
m.tWidth=g.Twidth';
m.tTail=g.Ttail';

m.noiseVar=g.bV;
m.fracNoise=[];
m.noiseCh=[];
m.noiseSpikes={};
m.noiseTimes={};
m.noiseMask={};
m.Channel=[];
m.Spikes={};
m.Times={};
m.shapeMask={};

%Spike Histograms, units in Hz
m.NullHisto=g.SpkHist*m.HzNorm;

%% get channels with max amplitude, segment into possible subclusters

m.Nsmth=[m.Namp+2*nsah-1 m.Nwidth+2*nswh m.Ntail+2*nsth];

Ind=0;
noiseInd=0;
m.subTimes={};
BorderMask=zeros(m.Ntail,m.Nwidth,m.Namp,'logical');
%BorderMask(1:end-1,2:end,:)=true;
BorderMask(2:end-1,2:end-1,2:end-1)=true;
BorderKernel=permute([0.5;1;0.5],[1 2 3]).*permute([0.5;1;0.5],[2 1 3]).*permute([0.5;1;0.5],[2 3 1]);
for i=1:m.Nch
    sTI=(spikeTimes{i}(:,1)>0);
    sT=max(spikeTimes{i}(sTI,2),0);
    SpkSH=spikeTimes{i}(sTI,1)';
    SHqp=reshape(histcounts(SpkSH(1,spikeMask(floor(sT/m.noisedt)+1,i)<m.noiseThr),1:m.Namp*m.Nwidth*m.Ntail+1),m.Ntail,m.Nwidth,m.Namp)*m.HzNorm;
    %SHqp=squeeze(g.SpkHist0(:,:,:,i))*m.HzNorm;
    %SH=reshape(histcounts(SpkSH,1:m.Namp*m.Nwidth*m.Ntail+1),m.Ntail,m.Nwidth,m.Namp)*m.HzNorm;
    SH=squeeze(g.SpkHist(:,:,:,i))*m.HzNorm;
    %split with watershed
    [L,Lpwt,Spwt,Lspk]=NspikeRad(SHqp,SH,m.Ntail,m.Nwidth,m.Namp);
    %Lmean=sum(SH.*AmpMat,'all')/sum(SH,'all');
    %minimum firing rate, no boundary cluster.
    Lmsk=find((Lspk>0.05).*(Lpwt(:,1)-0.5*Spwt(:,1)-Lpwt(:,2)/4>3).*...
        (Lpwt(:,2)-0.5*Spwt(:,2)>1).*(Lpwt(:,2)+0.5*Spwt(:,2)<m.Nwidth-1).*...
        (Lpwt(:,3)-0.5*Spwt(:,3)>1).*(Lpwt(:,3)+0.5*Spwt(:,3)<m.Ntail));
    %compute relative density at cluster boundaries
    Ncx=length(Lmsk);
    fracBorder=zeros(Ncx,1);
    Lunits=ones(m.Ntail,m.Nwidth,m.Namp,'logical');
    for q=1:Ncx
        Border=(convn(L==Lmsk(q),BorderKernel,'same')<6.1).*BorderMask.*(L==Lmsk(q));
        Nborder=sum(SH.*Border,'all')/sum(Border,'all');
        Nq=sum(SH.*(L==Lmsk(q)).*BorderMask,'all')/sum((L==Lmsk(q)).*BorderMask,'all');
        fracBorder(q,1)=Nborder/Nq;
        if fracBorder(q,1)<2/3
            Lunits=Lunits & (L~=Lmsk(q));
        end
    end
    %find spikes that match template clusters
    Ncx=length(Lmsk);
    if Ncx>0
        Msk=zeros(m.Ntail*m.Nwidth*m.Namp,Ncx);
        for j=1:Ncx
            Msk(:,j)=reshape(convn(1*(L==Lmsk(j)),Cks,'same'),[],1);
        end
        [SpkMx,SpkSorted]=max(Msk(min(max(SpkSH(1,:),1),m.Ntail*m.Nwidth*m.Namp),:),[],2);
        %assign clusters to spikes (0: no cluster found)
        SpkClust=(1*(SpkMx>0)).*SpkSorted;
    else
        %Msk=zeros(m.Ntail*m.Nwidth*m.Namp,1);
        SpkClust=zeros(1,size(SpkSH,2));
    end
    Lunits=reshape(Lunits,[],1);
    SpkHash=Lunits(min(max(SpkSH(1,:),1),m.Ntail*m.Nwidth*m.Namp),1);
    %new clusters
    for k=1:Ncx
        h=sum(spikeMask(floor(sT(SpkClust==k)/g.dt0)+1,i)>m.noiseThr)/...
            size(sT(SpkClust==k),1);
        if h<0.5
            Ind=Ind+1;
            m.fracNoise(Ind,1)=h;
            m.Channel(Ind,1)=i;
            [~,Tind]=sort(sT(SpkClust==k));
            Spk=SpkSH(1,SpkClust==k)';
            m.Spikes{Ind,1}=Spk(Tind);
            Tms=sT(SpkClust==k);
            m.Times{Ind,1}=Tms(Tind);
            m.fracBorder(Ind,1)=fracBorder(k,1);
            m.shapeMask{Ind,1}=(L==Lmsk(k));
            m.Nchan(Ind,1)=1;
            m.ChSorted(Ind,1)=1;
            SH1=histcounts(m.Spikes{Ind,1},1:m.Namp*m.Nwidth*m.Ntail+1);
            SH1=reshape(SH1,m.Ntail,m.Nwidth,m.Namp);
            m.Hist{Ind,1}=SH1;
            %SHbay(m.Nwidth*m.Ntail+1:end,k)=reshape(convn(SH1(:,:,2:end),CkernelXL,'same'),[],1);
            SH0=convn(SH1(:,:,2:end),Ckernel,'full')/Cks;
            %sValid(Ind,1)=sum(SH0,'all')/length(m.Spikes{Ind,1});no
            %need?
            m.SH{Ind,1}=SH0;
            m.nSH(Ind,1)=max(sum(SH0,'all'),1);
            m.rSH{Ind,1}=SH0/max(sum(SH0,'all'),Bsc);
            %m.rSHm{Ind,1}=max(SH0,Bsc)/sum(max(SH0,Bsc));
            m.pwt(Ind,:)=Lpwt(Lmsk(k),:);
            m.Spwt(Ind,:)=Spwt(Lmsk(k),:);
            m.nSpk(Ind,1)=Lspk(Lmsk(k));
            %compute ISI histograms (approx using upto 10 intervals)
            ISIh=histcounts([diff(m.Times{Ind,1})' m.Times{Ind,1}(3:end)'-m.Times{Ind,1}(1:end-2)'],1:2:301);
            %smoothen with a hanning window, reflective boundaries
            ISI0=conv([ISIh(1,nsih:-1:1) ISIh ISIh(1,end:-1:end-nsih+1)]',hwi,'valid');
            m.nISI(Ind,1)=sum(ISIh);
            m.nISI(Ind,2)=max(length(m.Times{Ind,1}),1);
            m.rISI(Ind,:)=ISI0/m.LenRec;
            m.rISIm(Ind,:)=max(ISI0,Bsc)/sum(max(ISI0,Bsc));
        else
            noiseInd=noiseInd+1;
            m.noiseCh(noiseInd)=i;
            [~,Tind]=sort(sT(SpkClust==k));
            Spk=SpkSH(1,SpkClust==k)';
            m.noiseSpikes{noiseInd,1}=Spk(Tind);
            Tms=sT(SpkClust==k);
            m.noiseTimes{noiseInd,1}=Tms(Tind);
            m.noiseMask{noiseInd,1}=(L==Lmsk(k));
            m.noisepwt(noiseInd,:)=Lpwt(Lmsk(k),:);
        end
    end
    %deal with false alarm spikes
    m.hashSpikes{i,1}=SpkSH(1,SpkHash)';
    m.hashTimes{i,1}=sT(SpkHash,1);
end
m.nUnits=size(m.Channel,1);

%plotting
if strcmp(subject,'Jo')
    m.Nch=57;
    ChMask=1:m.Nch;
    xMap=[4 4 4 4 3 3 3 3 1 1 1 1 2 3 3 3 3 4 4 4 4 6 6 6 6 5 5 5 5];
    xMap2=[3 3 3 3 1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 6 6 6 6 5 5 5 5];
    xMap=[xMap xMap2+4];
    yMap=[12 10 8 6 7 5 3 1 1 3 5 7 12 11 13 15 17 16 18 20 22 22 20 18 16 17 15 13 11];
    yMap2=[7 5 3 1 1 3 5 7 6 8 10 12 11 13 15 17 16 18 20 22 22 20 18 16 17 15 13 11];
    yMap=[yMap yMap2];
else
    m.Nch=64;
    ChMask=1:m.Nch;
    xMap=[3 3 3 3 4 4 4 4 6 6 6 6 5 5 5 5 4 4 4 4 3 3 3 3 1 1 1 1 2 2 2 2];
    xMap=[xMap+4 xMap];
    %yMap=[12 10 8 6 7 5 3 1 1 3 5 7 6 8 10 12 11 13 15 17 16 18 20 22 22 20 18 16 17 15 13 11];
    yMap=[11 13 15 17 16 18 20 22 22 20 18 16 17 15 13 11 12 10 8 6 7 5 3 1 1 3 5 7 6 8 10 12];
    yMap=[yMap yMap];
end

%%

xPos=(0:9)*0.1+0.005;
yPos=(21:-1:0)*0.042+0.01;
fig1=figure('Position',[2560 440 1600 1000]);
%CM=colormap('parula');
%width
for i=1:m.Nch
    ax1 = axes('Position',[xPos(xMap(ChMask(i))) yPos(yMap(ChMask(i))) 0.05 0.082]);
    imagesc(ax1,log10(squeeze(sum(m.NullHisto(2:end-1,:,:,i),1)))');
    colormap(ax1,'gray')
    hold(ax1,'on')
    for k=find(m.noiseCh==i)'
        plot(ax1,m.noisepwt(k,2),m.noisepwt(k,1),'wx')
    end
    for k=find(m.Channel==i)'
        if m.fracBorder(k,1)<0.6
            plot(ax1,m.pwt(k,2),m.pwt(k,1),'mo')
        else
            plot(ax1,m.pwt(k,2),m.pwt(k,1),'co')
        end
    end
    ax1.CLim=[-4 2];
    xticks(ax1,[])
    yticks(ax1,[])
    ax1.YDir='normal';
end
%tail
for i=1:m.Nch
    ax1 = axes('Position',[xPos(xMap(ChMask(i)))+0.05 yPos(yMap(ChMask(i))) 0.04 0.082]);
    imagesc(ax1,log10(squeeze(sum(m.NullHisto(:,2:end-1,:,i),2)))');
    colormap(ax1,'copper')
    hold(ax1,'on')
    for k=find(m.noiseCh==i)'
        plot(ax1,m.noisepwt(k,3),m.noisepwt(k,1),'wx')
    end
    for k=find(m.Channel==i)'
        if m.fracBorder(k,1)<0.6
            plot(ax1,m.pwt(k,3),m.pwt(k,1),'mo')
        else
            plot(ax1,m.pwt(k,3),m.pwt(k,1),'co')
        end
    end
    %plot(ax1,MaxAWT((ChannelOld+pkChOld-1)==i,3),MaxAWT((ChannelOld+pkChOld-1)==i,1),'bo')
    %plot(ax1,MaxAWTc((Channel0+peakCh0(iy)-1)==i,3),MaxAWTc((Channel0+peakCh0(iy)-1)==i,1),'ro')
    ax1.CLim=[-4 2];
    xticks(ax1,[])
    yticks(ax1,[])
    ax1.YDir='normal';
end
ax1 = axes('OuterPosition',[xPos(8) yPos(5)+0.05 0.28 0.21]);
imagesc(ax1,log10(squeeze(sum(sum(m.NullHisto(2:end-1,:,:,:),4),1)/m.Nch))');
colormap(ax1,'gray')
hold(ax1,'on')
plot(ax1,m.pwt(m.fracBorder(:,1)>=0.6,2),m.pwt(m.fracBorder(:,1)>=0.6,1),'co')
plot(ax1,m.pwt(m.fracBorder(:,1)<0.6,2),m.pwt(m.fracBorder(:,1)<0.6,1),'mo')
ax1.CLim=[-4 2];
xticks(ax1,interp1(m.tWidth,1:m.Nwidth,[0.1 0.2 0.5 1]))
xticklabels(ax1,[0.1 0.2 0.5 1])
xlabel(ax1,'width/ms')
yticks(ax1,[0.5 32.5 64.5])
yticklabels(ax1,m.xPwr(1:32:65))
ylabel(ax1,'power/stdev')
ax1.YDir='normal';
hC=colorbar(ax1,'Ticks',[-2 0 2],...
    'TickLabels',{'0.01','1','100'});
hC.Label.String = 'Hz';
%tail
ax1 = axes('OuterPosition',[xPos(1)-0.005 yPos(22)-0.01 0.28 0.21]);
imagesc(ax1,log10(squeeze(sum(sum(m.NullHisto(:,2:end-1,:,:),4),2)/m.Nch))');
colormap(ax1,'copper')
hold(ax1,'on')
plot(ax1,m.pwt(m.fracBorder(:,1)>=0.6,3),m.pwt(m.fracBorder(:,1)>=0.6,1),'co')
plot(ax1,m.pwt(m.fracBorder(:,1)<0.6,3),m.pwt(m.fracBorder(:,1)<0.6,1),'mo')
ax1.CLim=[-4 2];
xticks(ax1,interp1(m.tTail,1:m.Ntail,[1/3 1/2 2/3 8/9]))
xticklabels(ax1,{'2:1' '1:1' '1:2' '1:8'})
xlabel(ax1,'symmetry')
%xticks(ax1,1:5:16)
%xticklabels(ax1,round(m.tTail(1:5:16)*100)/100)
%xlabel(ax1,'frac. of tail')
yticks(ax1,[0.5 32.5 64.5])
yticklabels(ax1,m.xPwr(1:32:65))
ylabel(ax1,'power/stdev')
ax1.YDir='normal';
hC=colorbar(ax1,'Ticks',[-2 0 2],...
    'TickLabels',{'0.01','1','100'});
hC.Label.String = 'Hz';
if ~isfolder([plotBase 'MaxLoc/'])
    mkdir([plotBase 'MaxLoc/'])
end
tic; pause(2); toc;
saveas(fig1,[Clust.plotFolderE filesep Clust.plotFileE])
close(fig1)

%save to file.
%Need spike times,units, histograms (absolute spike numbers), channels,
%weights, some kind of cluster amplitude.
save(Target,'m','-v7.3');%(m.Namp+2*nsah-1)*(m.Nwidth+2*nswh)*(m.Ntail+2*nsth)


nRows=12;
dy=0.98/nRows;
fig1=figure('Position',[2560 440 1200 nRows*80]);
xPos=0.011+(0:25)*0.1;
yPos=0.01+((nRows-1):-1:0)*dy;

for i=1:m.nUnits
    ax1 = axes('Position',[xPos(floor((i-1)/nRows)+1)...
        yPos(mod(i-1,nRows)+1)+(dy-0.024)/4 0.035 (dy-0.024)]);
    imagesc(ax1,log10(squeeze(sum(m.Hist{i,1}(2:end-1,2:end-1,2:end-1),1))*m.Sampling/m.LenRec)');
    ax1.CLim=[-4 0];
    ylabel(ax1,['u' num2str(i) ',ch' num2str(m.Channel(i))])
    xticks(ax1,[])
    yticks(ax1,[])
    ax1.YDir='normal';
%     title(ax1,num2str(m.unitBgRatio(i,1)/100,2))
    ax1 = axes('Position',[xPos(floor((i-1)/nRows)+1)+0.035...
        yPos(mod(i-1,nRows)+1)+(dy-0.024)/4 0.025 (dy-0.024)]);
    %rSHtemp=reshape(m.rSH{i,1},m.Ntail,m.Nwidth,m.Namp);
    imagesc(ax1,log10(squeeze(sum(m.Hist{i,1}(2:end-1,2:end-1,2:end-1),2))*m.Sampling/m.LenRec)');
    colormap(ax1,'hot')
    ax1.CLim=[-4 0];
    xticks(ax1,[])
    yticks(ax1,[])
    ax1.YDir='normal';
%     title(ax1,num2str(m.unitBgAvg(i,1)*100,2))
%         %ISI dist.
%         ax1 = axes('Position',[xPos(floor((i-1)/nRows)+1)+0.06...
%             yPos(mod(i-1,nRows)+1) 0.029 (dy-0.024)]);
%         hold(ax1,'on')
%         for k=find(ClustId==iy(i))
%             plot(1:150,rISI(k,:),'k-')
%         end
%         plot(1:150,yISI(iy(i),:),'r-','LineWidth',2)
%         xticks(ax1,[])
%         yticks(ax1,[])
end
tic; pause(2); toc;
saveas(fig1,[Clust.plotFolderC filesep Clust.plotFileC])
close(fig1)

%save data
end