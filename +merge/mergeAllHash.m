function mergeAllHash(obj,IndList)
%only try to merge clusters if    plotBase,,plotMerge,plotMax
nSessions=length(IndList);
m0=load([obj.singleSessionFolder filesep obj.singleSessionFiles{IndList(1)}]);
m=m0.m;
for ch=1:m.Nch
    z=struct();
    z.RecId={};
    z.RawFile={};
    for ii=1:1
        m0=load([obj.singleSessionFolder filesep obj.singleSessionFiles{IndList(ii)}]);
        m=m0.m;
        z.Nch=m.Nch;
        z.Channel=ch;
        z.hashSpikes=cell(nSessions,1);
        z.hashTimes=cell(nSessions,1);
        z.RecId{1}=m.RecId;
        z.RawFile{1}=m.RawFile;
        z.LenRec(1)=m.LenRec;
        %general stuff
        z.Sampling=m.Sampling;
        z.ignore.Bsc=m.Bsc;
        z.ignore.nsTWA=m.nsTWA;
        z.ignore.nsTWAh=m.nsTWAh;
        z.Jthreshold=m.Jthreshold;
        z.ignore.Ckernel=m.Ckernel;
        z.ignore.Cks=m.Cks;
        z.HzNorm=m.HzNorm;
        z.xPwr=m.xPwr;
        z.ignore.rPwr=m.rPwr;
        z.ignore.aPwr=m.aPwr;
        z.tWidth=m.tWidth;
        z.tTail=m.tTail;
        z.Namp=m.Namp;
        z.Nwidth=m.Nwidth;
        z.Ntail=m.Ntail;
        z.ignore.Nsmth=m.Nsmth;
        z.hashSpikes{ii}=m.hashSpikes{ch,1};
        z.hashTimes{ii}=m.hashTimes{ch,1};
    end
    for ii=2:nSessions
        %ii0=ii-nMerge+1;
        disp(ii)
        %append another session to the data
        load([obj.singleSessionFolder filesep obj.singleSessionFiles{IndList(ii)}])
        %iSess=mod(ii-1,nMerge)+1;
        m=m0.m;
        z.RecId{ii}=m.RecId;
        z.RawFile{ii}=m.RawFile;
        z.hashSpikes{ii}=m.hashSpikes{ch,1};
        z.hashTimes{ii}=m.hashTimes{ch,1};
        z.LenRec(ii)=m.LenRec;
    end
    save([obj.hashFolder filesep obj.subject '_ch' num2str(ch) '_hash.mat'],'z','-v7.3')
end
end