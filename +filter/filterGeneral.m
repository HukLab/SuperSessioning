function [Filter,recParameter]=filterGeneral(RawFile, OutFile, Filter, recParameter)
    %Nch=length(recParameter.ChMask);
    Nchx=sum(recParameter.ChMask,2);
    ChMask=recParameter.ChMask;
    ChMap=recParameter.ChMap;
    
    if recParameter.readfromKWIK
        %iGroup=1;
        bitVolt=h5read(RawFile,['/recordings/' num2str(recParameter.iGroup) '/application_data/channel_bit_volts']);
        sRate=h5read(RawFile,['/recordings/' num2str(recParameter.iGroup) '/application_data/channel_sample_rates']);
        
        recParameter.bitVolt=median(double(bitVolt(ChMap(ChMask),1)));
        recParameter.sRate=median(double(sRate(ChMap(ChMask),1)));
        recParameter.bitVoltCh=double(bitVolt(ChMap(ChMask),1));
        recParameter.sRateCh=double(sRate(ChMap(ChMask),1));
        if ~isfield(recParameter,'HdfRawDataPath')
            recParameter.HdfRawDataPath=['/recordings/' num2str(recParameter.iGroup) '/data'];
        end
    else
        recParameter.HdfRawDataPath='/data';
    end
    a=h5info(RawFile,recParameter.HdfRawDataPath);
    DataSize=a.Dataspace.Size;
    if recParameter.tEnd==-1
        recParameter.nEnd=double(floor(DataSize(1,2)));
    else
        recParameter.nEnd=double(floor(recParameter.tEnd*60*sRate(1,1)));
    end
    if recParameter.tStart==0
        recParameter.nStart=1;
    else
        recParameter.nStart=double(floor(recParameter.tStart*60*sRate(1,1))+1);
    end
    
    h5create(OutFile,'/data',[Nchx Inf],'ChunkSize',[1 2048],'Datatype','int16','FillValue',int16(0))
    
    %Parameters
    recParameter.LenRec=recParameter.nEnd;%double(floor(DataSize(1,2)));
    recParameter.Nch=Nchx;
    
    F = filter.filterBase(recParameter,3e6);
    F=F.filterAll(RawFile, OutFile);
    F.templateSine=[];
    Filter.general=F;
    F.plotResults(Filter);
end