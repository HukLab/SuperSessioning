%% Extract recording parameter
%only need if no 60 Hz removal
function recParameter=ExtractRecParams(RawFile, recParameter)
ChMask=recParameter.ChMask;
ChMap=recParameter.ChMap;

if recParameter.readfromKWIK
    iGroup=0;
    bitVolt=h5read(RawFile,['/recordings/' num2str(iGroup) '/application_data/channel_bit_volts']);
    sRate=h5read(RawFile,['/recordings/' num2str(iGroup) '/application_data/channel_sample_rates']);
    
    recParameter.bitVolt=median(double(bitVolt(ChMap(ChMask),1)));
    recParameter.sRate=median(double(sRate(ChMap(ChMask),1)));
    recParameter.bitVoltCh=double(bitVolt(ChMap(ChMask),1));
    recParameter.sRateCh=double(sRate(ChMap(ChMask),1));
    HdfDataPath=['/recordings/' num2str(iGroup) '/data'];
else
    HdfDataPath='/data';
end
a=h5info(RawFile,HdfDataPath);
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
%Parameters
recParameter.LenRec=double(floor(DataSize(1,2)));
recParameter.Nch=sum(recParameter.ChMask,2);
end
