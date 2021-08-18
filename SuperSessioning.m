classdef SuperSessioning
    %class for file handling, shared variables and collection of
    %applied parameters.
    properties
        Sessions
        TimeStamps
        nRec
        %folders and files
        BaseFolder
        RawFiles
        RawFolder
        FiltFiles
        FiltFolder
        locMaxFiles
        locMaxFolder
        singleSessionFiles
        singleSessionFolder
        mergedList
        mergedFolder
        mergedFile
        hashFolder
        plotFolder
        %defaults
        defaultRecParameter
        defaultFilter
        defaultbTM
        defaultClust
        nC%naming conventions
        subject
        %session dependent
        recParameter
        Noise
        Filter
        bTM
        Clust
    end
    methods
        % create a file structure for analysis
        function obj = SuperSessioning(BaseFolder,RawFolder,subject)
            obj.BaseFolder=BaseFolder;
            %assert that BaseFolder exists
            assert(isdir(BaseFolder),'Folder does not exist')
            if nargin >1
                assert(isdir(RawFolder),'Raw data folder does not exist')
                obj.RawFolder=RawFolder;
            else
                obj.RawFolder=[BaseFolder filesep 'Raw'];
                if ~isfolder(obj.RawFolder)
                    mkdir(obj.RawFolder)
                end
            end
            %create folder structure 
            %(may need to change these later, as files are moved around)
            obj.FiltFolder='filtered';
            if ~isfolder(obj.XFolder(obj.FiltFolder))
                mkdir(obj.XFolder(obj.FiltFolder));
            end
            obj.locMaxFolder='locMax';
            if ~isfolder(obj.XFolder(obj.locMaxFolder))
                mkdir(obj.XFolder(obj.locMaxFolder));
            end
            obj.singleSessionFolder='singleSession';
            if ~isfolder(obj.XFolder(obj.singleSessionFolder))
                mkdir(obj.XFolder(obj.singleSessionFolder));
            end
            obj.mergedFolder='merged';
            if ~isfolder(obj.XFolder(obj.mergedFolder))
                mkdir(obj.XFolder(obj.mergedFolder));
            end
            obj.plotFolder='Plots';
            if ~isfolder(obj.XFolder(obj.plotFolder))
                mkdir(obj.XFolder(obj.plotFolder));
            end
            % initialize variables
            obj.TimeStamps=[];%allow for anything here, just make sure MATLAB can sort it.
            obj.mergedList=[];
            obj.RawFiles={};
            obj.FiltFiles={};
            obj.locMaxFiles={};
            obj.singleSessionFiles={};
            obj.nRec=0;
            %settings for file names
            obj.nC={};
            obj.nC.nExtRaw=7;
            obj.nC.Ext60Hz='_60Hz.kwd';
            obj.nC.nExt60Hz=4;
            obj.nC.ExtCAR='_CAR.kwd';
            obj.nC.nExtCAR=4;
            obj.nC.ExtLocMax='_LM.mat';
            obj.nC.nExtLocMax=4;
            obj.nC.plotExt='.png';
            %subject identifier/name
            obj.subject=subject;
            
            %default settings for filtering
            obj.defaultFilter=struct();
            %Origin: load from Raw data ('raw', default) or temporary directory ('temp')
            %Target: save in a temporary directory ('temp', default) or as filtered data ('filtered')
            %recording Parameter
            obj.defaultRecParameter.readfromKWIK=true;
            obj.defaultRecParameter.RawFormat='rhd';
            
            %can globally change the order of electrodes here if desired
            %(better use for individual recordings, i.e. if connectors accidentally switched during a recording)
            obj.defaultRecParameter.ChMap=1:64;
            %mask for excluding channels (after mapping)
            obj.defaultRecParameter.ChMask=ones(1,64,'logical');
            %sampling rate, amplifier gain (muV per ADCcount)
            obj.defaultRecParameter.sRate=30000;
            obj.defaultRecParameter.bitVolt=0.1950;
            %restricting temporal interval/recording electrodes 
            obj.defaultRecParameter.tStart=0;%change to exclude a segment at the beginning of the recording
            obj.defaultRecParameter.tEnd=-1;%change to exclude a segment at the end of the recording
            obj.defaultRecParameter.iGroup=0;
            %parameters for 60Hz filter
            obj.defaultFilter.pwl.origin='raw';
            obj.defaultFilter.pwl.target='filt';
            obj.defaultFilter.pwl.plotResults=true;
            obj.defaultFilter.pwl.plotFolder='Plot_60Hz';
            obj.defaultFilter.pwl.plotName='Powerline';
            obj.defaultFilter.pwl.returnPlotData=false;
            %parameters for common average referencing
            obj.defaultFilter.car.origin='filt';
            obj.defaultFilter.car.target='filt';
            obj.defaultFilter.car.plotResults=true;
            obj.defaultFilter.car.plotFolder='Plot_CAR';
            obj.defaultFilter.car.plotName='CommonAvg';
            %parameter for combined filter (as alternative for now)
            obj.defaultFilter.origin='raw';
            obj.defaultFilter.target='filt';
            obj.defaultFilter.plotResults=true;
            obj.defaultFilter.plotFolder='Plot_Filter';
            obj.defaultFilter.plotName='Filter';
            obj.defaultFilter.returnPlotData=true;
            %default settings for local minima
            obj.defaultbTM.plotResults=true;
            obj.defaultbTM.plotFolder='Plot_LocMax';
            obj.defaultbTM.plotName='LocMax';
            obj.defaultbTM.origin='filt';
            obj.defaultbTM.ChMap=obj.defaultRecParameter.ChMap(obj.defaultRecParameter.ChMask);
            %default settings for clustering
            obj.defaultClust.plotResults=true;
            obj.defaultClust.plotFolderE='Plot_Electrodes';
            obj.defaultClust.plotNameE='Electrodes';
            obj.defaultClust.plotFolderC='Plot_Cluster';
            obj.defaultClust.plotNameC='Clusters';
            %default settings for merging
            %obj.defaultMerge.plotResults=true;
            %obj.defaultMerge.plotFolder=[plotFolder filesep 'Plot_Merged'];
            %obj.defaultMerge.plotFolder=[plotFolder filesep 'Plot_Merged_cluster'];
        end
        function [obj,Ind] = addRaw(obj,FileName,TimeStamp)
            %add the path, check whether exists
            assert(isfile([obj.RawFolder filesep FileName]),['File not found, should be in ' obj.RawFolder])
            %test whether in chronological order
            if all(obj.TimeStamps<TimeStamp)
                obj.nRec=obj.nRec+1;
                obj.TimeStamps(obj.nRec,1)=TimeStamp;
                obj.RawFiles{obj.nRec,1}=FileName;
                obj.FiltFiles{obj.nRec,1}='';
                obj.locMaxFiles{obj.nRec,1}='';
                obj.Filter{obj.nRec,1}=obj.defaultFilter;
                obj.bTM{obj.nRec,1}=obj.defaultbTM;
                obj.Clust{obj.nRec,1}=obj.defaultClust;
                obj.recParameter{obj.nRec,1}=obj.defaultRecParameter;
                Ind=obj.nRec;
            elseif any(obj.TimeStamps==TimeStamp)
                %assume adding a Raw file to the structure later
                Ind=find(obj.TimeStamps==TimeStamp);
                assert(length(Ind)==1);
                assert(isempty(obj.RawFiles{Ind}),['existing file with same timestamp: ' obj.RawFolder filesep obj.RawFiles{Ind}])
                obj.RawFiles{Ind}=FileName;
                if ~isstruct(obj.Filter{Ind,1})
                    obj.Filter{Ind,1}=obj.defaultFilter;
                end
                if ~isstruct(obj.recParameter{Ind,1})
                    obj.recParameter{Ind,1}=obj.defaultRecParameter;
                end
                if ~isstruct(obj.bTM{Ind,1})
                    obj.bTM{Ind,1}=obj.defaultbTM;
                end
                if ~isstruct(obj.Clust{Ind,1})
                    obj.Clust{Ind,1}=obj.defaultClust;
                end
            else
                 %add an intermediate Raw file (allow adding files later),
                 %expand structures
                 Ind=sum(obj.TimeStamps<TimeStamp);
                 obj.TimeStamps=[obj.TimeStamps(1:Ind,1); TimeStamp; obj.TimeStamps(Ind+1:obj.nRec,1)];
                 obj.RawFiles=[{obj.RawFiles(1:Ind,1)}; {[obj.RawFolder filesep FileName]}; {obj.RawFiles(Ind+1:obj.nRec,1)}];
                 obj.FiltFiles=[{obj.FiltFiles(1:Ind,1)}; {}; {obj.FiltFiles(Ind+1:obj.nRec,1)}];
                 obj.locMaxFiles=[{obj.locMaxFiles(1:Ind,1)}; {}; {obj.locMaxFiles(Ind+1:obj.nRec,1)}];
                 obj.mergedList=[obj.mergedList(1:Ind,1); 0; obj.mergedList(Ind+1:obj.nRec,1)];
                 obj.recParameter=[{obj.recParameter(1:Ind,1)}; obj.defaultRecParameter; {obj.recParameter(Ind+1:obj.nRec,1)}];
                 obj.Filter=[{obj.Filter(1:Ind,1)}; obj.defaultFilter; {obj.Filter(Ind+1:obj.nRec,1)}];
                 obj.bTM=[{obj.bTM(1:Ind,1)}; obj.defaultbTM; {obj.bTM(Ind+1:obj.nRec,1)}];
                 obj.Clust=[{obj.Clust(1:Ind,1)}; obj.defaultClust; {obj.Clust(Ind+1:obj.nRec,1)}];
                 obj.nRec=obj.nRec+1;
                 Ind=Ind+1;
            end
        end
        %filter 60 Hz component
        function obj=filter60Hz(obj,Ind)
            switch lower(obj.Filter{Ind}.pwl.origin)
                case 'raw'
                    Origin=[obj.RawFolder filesep obj.RawFiles{Ind}];
                case 'filt'
                    Origin=[obj.XFolder(obj.FiltFolder) filesep obj.FiltFiles{Ind}];
                otherwise
                    Origin=[obj.RawFolder filesep obj.RawFiles{Ind}];
            end
            assert(isfile(Origin),['File not found, should be in ' Origin])
            obj.FiltFiles{Ind}=[obj.RawFiles{Ind}(1:end-obj.nC.nExtRaw-1) obj.nC.Ext60Hz];
            Target=[obj.XFolder(obj.FiltFolder) filesep obj.FiltFiles{Ind}];
            if obj.Filter{Ind}.pwl.plotResults
                if ~isfolder(obj.PFolder(obj.Filter{Ind}.pwl.plotFolder))
                    mkdir(obj.PFolder(obj.Filter{Ind}.pwl.plotFolder));
                end
            end
            %unique names for plots
            obj.Filter{Ind}.pwl.plotName=[obj.RawFiles{Ind}(1:end-obj.nC.nExtRaw-1) '_' ...
                obj.defaultFilter.pwl.plotName obj.nC.plotExt];
            %obj.Filter{Ind}.pwl.BaseFolder=obj.BaseFolder;
            obj.Filter{Ind}.pwl.PlotBase=obj.XFolder(obj.plotFolder);
            [obj.Filter{Ind},obj.recParameter{Ind}]=filter.Filter60Hz(Origin, Target, obj.Filter{Ind},obj.recParameter{Ind});
            %may not want a temporary folder! Just save as filtered output
        end
        %common average referencing
        function obj=filterCAR(obj,Ind)
            switch lower(obj.Filter{Ind}.car.origin)
                case 'raw'
                    Origin=[obj.RawFolder filesep obj.RawFiles{Ind}];
                case 'filt'
                    Origin=[obj.XFolder(obj.FiltFolder) filesep obj.FiltFiles{Ind}];
                otherwise
                    Origin=[obj.RawFolder filesep obj.FiltFiles{Ind}];
            end
            assert(isfile(Origin),['File not found, should be in ' Origin])
            obj.FiltFiles{Ind}=[obj.RawFiles{Ind}(1:end-obj.nC.nExtRaw-1) obj.nC.ExtCAR];
            Target=[obj.XFolder(obj.FiltFolder) filesep obj.FiltFiles{Ind}];
            if obj.Filter{Ind}.car.plotResults
                if ~isfolder(obj.PFolder(obj.Filter{Ind}.car.plotFolder))
                    mkdir(obj.PFolder(obj.Filter{Ind}.car.plotFolder));
                end
            end
            obj.Filter{Ind}.car.plotName=[obj.RawFiles{Ind}(1:end-obj.nC.nExtRaw-1) '_' ...
                obj.defaultFilter.car.plotName obj.nC.plotExt];
            obj.Filter{Ind}.car.PlotBase=obj.XFolder(obj.plotFolder);
            [obj.Filter{Ind},obj.recParameter{Ind}]=filter.FilterCAR(Origin,Target, obj.Filter{Ind},obj.recParameter{Ind});
        end
        function obj=filterGeneral(obj,Ind)
            switch lower(obj.Filter{Ind}.origin)
                case 'raw'
                    Origin=[obj.RawFolder filesep obj.RawFiles{Ind}];
                case 'filt'
                    Origin=[obj.XFolder(obj.FiltFolder) filesep obj.FiltFiles{Ind}];
                otherwise
                    Origin=[obj.RawFolder filesep obj.RawFiles{Ind}];
            end
            assert(isfile(Origin),['File not found, should be in ' Origin])
            obj.FiltFiles{Ind}=[obj.RawFiles{Ind}(1:end-obj.nC.nExtRaw-1) obj.nC.ExtCAR];
            Target=[obj.XFolder(obj.FiltFolder) filesep obj.FiltFiles{Ind}];
            if obj.Filter{Ind}.pwl.plotResults
                if ~isfolder(obj.PFolder(obj.Filter{Ind}.plotFolder))
                    mkdir(obj.PFolder(obj.Filter{Ind}.plotFolder));
                end
            end
            %unique names for plots
            obj.Filter{Ind}.plotName=[obj.RawFiles{Ind}(1:end-obj.nC.nExtRaw-1) '_' ...
                obj.defaultFilter.plotName obj.nC.plotExt];
            %obj.Filter{Ind}.pwl.BaseFolder=obj.BaseFolder;
            obj.Filter{Ind}.PlotBase=obj.XFolder(obj.plotFolder);
            [obj.Filter{Ind},obj.recParameter{Ind}]=filter.filterGeneral(Origin, Target, obj.Filter{Ind},obj.recParameter{Ind});
            %may not want a temporary folder! Just save as filtered output
        end
        %only extract parameter if filtered file available
        function obj=filterParams(obj,Ind)
            Origin=[obj.RawFolder filesep obj.RawFiles{Ind}];
            assert(isfile(Origin),['File not found, should be in ' Origin])
            obj.FiltFiles{Ind}=[obj.RawFiles{Ind}(1:end-obj.nC.nExtRaw-1) obj.nC.ExtCAR];
            obj.recParameter{Ind}=filter.ExtractRecParams(Origin, obj.recParameter{Ind});
        end
        %need another step here to estimate high-pass filtered variances
        function obj=estimateStd(obj,Ind)
            %see whether filtering was done
            Origin=[obj.XFolder(obj.FiltFolder) filesep obj.FiltFiles{Ind}];
            assert(isfile(Origin),['File not found, should be in ' Origin])
            %read data, high-pass filter and estimate variance
            %make histograms for (overlapping) chunks of data (~3s), determine
            %mode of histogram -- problem: outliers?
            %do electrode by electrode
            [b, a] = butter(4, 2*300/obj.recParameter{Ind}.sRate, 'high');%300 Hz highpass filter
            obj.Noise{Ind}.sigmaADC=zeros(obj.recParameter{Ind}.Nch,1);
            obj.Noise{Ind}.sigma=zeros(obj.recParameter{Ind}.Nch,1);
            obj.Noise{Ind}.ACF=zeros(obj.recParameter{Ind}.Nch,31);
            obj.Noise{Ind}.lowVarDT=500;
            N=floor(obj.recParameter{Ind}.LenRec/obj.Noise{Ind}.lowVarDT);
            obj.Noise{Ind}.tVar=zeros(1,N);
            for j=1:sum(obj.recParameter{Ind}.ChMask)
                x=double(h5read([obj.XFolder(obj.FiltFolder) filesep obj.FiltFiles{Ind}],...
                    '/data',[j,1],[1,obj.recParameter{Ind}.LenRec])');
                y=filtfilt(b,a,x);
                ystd=std(y);
                obj.Noise{Ind}.tVar=obj.Noise{Ind}.tVar+...
                    std(reshape(y(1:obj.Noise{Ind}.lowVarDT*N),obj.Noise{Ind}.lowVarDT,N),1);
                %mean(abs(reshape(y(1:500*N),500,N)),1)];
                %also want ACF
                obj.Noise{Ind}.ACF(j,:) = autocorr(y,'NumLags',30);
                %clip everything larger than 5 standard deviations (these
                %may be spikes and therefore activity-dependent)
                x=y((-5*ystd<y) & (y<5*ystd));
                %additionally take 30th percentile of 1000 frames cutouts!
                %(TODO)
                obj.Noise{Ind}.sigmaADC(j,1)=std(x);
                obj.Noise{Ind}.sigma(j,1)=std(x)*obj.recParameter{Ind}.bitVolt;
            end
            obj.Noise{Ind}.tVar=[obj.Noise{Ind}.tVar obj.Noise{Ind}.tVar(:,end)];
            a=conv([obj.Noise{Ind}.tVar(1,7:-1:1) obj.Noise{Ind}.tVar ...
                obj.Noise{Ind}.tVar(1,end-1:-1:end-7)],hanning(15)'/8,'valid');
            b=prctile(a,(10:1:90));
            c=sort(b(1:end-1)+diff(b).*(165:-1:86));
            d=sort(b(1:end-1)+diff(b).*(115:-1:36));
            obj.Noise{Ind}.lowVar=(a<c(7)) | ...
                (conv([obj.Noise{Ind}.tVar(1,2:-1:1) obj.Noise{Ind}.tVar ...
                obj.Noise{Ind}.tVar(1,end-1:-1:end-2)],hanning(5)'/3,'valid') < d(7));
        end
        %find local minima
        function obj=blindTemplateMatching(obj,Ind)
            switch lower(obj.bTM{Ind}.origin)
                case 'raw'
                    Origin=[obj.RawFolder filesep obj.RawFiles{Ind}];
                case 'filt'
                    Origin=[obj.XFolder(obj.FiltFolder) filesep obj.FiltFiles{Ind}];
                otherwise
                    Origin=[obj.RawFolder filesep obj.FiltFiles{Ind}];
            end
            obj.locMaxFiles{Ind}=[obj.RawFiles{Ind}(1:end-obj.nC.nExtRaw-1) obj.nC.ExtLocMax];
            Target=[obj.XFolder(obj.locMaxFolder) filesep obj.locMaxFiles{Ind}];
            if obj.bTM{Ind}.plotResults
                if ~isfolder(obj.PFolder(obj.bTM{Ind}.plotFolder))
                    mkdir(obj.PFolder(obj.bTM{Ind}.plotFolder));
                end
            end
            obj.bTM{Ind}.plotFile=[obj.RawFiles{Ind}(1:end-obj.nC.nExtRaw-1) '_' obj.defaultbTM.plotName obj.nC.plotExt];
            obj.bTM{Ind}.PlotBase=obj.XFolder(obj.plotFolder);
            obj.bTM{Ind}.Noise=obj.Noise{Ind};
            blindTM.blindTemplateMatchingGPU(Origin,Target,obj.bTM{Ind},obj.recParameter{Ind})% probably want to add a few inputs here, and some output?
        end
        %cluster local minima from a session
        function obj=sortSession(obj,Ind)
            Origin=[obj.XFolder(obj.locMaxFolder) filesep obj.locMaxFiles{Ind}];
            obj.singleSessionFiles{Ind}=[obj.RawFiles{Ind}(1:end-obj.nC.nExtRaw-1) obj.nC.ExtLocMax];
            Target=[obj.XFolder(obj.singleSessionFolder) filesep obj.singleSessionFiles{Ind}];
            if obj.Clust{Ind}.plotResults
                if ~isfolder(obj.PFolder(obj.Clust{Ind}.plotFolderE))
                    mkdir(obj.PFolder(obj.Clust{Ind}.plotFolderE));
                end
                if ~isfolder(obj.PFolder(obj.Clust{Ind}.plotFolderC))
                    mkdir(obj.PFolder(obj.Clust{Ind}.plotFolderC));
                end
            end
            obj.Clust{Ind}.RecId=obj.RawFiles{Ind}(1:end-obj.nC.nExtRaw-1);
            obj.Clust{Ind}.RawFile=[obj.RawFolder filesep obj.RawFiles{Ind}];
            obj.Clust{Ind}.plotNameE=[obj.RawFiles{Ind}(1:end-obj.nC.nExtRaw-1) '_' obj.defaultClust.plotNameE obj.nC.plotExt];
            obj.Clust{Ind}.plotNameC=[obj.RawFiles{Ind}(1:end-obj.nC.nExtRaw-1) '_' obj.defaultClust.plotNameC obj.nC.plotExt];
            obj.Clust{Ind}.PlotBase=obj.XFolder(obj.plotFolder);
            obj.Clust{Ind}=merge.clusterSession(Origin, Target, obj.Clust{Ind});
        end
        %the following functions can be run using a subset of
        %recordings--------------------------------------------------------
        %do I need to have a clear previous results function? -- looks like
        %merge list of recordings
        function obj=mergeAll(obj,FileName,IndList)
            %need to check if there is already a sorted file
            assert(~isfile([obj.XFolder(obj.mergedFolder) filesep FileName]))
            merge.mergeAllLocMax(obj,[obj.XFolder(obj.mergedFolder) filesep FileName],IndList);
            obj.mergedList=IndList;
            obj.mergedFile=FileName;
            %obj.TlastMerged=obj.TimeStamp(obj.mergedList(end));%better use a list of indices than Mask!!
        end
        %add next recording to merge
        function obj=mergeNext(obj,Ind)
            %make sure this session is temporally after sessions that were already merged
            %assert(obj.TimeStamps(Ind)>obj.TimeStamps(obj.mergedList(end)),...
            %    'Session in between already merged ones.')
            copyfile([obj.XFolder(obj.mergedFolder) filesep obj.mergedFile],[obj.XFolder(obj.mergedFolder) filesep obj.mergedFile(1:end-4) '.bak']);
            merge.mergeIncLocMax(obj,Ind)
            obj.mergedList=[obj.mergedList Ind];
        end
        %merge hash (no incremental version here)
        function obj=mergeHash(obj,IndList)
            obj.hashFolder='hash';
            if ~isfolder(obj.XFolder(obj.hashFolder))
                mkdir(obj.XFolder(obj.hashFolder));
            end
            merge.mergeAllHash(obj,IndList);
        end
        %some file management stuff, helper functions-------------------------------
        %main folders
        function FolderName=XFolder(obj,subFolder)
            FolderName=[obj.BaseFolder filesep subFolder];
        end
        %subfolders for plotting
        function FolderName=PFolder(obj,subFolder)
            FolderName=[obj.BaseFolder filesep obj.plotFolder filesep subFolder];
        end
    end
end