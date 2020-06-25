classdef SuperSessioning
    properties
        Sessions
        RawFiles
        RawFolder
        FiltFiles
        FiltFolder
        LocMaxFiles
        locMaxFolder
        MergedList
        mergedFolder
        TimeStamps
        defaultFilter
        nC
        nRec
        subject
        
    end
    methods
        % create a file structure for analysis
        function obj = SuperSessioning(BaseFolder,RawFolder)
            %assert that BaseFolder exists
            assert(isDir(BaseFolder),'Folder does not exist')
            if nargin >1
                assert(isDir(RawFolder),'Raw data folder does not exist')
                obj.RawFolder=RawFolder;
            else
                obj.RawFolder=[BaseFolder filesep 'Raw'];
                if ~isfolder(obj.RawFolder)
                    mkdir(obj.RawFolder)
                end
            end
            %create folder structure 
            %(may need to change these later, as files are moved around)
            obj.FiltFolder=[BaseFolder filesep 'filtered'];
            if ~isfolder(obj.FiltFolder)
                mkdir(obj.FiltFolder)
            end
            obj.locMaxFolder=[BaseFolder filesep 'locMax'];
            if ~isfolder(obj.locMaxFolder)
                mkdir(obj.locMaxFolder);
            end
            obj.singleSessionFolder=[BaseFolder filesep 'singleSession'];
            if ~isfolder(obj.singleSessionFolder)
                mkdir(obj.singleSessionFolder);
            end
            obj.mergedFolder=[BaseFolder filesep 'merged'];
            if ~isfolder(obj.mergedFolder)
                mkdir(obj.mergedFolder);
            end
            obj.plotFolder=[BaseFolder filesep 'Plots'];
            if ~isfolder(obj.plotFolder)
                mkdir(obj.plotFolder);
            end
            % initialize variables
            obj.TimeStamps=[];%allow for anything here, just make sure MATLAB can sort it.
            obj.RawFiles={};
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
            obj.subject='';
            
            %default settings for filtering
            obj.defaultFilter=struct();
            %Origin: load from Raw data ('raw', default) or temporary directory ('temp')
            %Target: save in a temporary directory ('temp', default) or as filtered data ('filtered')
            %parameters for 60Hz filter
            obj.defaultFilter.pwl.origin='raw';
            obj.defaultFilter.pwl.target='filt';
            obj.defaultFilter.pwl.plotResults=true;
            obj.defaultFilter.pwl.plotFolder=[plotFolder filesep 'Plot_60Hz'];
            obj.defaultFilter.pwl.plotName='Powerline';
            obj.defaultFilter.pwl.returnPlotData=false;
            %parameters for common average referencing
            obj.defaultFilter.car.origin='filt';
            obj.defaultFilter.car.target='filt';
            obj.defaultFilter.car.plotResults=true;
            obj.defaultFilter.car.plotFolder=[plotFolder filesep 'Plot_CAR'];
            obj.defaultFilter.car.plotName='CommonAvg';
            %restricting temporal interval/recording electrodes 
            obj.defaultFilter.tStart=0;%change to exclude a segment at the beginning of the recording
            obj.defaultFilter.tEnd=-1;%change to exclude a segment at the end of the recording
            %can globally change the order of electrodes here if desired
            %(better use for individual recordings, i.e. if connectors accidentally switched during a recording)
            obj.defaultFilter.ChMap=1:64;
            %mask for excluding channels (after mapping)
            obj.defaultFilter.ChMask=ones(1,64,'logical');
            
            %default settings for local minima
            obj.defaultbTM.plotResults=true;
            obj.defaultbTM.plotFolder=[plotFolder filesep 'Plot_LocMax'];
            obj.defaultbTM.plotName='LocMax';
            %default settings for clustering
            obj.defaultClust.plotResults=true;
            obj.defaultClust.plotFolderE=[plotFolder filesep 'Plot_Electrodes'];
            obj.defaultClust.plotNameE='Electrodes';
            obj.defaultClust.plotFolderC=[plotFolder filesep 'Plot_Cluster'];
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
                obj.mergedList(obj.nRec,1)=0;
                obj.Filter{obj.nRec,1}=obj.defaultFilter;
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
            else
                 %add an intermediate Raw file (allow adding files later),
                 %expand structures
                 Ind=sum(obj.TimeStamps<TimeStamp);
                 obj.TimeStamps=[obj.TimeStamps(1:Ind,1); TimeStamp; obj.TimeStamps(Ind+1:obj.nRec,1)];
                 obj.RawFiles=[{obj.RawFiles(1:Ind,1)}; {[obj.RawFolder filesep FileName]}; {obj.RawFiles(Ind+1:obj.nRec,1)}];
                 obj.FiltFiles=[{obj.FiltFiles(1:Ind,1)}; {}; {obj.FiltFiles(Ind+1:obj.nRec,1)}];
                 obj.locMaxFiles=[{obj.locMaxFiles(1:Ind,1)}; {}; {obj.locMaxFiles(Ind+1:obj.nRec,1)}];
                 obj.mergedList=[obj.mergedList(1:Ind,1); 0; obj.mergedList(Ind+1:obj.nRec,1)];
                 obj.Filter=[{obj.Filter(1:Ind,1)}; obj.defaultFilter; {obj.Filter(Ind+1:obj.nRec,1)}];
                 obj.nRec=obj.nRec+1;
                 Ind=Ind+1;
            end
        end
        %filter 60 Hz component
        function obj=filter60Hz(obj,Ind)
            switch lower(obj.Filter{Ind}.pwl.origin)
                case 'raw'
                    Origin=[obj.RawFolder filesep obj.RawFile{Ind}];
                case 'filt'
                    Origin=[obj.FiltFolder filesep obj.FiltFile{Ind}];
                otherwise
                    Origin=[obj.RawFolder filesep obj.RawFile{Ind}];
            end
            assert(isfile(Origin),['File not found, should be in ' Origin])
            obj.FiltFile=[obj.RawFile{Ind}(1:end-obj.nC.nExtRaw-1) obj.nC.Ext60Hz];
            Target=[obj.FiltFolder filesep obj.FiltFile];
            if obj.Filter.pwl.plotResults
                if ~isfolder(obj.Filter{Ind}.pwl.plotFolder)
                    mkdir(obj.Filter{Ind}.pwl.plotFolder);
                end
            end
            %unique names for plots
            obj.Filter{Ind}.pwl.plotName=[obj.RawFile{Ind}(1:end-obj.nC.nExtRaw-1) '_' obj.Filter.pwl.plotName obj.nC.plotExt];
            obj.Filter{Ind}=Filter60Hz(Origin, Target, obj.Filter{Ind});
            %may not want a temporary folder! Just save as filtered output
        end
        %common average referencing
        function obj=filterCAR(obj,Ind)
            switch lower(obj.Filter{Ind}.car.origin)
                case 'raw'
                    Origin=[obj.RawFolder filesep obj.RawFile{Ind}];
                case 'filt'
                    Origin=[obj.FiltFolder filesep obj.FiltFile{Ind}];
                otherwise
                    Origin=[obj.RawFolder filesep obj.Filtile{Ind}];
            end
            assert(isfile(Origin),['File not found, should be in ' Origin])
            obj.FiltFile=[obj.RawFile{Ind}(1:end-obj.nC.nExtRaw-1) obj.nC.ExtCAR];
            Target=[obj.FiltFolder filesep obj.FiltFile];
            if obj.Filter.car.plotResults
                if ~isfolder(obj.Filter{Ind}.car.plotFolder)
                    mkdir(obj.Filter{Ind}.car.plotFolder);
                end
            end
            obj.Filter{Ind}.car.plotName=[obj.RawFile{Ind}(1:end-obj.nC.nExtRaw-1) '_' obj.Filter.car.plotName obj.nC.plotExt];
            obj.Filter{Ind}=FilterCAR(Origin,Target, obj.Filter{Ind});
        end
        %need another step here to estimate high-pass filtered variances
        function obj=estimateStd(obj,Ind)
            %see whether filtering was done
            
            %read data, high-pass filter and estimate variance
            %make histograms for (overlapping) chunks of data, determine
            %mode of histogram -- problem: outliers?
        end
        %find local minima
        function obj=blindTemplateMatching(obj,Ind)
            switch lower(obj.bTM{Ind}.origin)
                case 'raw'
                    Origin=[obj.RawFolder filesep obj.RawFile{Ind}];
                case 'filt'
                    Origin=[obj.FiltFolder filesep obj.FiltFile{Ind}];
                otherwise
                    Origin=[obj.RawFolder filesep obj.FiltFile{Ind}];
            end
            obj.LocMaxFile=[obj.RawFile{Ind}(1:end-obj.nC.nExtRaw-1) obj.nC.ExtLocMax];
            Target=[obj.LocMaxFolder filesep obj.LocMaxFile];
            if obj.bTM.plotResults
                if ~isfolder(obj.bTM{Ind}.plotFolder)
                    mkdir(obj.bTM{Ind}.plotFolder);
                end
            end
            obj.bTM{Ind}.plotName=[obj.RawFile{Ind}(1:end-obj.nC.nExtRaw-1) '_' obj.bTM.plotName obj.nC.plotExt];
            blindTemplateMatchingGPU(Origin,Target,obj.bTM{Ind})% probably want to add a few inputs here, and some output?
        end
        %cluster local minima from a session
        function obj=sortSession(obj,Ind)
            Origin=[obj.LocMaxFolder filesep obj.LocMaxFile{Ind}];
            obj.singleSessionFile=[obj.RawFile{Ind}(1:end-obj.nC.nExtRaw-1) obj.nC.ExtLocMax];
            Target=[obj.singleSessionFolder filesep obj.singleSessionFile];
            if obj.Clust.plotResults
                if ~isfolder(obj.Clust{Ind}.plotFolderE)
                    mkdir(obj.Clust{Ind}.plotFolderE);
                end
                if ~isfolder(obj.Clust{Ind}.plotFolderC)
                    mkdir(obj.Clust{Ind}.plotFolderC);
                end
            end
            obj.Clust{Ind}.plotNameE=[obj.RawFile{Ind}(1:end-obj.nC.nExtRaw-1) '_' obj.Clust.plotNameE obj.nC.plotExt];
            obj.Clust{Ind}.plotNameC=[obj.RawFile{Ind}(1:end-obj.nC.nExtRaw-1) '_' obj.Clust.plotNameC obj.nC.plotExt];
            obj.Clust{Ind}=clusterSession(Origin, Target, obj.Clust{Ind});
        end
        %merge list of recordings
        function obj=mergeAll(obj,Mask,FileName)
            %need to check if there is already a sorted file
            assert ~isfile([obj.mergedFolder filesep FileName],'file')
            obj=mergeAllLocMax(obj,[obj.mergedFolder filesep FileName],Mask);
        end
        %add next recording to merge
        function obj=mergeNext(obj,Ind)
            mergeIncLocMax(mBase,outBase,OutBaseCAR,Files(i),subject)
        end
        %merge hash (no incremental version here)
        function obj=mergeHash(obj,Ind)
            mergeAllHash(mBase,outBase,Files(2:4),subject)
        end
    end
end