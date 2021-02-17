%% Example master script 
%This script can be used (and modified, renamed, etc.) to control the sorting process

%define base folder (to load instance of the sorter)
BaseFolder='/home/huklab/XXX';

%subject identifier
subject='test';
%where the raw data can be found
RawFolder='/home/huklab/YYY';

%% create an instance of the spike sorter
%needs to be evaluated once, doesn't take any data, but sets default values
%for the sorting, and creates folders for sorting.
%-needs: 1. Folder where the sorting metadata (including this instance of the sorter) will be stored
%        2. Folder where the raw data can be found
%
% comment this section once you have created an instance of the sorter!

X=SuperSessioning(BaseFolder,RawFolder,subject);

%can optionally change default values here

%save X if not existent already.
assert(~isfile([BaseFolder filesep 'spikeSorter.mat']),'An instance of the sorter already exists.')
save([BaseFolder filesep 'spikeSorter.mat'],'X','-v7.3');

%% Load the sorter and add data or perform sorting
%(can do this after each recording session, and add data as they are generated)

load([BaseFolder filesep 'spikeSorter.mat'])

%% add data
%can copy and paste from file browser (would be nice to have a timestamp in
%the file name, e.g. test_2017-11-13_13-12-21.raw.kwd, otherwise need to
%adapt method to create a sortable timestamp.

File='/home/huklab/YYY/test_2017-11-13_13-12-21.raw.kwd';
%%
[~,NAME,EXT] = fileparts(File);
FileName=[NAME EXT];
h=regexp(NAME,'(\d{4})-(\d{2})-(\d{2})_(\d{2})-(\d{2})-(\d{2})','tokens');
h=h{:};
TimeStamp=datenum(str2double(h{1}),str2double(h{2}),str2double(h{3}),str2double(h{4}),...
    str2double(h{5}),str2double(h{6}));

[X,Ind]=X.addRaw(FileName,TimeStamp);
%returns current index of the raw file, can change some session specific
%parameters here:
%% If needed: restrict analysis window
%take only first 5 min of data
X.recParameter{Ind}.tEnd=5;

%% If needed: change relevant dataset (KWIK format)
%when recording has initially been stopped and then continued. Not
%considering concatenating those data, as usually the first/second part was not
%recorded intentionally and likely has a lot of noise.
X.recParameter{Ind}.iGroup=1;%(0 default, i.e. first recording bout)

%% filter powerline noise and common average referencing
%filter for old version (to be deprecated)
%X=X.filter60Hz(Ind);
%X=X.filterCAR(Ind);
%new version
X=filterGeneral(X,Ind);
%% estimate standard deviation of high-pass filtered signal
X=X.estimateStd(Ind);

%save temporary backup of current sorter
save([BaseFolder filesep 'spikeSorter_temp.mat'],'X','-v7.3');

%% detect local matches with predefined templates
X=X.blindTemplateMatching(Ind);
%X=X.blindTemplateMatching(1);
%save temporary backup of current sorter
save([BaseFolder filesep 'spikeSorter_temp.mat'],'X','-v7.3');

%% cluster shape histograms
X=X.sortSession(Ind);

%save tempordilary backup of current sorter
save([BaseFolder filesep 'spikeSorter_temp.mat'],'X','-v7.3');

%% merge across sessions
%(this can be used to merge a list of sessions at once. When already merged
%sessions are available, comment this section and use the incremental
%version)
SessInd=[1 2];
FileName=[X.subject '_All_cat.mat'];
X=X.mergeAll(FileName,SessInd);

%save temporary backup of current sorter
save([BaseFolder filesep 'spikeSorter_temp.mat'],'X','-v7.3');

%% add another session to merge
X=X.mergeNext(Ind);

%% accumulate multi-unit Hash data
%need to do once, after recording a bunch of sessions (?)
%should run on a per electrode basis.
X=X.mergeHash(X.mergedList);

%% after each session (?!)
%backup old X, save new X
copyfile([BaseFolder filesep 'spikeSorter.mat'], [BaseFolder filesep 'spikeSorter.mat.bak']);
save([BaseFolder filesep 'spikeSorter.mat'],'X','-v7.3');