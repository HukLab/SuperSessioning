%% Example master script 
%This script can be used (and modified, renamed, etc.) to control the sorting process

%define folders 
BaseFolder='';
RawFolder='';

%% create an instance of the spike sorter
%needs to be evaluated once, doesn't take any data, but sets default values
%for the sorting, and creates folders for sorting.
%-needs: 1. Folder where the sorting metadata (including this instance of the sorter) will be stored
%        2. Folder where the raw data can be found
%
% comment this section once you have created an instance of the sorter!

X=SuperSessioning(BaseFolder,RawFolder);

%can optionally change default values here

%save X if not existent already.
assert(~isfile([BaseFolder filesep 'spikeSorter.mat']),'An instance of the sorter already exists.')
save([BaseFolder filesep 'spikeSorter.mat'],'X','-v7.3');

%% Load the sorter and add data or perform sorting
%(can do this after each recording session, and add data as they are generated)

load([BaseFolder filesep 'spikeSorter.mat'])

%% add data
%can copy and paste from file browser (would be nice to have a timestamp in
%the file name, e.g. XXX_2017-11-13_13-12-21.raw.kwd, otherwise need to
%adapt method to create a sortable timestamp.
File='';
[~,NAME,EXT] = fileparts(File);
Filename=[NAME EXT];
TimeStamp=second(datestr(datenum(regexp(NAME,...
    '^(\d{4})-(\d{2})-(\d{2})_(\d{2})-(\d{2})-(\d{2})','tokens'))));

[X,Ind]=X.addRaw(FileName,TimeStamp);
%returns current index of the raw file, can change some session specific
%parameters here

%% filter powerline noise and common average referencing
X=X.filter60Hz(Ind);
X=X.filterCAR(Ind);
% estimate standard deviation of high-pass filtered signal
X=X.estimateStd(Ind);

%save temporary backup of current sorter
save([BaseFolder filesep 'spikeSorter_temp.mat'],'X','-v7.3');

%% detect local matches with predefined templates
X=X.blindTemplateMatching(Ind);

%save temporary backup of current sorter
save([BaseFolder filesep 'spikeSorter_temp.mat'],'X','-v7.3');

%% cluster shape histograms
X=X.sortSession(Ind);

%save temporary backup of current sorter
save([BaseFolder filesep 'spikeSorter_temp.mat'],'X','-v7.3');

%% merge across sessions
%(this can be used to merge a list of sessions at once. When already merged
%sessions are available, comment this section and use the incremental
%version)
X=X.mergeAll(Mask,FileName);


%save temporary backup of current sorter
save([BaseFolder filesep 'spikeSorter_temp.mat'],'X','-v7.3');

%% add another session to merge
%make sure this session is temporally after sessions that were already merged
assert(X.TimeStamp(Ind)>X.TlastMerged,'Session in between already merged ones.')
X=X.mergeNext(Ind);

%% accumulate multi-unit Hash data
%need to do once, after recording a bunch of sessions (?)
%should run on a per electrode basis.
X=X.mergeHash(Ind);

%% after each session (?!)
%backup old X, save new X
copyfile([BaseFolder filesep 'spikeSorter.mat'], [BaseFolder filesep 'spikeSorter.mat.bak']);
save([BaseFolder filesep 'spikeSorter.mat'],'X','-v7.3');