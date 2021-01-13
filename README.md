#SuperSessioning
a toolbox for cross session sorting of neural data

##Scope
####Technical
- Characterize local minima in the voltage traces across the array, and obtain an estimate of overall recording quality.
- Establish a metric for individual spikes for cross session comparisons, which is robust to additive noise.
- Use a minimal temporal windows for sorting and reduce overlap between spikes

####Experimental
- Chronically implanted electrode arrays in awake behaving animals, with no or little electrode drift across recording sessions.
- Sorting on individual electrodes, assuming electrode distances >100 µm
- Perform clustering in a noise dependent way, and ignore bouts of high noise across the array in the clustering step.
- visualize 

##Overview of the toolbox
- roughly split into single session and cross session analysis
- use one .mat file to keep track of parameters, files etc. across sessions
- same file keeps track of which sessions have been merged

####Operating modes
- (TODO) function that (currently) does detection and sorting for single sessions only, plus RecordingExtractors and SortingExtractors, for integration with [SpikeInterface](https://github.com/SpikeInterface).
- (TODO) function that does the sorting for one recording session, without merging.
- script for stepwise analysis, and session merges: SortMaster.m

## Detection and sorting pipeline
1. **extract Data and recording parameter**
 	- extractDataFromKwik.m
	- manuallyAnnotatedHdf
[//]: # (1. (**test for potential amplifier recalibrations**))
[//]: # (	- StepTest (TODO))
1. **remove 60 Hz component**
	- filter60Hz
1. **common average referencing**
	- filterCAR
1. **determine signal variance**
and determine segments of low variance
	- estimateStd
1. **match templates**
and obtain 3D densities of local minima
	- blindTemplateMatching
1. **find clusters**
use only segments with 'low' variance
	- clusterSession

### Matching across sessions
1. **match across sessions**
	- mergeAll
	-mergeNext
	-mergeHash

##Structure of sorted data
####SpikeSorter instance X
####single session cluster
####cross session cluster
Output structure z, saved as 'All_cat.mat', with variables

variables for individual spikes | description
-------------|-----------------
z.Spikes\{id\} | cluster indices of spikes for recording session 'id'
z.Times\{id\} | time stamps  of spikes for recording session 'id'
[//]:(TODO (?) | spike amplitudes)
[//]:(TODO (?) | spike widths)
[//]:(TODO (?) | spike symmetry)

variables for units | description
------|--------------------
z.Channel(u,1) | Electrode site on which cluster 'u' was recorded
z.nSpk(id,u) |  firing rate of cluster 'u' for recording session 'id' \[Hz\]
z.fracBorder(id,u) | Unit isolation: relative density of spikes in border voxels
z.pwt(id,u,1) | average amplitude of spikes in cluster 'u' for recording session 'id' \[bins in histogram\]
z.pwt(id,u,2) | average width of spikes in cluster 'u' for recording session 'id'  \[bins in histogram\]
z.pwt(id,u,3) | average symmetry of spikes in cluster 'u' for recording session 'id'  \[bins in histogram\]
z.Shapes(id,u,:) | median spike shapes ofcluster 'u' for recording session 'id'
z.ShapesD(id,u,:) | first and third quartile of spike shapes of cluster 'u' for recording session 'id'
z.Hist\{id\}\{u\}| 3D Histogram of spike shapes in cluster 'u' for recording session 'id'
z.Sessions(u,id) | boolean matrix for which units found in which sessions

Scaling | description
------|--------------------
z.xPwr | bin edges for amplitudes
z.tWidth | bin centers for spike FWHM
z.tTail | bin centers for symmetry

Recording sessions | description
------|--------------------
z.RecId\{id\} | Identifier of recording session 'id' (yyyy-mm-dd_hh-mm-ss)
z.LenRec(id) | Duration of recording session 'id'

####cross session hash
Output structure z, saved as 'subject\_chXX\_hash.mat'

variables for individual events | description
-------------|-----------------
z.hashSpikes\{id\} | 3D voxel locations spikes for recording session 'id' (reshape histogram as \[z.Ntail x z.Nwidth x z.Namp\])
z.hashTimes\{id\} | time stamps  of spikes for recording session 'id'

## Parameter for sorting

#### sorting



####merging clusters
parameter description | name | default value
---------------|---------|-----------------
threshold for JS divergence | JSthreshold | 0.2*log(2)
number of consecutive sessions to consider for merges | nMerge | 3
acceptable range of (penalized) cross-session shifts in amplitude \[bins\] | Nampshift | 5
acceptable range of (penalized) cross-session  shifts in width \[bins\] | Nwshift | 1
acceptable range of (penalized) cross-session  shifts in symmetry \[bins\] | Ntshift | 1
peak time of spike in shape cutouts \[frames\] | NcutPre | 20
length of spike shape cutouts \[frames\] | Ncut | 60