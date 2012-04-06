# emacs =  -*- mode =  python; py-indent-offset =  4; indent-tabs-mode =  nil -*-
# vi =  set ft=python sts=4 ts=4 sw=4 et : 

""" 
**************** Distributed Processing Settings ****************
	run_on_grid : options True/False or 0/1
	qsub_args : For use only when run_on_grid = True.

	-q QUEUE.q : The queue of nodes to use on the cluster/grid, default is set to the queue of all nodes on the cluster(all.q).

	
	-pe A-B : Specify the range of the number of slots to use on each node. A is the lower bound on the range and B is the higher bound on the number of slots available

""" 
run_on_grid = False
qsub_args = '-q all.q'
num_cores = 8

"""
******************************************************************
""" 

""" 
******************** Directory Setup *******************
""" 

"""
-------------------- Nipype Specifications ---------------
"""

"""
Nipype Cache/Working/Debug Directory

"""
workingDirectory = '/home/ssikka/nki_nyu_pipeline/working_dir/'

"""
	Crash Log Directory
"""
crashLogDirectory = '/home/ssikka/nki_nyu_pipeline/crash'










"""
-------------------- Data Specifications ---------------
"""

"""
	Subject Data Directory and List
"""
subjectDirectory = '/home/ssikka/nki_nyu_pipeline/data'
subjectList = '/home/ssikka/nki_nyu_pipeline/data/subjects.txt'


"""
Anatomical File Name and Location within Subject Directory
Note: anatomicalFileName substituted into anatomicalFilePath 
"""
anatomicalFileName = 'mprage'
anatomicalFilePath = '%s/*/*/%s.nii.gz'


"""
	anatLogFile: In case of multiple anatomical scans per subject. Specify the name of the log file which indicates which single anatomical scan to process.

	anatLogFilePath: Specify the location of the log file relative to the subject directory location. The pipeline substitutes 1st %s with the subject number/name and the 2nd %s is the name of the log file. 
The User needs to specify only the structure as per his analysis. 


	Ex: /sam/wave1/sub8000/
						anat_1/mprage.nii.gz
						anat_2/mprage.nii.gz
						logs/log.txt
	anatLogFile = ‘log.txt’
	anatLogFilePath = ‘%s/*/%s’
	
"""
anatLogFile = ‘log.txt’
anatLogFilePath = '%s/*/*/%s.nii.gz'



"""
Functional File Name and Location within Subject Directory
Note: functionalFileName substituted into functionalDirectorySetup template
"""
functionalFileName = 'rest'
functionalFilePath = '%s/*/%s.nii.gz'




"""
	Output Target Directory 
"""
sinkDirectory = '/home/ssikka/nki_nyu_pipeline/data'


"""
-------------- Miscellaneous Specifications -------------------
"""


"""
	Set FSL Directory (For Purposes of Template Selection)
"""
FSLDIR = '/usr/share/fsl/'


"""
	Tissue Priors Directory
"""
priorDirectory = '/home/ssikka/nki_nyu_pipeline/tissuepriors'

"""
*********************************************************
"""

"""
******* Optional Header and Timeseries Overrides ********
n_vols: # of volumes in functional timeseries (defaults to image header specification)
TR = time of repetition (defaults to image header specification)
	start_idx : starting time point(defaults to 0)
	stop_idx : stop time point (defaults to [n_vols – 1])
"""

start_idx = 0
stop_idx = 119
n_vols = 120
TR = 2.5



"""
**************************************************************
"""

""" 
******************** Image Processing Specs ******************
NOTE: Specification of multiple parameters for any step indicates
execution of a unique processing path for each parameter from 
current step onward (e.g., fwhm = [4.5 6] will produce two unique output directories – one for 4.5 spatial smoothing and one for 6mm
""" 



""" 
MNI Template Resolution Specification For Registration
""" 
standardResolution = '2mm'

""" 
------------------------- Spatial Filter --------------------
Lowpass 3D Gaussian Filter Kernel Specification in mm [FSL: fslmaths]
note: Set to [] or 0 to skip spatial filtering
"""
fwhm = [6, 4.5, 9, 0 ]



""" 
------------------- Nuisance Signal Correction --------------
"""

"""
	Select Desired Approach:
-	Global mean signal regression
-	Compcor 
-	White Matter & CSF regression (WMCSF)
-	First principal component regression
-	Medial Angle Correction
-	Motion Regression (run in combination with any approaches specified, or alone if no other approach selected)

"""
#global_signal, compcor, wmcsf, firstprinc, median angle, motion
Corrections = [
[True, False, False, False, False, False],
[False, True, False, False, False, True ] 
    ]


""" 
	For Compcor Use Only: Number of components to regress 
"""
ncomponents = [5, 6]


"""
	For Median Angle Correction Only: Target Angle (degrees)
"""
target_angle_deg = [90, 60]


""" 
***************************************************************
"""



 
"""
------------------------- Temporal Filtering -----------------
Temporal Filter Specification in Hz performed after Nuisance Signal Correction
note: Set to [] or 0 to skip high- or low-pass filter
"""
highPassFreqNuisance = [0.005, 0.01]
lowPassFreqNuisance = 0.1


""" 
******************** Derivative Specification *****************
""" 

""" 
Scrub data prior to derivate generation: In accord with Power et al. (2012); forking not enable yet for this step (next version).
Default value True/False
""" 
scrubData = False


""" 
Select Derivates:
# alff, sca, vmhc, reho, unit (e.g., parcellation unit, mask) timeseries extraction, voxel timseries extraction, vertices extraction
-ALFF/fALFF
-SCA(seed based correlation analysis)
-VMHC
-REHO
-Unit TimeSeries Extraction
-Vertices Extraction
""" 
derivates: = [True, True, True, True, True, True, True, True]


""" 
------------------------ ALFF/fALFF Options ---------------------------
""" 

"""
	For ALFF/fALFF only
Notes: 1) this derivative is allergic to scrubbed data and thus will never use them. 
"""
highPassFreqALFF = [0.005]
lowPassFreqALFF = 0.1



""" 
---------------- Seed-Based Correlation Analysis Options -----------------
""" 

""" 
	For SCA use only
"""
seed_file = '/home/ssikka/nki_nyu_pipeline/seed_list.txt'


""" 
------------- Timeseries Extraction Options ---------------------------
""" 

"""
	For Unit Timeseries Extraction Only
Note: Definitions Directory should contain one subdirectory for each set of units to be generated (e.g., Harvard-Oxford Atlas, AAL, Craddock, Dosenbach-160); one output file / set define   
""" 
unitDefinitionsDirectory =’/home/ssikka/unitDefSets’

# Output type: .csv, numPy
unitTSOutputs = [true, true]

""" 
	For Voxel Timeseries Extraction Only
Note: Definitions Directory should contain one subdirectory for each mask/mask set to be used to select voxels to be output; one output file / mask 
""" 
voxelMasksDirectory =’/home/ssikka/unitDefSets’


# Output type: .csv, numPy
voxelTSOutputs = [False, True]

""" 
	For Vertices Timeseries Extraction Only
""" 
# Output type: .csv, numPy
verticesTSOutputs = [False, True]


""" 
**************************************************************
""" 

""" 
********************* FSL Group Analysis *********************
Notes: 
- Separate group analysis conducted for each derivative
- Not applicable to time series extraction derivatives
"""

""" 
Specify Model List File that contains one or more models to be executed per derivate
""" 
modelsList=’home/ssikka/myModels.txt’

z_threshold=2.3
p_threshold=0.05
f_test=’yes’

“””
********************** THAT’S ALL FOLKS! **********************
“””

