# emacs =  -*- mode =  python; py-indent-offset =  4; indent-tabs-mode =  nil -*-
# vi =  set ft=python sts=4 ts=4 sw=4 et = 

"""
    Set which FSL to use
"""
FSLDIR = '/usr/share/fsl/4.1/'

"""
    Point to directory where your subjects reside
"""
subj_dir = '/home/ssikka/nki_nyu_pipeline/data'

"""
    Point to directory where pipeline can store results 
"""
sink_dir = '/home/ssikka/nki_nyu_pipeline/data'

"""
    Set Temporary Directory where Nipype can store
    temporary results
"""
working_dir = '/home/ssikka/nki_nyu_pipeline/working_dir/'

"""
    Set the location where 3mm and 2mm Tissue Priors
    located.
"""
prior_dir = '/home/ssikka/nki_nyu_pipeline/tissuepriors'

"""
    Point to directory where pipeline can store crash .npz files
    if it crashes
"""
crash_dir = '/home/ssikka/nki_nyu_pipeline/crash'

"""
    Functional volumes to keep
    start_idx : starting volume
    stop_idx : Last volume
"""
start_idx = 0
stop_idx = 119
n_vols = 120
TR = 2.5

"""
    Seed file

    Each line contains full path to seed file

"""
seed_file = '/home/ssikka/nki_nyu_pipeline/seed_list.txt'

"""
    subj_file : '/data/ADHD200/docs/subjects.txt'
    Each line in the file contains subjects
    0010001
    0010002
    etc.
"""
func_session_file = '/home/ssikka/nki_nyu_pipeline/sessions.txt'
anat_session_file = '/home/data/Projects/nuisance_reliability_paper/work_lists/session_list.txt'
subj_file = '/home/ssikka/nki_nyu_pipeline/data/subjects.txt'
log_file = None
standard_res = '2mm'
fwhm = [6]

"""
	Set various thresholds for cerebral spinal fluid , white matter 
	and gray matter mask generation during segmentation
"""
cerebralSpinalFluidThreshold = [0.4]
whiteMatterThreshold = [0.66]
grayMatterThreshold = [0.2]


"""
    Value : 1 or 0
"""
nuisanceHighPassFilter = 1
nuisanceLowPassFilter = 1

"""
    When nuisanceHighPassFilter : 1
    Set nuisanceHighPassLowCutOff to a decimal value

    When nuisanceHighPassFilter : 0
    Set nuisanceHighPassLowCutOff to None

    Same holds for nuisanceLowPassFilter & nuisanceLowPassHighCutOff

"""
nuisanceHighPassLowCutOff = [0.01]
nuisanceLowPassHighCutOff = [0.1]


""" 
------------------------ ALFF/fALFF Options ---------------------------
"""

"""
    For ALFF/fALFF only
    Notes: 1) this derivative is allergic to scrubbed data and thus will never use them.
           2) You need to specify both highPassFreqALFF and lowPassFreqALFF if you intend
              to use this derivative. The Default values are set below.  
"""

highPassFreqALFF = [0.005]
lowPassFreqALFF = 0.1

""" 
    Scrub data prior to derivate generation: In accord with Power et al. (2012); forking not enable yet for this step (next version).
    Default value 1/0
"""
scrubData = [1, 0]



"""
    Number of components to regress with compcor
"""
ncomponents = [5, 6]

"""
    Target Angle in Degree for Median Angle Correction
"""

target_angle_deg = [90, 60]


"""
    Which Signals do you which to regress out
"""
#['global', 'compcor', 'wm', 'csf', 'gm', 'firstprinc', 'motion']
#regressors = [0, 0, 0, 0, 0, 1]
Corrections = [
                [1, 0, 0, 0, 0, 0, 1],
                [0, 0, 1, 1, 1, 1, 1 ]
              ]


"""
    Mandatory
    where are your anatomical scans located relative to your Subjects Directory
"""
anat_template = '%s/anat/%s.nii.gz'
anat_template_list = ['subject', 'mprage']


"""
    Mandatory
    where are your functional scans located relative to your Subjects Directory
"""

func_template = '%s/%s/%s.nii.gz'
func_template_list = ['subject', 'session', 'rest']


"""
    SET ONLY when analysis is not set to 'all' and u need to run alff

    where are  rest_res, rest_mask, rest_mask2standard  scans located
    relative to your Subjects Directory
"""
alff_template = '%s/*/%s.nii.gz'

"""
    SET ONLY when analysis is not set to 'all' and u need to run alff

    where are  example_func2highres.mat, highres2standard_warp.nii.gz scans located
    relative to your Subjects Directory
"""
alff_warp_template = '%s/*/*/%s'

"""
    SET ONLY when analysis is not set to 'all' and u need to run 
    generate functional connectivity maps. 

    where are  rest_res2standard.nii.gz, rest_res_filt.nii.gz,
    rest_mask2standard.nii.gz, example_func.nii.gz scans located
    relative to your Subjects Directory
"""
ifc_template = '%s/*/%s.nii.gz'

"""
    SET ONLY when analysis is not set to 'all' and u need to run 
    generate functional connectivity maps. 

    where are example_func2highres.mat,
              highres2standard_warp.nii.gz,
              stand2highres_warp.nii.gz
              highres2example_func.mat 
    scans located
    relative to your Subjects Directory
"""
ifc_warp_template = '%s/*/*/%s'

vmhc_rest_res_template = '%s/*/%s.nii.gz'
vmhc_anat_reorient_template = '%s/*/%s.nii.gz'
vmhc_example_func2highres_mat_template = '%s/*/*/%s'

""" 
------------- Timeseries Extraction Options ---------------------------
"""

"""
    For Unit Timeseries Extraction Only
Note: Definitions Directory should contain one subdirectory for each set of units to be generated (e.g., Harvard-Oxford Atlas, AAL, Craddock, Dosenbach-160); one output file / set define   
"""
unitDefinitionsDirectory = '/home/ssikka/nki_nyu_pipeline/tsdata'

# Output type: .csv, numPy
unitTSOutputs = [1, 1]

""" 
    For Voxel Timeseries Extraction Only
Note: Definitions Directory should contain one subdirectory for each mask/mask set to be used to select voxels to be output; one output file / mask 
"""
voxelMasksDirectory = '/home/ssikka/nki_nyu_pipeline/tsdata'


# Output type: .csv, numPy
voxelTSOutputs = [0, 1]

""" 
    For Vertices Timeseries Extraction Only
"""
# Output type: .csv, numPy
verticesTSOutputs = [0, 1]
runSurfaceRegistraion = [0]
reconSubjectsDirectory = '/home/ssikka/nki_nyu_pipeline/recon_subjects'
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
    List of one or more derivatives on which the group anlaysis is be run
    As a default, some of derivatives which are generated from the pipleine have
    been specified
"""
derivativeFile = '/Users/ranjeet.khanuja/Documents/workspace/GroupAnalysis/seed_list.txt'

""" 
Specify Model List File that contains one or more models to be executed per derivate
"""
modelsFile = '/Users/ranjeet.khanuja/Documents/workspace/GroupAnalysis/model_list.txt'
"""
     Templates for Group Analysis
     Design Files
     .con -> list of contrasts requested.
     .fts -> list of F-tests requested.
     .mat -> the actual values in the design matrix   
     
     model_name in the template_list will be fetch from  modesList above
"""
mat = '/Users/ranjeet.khanuja/Desktop/data2/models/%s/%s.mat'
con = '/Users/ranjeet.khanuja/Desktop/data2/models/%s/%s.con'
fts = '/Users/ranjeet.khanuja/Desktop/data2/models/%s/%s.fts'
grp = '/Users/ranjeet.khanuja/Desktop/data2/models/%s/%s.grp'

mat_template_list = ['model_name', 'model_name']
con_template_list = ['model_name', 'model_name']
fts_template_list = ['model_name', 'model_name']
grp_template_list = ['model_name', 'model_name']

"""
    Derivative Template
    derivative in the template_list will be fetched from the derivative file 
    defined above
"""
derv_template = '/Users/ranjeet.khanuja/Desktop/data2/subjects/*/%s.nii.gz'
derv_template_list = ['derivative']

z_threshold = 2.3
p_threshold = 0.05
f_test = 1

# all, basic, scrubbing, nuisance, alff, ifc, vmhc, reho, group_analysis
analysis = [0, 1, 0, 0, 0, 0, 0, 0, 0 ]
run_on_grid = 0
qsub_args = '-q all.q'
num_cores = 10
