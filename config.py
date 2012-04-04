# emacs =  -*- mode =  python; py-indent-offset =  4; indent-tabs-mode =  nil -*-
# vi =  set ft=python sts=4 ts=4 sw=4 et = 

"""
	Set which FSL to use
"""
FSLDIR = '/usr/share/fsl/'

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
subj_file = '/home/ssikka/nki_nyu_pipeline/data/subjects.txt'
log_file = None
rest_name = 'rest'
anat_name = 'mprage'
standard_res = '3mm'
fwhm = [6]
rest_hp = [0.005, 0.01]
rest_lp = 0.1
alff_HP = [0.005]
alff_LP = 0.1


"""
	Mandatory
	where are your anatomical scans located relative to your Subjects Directory
"""

anat_template = '%s/*/*/%s.nii.gz'

"""
	Mandatory
	where are your functional scans located relative to your Subjects Directory
"""
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
rsfc_template = '%s/*/*/%s.nii.gz'

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
rsfc_warp_template = '%s/*/*/*/%s'


# all, basic, scrubbing, nuisance, alff, rsfc, vmhc, reho, group_analysis
analysis = [False, False, False, False, False, True, False, False, False ]
run_on_grid = False
qsub_args = '-q all.q'
num_cores = 8
