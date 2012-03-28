# emacs =  -*- mode =  python; py-indent-offset =  4; indent-tabs-mode =  nil -*-
# vi =  set ft=python sts=4 ts=4 sw=4 et = 

FSLDIR = '/frodo/shared/RH5_fsl/'
scripts_dir = '/home/sharad/nki_nyu_pipeline'
working_dir = '/home/sharad/nki_nyu_pipeline/working_dir/'
prior_dir = '/home/sharad/nki_nyu_pipeline/tissuepriors'
seed_file = 'seed_list.txt'
batch_file = 'batch_list.txt'
rest_name = 'rest'
anat_name = 'mprage'
standard_res = '3mm'
fwhm = [6.0, 5.0]
rest_hp = [0.005, 0.01]
rest_lp = [0.1]
alff_HP = [0.01]
alff_LP[0.1]

#defaulf,no_global,median_angle,compcor
which_regression = [['compcor', '/home/sharad/nki_nyu_pipeline/templates/nuisance_cc.fsf']]

ncomponents = []

#all, basic, nuisance,nuisance+alff,nuisance+alff+rsfc,nuisance+rsfc,rssc+alff,basic+alff,basic+rsfc
analysis = ['basic']
