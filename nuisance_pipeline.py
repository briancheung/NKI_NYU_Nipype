#!/frodo/shared/epd/bin/python
import e_afni
import nki_nodes
import sys
import os
import commands
import nipype.pipeline.engine as pe
import nipype.algorithms.rapidart as ra
import nipype.interfaces.afni as afni
import nipype.interfaces.fsl as fsl
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
from nipype.interfaces.afni.base import Info, AFNITraitedSpec, AFNICommand
from nipype.interfaces.base import Bunch, TraitedSpec, File, Directory, InputMultiPath
from nipype.utils.filemanip import fname_presuffix, list_to_filename, split_filename
from nipype.utils.docparse import get_doc
from nipype.interfaces.fsl import MultiImageMaths
from nipype.interfaces.fsl import ExtractROI
from nipype.interfaces.fsl.maths import DilateImage
from multiprocessing import Process
from multiprocessing import Pool
from ConfigParser import SafeConfigParser

#*******************************************************************************************
dir_file = None
subjectList = None
FWHM = None
FSLDIR = ''
scripts_dir = ''
prior_dir = ''
anat_name = ''
anat_dir = ''
rest_name = ''
rest_dir = ''
working_dir = ''
seed_list= ''
standard_brain = ''
batch_list = ''
standard_res = ''
logFile = ''
nuisance_template = ''
nuisance_template_cc = ''
nuisance_template_noglobal = ''
recon_subjects = ''
seeds_dir = ''
which_regression = ''
infosource = None
datasource = None
datasource_warp = None
workflow = None
renamer = None
HP = None
LP = None
hp = None
lp = None
#*******************************************************************************************



def signalExtractionEfficient():

	import sys
	global workflow
	global scripts_dir
	global standard_brain
	global rest_name
	global FSLDIR
	global standard_res
	global infosource
	global datasource_warp
	global datasource
	global which_regression
	
	workflow.connect(datasource,'rest_pp',nki_nodes.lifeSaverNuisancepp,'in_file')
	workflow.connect(datasource,'rest_pp',nki_nodes.adjustPath,'rest_path')
	workflow.connect(infosource,'strategy',nki_nodes.adjustPath,'strategy')
	workflow.connect(infosource,'sublist',nki_nodes.adjustPath,'sublist')
	workflow.connect(nki_nodes.adjustPath,'result_path',nki_nodes.lifeSaverNuisancepp,'format_string')

	workflow.connect(datasource,'rest_mask',nki_nodes.lifeSaverNuisanceMask,'in_file')
	workflow.connect(datasource,'rest_mask',nki_nodes.adjustPath1,'rest_path')
	workflow.connect(infosource,'strategy',nki_nodes.adjustPath1,'strategy')
	workflow.connect(infosource,'sublist',nki_nodes.adjustPath1,'sublist')
	workflow.connect(nki_nodes.adjustPath1,'result_path',nki_nodes.lifeSaverNuisanceMask,'format_string')


	workflow.connect(nki_nodes.lifeSaverNuisancepp,'out_file',nki_nodes.TR,'in_files')
	workflow.connect(nki_nodes.lifeSaverNuisancepp,'out_file',nki_nodes.NVOLS,'in_files')
	workflow.connect(nki_nodes.NVOLS,'nvols',nki_nodes.nuisance_poly,'nvols')	
	if(which_regression.lower() == 'default'):
		#make initial connections to facilitate the default regression : global, wm, csf
		workflow.connect(nki_nodes.lifeSaverNuisancepp,'out_file',nki_nodes.nuisance_globalE,'in_file')
		workflow.connect(datasource_warp,'global_mask',nki_nodes.nuisance_globalE,'mask')
		workflow.connect(nki_nodes.lifeSaverNuisancepp,'out_file',nki_nodes.nuisance_csf,'in_file')
		workflow.connect(datasource_warp,'csf_mask',nki_nodes.nuisance_csf,'mask')
		workflow.connect(nki_nodes.lifeSaverNuisancepp,'out_file',nki_nodes.nuisance_wm,'in_file')
		workflow.connect(datasource_warp,'wm_mask',nki_nodes.nuisance_wm,'mask')

		workflow.connect(nki_nodes.lifeSaverNuisancepp,'out_file', nki_nodes.MC,'pp')
		workflow.connect(infosource,'nuisance_template',nki_nodes.FSF,'nuisance_template')
		workflow.connect(nki_nodes.lifeSaverNuisancepp,'out_file',nki_nodes.FSF,'rest_pp')
		workflow.connect(nki_nodes.TR,'TR', nki_nodes.FSF,'TR')
		workflow.connect(nki_nodes.NVOLS,'nvols',nki_nodes.FSF,'n_vols')
		workflow.connect(nki_nodes.MC,'EV_Lists',nki_nodes.copyS,'EV_Lists')
		workflow.connect(nki_nodes.FSF,'nuisance_files',nki_nodes.copyS,'nuisance_files')
		workflow.connect(nki_nodes.nuisance_globalE,'out_file',nki_nodes.copyS,'global1Ds')
		workflow.connect(nki_nodes.nuisance_csf,'out_file',nki_nodes.copyS,'csf1Ds')
		workflow.connect(nki_nodes.nuisance_wm,'out_file',nki_nodes.copyS,'wm1Ds')
		
	
	elif(which_regression.lower() == "no_global_signal"):

		workflow.connect(nki_nodes.lifeSaverNuisancepp,'out_file',nki_nodes.nuisance_csf,'in_file')
		workflow.connect(datasource_warp,'csf_mask',nki_nodes.nuisance_csf,'mask')
		workflow.connect(nki_nodes.lifeSaverNuisancepp,'out_file',nki_nodes.nuisance_wm,'in_file')
		workflow.connect(datasource_warp,'wm_mask',nki_nodes.nuisance_wm,'mask')
		
		workflow.connect(nki_nodes.lifeSaverNuisancepp,'out_file', nki_nodes.MC,'pp')
		workflow.connect(infosource,'nuisance_template_noglobal',nki_nodes.FSF,'nuisance_template')
		workflow.connect(nki_nodes.lifeSaverNuisancepp,'out_file',nki_nodes.FSF,'rest_pp')
		workflow.connect(nki_nodes.TR,'TR', nki_nodes.FSF,'TR')
		workflow.connect(nki_nodes.NVOLS,'nvols',nki_nodes.FSF,'n_vols')
		workflow.connect(nki_nodes.MC,'EV_Lists',nki_nodes.copyS,'EV_Lists')
		workflow.connect(nki_nodes.FSF,'nuisance_files',nki_nodes.copyS,'nuisance_files')
		workflow.connect(infosource,'dummy',nki_nodes.copyS,'global1Ds')
		workflow.connect(nki_nodes.nuisance_csf,'out_file',nki_nodes.copyS,'csf1Ds')
		workflow.connect(nki_nodes.nuisance_wm,'out_file',nki_nodes.copyS,'wm1Ds')

	workflow.connect(nki_nodes.nuisance_poly,'regressors',nki_nodes.copyS,'regressors')
	workflow.connect(datasource,'mc_oneD',nki_nodes.MC,'in_file')


def nuisance():

	

	global workflow
	global scripts_dir
	global standard_brain
	global rest_name
	global FSLDIR
	global standard_res
	global infosource
	global which_regression
	global datasource_warp
	global datasource

	signalExtractionEfficient()
	
	workflow.connect(nki_nodes.FSF,'nuisance_files',nki_nodes.nuisance_featM,'fsf_file')
	workflow.connect(nki_nodes.copyS,'EV_final_lists_set',nki_nodes.nuisance_featM,'ev_files')

	#which way to go
	if(which_regression == 'no_global_signal'):
		workflow.connect(nki_nodes.lifeSaverNuisancepp,'out_file',nki_nodes.nuisance_brick,'in_file')
		workflow.connect(nki_nodes.lifeSaverNuisancepp,'out_file',nki_nodes.nuisance_fgls,'in_file')
	
	else:
		workflow.connect(nki_nodes.lifeSaverNuisancepp,'out_file',nki_nodes.nuisance_brick,'in_file')
		workflow.connect(nki_nodes.lifeSaverNuisancepp,'out_file',nki_nodes.nuisance_fgls,'in_file')

	workflow.connect(datasource,'rest_pp_mask',nki_nodes.nuisance_brick,'mask')
	workflow.connect(nki_nodes.nuisance_featM, 'design_file',nki_nodes.nuisance_fgls,'design_file')
	workflow.connect(nki_nodes.lifeSaverNuisancepp, ('out_file',nki_nodes.getStatsDir),nki_nodes.nuisance_fgls, 'results_dir' )
	workflow.connect(nki_nodes.nuisance_brick,'min_val',nki_nodes.nuisance_fgls,'threshold')
	workflow.connect(nki_nodes.nuisance_fgls,'residual4d',nki_nodes.nuisance_stat,'in_file')
	workflow.connect(nki_nodes.nuisance_fgls,'residual4d',nki_nodes.nuisance_calc,'infile_a')
	workflow.connect(nki_nodes.nuisance_stat,'out_file', nki_nodes.nuisance_calc,'infile_b')
	workflow.connect(nki_nodes.nuisance_calc,'out_file',nki_nodes.func_filter,'in_file')
	workflow.connect( infosource,'hp',nki_nodes.func_filter,'highpass')
	workflow.connect( infosource,'lp',nki_nodes.func_filter,'lowpass')
	workflow.connect(nki_nodes.func_filter,'out_file', nki_nodes.nuisance_warp,'in_file')
	workflow.connect(infosource,'standard',nki_nodes.nuisance_warp,'ref_file')
	workflow.connect(datasource_warp, 'fieldcoeff_file',nki_nodes.nuisance_warp,'field_file')
	workflow.connect(datasource_warp,'premat',nki_nodes.nuisance_warp,'premat')
	workflow.connect(nki_nodes.lifeSaverNuisanceMask, 'out_file', nki_nodes.nuisance_warp_1,'in_file')
	workflow.connect(infosource,'standard',nki_nodes.nuisance_warp_1,'ref_file')
	workflow.connect(datasource_warp, 'fieldcoeff_file',nki_nodes.nuisance_warp_1,'field_file')
	workflow.connect(datasource_warp,'premat',nki_nodes.nuisance_warp_1,'premat')



def extract_seed_name(seed):

	fname = os.path.basename(seed)

	seed_name = (fname.split('.'))[0]

	return seed_name

def printToFile(time_series,seed_ts_dir,seed_name):

	print seed_ts_dir + '/%s.1D'%(seed_name)	
	print time_series
	f = open(seed_ts_dir + '/%s.1D'%(seed_name),'w')

	for tr in time_series:
		print >>f,"%f" %(tr)
	
	f.close()




def readSeedList():

	global seed_list
	global seeds_dir

	f= open(seed_list,'r')

	lines = f.readlines()

	f.close()

	seeds = []
	for line in lines:
		line = line.rstrip('\r\n')
		seeds.append(os.path.basename(line))
		seeds_dir = os.path.dirname(line)
	
	return seeds

def gatherData(sublist, analysisdirectory):

	global datasource
	global datasource_warp
	global working_dir
	global anat_name
	global anat_dir
	global rest_name
	global rest_dir
	global seed_list


	datasource = pe.Node( interface = nio.DataGrabber(infields=['subject_id'], outfields = [ 'mc_oneD', 'rest_pp', 'rest_pp_mask', 'rest_mask' ]) , name= 'datasource')
	datasource.inputs.base_directory = analysisdirectory
	#datasource.inputs.template = '%s/*/%s.nii.gz'
	datasource.inputs.template = '%s/*/*/%s'
	datasource.inputs.template_args = dict(mc_oneD  = [['subject_id',rest_name + '_mc.1D']] , rest_pp = [['subject_id',rest_name + '_pp.nii.gz']],
		rest_pp_mask = [['subject_id',rest_name + '_pp_mask.nii.gz']], rest_mask = [['subject_id',rest_name + '_mask.nii.gz']])
	datasource.iterables = ('subject_id', sublist)


	datasource_warp = pe.Node( interface = nio.DataGrabber(infields=['subject_id'], outfields = [ 'premat', 'fieldcoeff_file', 'global_mask','csf_mask', 'wm_mask' ]) , name= 'datasource_warp')
	datasource_warp.inputs.base_directory = analysisdirectory
	#datasource_warp.inputs.template = '%s/*/%s.nii.gz'
	datasource_warp.inputs.template = '%s/*/*/*/%s'
	datasource_warp.inputs.template_args = dict( premat = [['subject_id','example_func2highres.mat' ]] , fieldcoeff_file = [['subject_id','highres2standard_warp.nii.gz']], 
		global_mask = [['subject_id','global_mask.nii.gz' ]], csf_mask = [['subject_id','csf_mask.nii.gz' ]], wm_mask = [['subject_id','wm_mask.nii.gz' ]])
	datasource_warp.iterables = ('subject_id', sublist)


	nki_nodes.setSeedList(seed_list)

	#datasource_seeds = pe.Node( interface = nio.DataGrabber(infields=['seeds'],outfields = ['seeds']), name = 'datasource_seeds')
	#datasource_seeds.inputs.base_directory = seeds_dir
	#datasource_seeds.inputs.template = '%s'
	#datasource_seeds.inputs.seeds = seeds



def getFields():

	fields = ['strategy','sublist','brain_symmetric','symm_standard','standard_res_brain','standard','standard_brain_mask_dil','twomm_brain_mask_dil','config_file_twomm','config_file','PRIOR_CSF','PRIOR_WHITE','nuisance_template','nuisance_template_cc','nuisance_template_noglobal','HP','LP','hp','lp','which_regression','dummy']
	return fields

def getInfoSource(sublist, analysisdirectory):
	
	global infosource
	global FSLDIR
	global standard_res
	global prior_dir
	global nuisance_template
	global nuisance_template_cc
	global which_regression
	global HP
	global LP
	global hp
	global lp


	formula  = getFields()
	infosource = pe.Node(interface=util.IdentityInterface(fields= formula), name="infosource")

	infosource.inputs.standard_res_brain = os.path.abspath(FSLDIR + '/data/standard/MNI152_T1_%s_brain.nii.gz' %(standard_res))

	infosource.inputs.standard = os.path.abspath(FSLDIR + '/data/standard/MNI152_T1_%s.nii.gz' %(standard_res))
	infosource.inputs.standard_brain_mask_dil = os.path.abspath(FSLDIR + '/data/standard/MNI152_T1_%s_brain_mask_dil.nii.gz' %(standard_res))
	infosource.inputs.config_file = os.path.abspath(FSLDIR+'/etc/flirtsch/T1_2_MNI152_%s.cnf'%(standard_res))
	infosource.inputs.config_file_twomm = os.path.abspath(FSLDIR + '/etc/flirtsch/T1_2_MNI152_2mm.cnf')
	infosource.inputs.twomm_brain_mask_dil = os.path.abspath(FSLDIR+ '/data/standard/MNI152_T1_2mm_brain_mask_symmetric_dil.nii.gz')
	infosource.inputs.PRIOR_CSF = os.path.abspath(prior_dir + '/avg152T1_csf_bin.nii.gz')
	infosource.inputs.PRIOR_WHITE = os.path.abspath(prior_dir + '/avg152T1_white_bin.nii.gz')
	infosource.inputs.HP = float(HP)
	infosource.inputs.LP = float(LP)
	infosource.inputs.hp = float(hp)
	infosource.inputs.lp = float(lp)
	infosource.inputs.which_regression = which_regression
	infosource.inputs.strategy = strategy
	infosource.inputs.sublist = subjectList

	infosource.inputs.brain_symmetric = os.path.abspath(FSLDIR + '/data/standard/MNI152_T1_2mm_brain_symmetric.nii.gz')
	infosource.inputs.symm_standard = os.path.abspath(FSLDIR + '/data/standard/MNI152_T1_2mm_symmetric.nii.gz')

	infosource.inputs.nuisance_template = os.path.abspath(nuisance_template)
	infosource.inputs.nuisance_template_cc = os.path.abspath(nuisance_template_cc)
	infosource.inputs.nuisance_template_noglobal = os.path.abspath(nuisance_template_noglobal)
	infosource.inputs.dummy = os.path.abspath(FSLDIR + '/data/standard/MNI152_T1_%s_symmetric.nii.gz' %(standard_res))




def getSinkOther():

	datasink = None
	
	datasink = pe.Node(nio.DataSink(), name = 'sinker-other')
	datasink.inputs.base_directory = os.path.abspath('/')
	datasink.inputs.container = os.path.abspath('/')
	datasink.inputs.regexp_substitutions = [(r"_subject_id_(\w|\d)+", ''),(r"-",'/'),(r"(.)+_func_(\w|\d)+",''),(r"(.)+_reg_(\w|\d)+",''),(r"(.)+_seg_(\w|\d)+",''),(r"(.)+_nuisance_(\w|\d)+",''),(r"(.)+_alff_(\w|\d)+/",''),(r"(.)+_falff(\w|\d)+/",''),(r"(.)+_RSFC_(\w|\d)+",''),(r"(.)+/_seed_(.)+seeds\.\./",'')]


	return datasink




def adjustPath(rest_path):

	import re
	result_path = []
	
	if(isinstance(rest_path,list)):
		for r_p in rest_path:
			result_path.append(re.sub(r'\/','-',r_p))

		print str(result_path)
	else:
		result_path.append(re.sub(r'\/','-',rest_path))
		print str(result_path)
	return result_path



def adjustNuisance(rest_path):

	import re
	import sys
	result_path = []

	if(isinstance(rest_path,list)):
		for r_p in rest_path:

			split_path = r_p.split('-')

			path = ""
			for index in range(0,len(split_path) - 1):
				path += split_path[index] + "-"

			path += "nuisance" + "-"
			path += split_path[len(split_path) - 1]
		
			print "\nprocessed path %s\n" %(path)
			result_path.append(path)
		

	else:
		split_path = rest_path.split('-')

		path = ""
		for index in range(0,len(split_path) - 1):
			path += split_path[index] + "-"

		path += "nuisance" + "-"
		path += split_path[len(split_path) - 1]

		print "\nprocessed path %s\n" %(path)
		result_path.append(path)
		
	return result_path




def makeOutputConnectionsNuisance(datasink):

	global which_regression


	if not (which_regression == "compcor" ) and not(which_regression == "median_angle") :
	## Connect nuisance nodes to datasink	
	
		if not which_regression == "no_global_signal":

			workflow.connect(nki_nodes.nuisance_globalE_r, 'out_file',nki_nodes.lifeSaver_nuisance_globalE,'in_file')
			workflow.connect(nki_nodes.nuisance_globalE_r, ('out_file',adjustNuisance),nki_nodes.lifeSaver_nuisance_globalE,'format_string')
			workflow.connect(nki_nodes.lifeSaver_nuisance_globalE,'out_file',datasink,'@nuisance_globalE')
			#workflow.connect(nki_nodes.nuisance_globalE,'out_file',datasink,'@nuisance_globalE')
	
		workflow.connect(nki_nodes.nuisance_csf_r, 'out_file',nki_nodes.lifeSaver_nuisance_csf,'in_file')
		workflow.connect(nki_nodes.nuisance_csf_r, ('out_file',adjustNuisance),nki_nodes.lifeSaver_nuisance_csf,'format_string')
		workflow.connect(nki_nodes.lifeSaver_nuisance_csf,'out_file',datasink,'@nuisance_csf')
		#workflow.connect(nki_nodes.nuisance_csf,'out_file',datasink,'@nuisance_csf')
	
		workflow.connect(nki_nodes.nuisance_wm_r, 'out_file',nki_nodes.lifeSaver_nuisance_wm,'in_file')
		workflow.connect(nki_nodes.nuisance_wm_r, ('out_file',adjustNuisance),nki_nodes.lifeSaver_nuisance_wm,'format_string')
		workflow.connect(nki_nodes.lifeSaver_nuisance_wm,'out_file',datasink,'@nuisance_wm')
		#workflow.connect(nki_nodes.nuisance_wm,'out_file',datasink,'@nuisance_wm')


	workflow.connect(nki_nodes.nuisance_featM_r, 'out_file',nki_nodes.Saver_nuisance_featM,'in_file')
	workflow.connect(nki_nodes.lifeSaverNuisancepp,'out_file',nki_nodes.Saver_nuisance_featM,'pp')
	workflow.connect(nki_nodes.nuisance_featM_r, 'out_file',nki_nodes.lifeSaver_nuisance_featM,'in_file')
	workflow.connect(nki_nodes.Saver_nuisance_featM, 'out_file',nki_nodes.lifeSaver_nuisance_featM,'format_string')
	workflow.connect(nki_nodes.lifeSaver_nuisance_featM,'out_file',datasink,'@nuisance_featM')
#	workflow.connect(nki_nodes.nuisance_featM,'design_file',datasink,'@nuisance_featM')

#	workflow.connect(nki_nodes.nuisance_fgls, 'results_dir',nki_nodes.Saver_nuisance_fglsd,'in_file')
#	workflow.connect(nki_nodes.func_detrendc_r,'out_file',nki_nodes.Saver_nuisance_fglsd,'pp')
#	workflow.connect(nki_nodes.nuisance_fgls, 'results_dir',nki_nodes.lifeSaver_nuisance_fglsd,'in_file')
#	workflow.connect(nki_nodes.Saver_nuisance_fglsd, 'out_file',nki_nodes.lifeSaver_nuisance_fglsd,'format_string')
#	workflow.connect(nki_nodes.nuisance_fgls,'results_dir',datasink,'@nuisance_fgls')


	workflow.connect(nki_nodes.nuisance_fgls_r, 'out_file',nki_nodes.Saver_nuisance_fgls,'in_file')
	workflow.connect(nki_nodes.lifeSaverNuisancepp,'out_file',nki_nodes.Saver_nuisance_fgls,'pp')
	workflow.connect(nki_nodes.nuisance_fgls_r, 'out_file',nki_nodes.lifeSaver_nuisance_fgls,'in_file')
	workflow.connect(nki_nodes.Saver_nuisance_fgls, 'out_file',nki_nodes.lifeSaver_nuisance_fgls,'format_string')
	workflow.connect(nki_nodes.lifeSaver_nuisance_fgls,'out_file',datasink,'@nuisance_fgls')
#	workflow.connect(nki_nodes.nuisance_fgls,'residual4d',datasink,'@nuisance_fgls_')

	workflow.connect(nki_nodes.nuisance_stat_r, 'out_file',nki_nodes.Saver_nuisance_stat,'in_file')
	workflow.connect(nki_nodes.lifeSaverNuisancepp,'out_file',nki_nodes.Saver_nuisance_stat,'pp')
	workflow.connect(nki_nodes.nuisance_stat_r, 'out_file',nki_nodes.lifeSaver_nuisance_stat,'in_file')
	workflow.connect(nki_nodes.Saver_nuisance_stat, 'out_file',nki_nodes.lifeSaver_nuisance_stat,'format_string')
	workflow.connect(nki_nodes.lifeSaver_nuisance_stat,'out_file',datasink,'@nuisance_stat')
#	workflow.connect(nki_nodes.nuisance_stat,'out_file',datasink,'@nuisance_stat')

	
	workflow.connect(nki_nodes.nuisance_calc_r, 'out_file',nki_nodes.Saver_nuisance_calc,'in_file')
	workflow.connect(nki_nodes.lifeSaverNuisancepp,'out_file',nki_nodes.Saver_nuisance_calc,'pp')
	workflow.connect(nki_nodes.nuisance_calc_r, 'out_file',nki_nodes.lifeSaver_nuisance_calc,'in_file')
	workflow.connect(nki_nodes.Saver_nuisance_calc, 'out_file',nki_nodes.lifeSaver_nuisance_calc,'format_string')
	workflow.connect(nki_nodes.lifeSaver_nuisance_calc,'out_file',datasink,'@nuisance_calc')
#	workflow.connect(nki_nodes.nuisance_calc,'out_file',datasink,'@nuisance_calc')

	workflow.connect(nki_nodes.nuisance_warp_r, 'out_file',nki_nodes.Saver_nuisance_warp,'in_file')
	workflow.connect(nki_nodes.lifeSaverNuisancepp,'out_file',nki_nodes.Saver_nuisance_warp,'pp')
	workflow.connect(nki_nodes.nuisance_warp_r, 'out_file',nki_nodes.lifeSaver_nuisance_warp,'in_file')
	workflow.connect(nki_nodes.Saver_nuisance_warp, 'out_file',nki_nodes.lifeSaver_nuisance_warp,'format_string')
	workflow.connect(nki_nodes.lifeSaver_nuisance_warp,'out_file',datasink,'@nuisance_warp')
#	workflow.connect(nki_nodes.nuisance_warp,'out_file',datasink,'@nuisance_warp')

	workflow.connect(nki_nodes.nuisance_warp_1_r,'out_file',datasink,'@nuisance_warp_1')


def makeOutputConnections(analysisdirectory):



	datasink = getSinkOther()

	makeOutputConnectionsNuisance(datasink)


def rename_nuisance_outputs():

	global which_regression

	if (not which_regression == "compcor") and (not which_regression == "median_angle"):

		if (not which_regression == "no_global_signal"):
			workflow.connect(nki_nodes.nuisance_globalE, 'out_file',nki_nodes.nuisance_globalE_o,'in_file')
			workflow.connect(renamer,'nuisance_globalE_out_file',nki_nodes.nuisance_globalE_o,'name')
			workflow.connect(nki_nodes.nuisance_globalE, 'out_file',nki_nodes.nuisance_globalE_r,'in_file')
			workflow.connect(nki_nodes.nuisance_globalE_o,'out_file',nki_nodes.nuisance_globalE_r,'format_string')
			#workflow.connect(nki_nodes.nuisance_globalE,'out_file',datasink,'nuisance.@nuisance_globalE')
	
		workflow.connect(nki_nodes.nuisance_csf, 'out_file',nki_nodes.nuisance_csf_o,'in_file')
		workflow.connect(renamer,'nuisance_csf_out_file',nki_nodes.nuisance_csf_o,'name')
		workflow.connect(nki_nodes.nuisance_csf, 'out_file',nki_nodes.nuisance_csf_r,'in_file')
		workflow.connect(nki_nodes.nuisance_csf_o,'out_file',nki_nodes.nuisance_csf_r,'format_string')
		#workflow.connect(nki_nodes.nuisance_csf,'out_file',datasink,'nuisance.@nuisance_csf')
	
		workflow.connect(nki_nodes.nuisance_wm, 'out_file',nki_nodes.nuisance_wm_o,'in_file')
		workflow.connect(renamer,'nuisance_wm_out_file',nki_nodes.nuisance_wm_o,'name')
		workflow.connect(nki_nodes.nuisance_wm, 'out_file',nki_nodes.nuisance_wm_r,'in_file')
		workflow.connect(nki_nodes.nuisance_wm_o,'out_file',nki_nodes.nuisance_wm_r,'format_string')
		#workflow.connect(nki_nodes.nuisance_wm,'out_file',datasink,'nuisance.@nuisance_wm')

	
	workflow.connect(nki_nodes.nuisance_featM, 'design_file',nki_nodes.nuisance_featM_o,'in_file')
	workflow.connect(renamer,'nuisance_featM_design_file',nki_nodes.nuisance_featM_o,'name')
	workflow.connect(nki_nodes.nuisance_featM, 'design_file',nki_nodes.nuisance_featM_r,'in_file')
	workflow.connect(nki_nodes.nuisance_featM_o,'out_file',nki_nodes.nuisance_featM_r,'format_string')
	#workflow.connect(nki_nodes.nuisance_featM,'design_file',datasink,'nuisance.@nuisance_featM')
	
	#workflow.connect(nki_nodes.nuisance_fgls,'results_dir',datasink,'nuisance.@nuisance_fgls')
	
	workflow.connect(nki_nodes.nuisance_fgls, 'residual4d',nki_nodes.nuisance_fgls_o,'in_file')
	workflow.connect(renamer,'nuisance_fgls_residual4d',nki_nodes.nuisance_fgls_o,'name')
	workflow.connect(nki_nodes.nuisance_fgls, 'residual4d',nki_nodes.nuisance_fgls_r,'in_file')
	workflow.connect(nki_nodes.nuisance_fgls_o,'out_file',nki_nodes.nuisance_fgls_r,'format_string')
	#workflow.connect(nki_nodes.nuisance_fgls,'residual4d',datasink,'nuisance.@nuisance_fgls_')
	
	workflow.connect(nki_nodes.nuisance_stat, 'out_file',nki_nodes.nuisance_stat_o,'in_file')
	workflow.connect(renamer,'nuisance_stat_out_file',nki_nodes.nuisance_stat_o,'name')
	workflow.connect(nki_nodes.nuisance_stat, 'out_file',nki_nodes.nuisance_stat_r,'in_file')
	workflow.connect(nki_nodes.nuisance_stat_o,'out_file',nki_nodes.nuisance_stat_r,'format_string')
	#workflow.connect(nki_nodes.nuisance_stat,'out_file',datasink,'nuisance.@nuisance_stat')
	
	workflow.connect(nki_nodes.nuisance_calc, 'out_file',nki_nodes.nuisance_calc_o,'in_file')
	workflow.connect(renamer,'nuisance_calc_out_file',nki_nodes.nuisance_calc_o,'name')
	workflow.connect(nki_nodes.nuisance_calc, 'out_file',nki_nodes.nuisance_calc_r,'in_file')
	workflow.connect(nki_nodes.nuisance_calc_o,'out_file',nki_nodes.nuisance_calc_r,'format_string')
	#workflow.connect(nki_nodes.nuisance_calc,'out_file',datasink,'nuisance.@nuisance_calc')
	
	workflow.connect(nki_nodes.nuisance_warp, 'out_file',nki_nodes.nuisance_warp_o,'in_file')
	workflow.connect(renamer,'nuisance_warp_out_file',nki_nodes.nuisance_warp_o,'name')
	workflow.connect(nki_nodes.nuisance_warp, 'out_file',nki_nodes.nuisance_warp_r,'in_file')
	workflow.connect(nki_nodes.nuisance_warp_o,'out_file',nki_nodes.nuisance_warp_r,'format_string')
	#workflow.connect(nki_nodes.nuisance_warp,'out_file',datasink,'nuisance.@nuisance_warp')

	workflow.connect(nki_nodes.nuisance_warp_1, 'out_file',nki_nodes.nuisance_warp_1_o,'in_file')
	workflow.connect(renamer,'nuisance_warp_1_out_file',nki_nodes.nuisance_warp_1_o,'name')
	workflow.connect(nki_nodes.nuisance_warp_1, 'out_file',nki_nodes.nuisance_warp_1_r,'in_file')
	workflow.connect(nki_nodes.nuisance_warp_1_o,'out_file',nki_nodes.nuisance_warp_1_r,'format_string')


def renameOutputs():

	rename_nuisance_outputs()

def getNodeOutputNames():

	names = ['nuisance_globalE_out_file','nuisance_csf_out_file','nuisance_wm_out_file','nuisance_featM_design_file','nuisance_fgls_residual4d','nuisance_stat_out_file','nuisance_calc_out_file','nuisance_warp_out_file','nuisance_erosion_csf1_out_file','nuisance_erosion_wm1_out_file','nuisance_compcor_out_file','nuisance_MedianAngle_out_file','nuisance_warp_1_out_file']


	return names

def getRenamer():

	global renamer
	global rest_name
	global anat_name

	nodeOutputNames = getNodeOutputNames()
	renamer = pe.Node(interface=util.IdentityInterface(fields= nodeOutputNames), name="renamer")

	renamer.inputs.func_filter_out_file = rest_name + '_res_filt.nii.gz'
	#renamer for nuisance
	renamer.inputs.nuisance_globalE_out_file = 'global.1D'
	renamer.inputs.nuisance_csf_out_file = 'csf.1D'
	renamer.inputs.nuisance_wm_out_file = 'wm.1D'
	renamer.inputs.nuisance_featM_design_file = 'nuisance.mat'
	renamer.inputs.nuisance_fgls_residual4d = 'res4d.nii.gz'
	renamer.inputs.nuisance_stat_out_file = 'res4d_mean.nii.gz'
	renamer.inputs.nuisance_calc_out_file = rest_name + '_res.nii.gz'
	renamer.inputs.nuisance_warp_out_file = rest_name + '_res2standard.nii.gz'
	renamer.inputs.nuisance_warp_1_out_file = rest_name + '_mask2standard.nii.gz'
	renamer.inputs.nuisance_erosion_csf1_out_file = "csf_mask_erosion.nii.gz"
	renamer.inputs.nuisance_erosion_wm1_out_file = "wm_mask_erosion.nii.gz"
	renamer.inputs.nuisance_compcor_out_file = rest_name + '_pp_cc.nii.gz'
	renamer.inputs.nuisance_MedianAngle_out_file = rest_name + '_pp_mac.nii.gz'
	


def setParameters(sublist):

	global FWHM

	global subjectList
	subjectList = sublist
	nki_nodes.setParameters(FWHM)


def processS(sublist, analysisdirectory):

	from nipype.utils.config import config
	global infosource
	global datasource
	global workflow
	global strategy

	setParameters(sublist)
	
	numCores = int(sys.argv[1])
	
	getInfoSource(sublist, analysisdirectory)
	getRenamer()
	gatherData(sublist, analysisdirectory)
	wfname =  'nuisance_' + strategy
	workflow = pe.Workflow(name=wfname)
	workflow.base_dir = working_dir
	nuisance()
	renameOutputs()
	makeOutputConnections(analysisdirectory)

	workflow.run(plugin='MultiProc', plugin_args={'n_procs' : numCores})
	#workflow.write_graph()


def processSubjects(analysisdirectory,subject_list):

	#global subject_list

	
	try:
		fp = open(subject_list,'r')

		lines = fp.readlines()

		fp.close()
		line1 = []
		cnt = 0
		for line in lines:
			line = line.rstrip('\r\n')

			flag = 0

			if (flag == 0):
				line1.append(line)
				
				cnt += 1


		if not (cnt == 0):
			processS(line1,analysisdirectory)

	except:
		raise
def readSubjects():


	fp = open(batch_list,'r')

	lines = fp.readlines()
	line1 = []

	for line in lines:
		line = line.rstrip('\r\n')
		line1.append(line)
	
	fp.close()

	processes = [ Process(target=processSubjects, args=((l.split(' '))[0],(l.split(' '))[1])) for l in line1 ]
	processes[0].start()	

def readDirSetup():

	global FSLDIR
	global scripts_dir
	global working_dir
	global seed_list
	global batch_list
	global rest_name
	global rest_dir
	global anat_name
	global anat_dir
	global standard_res
	global standard_brain
	global logFile
	global prior_dir
	global nuisance_template
	global nuisance_template_cc
	global nuisance_template_noglobal
	global recon_subjects
	global HP
	global LP
	global hp
	global lp
	global which_regression
	global FWHM
	global dir_file
	global strategy
	parsermap = {}
	
	if(not(os.path.exists(sys.argv[2]))):
		sys.stderr.write('\ndefinition file expected: '+sys.argv[2])

	dir_file = os.path.abspath(sys.argv[2])
	strategy = sys.argv[3]
	parser = SafeConfigParser()
	parser.read('dir_setup.ini')

	for section in parser.sections():

		for variable , value in parser.items(section):
			print variable + ' ' + value
			parsermap[variable] = value
	
	FSLDIR = parsermap['fsldir']
	scripts_dir = parsermap['scripts_dir']
	prior = parsermap['prior_dir']
	working_dir = parsermap['working_dir']
	seed_list = scripts_dir + '/' + parsermap['seed_file']
	batch_list = scripts_dir + '/' + parsermap['batch_file']
	rest_name = parsermap['rest_name']
	rest_dir = parsermap['rest_dir']
	anat_name = parsermap['anat_name']
	anat_dir = parsermap['anat_dir']
	standard_res = parsermap['standard_res']
	nuisance_template = parsermap['nuisance_template']
	nuisance_template_cc = parsermap['nuisance_template_cc']
	nuisance_template_noglobal = parsermap['nuisance_template_noglobal']
	prior_dir = prior + '/%s' %(standard_res)
	standard_brain = FSLDIR + '/data/standard/MNI152_T1_%s_brain.nii.gz' %(standard_res)
	logFile = parsermap['logfile']
	recon_subjects = parsermap['recon_subjects']
	HP = parsermap['alff_hp']
	LP = parsermap['alff_lp']
	hp = parsermap['rest_hp']
	lp = parsermap['rest_lp']
	FWHM = float(parsermap['fwhm'])
	which_regression = parsermap['which_regression']

def main():

	if ( len(sys.argv) < 4 ):
		sys.stderr.write("./alff_pipeline.py <Number_of_cores> </path/to/dir_setup_X.ini> <strategy-label>")
	else:
	
		readDirSetup()
		readSubjects()

	

if __name__ == "__main__":

	sys.exit(main())
