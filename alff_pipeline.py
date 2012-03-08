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
reho_register = 'flirt'
infosource = None
datasource = None
datasource_warp = None
workflow = None
renamer = None
HP = None
LP = None
hp = None
lp = None
zeroSM = ''
#*******************************************************************************************



	


tolist = lambda x: [x]

def takemod(nvols):

	print '>>>>' + str(nvols)
	decisions = []
	for vol in nvols:
		mod = int (int(vol) % 2)

		if mod == 1:
			decisions.append([0])
			#return [0]
		else:
			decisions.append([1])
			#return [1]

	return decisions
	

def set_op_str(n2):

	strs = []
	for n in n2:	
		str = "-Tmean -mul %f" %(n)
		strs.append(str)
	return strs


def set_op1_str(nvols):

	strs = []
	for vol in nvols:
		str = '-Tmean -mul %d -div 2' %(int(vol))
		strs.append(str)
	
	return strs

def FALFF_NORM_REGISTER_WF():

	global workflow
	global scripts_dir
	global standard_brain
	global rest_name
	global FSLDIR
	global standard_res
	global infosource


	workflow.connect(nki_nodes.alff_sqrt,'out_file',nki_nodes.alff_falff,'in_file')
	workflow.connect(nki_nodes.NVOLS,('nvols',set_op1_str),nki_nodes.alff_falff,'op_string')
	workflow.connect(nki_nodes.alff_sum,'out_file',nki_nodes.alff_falff1,'in_file')
	workflow.connect(nki_nodes.alff_falff,'out_file',nki_nodes.alff_falff1,'operand_files')

	workflow.connect(nki_nodes.alff_sum,'out_file',nki_nodes.alff_normM,'in_file')
	workflow.connect(datasource,'rest_mask',nki_nodes.alff_normM,'mask_file')
	workflow.connect(nki_nodes.alff_sum,'out_file',nki_nodes.alff_normS,'in_file')
	workflow.connect(datasource,'rest_mask',nki_nodes.alff_normS,'mask_file')
	workflow.connect(nki_nodes.alff_falff1,'out_file',nki_nodes.alff_normM1,'in_file')
	workflow.connect(datasource,'rest_mask',nki_nodes.alff_normM1,'mask_file')
	workflow.connect(nki_nodes.alff_falff1,'out_file',nki_nodes.alff_normS1,'in_file')
	workflow.connect(datasource,'rest_mask',nki_nodes.alff_normS1,'mask_file')

	workflow.connect(nki_nodes.alff_normM,'out_stat', nki_nodes.alff_op_string,'mean')
	workflow.connect(nki_nodes.alff_normS,'out_stat', nki_nodes.alff_op_string,'std_dev')
	workflow.connect(nki_nodes.alff_op_string,'op_string', nki_nodes.alff_Z_alff, 'op_string')
	workflow.connect(nki_nodes.alff_sum,'out_file',nki_nodes.alff_Z_alff,'in_file')
	workflow.connect(datasource,'rest_mask',nki_nodes.alff_Z_alff,'operand_files')

	workflow.connect(nki_nodes.alff_normM1,'out_stat', nki_nodes.alff_op_string1,'mean')
	workflow.connect(nki_nodes.alff_normS1,'out_stat', nki_nodes.alff_op_string1,'std_dev')
	workflow.connect(nki_nodes.alff_op_string1,'op_string', nki_nodes.alff_Z_falff, 'op_string')
	workflow.connect(nki_nodes.alff_falff1,'out_file',nki_nodes.alff_Z_falff,'in_file')
	workflow.connect(datasource,'rest_mask',nki_nodes.alff_Z_falff,'operand_files')

	workflow.connect(infosource,'standard', nki_nodes.alff_warp_alff,'ref_file')
	workflow.connect(nki_nodes.alff_Z_alff,'out_file',nki_nodes.alff_warp_alff,'in_file')
	workflow.connect(datasource_warp, 'fieldcoeff_file',nki_nodes.alff_warp_alff,'field_file')
	workflow.connect(datasource_warp,'premat',nki_nodes.alff_warp_alff,'premat')

	workflow.connect(infosource,'standard', nki_nodes.alff_warp_falff,'ref_file')
	workflow.connect(nki_nodes.alff_Z_falff,'out_file',nki_nodes.alff_warp_falff,'in_file')
	workflow.connect(datasource_warp, 'fieldcoeff_file',nki_nodes.alff_warp_falff,'field_file')
	workflow.connect(datasource_warp,'premat',nki_nodes.alff_warp_falff,'premat')

	workflow.connect(nki_nodes.alff_warp_alff,'out_file', nki_nodes.alff_smooth,'in_file')
	workflow.connect(datasource,'rest_mask2standard', nki_nodes.alff_smooth,'operand_files')

	workflow.connect(nki_nodes.alff_warp_falff,'out_file', nki_nodes.falff_smooth,'in_file')
	workflow.connect(datasource,'rest_mask2standard', nki_nodes.falff_smooth,'operand_files')




def falff():

	global workflow
	global scripts_dir
	global standard_brain
	global rest_name
	global FSLDIR
	global standard_res
	global infosource
	global strategy

	
	workflow.connect(datasource,'rest_res',nki_nodes.lifeSaverALFF,'in_file')
	workflow.connect(datasource,'rest_res',nki_nodes.adjustPath,'rest_path')
	workflow.connect(infosource,'strategy',nki_nodes.adjustPath,'strategy')
	workflow.connect(infosource,'sublist',nki_nodes.adjustPath,'sublist')
	workflow.connect(nki_nodes.adjustPath,'result_path',nki_nodes.lifeSaverALFF,'format_string')

	workflow.connect(nki_nodes.lifeSaverALFF, 'out_file',nki_nodes.TR,'in_files')
	workflow.connect(nki_nodes.lifeSaverALFF, 'out_file',nki_nodes.NVOLS,'in_files')
	workflow.connect(nki_nodes.lifeSaverALFF,'out_file',nki_nodes.alff_mean,'in_file')
	workflow.connect(nki_nodes.lifeSaverALFF,'out_file',nki_nodes.alff_roi,'in_file')
	workflow.connect(nki_nodes.NVOLS,'nvols',nki_nodes.alff_roi,'t_size')
	workflow.connect(nki_nodes.lifeSaverALFF,'out_file',nki_nodes.alff_cp1,'in_file')

	workflow.connect(nki_nodes.alff_roi,'roi_file',nki_nodes.alff_concatnode,'in1')
	workflow.connect(nki_nodes.alff_cp1,'out_file',nki_nodes.alff_concatnode,'in2')
	workflow.connect(nki_nodes.alff_concatnode, 'out', nki_nodes.alff_selectnode, 'inlist')
	workflow.connect(nki_nodes.NVOLS,('nvols',takemod), nki_nodes.alff_selectnode, 'index')
	workflow.connect(nki_nodes.alff_selectnode,'out', nki_nodes.alff_pspec,'in_file')
	workflow.connect(nki_nodes.alff_pspec,'out_file',nki_nodes.alff_sqrt,'in_file')

	workflow.connect(nki_nodes.NVOLS,'nvols',nki_nodes.calcN1,'nvols')
	workflow.connect(nki_nodes.TR,'TR',nki_nodes.calcN1,'TR')
	workflow.connect(infosource,'HP',nki_nodes.calcN1,'HP')

	workflow.connect(nki_nodes.NVOLS,'nvols',nki_nodes.calcN2,'nvols')
	workflow.connect(nki_nodes.TR,'TR',nki_nodes.calcN2,'TR')
	workflow.connect(infosource,'LP',nki_nodes.calcN2,'LP')
	workflow.connect(infosource,'HP',nki_nodes.calcN2,'HP')

	workflow.connect(nki_nodes.alff_sqrt,'out_file',nki_nodes.alff_roi1,'in_file')
	workflow.connect(nki_nodes.calcN1,'n1',nki_nodes.alff_roi1,'t_min')
	workflow.connect(nki_nodes.calcN2,'n2',nki_nodes.alff_roi1,'t_size')
	workflow.connect(nki_nodes.alff_roi1,'roi_file',nki_nodes.alff_sum,'in_file')
	workflow.connect(nki_nodes.calcN2,('n2',set_op_str),nki_nodes.alff_sum,'op_string')

	FALFF_NORM_REGISTER_WF()

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


	datasource = pe.Node( interface = nio.DataGrabber(infields=['subject_id'], outfields = [ 'rest_res', 'rest_mask', 'rest_mask2standard' ]) , name= 'datasource')
	datasource.inputs.base_directory = analysisdirectory
	#datasource.inputs.template = '%s/*/%s.nii.gz'
	datasource.inputs.template = '%s/*/*/%s.nii.gz'
	datasource.inputs.template_args = dict( rest_res = [['subject_id',rest_name +'_res' ]] , rest_mask = [['subject_id',rest_name + '_mask']], rest_mask2standard = [['subject_id',rest_name + '_mask2standard']])
	datasource.iterables = ('subject_id', sublist)

	datasource_warp = pe.Node( interface = nio.DataGrabber(infields=['subject_id'], outfields = [ 'premat', 'fieldcoeff_file' ]) , name= 'datasource_warp')
	datasource_warp.inputs.base_directory = analysisdirectory
	#datasource_warp.inputs.template = '%s/*/%s.nii.gz'
	datasource_warp.inputs.template = '%s/*/*/*/%s'
	datasource_warp.inputs.template_args = dict( premat = [['subject_id','example_func2highres.mat' ]] , fieldcoeff_file = [['subject_id','highres2standard_warp.nii.gz']])
	datasource_warp.iterables = ('subject_id', sublist)



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
	global strategy
	global subjectList

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



def getSink(analysisdirectory):

	global anat_name
	global rest_name
	sink = pe.Node(nio.DataSink(), name = 'sinker-anatomical')
	sink.inputs.base_directory = os.path.abspath(analysisdirectory)
	sink.inputs.regexp_substitutions = [(r"_subject_id_(\w|\d)+", ''),(r"%s.nii.gz/" %(anat_name),''),(r"_anat_(\w|\d)+",'')]

	return sink

def getSinkOther():

	datasink = None
	
	datasink = pe.Node(nio.DataSink(), name = 'sinker-other')
	datasink.inputs.base_directory = os.path.abspath('/')
	datasink.inputs.container = os.path.abspath('/')
	datasink.inputs.regexp_substitutions = [(r"-",'/'),(r"(.)+_func_(\w|\d)+",''),(r"(.)+_reg_(\w|\d)+",''),(r"(.)+_seg_(\w|\d)+",''),(r"(.)+_nuisance_(\w|\d)+",''),(r"(.)+_alff_(\w|\d)+/",''),(r"(.)+_falff(\w|\d)+/",''),(r"(.)+_RSFC_(\w|\d)+",''),(r"(.)+/_seed_(.)+seeds\.\./",'')]


	return datasink


def adjustALFF(rest_path):

	import re
	import sys
	result_path = []

	if(isinstance(rest_path,list)):
		for r_p in rest_path:

			split_path = r_p.split('-')

			path = ""
			for index in range(0,len(split_path) - 1):
				path += split_path[index] + "-"

			path += "ALFF" + "-"
			path += split_path[len(split_path) - 1]

			print "\nprocessed path %s\n" %(path)
			result_path.append(path)
		
	else:
		split_path = rest_path.split('-')

		path = ""
		for index in range(0,len(split_path) - 1):
			path += split_path[index] + "-"

		path += "ALFF" + "-"
		path += split_path[len(split_path) - 1]

		print "\nprocessed path %s\n" %(path)
		result_path.append(path)
		

	return result_path


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


def makeOutputConnectionsALFF(datasink):

#	## Connect falff nodes to datasink	
	
	
	workflow.connect(nki_nodes.alff_mean_r, 'out_file',nki_nodes.lifeSaver_alff_mean,'in_file')
	workflow.connect(nki_nodes.alff_mean_r, ('out_file',adjustALFF),nki_nodes.lifeSaver_alff_mean,'format_string')
	workflow.connect(nki_nodes.lifeSaver_alff_mean,'out_file',datasink,'@alff_mean')
	
	workflow.connect(nki_nodes.alff_pspec_r, 'out_file',nki_nodes.lifeSaver_alff_pspec,'in_file')
	workflow.connect(nki_nodes.alff_pspec_r, ('out_file',adjustALFF),nki_nodes.lifeSaver_alff_pspec,'format_string')
	workflow.connect(nki_nodes.lifeSaver_alff_pspec,'out_file',datasink,'@alff_pspec')
	
	workflow.connect(nki_nodes.alff_sum_r, 'out_file',nki_nodes.lifeSaver_alff_sum,'in_file')
	workflow.connect(nki_nodes.alff_sum_r, ('out_file',adjustALFF),nki_nodes.lifeSaver_alff_sum,'format_string')
	workflow.connect(nki_nodes.lifeSaver_alff_sum,'out_file',datasink,'@alff_sum')
	
	workflow.connect(nki_nodes.alff_falff1_r, 'out_file',nki_nodes.lifeSaver_alff_falff1,'in_file')
	workflow.connect(nki_nodes.alff_falff1_r, ('out_file',adjustALFF),nki_nodes.lifeSaver_alff_falff1,'format_string')
	workflow.connect(nki_nodes.lifeSaver_alff_falff1,'out_file',datasink,'@alff_falff1')
	
	workflow.connect(nki_nodes.alff_Z_falff_r, 'out_file',nki_nodes.lifeSaver_alff_Z_falff,'in_file')
	workflow.connect(nki_nodes.alff_Z_falff_r, ('out_file',adjustALFF),nki_nodes.lifeSaver_alff_Z_falff,'format_string')
	workflow.connect(nki_nodes.lifeSaver_alff_Z_falff,'out_file',datasink,'@alff_Z_falff')
	
	workflow.connect(nki_nodes.alff_Z_alff_r, 'out_file',nki_nodes.lifeSaver_alff_Z_alff,'in_file')
	workflow.connect(nki_nodes.alff_Z_alff_r, ('out_file',adjustALFF),nki_nodes.lifeSaver_alff_Z_alff,'format_string')
	workflow.connect(nki_nodes.lifeSaver_alff_Z_alff,'out_file',datasink,'@alff_Z_alff')
	
	workflow.connect(nki_nodes.alff_warp_alff_r, 'out_file',nki_nodes.lifeSaver_alff_warp_alff,'in_file')
	workflow.connect(nki_nodes.alff_warp_alff_r, ('out_file',adjustALFF),nki_nodes.lifeSaver_alff_warp_alff,'format_string')
	workflow.connect(nki_nodes.lifeSaver_alff_warp_alff,'out_file',datasink,'@alff_warp_alff')
	
	workflow.connect(nki_nodes.alff_warp_falff_r, 'out_file',nki_nodes.lifeSaver_alff_warp_falff,'in_file')
	workflow.connect(nki_nodes.alff_warp_falff_r, ('out_file',adjustALFF),nki_nodes.lifeSaver_alff_warp_falff,'format_string')
	workflow.connect(nki_nodes.lifeSaver_alff_warp_falff,'out_file',datasink,'@alff_warp_falff')

	workflow.connect(nki_nodes.alff_smooth_r, 'out_file',nki_nodes.lifeSaver_alff_smooth,'in_file')
	workflow.connect(nki_nodes.alff_smooth_r, ('out_file',adjustALFF),nki_nodes.lifeSaver_alff_smooth,'format_string')
	workflow.connect(nki_nodes.lifeSaver_alff_smooth,'out_file',datasink,'@alff_smooth')

	workflow.connect(nki_nodes.falff_smooth_r, 'out_file',nki_nodes.lifeSaver_falff_smooth,'in_file')
	workflow.connect(nki_nodes.falff_smooth_r, ('out_file',adjustALFF),nki_nodes.lifeSaver_falff_smooth,'format_string')
	workflow.connect(nki_nodes.lifeSaver_falff_smooth,'out_file',datasink,'@falff_smooth')



def makeOutputConnections(analysisdirectory):



	datasink = getSinkOther()
	makeOutputConnectionsALFF(datasink)



def rename_alff_outputs():

	


	workflow.connect(nki_nodes.alff_mean, 'out_file',nki_nodes.alff_mean_o,'in_file')
	workflow.connect(renamer,'alff_mean_out_file',nki_nodes.alff_mean_o,'name')
	workflow.connect(nki_nodes.alff_mean, 'out_file',nki_nodes.alff_mean_r,'in_file')
	workflow.connect(nki_nodes.alff_mean_o,'out_file',nki_nodes.alff_mean_r,'format_string')
#	workflow.connect(nki_nodes.alff_mean,'out_file',datasink,'ALFF.@alff_mean')

	workflow.connect(nki_nodes.alff_pspec, 'out_file',nki_nodes.alff_pspec_o,'in_file')
	workflow.connect(renamer,'alff_pspec_out_file',nki_nodes.alff_pspec_o,'name')
	workflow.connect(nki_nodes.alff_pspec, 'out_file',nki_nodes.alff_pspec_r,'in_file')
	workflow.connect(nki_nodes.alff_pspec_o,'out_file',nki_nodes.alff_pspec_r,'format_string')
#	workflow.connect(nki_nodes.alff_pspec,'out_file',datasink,'ALFF.@alff_pspec')

	workflow.connect(nki_nodes.alff_sum, 'out_file',nki_nodes.alff_sum_o,'in_file')
	workflow.connect(renamer,'alff_sum_out_file',nki_nodes.alff_sum_o,'name')
	workflow.connect(nki_nodes.alff_sum, 'out_file',nki_nodes.alff_sum_r,'in_file')
	workflow.connect(nki_nodes.alff_sum_o,'out_file',nki_nodes.alff_sum_r,'format_string')
#	workflow.connect(nki_nodes.alff_sum,'out_file',datasink,'ALFF.@alff_sum')

	workflow.connect(nki_nodes.alff_falff1, 'out_file',nki_nodes.alff_falff1_o,'in_file')
	workflow.connect(renamer,'alff_falff1_out_file',nki_nodes.alff_falff1_o,'name')
	workflow.connect(nki_nodes.alff_falff1, 'out_file',nki_nodes.alff_falff1_r,'in_file')
	workflow.connect(nki_nodes.alff_falff1_o,'out_file',nki_nodes.alff_falff1_r,'format_string')
#	workflow.connect(nki_nodes.alff_falff1,'out_file',datasink,'ALFF.@alff_falff1')

	workflow.connect(nki_nodes.alff_Z_falff, 'out_file',nki_nodes.alff_Z_falff_o,'in_file')
	workflow.connect(renamer,'alff_Z_falff_out_file',nki_nodes.alff_Z_falff_o,'name')
	workflow.connect(nki_nodes.alff_Z_falff, 'out_file',nki_nodes.alff_Z_falff_r,'in_file')
	workflow.connect(nki_nodes.alff_Z_falff_o,'out_file',nki_nodes.alff_Z_falff_r,'format_string')

	workflow.connect(nki_nodes.alff_Z_alff, 'out_file',nki_nodes.alff_Z_alff_o,'in_file')
	workflow.connect(renamer,'alff_Z_alff_out_file',nki_nodes.alff_Z_alff_o,'name')
	workflow.connect(nki_nodes.alff_Z_falff, 'out_file',nki_nodes.alff_Z_alff_r,'in_file')
	workflow.connect(nki_nodes.alff_Z_alff_o,'out_file',nki_nodes.alff_Z_alff_r,'format_string')

	workflow.connect(nki_nodes.alff_warp_alff, 'out_file',nki_nodes.alff_warp_alff_o,'in_file')
	workflow.connect(renamer,'alff_warp_alff_out_file',nki_nodes.alff_warp_alff_o,'name')
	workflow.connect(nki_nodes.alff_warp_alff, 'out_file',nki_nodes.alff_warp_alff_r,'in_file')
	workflow.connect(nki_nodes.alff_warp_alff_o,'out_file',nki_nodes.alff_warp_alff_r,'format_string')

	workflow.connect(nki_nodes.alff_warp_falff, 'out_file',nki_nodes.alff_warp_falff_o,'in_file')
	workflow.connect(renamer,'alff_warp_falff_out_file',nki_nodes.alff_warp_falff_o,'name')
	workflow.connect(nki_nodes.alff_warp_falff, 'out_file',nki_nodes.alff_warp_falff_r,'in_file')
	workflow.connect(nki_nodes.alff_warp_falff_o,'out_file',nki_nodes.alff_warp_falff_r,'format_string')

	workflow.connect(nki_nodes.alff_smooth, 'out_file',nki_nodes.alff_smooth_o,'in_file')
	workflow.connect(renamer,'alff_smooth_out_file',nki_nodes.alff_smooth_o,'name')
	workflow.connect(nki_nodes.alff_smooth, 'out_file',nki_nodes.alff_smooth_r,'in_file')
	workflow.connect(nki_nodes.alff_smooth_o,'out_file',nki_nodes.alff_smooth_r,'format_string')

	workflow.connect(nki_nodes.falff_smooth, 'out_file',nki_nodes.falff_smooth_o,'in_file')
	workflow.connect(renamer,'falff_smooth_out_file',nki_nodes.falff_smooth_o,'name')
	workflow.connect(nki_nodes.falff_smooth, 'out_file',nki_nodes.falff_smooth_r,'in_file')
	workflow.connect(nki_nodes.falff_smooth_o,'out_file',nki_nodes.falff_smooth_r,'format_string')



def renameOutputs():

	rename_alff_outputs()

def getNodeOutputNames():

	names = ['alff_smooth_out_file','falff_smooth_out_file','reg_warp_out_file','reg_fnt_jacobian_file','reg_fnt_fieldcoeff_file','reg_fnt_warped_file','reg_flirt1_out_file','reg_flirt1_out_matrix_file','reg_xfm2_out_file','anat_refit_out_file','anat_reorient_out_file','anat_skullstrip_out_file','anat_calc_out_file','func_calc_out_file','func_refit_out_file','func_reorient_out_file','func_tstat_out_file','func_volreg_oned_file','func_volreg_out_file','func_automask_out_file','func_calcR_out_file','func_calcI_out_file','func_despike_out_file','func_smooth_out_file','func_scale_out_file','func_filter_out_file','func_detrenda_out_file','func_detrendb_out_file','func_detrendc_out_file','func_mask_out_file','reg_flirt_out_file','reg_flirt_out_matrix_file','reg_xfm1_out_file','reg_xfm3_out_file','reg_flirt2_out_file','reg_xfm4_out_file','seg_flirt_out_file','seg_smooth_out_file','seg_flirt1_out_file','seg_smooth1_out_file','seg_flirt2_out_file','seg_thresh_out_file','seg_mask_out_file','seg_copy_out_file','seg_flirt3_out_file','seg_smooth2_out_file','seg_flirt4_out_file','seg_prior1_out_file','seg_flirt5_out_file','seg_thresh1_out_file','seg_mask1_out_file','nuisance_globalE_out_file','nuisance_csf_out_file','nuisance_wm_out_file','nuisance_featM_design_file','nuisance_fgls_residual4d','nuisance_stat_out_file','nuisance_calc_out_file','nuisance_warp_out_file','alff_cp_out_file','alff_mean_out_file','alff_pspec_out_file','alff_sum_out_file','alff_falff1_out_file','alff_Z_falff_out_file','alff_Z_alff_out_file','alff_warp_alff_out_file','alff_warp_falff_out_file','corrs','z_trans','register','nuisance_erosion_csf1_out_file','nuisance_erosion_wm1_out_file','nuisance_compcor_out_file','nuisance_MedianAngle_out_file']


	return names

def getRenamer():

	global renamer
	global rest_name
	global anat_name
	global FWHM
	nodeOutputNames = getNodeOutputNames()
	renamer = pe.Node(interface=util.IdentityInterface(fields= nodeOutputNames), name="renamer")

	
	#renamer for alff/falff
	renamer.inputs.alff_cp_out_file = rest_name + '_alff_pp.nii.gz'
	renamer.inputs.alff_mean_out_file = 'mean_'+ rest_name + '_alff_pp.nii.gz'
	renamer.inputs.alff_pspec_out_file = 'power_spectrum_distribution.nii.gz'
	renamer.inputs.alff_sum_out_file = 'ALFF.nii.gz'
	renamer.inputs.alff_falff1_out_file = 'fALFF.nii.gz'
	renamer.inputs.alff_Z_alff_out_file = 'ALFF_Z.nii.gz'
	renamer.inputs.alff_Z_falff_out_file = 'fALFF_Z.nii.gz'
	renamer.inputs.alff_warp_alff_out_file = 'ALFF_Z_2standard.nii.gz'
	renamer.inputs.alff_warp_falff_out_file = 'fALFF_Z_2standard.nii.gz'
	renamer.inputs.alff_smooth_out_file = 'ALFF_Z_2standard_FWHM' + str(FWHM) + '.nii.gz'
	renamer.inputs.falff_smooth_out_file = 'fALFF_Z_2standard_FWHM' + str(FWHM) + '.nii.gz'

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

	numCores = int(sys.argv[1])
	
	setParameters(sublist)	
	getInfoSource(sublist, analysisdirectory)
	getRenamer()
	gatherData(sublist, analysisdirectory)
	wfname =  'alff_'+strategy
	workflow = pe.Workflow(name=wfname)
	workflow.base_dir = working_dir
	falff()
	renameOutputs()
	makeOutputConnections(analysisdirectory)

	workflow.run(plugin='MultiProc', plugin_args={'n_procs' : numCores})
	#workflow.run()
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
			print 'inside processSubjects ; before check processed subject ',line
			#flag = checkProcessedSubject(analysisdirectory, line)

			print 'returned %d' %(flag)
			if (flag == 0):
				line1.append(line)
				
				cnt += 1


		if not (cnt == 0):
			#process =  pool.map(processS,line1)
			processS(line1,analysisdirectory)
			#pool.close()
			#pool.join()

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
	global zeroSM
	global reho_register
	global FWHM
	global dir_file
	global strategy

	parsermap = {}

	if(not(os.path.exists(sys.argv[2]))):
		sys.stderr.write('\ndefinition file expected: '+sys.argv[2])

	dir_file = os.path.abspath(sys.argv[2])
	strategy = sys.argv[3]



	parser = SafeConfigParser()
	parser.read(dir_file)

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
	anat_name = parsermap['anat_name']
	standard_res = parsermap['standard_res']
	prior_dir = prior + '/%s' %(standard_res)
	standard_brain = FSLDIR + '/data/standard/MNI152_T1_%s_brain.nii.gz' %(standard_res)
	logFile = parsermap['logfile']
	FWHM = float(parsermap['fwhm'])
	HP = parsermap['alff_hp']
	LP = parsermap['alff_lp']
	hp = parsermap['rest_hp']
	lp = parsermap['rest_lp']


def main():

	if ( len(sys.argv) < 4 ):
		sys.stderr.write("./alff_pipeline.py <Number_of_cores> </path/to/dir_setup_X.ini> <strategy-label>")
	else:
	
		readDirSetup()
		readSubjects()

	

if __name__ == "__main__":

	sys.exit(main())
