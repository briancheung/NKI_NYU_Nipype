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
subjectList = None
dir_file = None
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



	
def last_vol(vols):

	v = []
	for vol in vols:
		v.append(int(vol) - 1)

	return v




tolist = lambda x: [x]





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


def RSFC():
	
	global workflow
	global scripts_dir
	global standard_brain
	global rest_name
	global FSLDIR
	global standard_res
	global infosource
	#global datasource_seeds

	workflow.connect(datasource,'ref',nki_nodes.RSFC_warp,'ref_file')
	workflow.connect(datasource_warp,'warp',nki_nodes.RSFC_warp,'field_file')
	workflow.connect(datasource_warp,'postmat',nki_nodes.RSFC_warp,'postmat')
#	workflow.connect(nki_nodes.RSFC_warp,'out_file',nki_nodes.lifeSaverRSFC,'in_file')
#	workflow.connect(datasource,'rest_res2standard',nki_nodes.adjustPath,'rest_path')

#	workflow.connect(datasource,'rest_res2standard',nki_nodes.lifeSaverRSFC,'in_file')
#	workflow.connect(datasource,'rest_res2standard',nki_nodes.adjustPath,'rest_path')
#	workflow.connect(infosource,'strategy',nki_nodes.adjustPath,'strategy')
#	workflow.connect(infosource,'sublist',nki_nodes.adjustPath,'sublist')
#	workflow.connect(nki_nodes.adjustPath,'result_path',nki_nodes.lifeSaverRSFC,'format_string')
#	workflow.connect(nki_nodes.lifeSaverRSFC,'out_file',nki_nodes.RSFC_time_series,'in_file')
#	workflow.connect(nki_nodes.RSFC_time_series,'stats',nki_nodes.RSFC_printToFile,'time_series')
#	workflow.connect(nki_nodes.RSFC_printToFile,'ts_oneD',nki_nodes.RSFC_corr,'ideal_file')
#	workflow.connect(datasource,'rest_res_filt',nki_nodes.lifeSaverRSFC1,'in_file')
#	workflow.connect(datasource,'rest_res_filt',nki_nodes.adjustPath1,'rest_path')
#	workflow.connect(infosource,'strategy',nki_nodes.adjustPath1,'strategy')
#	workflow.connect(infosource,'sublist',nki_nodes.adjustPath1,'sublist')
#	workflow.connect(nki_nodes.adjustPath1,'result_path',nki_nodes.lifeSaverRSFC1,'format_string')
#	workflow.connect(nki_nodes.lifeSaverRSFC1,'out_file',nki_nodes.RSFC_corr,'in_file')
#	workflow.connect(nki_nodes.RSFC_corr,'out_file',nki_nodes.RSFC_z_trans,'infile_a')
#	workflow.connect(nki_nodes.RSFC_z_trans,'out_file',nki_nodes.RSFC_register,'in_file')
#	workflow.connect(infosource,'standard',nki_nodes.RSFC_register,'ref_file')
#	workflow.connect(datasource_warp, 'fieldcoeff_file',nki_nodes.RSFC_register,'field_file')
#	workflow.connect(datasource_warp,'premat',nki_nodes.RSFC_register,'premat')
#	workflow.connect(nki_nodes.RSFC_register,'out_file',nki_nodes.RSFC_smooth,'in_file')
#	workflow.connect(datasource,'rest_mask2standard',nki_nodes.RSFC_smooth,'operand_files')


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


	datasource = pe.Node( interface = nio.DataGrabber(infields=['subject_id'], outfields = [ 'rest_res2standard', 'rest_res_filt', 'rest_mask2standard', 'ref' ]) , name= 'datasource')
	datasource.inputs.base_directory = analysisdirectory
	#datasource.inputs.template = '%s/*/%s.nii.gz'
	datasource.inputs.template = '%s/*/*/%s.nii.gz'
	datasource.inputs.template_args = dict( rest_res2standard = [['subject_id',rest_name + '_res2standard']] , rest_res_filt = [['subject_id',rest_name + '_res_filt']], rest_mask2standard = [['subject_id',rest_name + '_mask2standard']] , ref = [['subject_id','example_func' ]])
	datasource.iterables = ('subject_id', sublist)


	datasource_warp = pe.Node( interface = nio.DataGrabber(infields=['subject_id'], outfields = [ 'premat', 'fieldcoeff_file', 'warp', 'postmat']) , name= 'datasource_warp')
	datasource_warp.inputs.base_directory = analysisdirectory
	#datasource_warp.inputs.template = '%s/*/%s.nii.gz'
	datasource_warp.inputs.template = '%s/*/*/*/%s'
	datasource_warp.inputs.template_args = dict( premat = [['subject_id','example_func2highres.mat' ]] , fieldcoeff_file = [['subject_id','highres2standard_warp.nii.gz']],warp = [['subject_id','stand2highres_warp.nii.gz' ]], postmat = [['subject_id','highres2example_func.mat' ]])
	datasource_warp.iterables = ('subject_id', sublist)


	nki_nodes.setSeedList(seed_list)




def getFields():

	fields = ['strategy','sublist','brain_symmetric','symm_standard','standard_res_brain','standard','standard_brain_mask_dil','twomm_brain_mask_dil','config_file_twomm','config_file','PRIOR_CSF','PRIOR_WHITE']
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
	infosource.inputs.sublist = subjectList
	infosource.inputs.strategy = strategy

	infosource.inputs.brain_symmetric = os.path.abspath(FSLDIR + '/data/standard/MNI152_T1_2mm_brain_symmetric.nii.gz')
	infosource.inputs.symm_standard = os.path.abspath(FSLDIR + '/data/standard/MNI152_T1_2mm_symmetric.nii.gz')




def getSinkRSFC():

	datasinkRSFC = None
	
	datasinkRSFC = pe.Node(nio.DataSink(), name = 'sinker-RSFC')
	datasinkRSFC.inputs.base_directory = os.path.abspath('/')
	datasinkRSFC.inputs.container = os.path.abspath('/')
	datasinkRSFC.inputs.regexp_substitutions = [(r"-",'/'),(r"(.)+_func_(\w|\d)+",''),(r"(.)+_reg_(\w|\d)+",''),(r"(.)+_seg_(\w|\d)+",''),(r"(.)+_nuisance_(\w|\d)+",''),(r"(.)+_alff_(\w|\d)+/",''),(r"(.)+_RSFC_(\w|\d)+",''),(r"(.)+/_seed_(.)+seeds\.\./",'')]


	return datasinkRSFC


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



def makeOutputConnectionsRSFC(datasink):
	## connect RSFC nodes to datasink

	workflow.connect(nki_nodes.RSFC_printToFile, 'ts_oneD',nki_nodes.Saver_RSFC_printToFile,'in_file')
	workflow.connect(nki_nodes.lifeSaverRSFC,'out_file',nki_nodes.Saver_RSFC_printToFile,'pp')
	workflow.connect(nki_nodes.RSFC_printToFile, 'ts_oneD',nki_nodes.lifeSaver_RSFC_printToFile,'in_file')
	workflow.connect(nki_nodes.Saver_RSFC_printToFile, 'out_file',nki_nodes.lifeSaver_RSFC_printToFile,'format_string')
	workflow.connect(nki_nodes.lifeSaver_RSFC_printToFile,'out_file',datasink,'@RSFC_printToFile')
#	workflow.connect(nki_nodes.RSFC_printToFile,'ts_oneD',datasink,'RSFC.@ts_oneD')

	workflow.connect(nki_nodes.RSFC_corr_r, 'out_file',nki_nodes.Saver_RSFC_corr,'in_file')
	workflow.connect(nki_nodes.lifeSaverRSFC,'out_file',nki_nodes.Saver_RSFC_corr,'pp')
	workflow.connect(nki_nodes.RSFC_corr_r, 'out_file',nki_nodes.lifeSaver_RSFC_corr,'in_file')
	workflow.connect(nki_nodes.Saver_RSFC_corr, 'out_file',nki_nodes.lifeSaver_RSFC_corr,'format_string')
	workflow.connect(nki_nodes.lifeSaver_RSFC_corr,'out_file',datasink,'@RSFC_corr')
	#workflow.connect(nki_nodes.RSFC_corr,'out_file',datasink,'RSFC.@RSFC_corr')
	
	workflow.connect(nki_nodes.RSFC_z_trans_r, 'out_file',nki_nodes.Saver_RSFC_z_trans,'in_file')
	workflow.connect(nki_nodes.lifeSaverRSFC,'out_file',nki_nodes.Saver_RSFC_z_trans,'pp')
	workflow.connect(nki_nodes.RSFC_z_trans_r, 'out_file',nki_nodes.lifeSaver_RSFC_z_trans,'in_file')
	workflow.connect(nki_nodes.Saver_RSFC_z_trans, 'out_file',nki_nodes.lifeSaver_RSFC_z_trans,'format_string')
	workflow.connect(nki_nodes.lifeSaver_RSFC_z_trans,'out_file',datasink,'@RSFC_z_trans')
	#workflow.connect(nki_nodes.RSFC_z_trans,'out_file',datasink,'RSFC.@RSFC_z_trans')
	
	workflow.connect(nki_nodes.RSFC_register_r, 'out_file',nki_nodes.Saver_RSFC_register,'in_file')
	workflow.connect(nki_nodes.lifeSaverRSFC,'out_file',nki_nodes.Saver_RSFC_register,'pp')
	workflow.connect(nki_nodes.RSFC_register_r, 'out_file',nki_nodes.lifeSaver_RSFC_register,'in_file')
	workflow.connect(nki_nodes.Saver_RSFC_register, 'out_file',nki_nodes.lifeSaver_RSFC_register,'format_string')
	workflow.connect(nki_nodes.lifeSaver_RSFC_register,'out_file',datasink,'@RSFC_register')
	#workflow.connect(nki_nodes.RSFC_register,'out_file',datasink,'RSFC.@RSFC_register')

	workflow.connect(nki_nodes.RSFC_smooth_r, 'out_file',nki_nodes.Saver_RSFC_smooth,'in_file')
	workflow.connect(nki_nodes.lifeSaverRSFC,'out_file',nki_nodes.Saver_RSFC_smooth,'pp')
	workflow.connect(nki_nodes.RSFC_smooth_r, 'out_file',nki_nodes.lifeSaver_RSFC_smooth,'in_file')
	workflow.connect(nki_nodes.Saver_RSFC_smooth, 'out_file',nki_nodes.lifeSaver_RSFC_smooth,'format_string')
	workflow.connect(nki_nodes.lifeSaver_RSFC_smooth,'out_file',datasink,'@RSFC_smooth')

def makeOutputConnections(analysisdirectory):



	datasinkRSFC = getSinkRSFC()
	makeOutputConnectionsRSFC(datasinkRSFC)




def rename_RSFC_outputs():

	
	workflow.connect(nki_nodes.RSFC_corr, 'out_file',nki_nodes.RSFC_corr_o,'in_file')
	workflow.connect(nki_nodes.RSFC_printToFile,'ts_oneD',nki_nodes.RSFC_corr_o,'name')
	workflow.connect(renamer,'corrs',nki_nodes.RSFC_corr_o,'whichNode')
	workflow.connect(nki_nodes.RSFC_corr, 'out_file',nki_nodes.RSFC_corr_r,'in_file')
	workflow.connect(nki_nodes.RSFC_corr_o,'out_file',nki_nodes.RSFC_corr_r,'format_string')
	
	workflow.connect(nki_nodes.RSFC_z_trans, 'out_file',nki_nodes.RSFC_z_trans_o,'in_file')
	workflow.connect(nki_nodes.RSFC_printToFile,'ts_oneD',nki_nodes.RSFC_z_trans_o,'name')
	workflow.connect(renamer,'z_trans',nki_nodes.RSFC_z_trans_o,'whichNode')
	workflow.connect(nki_nodes.RSFC_z_trans, 'out_file',nki_nodes.RSFC_z_trans_r,'in_file')
	workflow.connect(nki_nodes.RSFC_z_trans_o,'out_file',nki_nodes.RSFC_z_trans_r,'format_string')
	
	workflow.connect(nki_nodes.RSFC_register, 'out_file',nki_nodes.RSFC_register_o,'in_file')
	workflow.connect(nki_nodes.RSFC_printToFile,'ts_oneD',nki_nodes.RSFC_register_o,'name')
	workflow.connect(renamer,'register',nki_nodes.RSFC_register_o,'whichNode')
	workflow.connect(nki_nodes.RSFC_register, 'out_file',nki_nodes.RSFC_register_r,'in_file')
	workflow.connect(nki_nodes.RSFC_register_o,'out_file',nki_nodes.RSFC_register_r,'format_string')

	workflow.connect(nki_nodes.RSFC_smooth, 'out_file',nki_nodes.RSFC_smooth_o,'in_file')
	workflow.connect(nki_nodes.RSFC_printToFile,'ts_oneD',nki_nodes.RSFC_smooth_o,'name')
	workflow.connect(renamer,'smooth',nki_nodes.RSFC_smooth_o,'whichNode')
	workflow.connect(nki_nodes.RSFC_smooth, 'out_file',nki_nodes.RSFC_smooth_r,'in_file')
	workflow.connect(nki_nodes.RSFC_smooth_o,'out_file',nki_nodes.RSFC_smooth_r,'format_string')




def renameOutputs():

	rename_RSFC_outputs()

def getNodeOutputNames():

	names = ['smooth','corrs','z_trans','register']


	return names

def getRenamer():

	global renamer
	global rest_name
	global anat_name
	global FWHM

	nodeOutputNames = getNodeOutputNames()
	renamer = pe.Node(interface=util.IdentityInterface(fields= nodeOutputNames), name="renamer")

	#renamer for RSFC
	renamer.inputs.corrs = 'corrs' 
	renamer.inputs.z_trans = 'z_trans'
	renamer.inputs.register = 'register'
	renamer.inputs.smooth = 'Z_2standard_FWHM' + str(FWHM)

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
	
	wfname =  'RSFC_'+strategy

	workflow = pe.Workflow(name=wfname)
	workflow.base_dir = working_dir
	RSFC()
	#renameOutputs()
	#makeOutputConnections(analysisdirectory)

	workflow.run(plugin='MultiProc', plugin_args={'n_procs' : numCores})
	#workflow.write_graph()


def processSubjects(analysisdirectory,subject_list):

	
	try:
		fp = open(subject_list,'r')

		lines = fp.readlines()

		fp.close()
		line1 = []
		cnt = 0
		for line in lines:
			line = line.rstrip('\r\n')

			flag = 0

			print 'returned %d' %(flag)
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
	prior_dir = prior + '/%s' %(standard_res)
	standard_brain = FSLDIR + '/data/standard/MNI152_T1_%s_brain.nii.gz' %(standard_res)
	logFile = parsermap['logfile']
	FWHM = float(parsermap['fwhm'])


def main():

	if ( len(sys.argv) < 4 ):
		sys.stderr.write("./rsfc_pipeline.py <Number_of_cores> </path/to/dir_setup_X.ini> <strategy-label>")
	else:
	
		readDirSetup()
		readSubjects()

	

if __name__ == "__main__":

	sys.exit(main())
