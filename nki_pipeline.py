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
from nipype.interfaces.afni.preprocess import ThreedSkullStrip
from nipype.interfaces.afni.preprocess import Threedcalc
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
datasource_seeds = None
workflow = None
renamer = None
HP = None
LP = None
hp = None
lp = None
zeroSM = ''
#*******************************************************************************************



def anatpreproc():

	global anat_name
	global infosource
	global datasource
	global workflow
	
	#workflow.connect( infosource,'subject_id',datasource, 'subject_id' )
	workflow.connect(datasource, 'anat',nki_nodes.anat_refit,'in_file')
	workflow.connect(nki_nodes.anat_refit,'out_file',nki_nodes.anat_reorient,'in_file')
	workflow.connect(nki_nodes.anat_reorient,'out_file',nki_nodes.anat_skullstrip,'in_file')
	workflow.connect(nki_nodes.anat_skullstrip,'out_file',nki_nodes.anat_calc,'infile_b')
	workflow.connect(nki_nodes.anat_reorient,'out_file',nki_nodes.anat_calc,'infile_a')
	
def last_vol(vols):

	v = []
	for vol in vols:
		v.append(int(vol) - 1)

	return v

def funcpreproc():

	global rest_name
	global infosource
	global datasource
	global workflow
	global zeroSM

	workflow.connect(datasource,'rest',nki_nodes.lifeSaver,'in_file')
	workflow.connect(datasource,('rest',adjustPath),nki_nodes.lifeSaver,'format_string')
	workflow.connect(nki_nodes.lifeSaver, 'out_file',nki_nodes.TR,'in_files')
	workflow.connect(nki_nodes.lifeSaver, 'out_file',nki_nodes.NVOLS,'in_files')
	#workflow.connect(datasource, 'rest',nki_nodes.TR,'in_files')
	#workflow.connect(datasource, 'rest',nki_nodes.NVOLS,'in_files')
	workflow.connect(nki_nodes.NVOLS,('nvols',last_vol),nki_nodes.func_calc,'stop_idx')
	workflow.connect(nki_nodes.lifeSaver, 'out_file', nki_nodes.func_calc,'infile_a')
	#workflow.connect(datasource, 'rest', nki_nodes.func_calc,'infile_a')
	workflow.connect( nki_nodes.func_calc , 'out_file', nki_nodes.func_refit , 'in_file' )
	workflow.connect( nki_nodes.func_refit , 'out_file' , nki_nodes.func_reorient , 'in_file' )
	workflow.connect( nki_nodes.func_reorient , 'out_file' , nki_nodes.func_tstat , 'in_file' )
	workflow.connect( nki_nodes.func_reorient , 'out_file' , nki_nodes.func_volreg , 'in_file' )
	workflow.connect( nki_nodes.func_tstat , 'out_file' , nki_nodes.func_volreg , 'basefile' )
	workflow.connect( nki_nodes.func_volreg,'out_file',nki_nodes.func_tstat_1,'in_file')
	workflow.connect( nki_nodes.func_reorient , 'out_file' , nki_nodes.func_volreg_1 , 'in_file' )
	workflow.connect( nki_nodes.func_tstat_1 , 'out_file' , nki_nodes.func_volreg_1 , 'basefile' )
	workflow.connect( nki_nodes.func_volreg_1 , 'out_file', nki_nodes.func_automask , 'in_file')
	workflow.connect( nki_nodes.func_volreg_1,'out_file', nki_nodes.func_calcR,'infile_a')
	workflow.connect( nki_nodes.func_automask,'out_file', nki_nodes.func_calcR,'infile_b')
	workflow.connect( nki_nodes.func_calcR,'out_file', nki_nodes.func_mean,'in_file')
	workflow.connect( nki_nodes.func_calcR,'out_file', nki_nodes.func_scale,'in_file')
	workflow.connect( nki_nodes.func_scale,'out_file',nki_nodes.func_mask,'in_file')



def registration():

	global standard_brain
	global rest_name
	global FSLDIR
	global standard_res
	global workflow
	global infosource

	
	workflow.connect(nki_nodes.func_mean, 'out_file', nki_nodes.reg_flirt, 'in_file')
	workflow.connect(nki_nodes.anat_calc,'out_file', nki_nodes.reg_flirt, 'reference')
	workflow.connect(nki_nodes.reg_flirt, 'out_matrix_file', nki_nodes.reg_xfm1, 'in_file')
	workflow.connect(nki_nodes.anat_calc,'out_file', nki_nodes.reg_flirt1, 'in_file')
	workflow.connect(infosource,'standard_res_brain', nki_nodes.reg_flirt1, 'reference')
	workflow.connect(nki_nodes.reg_flirt1,'out_matrix_file',nki_nodes.reg_xfm2, 'in_file')
	workflow.connect(nki_nodes.reg_flirt,'out_matrix_file',nki_nodes.reg_xfm3, 'in_file')
	workflow.connect(nki_nodes.reg_flirt1, 'out_matrix_file',nki_nodes.reg_xfm3, 'in_file2')
	workflow.connect(nki_nodes.reg_xfm3,'out_file',nki_nodes.reg_flirt2,'in_matrix_file')
	workflow.connect(infosource,'standard_res_brain', nki_nodes.reg_flirt2, 'reference')
	workflow.connect(nki_nodes.func_mean, 'out_file', nki_nodes.reg_flirt2, 'in_file')
	workflow.connect(nki_nodes.reg_xfm3, 'out_file', nki_nodes.reg_xfm4, 'in_file')
	workflow.connect(nki_nodes.anat_reorient, 'out_file', nki_nodes.reg_fnt , 'in_file')
	workflow.connect(nki_nodes.reg_flirt1,'out_matrix_file', nki_nodes.reg_fnt, 'affine_file')
	workflow.connect(infosource,'standard',nki_nodes.reg_fnt,'ref_file')
	workflow.connect(infosource,'standard_brain_mask_dil',nki_nodes.reg_fnt,'refmask_file' )
	workflow.connect(infosource,'config_file',nki_nodes.reg_fnt,'config_file')
	workflow.connect(nki_nodes.func_mean, 'out_file',nki_nodes.reg_warp,'in_file')
	workflow.connect(infosource,'standard',nki_nodes.reg_warp, 'ref_file')
	workflow.connect(nki_nodes.reg_fnt, 'fieldcoeff_file', nki_nodes.reg_warp,'field_file')
	workflow.connect(nki_nodes.reg_flirt,'out_matrix_file', nki_nodes.reg_warp,'premat')

def pick_wm_0(probability_maps):

	import sys
	import os
	if(isinstance(probability_maps,list)):
		for file in probability_maps:
			print file
			if file.endswith("prob_0.nii.gz"):
				
				return file
	return None


def pick_wm_1(probability_maps):
	
	import sys
	if(isinstance(probability_maps,list)):
		for file in probability_maps:
			print file
			if file.endswith("prob_2.nii.gz"):
			
				return file
	return None

def segment():
	
	global workflow
	global scripts_dir
	global standard_brain
	global rest_name
	global FSLDIR
	global standard_res
	global infosource
	global zeroSM

	workflow.connect(nki_nodes.anat_calc,'out_file',nki_nodes.seg_segment, 'in_files' )
	workflow.connect( nki_nodes.seg_segment, ('probability_maps',pick_wm_0), nki_nodes.seg_flirt , 'in_file' )
	workflow.connect(nki_nodes.func_mean,'out_file',nki_nodes.seg_flirt,'reference')
	workflow.connect(nki_nodes.reg_xfm1, 'out_file', nki_nodes.seg_flirt,'in_matrix_file')

	workflow.connect( nki_nodes.seg_flirt , 'out_file' , nki_nodes.seg_smooth , 'in_file' )
	workflow.connect( nki_nodes.seg_smooth, 'out_file', nki_nodes.seg_flirt1, 'in_file')

	workflow.connect(infosource,'standard_res_brain', nki_nodes.seg_flirt1, 'reference')
	workflow.connect(nki_nodes.reg_xfm3, 'out_file',nki_nodes.seg_flirt1, 'in_matrix_file')
	workflow.connect( nki_nodes.seg_flirt1, 'out_file', nki_nodes.seg_smooth1, 'in_file')
	workflow.connect( infosource, 'PRIOR_CSF', nki_nodes.seg_smooth1, 'operand_files')
	workflow.connect( nki_nodes.seg_smooth1,'out_file', nki_nodes.seg_flirt2,'in_file')
	workflow.connect( nki_nodes.func_mean,'out_file', nki_nodes.seg_flirt2,'reference')
	workflow.connect( nki_nodes.reg_xfm4, 'out_file', nki_nodes.seg_flirt2,'in_matrix_file')
	workflow.connect( nki_nodes.seg_flirt2, 'out_file', nki_nodes.seg_thresh, 'in_file')
	workflow.connect( nki_nodes.seg_thresh,'out_file', nki_nodes.seg_mask,'in_file')
	workflow.connect(nki_nodes.func_mask, 'out_file' ,nki_nodes.seg_copy, 'in_file' )
	workflow.connect( nki_nodes.seg_copy,'out_file',  nki_nodes.seg_mask,'operand_files')
	workflow.connect(  nki_nodes.seg_segment, ('probability_maps',pick_wm_1),  nki_nodes.seg_flirt3 , 'in_file' )
	workflow.connect( nki_nodes.func_mean,'out_file', nki_nodes.seg_flirt3 , 'reference' )
	workflow.connect( nki_nodes.reg_xfm1, 'out_file', nki_nodes.seg_flirt3 , 'in_matrix_file' )

	workflow.connect( nki_nodes.seg_flirt3, 'out_file',  nki_nodes.seg_smooth2, 'in_file')
	workflow.connect( nki_nodes.seg_smooth2,'out_file', nki_nodes.seg_flirt4, 'in_file')
	workflow.connect( infosource,'standard_res_brain', nki_nodes.seg_flirt4, 'reference')
	workflow.connect( nki_nodes.reg_xfm3, 'out_file', nki_nodes.seg_flirt4, 'in_matrix_file')
	workflow.connect( nki_nodes.seg_flirt4,'out_file', nki_nodes.seg_prior1, 'in_file')
	workflow.connect( infosource, 'PRIOR_WHITE', nki_nodes.seg_prior1, 'operand_files')
	workflow.connect( nki_nodes.seg_prior1,'out_file', nki_nodes.seg_flirt5, 'in_file')
	workflow.connect( nki_nodes.func_mean,'out_file', nki_nodes.seg_flirt5, 'reference')
	workflow.connect( nki_nodes.reg_xfm4,'out_file', nki_nodes.seg_flirt5, 'in_matrix_file')
	workflow.connect( nki_nodes.seg_flirt5,'out_file', nki_nodes.seg_thresh1,'in_file')
	workflow.connect( nki_nodes.seg_thresh1,'out_file', nki_nodes.seg_mask1,'in_file')
	workflow.connect( nki_nodes.seg_copy,'out_file', nki_nodes.seg_mask1,'operand_files')



def signalExtractionEfficient():

	import sys
	global workflow
	global scripts_dir
	global standard_brain
	global rest_name
	global FSLDIR
	global standard_res
	global infosource
	global which_regression
	

	workflow.connect(nki_nodes.NVOLS,'nvols',nki_nodes.nuisance_poly,'nvols')	
	if(which_regression.lower() == 'default'):
		#make initial connections to facilitate the default regression : global, wm, csf
		workflow.connect(nki_nodes.func_scale,'out_file',nki_nodes.nuisance_globalE,'in_file')
		workflow.connect(nki_nodes.seg_copy,'out_file',nki_nodes.nuisance_globalE,'mask')
		workflow.connect(nki_nodes.func_scale,'out_file',nki_nodes.nuisance_csf,'in_file')
		workflow.connect(nki_nodes.seg_mask,'out_file',nki_nodes.nuisance_csf,'mask')
		workflow.connect(nki_nodes.func_scale,'out_file',nki_nodes.nuisance_wm,'in_file')
		workflow.connect(nki_nodes.seg_mask1,'out_file',nki_nodes.nuisance_wm,'mask')

		workflow.connect(nki_nodes.func_scale,'out_file', nki_nodes.MC,'pp')
		workflow.connect(infosource,'nuisance_template',nki_nodes.FSF,'nuisance_template')
		workflow.connect(nki_nodes.func_scale,'out_file',nki_nodes.FSF,'rest_pp')
		workflow.connect(nki_nodes.TR,'TR', nki_nodes.FSF,'TR')
		workflow.connect(nki_nodes.NVOLS,'nvols',nki_nodes.FSF,'n_vols')
		workflow.connect(nki_nodes.MC,'EV_Lists',nki_nodes.copyS,'EV_Lists')
		workflow.connect(nki_nodes.FSF,'nuisance_files',nki_nodes.copyS,'nuisance_files')
		workflow.connect(nki_nodes.nuisance_globalE,'out_file',nki_nodes.copyS,'global1Ds')
		workflow.connect(nki_nodes.nuisance_csf,'out_file',nki_nodes.copyS,'csf1Ds')
		workflow.connect(nki_nodes.nuisance_wm,'out_file',nki_nodes.copyS,'wm1Ds')
		
	
	elif(which_regression.lower() == "no_global_signal"):

		workflow.connect(nki_nodes.func_scale,'out_file',nki_nodes.nuisance_csf,'in_file')
		workflow.connect(nki_nodes.seg_mask,'out_file',nki_nodes.nuisance_csf,'mask')
		workflow.connect(nki_nodes.func_scale,'out_file',nki_nodes.nuisance_wm,'in_file')
		workflow.connect(nki_nodes.seg_mask1,'out_file',nki_nodes.nuisance_wm,'mask')
		
		workflow.connect(nki_nodes.func_scale,'out_file', nki_nodes.MC,'pp')
		workflow.connect(infosource,'nuisance_template_noglobal',nki_nodes.FSF,'nuisance_template')
		workflow.connect(nki_nodes.func_scale,'out_file',nki_nodes.FSF,'rest_pp')
		workflow.connect(nki_nodes.TR,'TR', nki_nodes.FSF,'TR')
		workflow.connect(nki_nodes.NVOLS,'nvols',nki_nodes.FSF,'n_vols')
		workflow.connect(nki_nodes.MC,'EV_Lists',nki_nodes.copyS,'EV_Lists')
		workflow.connect(nki_nodes.FSF,'nuisance_files',nki_nodes.copyS,'nuisance_files')
		workflow.connect(infosource,'dummy',nki_nodes.copyS,'global1Ds')
		workflow.connect(nki_nodes.nuisance_csf,'out_file',nki_nodes.copyS,'csf1Ds')
		workflow.connect(nki_nodes.nuisance_wm,'out_file',nki_nodes.copyS,'wm1Ds')

	workflow.connect(nki_nodes.nuisance_poly,'regressors',nki_nodes.copyS,'regressors')
	workflow.connect(nki_nodes.func_volreg , 'oned_file',nki_nodes.MC,'in_file')


def nuisance():

	

	global workflow
	global scripts_dir
	global standard_brain
	global rest_name
	global FSLDIR
	global standard_res
	global infosource
	global which_regression


	signalExtractionEfficient()
	
	workflow.connect(nki_nodes.FSF,'nuisance_files',nki_nodes.nuisance_featM,'fsf_file')
	workflow.connect(nki_nodes.copyS,'EV_final_lists_set',nki_nodes.nuisance_featM,'ev_files')

	#which way to go
	if(which_regression == 'no_global_signal'):
		workflow.connect(nki_nodes.func_scale,'out_file',nki_nodes.nuisance_brick,'in_file')
		workflow.connect(nki_nodes.func_scale,'out_file',nki_nodes.nuisance_fgls,'in_file')
	
	else:
		workflow.connect(nki_nodes.func_scale,'out_file',nki_nodes.nuisance_brick,'in_file')
		workflow.connect(nki_nodes.func_scale,'out_file',nki_nodes.nuisance_fgls,'in_file')

#	workflow.connect(nki_nodes.nuisance_selectnode,'out',nki_nodes.nuisance_brick,'in_file')
#	workflow.connect(nki_nodes.func_detrendc,'out_file',nki_nodes.nuisance_brick,'in_file')
	workflow.connect(nki_nodes.func_mask, 'out_file' ,nki_nodes.nuisance_brick,'mask')
	workflow.connect(nki_nodes.nuisance_featM, 'design_file',nki_nodes.nuisance_fgls,'design_file')
	workflow.connect(nki_nodes.func_scale, ('out_file',nki_nodes.getStatsDir),nki_nodes.nuisance_fgls, 'results_dir' )
	#workflow.connect(nki_nodes.func_detrendc,'out_file',nki_nodes.nuisance_fgls,'in_file')
	workflow.connect(nki_nodes.nuisance_brick,'min_val',nki_nodes.nuisance_fgls,'threshold')
	workflow.connect(nki_nodes.nuisance_fgls,'residual4d',nki_nodes.nuisance_stat,'in_file')
	workflow.connect(nki_nodes.nuisance_fgls,'residual4d',nki_nodes.nuisance_calc,'infile_a')
	workflow.connect(nki_nodes.nuisance_stat,'out_file', nki_nodes.nuisance_calc,'infile_b')
	workflow.connect(nki_nodes.nuisance_calc,'out_file',nki_nodes.func_filter,'in_file')
	workflow.connect( infosource,'hp',nki_nodes.func_filter,'highpass')
	workflow.connect( infosource,'lp',nki_nodes.func_filter,'lowpass')
	workflow.connect(nki_nodes.func_filter,'out_file', nki_nodes.nuisance_warp,'in_file')
	workflow.connect(infosource,'standard',nki_nodes.nuisance_warp,'ref_file')
	workflow.connect(nki_nodes.reg_fnt, 'fieldcoeff_file',nki_nodes.nuisance_warp,'field_file')
	workflow.connect(nki_nodes.reg_flirt,'out_matrix_file',nki_nodes.nuisance_warp,'premat')
	workflow.connect(nki_nodes.func_automask,'out_file', nki_nodes.nuisance_warp_1,'in_file')
	workflow.connect(infosource,'standard',nki_nodes.nuisance_warp_1,'ref_file')
	workflow.connect(nki_nodes.reg_fnt, 'fieldcoeff_file',nki_nodes.nuisance_warp_1,'field_file')
	workflow.connect(nki_nodes.reg_flirt,'out_matrix_file',nki_nodes.nuisance_warp_1,'premat')

def despkingWF():
	
	global workflow
	global scripts_dir
	global standard_brain
	global rest_name
	global FSLDIR
	global standard_res
	global infosource
	global zeroSM

	workflow.connect(nki_nodes.func_despike, 'out_file',nki_nodes.alff_detrend,'in_file')

	if not(zeroSM.lower() == 'on'):
		workflow.connect(nki_nodes.alff_detrend,'out_file',nki_nodes.alff_smooth,'in_file')
		workflow.connect(nki_nodes.func_automask,'out_file',nki_nodes.alff_smooth,'operand_files')
		workflow.connect(nki_nodes.alff_smooth,'out_file',nki_nodes.alff_scale,'in_file')
	else:
		workflow.connect(nki_nodes.alff_detrend,'out_file',nki_nodes.alff_scale,'in_file')
	workflow.connect(nki_nodes.alff_scale,'out_file',nki_nodes.alff_cp,'in_file')
	workflow.connect(nki_nodes.alff_cp,'out_file',nki_nodes.alff_mean,'in_file')

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




def readRestingImageStats(subject,analysisdirectory):

	global rest_name
	global rest_dir

	img = analysisdirectory + '/' + subject + '/' + rest_dir + '/' + rest_name + '.nii.gz'
	cmd = "3dinfo %s | grep 'Time step = *' > %s/TRInfo_N_%s" %(img,scripts_dir,(subject.split('/'))[0])
	print cmd
	sys.stderr.write('\n' + commands.getoutput(cmd))
	
	file = open('%s/TRInfo_N_%s' %(scripts_dir,(subject.split('/'))[0]),'r')

	lines = file.readlines()

	TR = ""
	n_vols = ""
	for line in lines:
		line = line.rstrip('\r\n')
		
		fields = line.split(' ')
		
		i = fields.index('step')
		TR = (fields[i+2]).rstrip('s')

		i = fields.index('steps')
		n_vols = fields[i+2]
		
	file.close()
	
	cmd = "rm -f %s/TRInfo_N_%s" %(scripts_dir,(subject.split('/'))[0])
	print cmd
	sys.stderr.write('\n'+commands.getoutput(cmd))

	return TR+' '+n_vols

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
	global datasource_seeds
	global working_dir
	global anat_name
	global anat_dir
	global rest_name
	global rest_dir
	global seed_list


	datasource = pe.Node( interface = nio.DataGrabber(infields=['subject_id'], outfields = [ 'anat', 'rest' ]) , name= 'datasource')
	datasource.inputs.base_directory = analysisdirectory
	#datasource.inputs.template = '%s/*/%s.nii.gz'
	datasource.inputs.template = '%s/*/*/%s.nii.gz'
	datasource.inputs.template_args = dict( anat = [['subject_id',anat_name]] , rest = [['subject_id',rest_name]])
	datasource.iterables = ('subject_id', sublist)


	nki_nodes.setSeedList(seed_list)

	#datasource_seeds = pe.Node( interface = nio.DataGrabber(infields=['seeds'],outfields = ['seeds']), name = 'datasource_seeds')
	#datasource_seeds.inputs.base_directory = seeds_dir
	#datasource_seeds.inputs.template = '%s'
	#datasource_seeds.inputs.seeds = seeds



def getFields():

	fields = ['brain_symmetric','symm_standard','standard_res_brain','standard','standard_brain_mask_dil','twomm_brain_mask_dil','config_file_twomm','config_file','PRIOR_CSF','PRIOR_WHITE','nuisance_template','nuisance_template_cc','nuisance_template_noglobal','HP','LP','hp','lp','which_regression','dummy']
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
	datasink.inputs.regexp_substitutions = [(r"_subject_id_(\w|\d)+", ''),(r"-",'/'),(r"(.)+_func_(\w|\d)+",''),(r"(.)+_reg_(\w|\d)+",''),(r"(.)+_seg_(\w|\d)+",''),(r"(.)+_nuisance_(\w|\d)+",''),(r"(.)+_alff_(\w|\d)+/",''),(r"(.)+_falff(\w|\d)+/",''),(r"(.)+_RSFC_(\w|\d)+",''),(r"(.)+/_seed_(.)+seeds\.\./",'')]


	return datasink


def getSinkRSFC():

	datasinkRSFC = None
	
	datasinkRSFC = pe.Node(nio.DataSink(), name = 'sinker-RSFC')
	datasinkRSFC.inputs.base_directory = os.path.abspath('/')
	datasinkRSFC.inputs.container = os.path.abspath('/')
	datasinkRSFC.inputs.regexp_substitutions = [(r"_subject_id_(\w|\d)+", ''),(r"-",'/'),(r"(.)+_func_(\w|\d)+",''),(r"(.)+_reg_(\w|\d)+",''),(r"(.)+_seg_(\w|\d)+",''),(r"(.)+_nuisance_(\w|\d)+",''),(r"(.)+_alff_(\w|\d)+/",''),(r"(.)+_RSFC_(\w|\d)+",''),(r"(.)+/_seed_(.)+seeds\.\./",'')]


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


def adjustReg(rest_path):

	import re
	import sys
	result_path = []

	if(isinstance(rest_path,list)):
		for r_p in rest_path:

			split_path = r_p.split('-')
	
			path = ""
			for index in range(0,len(split_path) - 1):
				path += split_path[index] + "-"

			path += "reg" + "-"
			path += split_path[len(split_path) - 1]
	
			print "\nprocessed path %s\n" %(path)
			result_path.append(path)
	else:

		split_path = rest_path.split('-')
	
		path = ""
		for index in range(0,len(split_path) - 1):
			path += split_path[index] + "-"

		path += "reg" + "-"
		path += split_path[len(split_path) - 1]
	
		print "\nprocessed path %s\n" %(path)
		result_path.append(path)
		

	return result_path

def adjustSegment(rest_path):

	import re
	import sys
	result_path = []


	if(isinstance(rest_path,list)):
		for r_p in rest_path:

			split_path = r_p.split('-')
	
			path = ""
			for index in range(0,len(split_path) - 1):
				path += split_path[index] + "-"

			path += "segment" + "-"
			path += split_path[len(split_path) - 1]

			print "\nprocessed path %s\n" %(path)
			result_path.append(path)
	else:

		split_path = rest_path.split('-')

		path = ""
		for index in range(0,len(split_path) - 1):
			path += split_path[index] + "-"

		path += "segment" + "-"
		path += split_path[len(split_path) - 1]

		print "\nprocessed path %s\n" %(path)
		result_path.append(path)
		

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


def makeOutputConnectionsAnat(datasinkAnat):

	global workflow
	global infosource

	workflow.connect(datasource, 'anat', datasinkAnat, 'container')
	
	## Connect anatpreproc nodes to datasink
	workflow.connect(nki_nodes.anat_refit,'out_file',datasinkAnat,'@mprage')
	workflow.connect(nki_nodes.anat_reorient,'out_file',datasinkAnat,'@mprage_RPI')
	workflow.connect(nki_nodes.anat_skullstrip_r,'out_file', datasinkAnat,'@mprage_surf')
	workflow.connect(nki_nodes.anat_calc,'out_file',datasinkAnat,'@mprage_brain')
	

def makeOutputConnectionsFunc(datasink):
	## Connect funcpreproc nodes to datasink
#	workflow.connect(datasource, 'rest', datasink, 'container')
	workflow.connect( nki_nodes.func_calc_r , 'out_file', datasink , '@rest_dr' )
	workflow.connect( nki_nodes.func_refit_r , 'out_file' , datasink, '@rest_dr_1' )
	workflow.connect( nki_nodes.func_reorient_r , 'out_file' , datasink, '@rest_ro' )
	workflow.connect( nki_nodes.func_tstat_r , 'out_file' , datasink, '@rest_ro_mean' )
	workflow.connect( nki_nodes.func_volrego_r , 'out_file', datasink, '@rest_mc_1D')
	workflow.connect( nki_nodes.func_volreg_r , 'out_file', datasink, '@rest_mc')
	workflow.connect( nki_nodes.func_automask_r,'out_file', datasink,'@rest_mask')
	workflow.connect( nki_nodes.func_calcR_r,'out_file', datasink,'@rest_ss')

	
	workflow.connect( nki_nodes.func_mean_r, 'out_file', datasink, '@rest_mean')
	workflow.connect( nki_nodes.func_scale_r, 'out_file', datasink, '@rest_gms')
	workflow.connect( nki_nodes.func_mask_r, 'out_file' , datasink, '@rest_pp_mask')

def makeOutputConnectionsReg(datasink, datasinkAnat):

	## Connect registration nodes to datasink

	workflow.connect(nki_nodes.reg_flirt_r, 'out_file',nki_nodes.lifeSaver_reg_flirt,'in_file')
	workflow.connect(nki_nodes.reg_flirt_r,('out_file',adjustReg),nki_nodes.lifeSaver_reg_flirt,'format_string')
	workflow.connect(nki_nodes.lifeSaver_reg_flirt, 'out_file',datasink,'@example_func2highres')

	workflow.connect(nki_nodes.reg_flirto_r, 'out_file',nki_nodes.lifeSaver_func2highresmat,'in_file')
	workflow.connect(nki_nodes.reg_flirto_r,('out_file',adjustReg),nki_nodes.lifeSaver_func2highresmat,'format_string')
	workflow.connect(nki_nodes.lifeSaver_func2highresmat, 'out_file',datasink,'@example_func2highresmat')

	workflow.connect(nki_nodes.reg_xfm1_r, 'out_file',nki_nodes.lifeSaver_highres2example_funcmat,'in_file')
	workflow.connect(nki_nodes.reg_xfm1_r,('out_file',adjustReg),nki_nodes.lifeSaver_highres2example_funcmat,'format_string')
	workflow.connect(nki_nodes.lifeSaver_highres2example_funcmat, 'out_file',datasink,'@highres2example_funcmat')

	workflow.connect(nki_nodes.reg_flirt1o_r, 'out_file',datasinkAnat, 'reg.@highres2standard')
	
	workflow.connect(nki_nodes.reg_flirt1_r, 'out_file',datasinkAnat, 'reg.@highres2standardmat')
	workflow.connect(nki_nodes.reg_xfm2_r, 'out_file',datasinkAnat,'reg.@standard2highresmat')

	workflow.connect(nki_nodes.reg_xfm3_r, 'out_file',nki_nodes.lifeSaver_func2standardmat,'in_file')
	workflow.connect(nki_nodes.reg_xfm3_r, ('out_file',adjustReg),nki_nodes.lifeSaver_func2standardmat,'format_string')
	workflow.connect(nki_nodes.lifeSaver_func2standardmat,'out_file',datasink,'@example_func2standardmat')
#	workflow.connect(nki_nodes.reg_xfm3, 'out_file',datasink,'@example_func2standardmat')


	workflow.connect(nki_nodes.reg_flirt2_r, 'out_file',nki_nodes.lifeSaver_func2standard,'in_file')
	workflow.connect(nki_nodes.reg_flirt2_r, ('out_file',adjustReg),nki_nodes.lifeSaver_func2standard,'format_string')
	workflow.connect(nki_nodes.lifeSaver_func2standard,'out_file',datasink,'@example_func2standard')
#	workflow.connect(nki_nodes.reg_flirt2, 'out_file',datasink,'@example_func2standard')
	

	workflow.connect(nki_nodes.reg_xfm4_r, 'out_file',nki_nodes.lifeSaver_standard2example_funcmat,'in_file')
	workflow.connect(nki_nodes.reg_xfm4_r, ('out_file',adjustReg),nki_nodes.lifeSaver_standard2example_funcmat,'format_string')
	workflow.connect(nki_nodes.lifeSaver_standard2example_funcmat,'out_file',datasink,'@standard2example_funcmat')
	#workflow.connect(nki_nodes.reg_xfm4, 'out_file',datasink,'@standard2example_funcmat')

	workflow.connect(nki_nodes.reg_fnt_r , 'out_file',datasinkAnat,'reg.@highres2standard_NL')
	workflow.connect(nki_nodes.reg_fntj_r, 'out_file',datasinkAnat,'reg.@highres2standard_jac')
	workflow.connect(nki_nodes.reg_fntf_r, 'out_file',datasinkAnat,'reg.@highres2standard_warp')
	
	workflow.connect(nki_nodes.reg_warp_r, 'out_file',nki_nodes.lifeSaver_example_func2standard_NL,'in_file')
	workflow.connect(nki_nodes.reg_warp_r, ('out_file',adjustReg),nki_nodes.lifeSaver_example_func2standard_NL,'format_string')
	workflow.connect(nki_nodes.lifeSaver_example_func2standard_NL,'out_file',datasink,'@example_func2standard_NL')
	#workflow.connect(nki_nodes.reg_warp_r,'out_file',datasinkAnat,'reg.@example_func2standard_NL')


def makeOutputConnectionsSeg(datasink, datasinkAnat):
	## Connect segmentation nodes to datasink	
	global zeroSM

	workflow.connect( nki_nodes.seg_segment, 'probability_maps',datasinkAnat,'segment.@seg_segment' )
	
	workflow.connect(nki_nodes.seg_flirt_r, 'out_file',nki_nodes.Saver_seg_flirt,'in_file')
	workflow.connect(nki_nodes.func_mean_r,'out_file',nki_nodes.Saver_seg_flirt,'pp')
	workflow.connect(nki_nodes.seg_flirt_r, 'out_file',nki_nodes.lifeSaver_seg_flirt,'in_file')
	workflow.connect(nki_nodes.Saver_seg_flirt,'out_file',nki_nodes.lifeSaver_seg_flirt,'format_string')
	workflow.connect(nki_nodes.lifeSaver_seg_flirt,'out_file',datasink,'@seg_flirt')
#	workflow.connect( nki_nodes.seg_flirt , 'out_file', datasink, 'segment.@seg_flirt' )

	workflow.connect(nki_nodes.seg_smooth_r, 'out_file',nki_nodes.Saver_seg_smooth,'in_file')
	workflow.connect(nki_nodes.func_mean_r,'out_file',nki_nodes.Saver_seg_smooth,'pp')
	workflow.connect(nki_nodes.seg_smooth_r, 'out_file',nki_nodes.lifeSaver_seg_smooth,'in_file')
	workflow.connect(nki_nodes.Saver_seg_smooth, 'out_file',nki_nodes.lifeSaver_seg_smooth,'format_string')
	workflow.connect(nki_nodes.lifeSaver_seg_smooth,'out_file',datasink,'@seg_smooth')
#	workflow.connect( nki_nodes.seg_smooth, 'out_file', datasink, 'segment.@seg_smooth')


	workflow.connect(nki_nodes.seg_flirt1_r, 'out_file',nki_nodes.Saver_seg_flirt1,'in_file')
	workflow.connect(nki_nodes.func_mean_r,'out_file',nki_nodes.Saver_seg_flirt1,'pp')
	workflow.connect(nki_nodes.seg_flirt1_r, 'out_file',nki_nodes.lifeSaver_seg_flirt1,'in_file')
	workflow.connect(nki_nodes.Saver_seg_flirt1, 'out_file',nki_nodes.lifeSaver_seg_flirt1,'format_string')
	workflow.connect(nki_nodes.lifeSaver_seg_flirt1,'out_file',datasink,'@seg_flirt1')
#	workflow.connect( nki_nodes.seg_flirt1, 'out_file', datasink, 'segment.@seg_flirt1')

	
	workflow.connect(nki_nodes.seg_smooth1_r, 'out_file',nki_nodes.Saver_seg_smooth1,'in_file')
	workflow.connect(nki_nodes.func_mean_r,'out_file',nki_nodes.Saver_seg_smooth1,'pp')
	workflow.connect(nki_nodes.seg_smooth1_r, 'out_file',nki_nodes.lifeSaver_seg_smooth1,'in_file')
	workflow.connect(nki_nodes.Saver_seg_flirt1, 'out_file',nki_nodes.lifeSaver_seg_smooth1,'format_string')
	workflow.connect(nki_nodes.lifeSaver_seg_smooth1,'out_file',datasink,'@seg_smooth1')
#	workflow.connect( nki_nodes.seg_smooth1,'out_file', datasink,'segment.@seg_smooth1')
	
	workflow.connect(nki_nodes.seg_flirt2_r, 'out_file',nki_nodes.Saver_seg_flirt2,'in_file')
	workflow.connect(nki_nodes.func_mean_r,'out_file',nki_nodes.Saver_seg_flirt2,'pp')
	workflow.connect(nki_nodes.seg_flirt2_r, 'out_file',nki_nodes.lifeSaver_seg_flirt2,'in_file')
	workflow.connect(nki_nodes.Saver_seg_flirt2, 'out_file',nki_nodes.lifeSaver_seg_flirt2,'format_string')
	workflow.connect(nki_nodes.lifeSaver_seg_flirt2,'out_file',datasink,'@seg_flirt2')
#	workflow.connect( nki_nodes.seg_flirt2, 'out_file', datasink, 'segment.@seg_flirt2')
	
	workflow.connect(nki_nodes.seg_thresh_r, 'out_file',nki_nodes.Saver_seg_thresh,'in_file')
	workflow.connect(nki_nodes.func_mean_r,'out_file',nki_nodes.Saver_seg_thresh,'pp')
	workflow.connect(nki_nodes.seg_thresh_r, 'out_file',nki_nodes.lifeSaver_seg_thresh,'in_file')
	workflow.connect(nki_nodes.Saver_seg_thresh, 'out_file',nki_nodes.lifeSaver_seg_thresh,'format_string')
	workflow.connect(nki_nodes.lifeSaver_seg_thresh,'out_file',datasink,'@seg_thresh')
#	workflow.connect( nki_nodes.seg_thresh,'out_file', datasink,'segment.@seg_thresh')
	
	workflow.connect(nki_nodes.seg_mask_r, 'out_file',nki_nodes.Saver_seg_mask,'in_file')
	workflow.connect(nki_nodes.func_mean_r,'out_file',nki_nodes.Saver_seg_mask,'pp')
	workflow.connect(nki_nodes.seg_mask_r, 'out_file',nki_nodes.lifeSaver_seg_mask,'in_file')
	workflow.connect(nki_nodes.Saver_seg_mask, 'out_file',nki_nodes.lifeSaver_seg_mask,'format_string')
	workflow.connect(nki_nodes.lifeSaver_seg_mask,'out_file',datasink,'@seg_mask')
#	workflow.connect( nki_nodes.seg_mask,'out_file', datasink,'segment.@seg_mask')
	
	workflow.connect(nki_nodes.seg_copy_r, 'out_file',nki_nodes.lifeSaver_seg_copy,'in_file')
	workflow.connect(nki_nodes.seg_copy_r, ('out_file',adjustSegment),nki_nodes.lifeSaver_seg_copy,'format_string')
	workflow.connect(nki_nodes.lifeSaver_seg_copy,'out_file',datasink,'@seg_copy')
#	workflow.connect( nki_nodes.seg_copy,'out_file',  datasink,'segment.@seg_copy')
	


	workflow.connect(nki_nodes.seg_flirt3_r, 'out_file',nki_nodes.Saver_seg_flirt3,'in_file')
	workflow.connect(nki_nodes.func_mean_r,'out_file',nki_nodes.Saver_seg_flirt3,'pp')
	workflow.connect(nki_nodes.seg_flirt3_r, 'out_file',nki_nodes.lifeSaver_seg_flirt3,'in_file')
	workflow.connect(nki_nodes.Saver_seg_flirt3, 'out_file',nki_nodes.lifeSaver_seg_flirt3,'format_string')
	workflow.connect(nki_nodes.lifeSaver_seg_flirt3,'out_file',datasink,'@seg_flirt3')
#	workflow.connect( nki_nodes.seg_flirt3, 'out_file', datasink, 'segment.@seg_flirt3')

	workflow.connect(nki_nodes.seg_smooth2_r, 'out_file',nki_nodes.Saver_seg_smooth2,'in_file')
	workflow.connect(nki_nodes.func_mean_r,'out_file',nki_nodes.Saver_seg_smooth2,'pp')
	workflow.connect(nki_nodes.seg_smooth2_r, 'out_file',nki_nodes.lifeSaver_seg_smooth2,'in_file')
	workflow.connect(nki_nodes.Saver_seg_smooth2, 'out_file',nki_nodes.lifeSaver_seg_smooth2,'format_string')
	workflow.connect(nki_nodes.lifeSaver_seg_smooth2,'out_file',datasink,'@seg_smooth2')
#	workflow.connect( nki_nodes.seg_smooth2,'out_file', datasink, 'segment.@seg_smooth2')
	
	workflow.connect(nki_nodes.seg_flirt4_r, 'out_file',nki_nodes.Saver_seg_flirt4,'in_file')
	workflow.connect(nki_nodes.func_mean_r,'out_file',nki_nodes.Saver_seg_flirt4,'pp')
	workflow.connect(nki_nodes.seg_flirt4_r, 'out_file',nki_nodes.lifeSaver_seg_flirt4,'in_file')
	workflow.connect(nki_nodes.Saver_seg_flirt4, 'out_file',nki_nodes.lifeSaver_seg_flirt4,'format_string')
	workflow.connect(nki_nodes.lifeSaver_seg_flirt4,'out_file',datasink,'@seg_flirt4')
#	workflow.connect( nki_nodes.seg_flirt4,'out_file', datasink, 'segment.@seg_flirt4')
#	
	workflow.connect(nki_nodes.seg_prior1_r, 'out_file',nki_nodes.Saver_seg_prior1,'in_file')
	workflow.connect(nki_nodes.func_mean_r,'out_file',nki_nodes.Saver_seg_prior1,'pp')
	workflow.connect(nki_nodes.seg_prior1_r, 'out_file',nki_nodes.lifeSaver_seg_prior1,'in_file')
	workflow.connect(nki_nodes.Saver_seg_prior1, 'out_file',nki_nodes.lifeSaver_seg_prior1,'format_string')
	workflow.connect(nki_nodes.lifeSaver_seg_prior1,'out_file',datasink,'@seg_prior1')
#	workflow.connect( nki_nodes.seg_prior1,'out_file', datasink, 'segment.@seg_prior1')
#	
	workflow.connect(nki_nodes.seg_flirt5_r, 'out_file',nki_nodes.Saver_seg_flirt5,'in_file')
	workflow.connect(nki_nodes.func_mean_r,'out_file',nki_nodes.Saver_seg_flirt5,'pp')
	workflow.connect(nki_nodes.seg_flirt5_r, 'out_file',nki_nodes.lifeSaver_seg_flirt5,'in_file')
	workflow.connect(nki_nodes.Saver_seg_flirt5, 'out_file',nki_nodes.lifeSaver_seg_flirt5,'format_string')
	workflow.connect(nki_nodes.lifeSaver_seg_flirt5,'out_file',datasink,'@seg_flirt5')
#	workflow.connect( nki_nodes.seg_flirt5,'out_file', datasink,'segment.@seg_flirt5')
#	
	workflow.connect(nki_nodes.seg_thresh1_r, 'out_file',nki_nodes.Saver_seg_thresh1,'in_file')
	workflow.connect(nki_nodes.func_mean_r,'out_file',nki_nodes.Saver_seg_thresh1,'pp')
	workflow.connect(nki_nodes.seg_thresh1_r, 'out_file',nki_nodes.lifeSaver_seg_thresh1,'in_file')
	workflow.connect(nki_nodes.Saver_seg_thresh1, 'out_file',nki_nodes.lifeSaver_seg_thresh1,'format_string')
	workflow.connect(nki_nodes.lifeSaver_seg_thresh1,'out_file',datasink,'@seg_thresh1')
#	workflow.connect( nki_nodes.seg_thresh1,'out_file', datasink,'segment.@seg_thresh1')
#	
	workflow.connect(nki_nodes.seg_mask1_r, 'out_file',nki_nodes.Saver_seg_mask1,'in_file')
	workflow.connect(nki_nodes.func_mean_r,'out_file',nki_nodes.Saver_seg_mask1,'pp')
	workflow.connect(nki_nodes.seg_mask1_r, 'out_file',nki_nodes.lifeSaver_seg_mask1,'in_file')
	workflow.connect(nki_nodes.Saver_seg_mask1, 'out_file',nki_nodes.lifeSaver_seg_mask1,'format_string')
	workflow.connect(nki_nodes.lifeSaver_seg_mask1,'out_file',datasink,'@seg_mask1')
#	workflow.connect( nki_nodes.seg_mask1,'out_file', datasink,'segment.@seg_mask1')


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
	workflow.connect(nki_nodes.func_mean_r,'out_file',nki_nodes.Saver_nuisance_featM,'pp')
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
	workflow.connect(nki_nodes.func_mean_r,'out_file',nki_nodes.Saver_nuisance_fgls,'pp')
	workflow.connect(nki_nodes.nuisance_fgls_r, 'out_file',nki_nodes.lifeSaver_nuisance_fgls,'in_file')
	workflow.connect(nki_nodes.Saver_nuisance_fgls, 'out_file',nki_nodes.lifeSaver_nuisance_fgls,'format_string')
	workflow.connect(nki_nodes.lifeSaver_nuisance_fgls,'out_file',datasink,'@nuisance_fgls')
#	workflow.connect(nki_nodes.nuisance_fgls,'residual4d',datasink,'@nuisance_fgls_')

	workflow.connect(nki_nodes.nuisance_stat_r, 'out_file',nki_nodes.Saver_nuisance_stat,'in_file')
	workflow.connect(nki_nodes.func_mean_r,'out_file',nki_nodes.Saver_nuisance_stat,'pp')
	workflow.connect(nki_nodes.nuisance_stat_r, 'out_file',nki_nodes.lifeSaver_nuisance_stat,'in_file')
	workflow.connect(nki_nodes.Saver_nuisance_stat, 'out_file',nki_nodes.lifeSaver_nuisance_stat,'format_string')
	workflow.connect(nki_nodes.lifeSaver_nuisance_stat,'out_file',datasink,'@nuisance_stat')
#	workflow.connect(nki_nodes.nuisance_stat,'out_file',datasink,'@nuisance_stat')

	
	workflow.connect(nki_nodes.nuisance_calc_r, 'out_file',nki_nodes.Saver_nuisance_calc,'in_file')
	workflow.connect(nki_nodes.func_mean_r,'out_file',nki_nodes.Saver_nuisance_calc,'pp')
	workflow.connect(nki_nodes.nuisance_calc_r, 'out_file',nki_nodes.lifeSaver_nuisance_calc,'in_file')
	workflow.connect(nki_nodes.Saver_nuisance_calc, 'out_file',nki_nodes.lifeSaver_nuisance_calc,'format_string')
	workflow.connect(nki_nodes.lifeSaver_nuisance_calc,'out_file',datasink,'@nuisance_calc')
#	workflow.connect(nki_nodes.nuisance_calc,'out_file',datasink,'@nuisance_calc')

	workflow.connect(nki_nodes.nuisance_warp_r, 'out_file',nki_nodes.Saver_nuisance_warp,'in_file')
	workflow.connect(nki_nodes.func_mean_r,'out_file',nki_nodes.Saver_nuisance_warp,'pp')
	workflow.connect(nki_nodes.nuisance_warp_r, 'out_file',nki_nodes.lifeSaver_nuisance_warp,'in_file')
	workflow.connect(nki_nodes.Saver_nuisance_warp, 'out_file',nki_nodes.lifeSaver_nuisance_warp,'format_string')
	workflow.connect(nki_nodes.lifeSaver_nuisance_warp,'out_file',datasink,'@nuisance_warp')
#	workflow.connect(nki_nodes.nuisance_warp,'out_file',datasink,'@nuisance_warp')

	workflow.connect(nki_nodes.nuisance_warp_1_r, 'out_file',nki_nodes.Saver_nuisance_warp_1,'in_file')
	workflow.connect(nki_nodes.func_mean_r,'out_file',nki_nodes.Saver_nuisance_warp_1,'pp')
	workflow.connect(nki_nodes.nuisance_warp_1_r, 'out_file',nki_nodes.lifeSaver_nuisance_warp_1,'in_file')
	workflow.connect(nki_nodes.Saver_nuisance_warp_1, 'out_file',nki_nodes.lifeSaver_nuisance_warp_1,'format_string')
	workflow.connect(nki_nodes.lifeSaver_nuisance_warp_1,'out_file',datasink,'@nuisance_warp_1')


def makeOutputConnections(analysisdirectory):


	datasinkAnat = getSink(analysisdirectory)
	makeOutputConnectionsAnat(datasinkAnat)

	datasink = getSinkOther()
	makeOutputConnectionsFunc(datasink)

	makeOutputConnectionsReg(datasink,datasinkAnat)
	makeOutputConnectionsSeg(datasink,datasinkAnat)
	makeOutputConnectionsNuisance(datasink)

def rename_func_outputs():


	global renamer
	global zeroSM

	workflow.connect(nki_nodes.func_calc, 'out_file',nki_nodes.func_calc_o,'in_file')
	workflow.connect(renamer,'func_calc_out_file',nki_nodes.func_calc_o,'name')
	workflow.connect(nki_nodes.func_calc, 'out_file',nki_nodes.func_calc_r,'in_file')
	workflow.connect(nki_nodes.func_calc_o,'out_file',nki_nodes.func_calc_r,'format_string')
#	workflow.connect( nki_nodes.func_calc , 'out_file', datasink , '@rest_dr' )

	workflow.connect(nki_nodes.func_refit, 'out_file',nki_nodes.func_refit_o,'in_file')
	workflow.connect(renamer,'func_refit_out_file',nki_nodes.func_refit_o,'name')
	workflow.connect(nki_nodes.func_refit, 'out_file',nki_nodes.func_refit_r,'in_file')
	workflow.connect(nki_nodes.func_refit_o,'out_file',nki_nodes.func_refit_r,'format_string')
#	workflow.connect( nki_nodes.func_refit , 'out_file' , datasink, '@rest_dr_1' )

	workflow.connect(nki_nodes.func_reorient, 'out_file',nki_nodes.func_reorient_o,'in_file')
	workflow.connect(renamer,'func_reorient_out_file',nki_nodes.func_reorient_o,'name')
	workflow.connect(nki_nodes.func_reorient, 'out_file',nki_nodes.func_reorient_r,'in_file')
	workflow.connect(nki_nodes.func_reorient_o,'out_file',nki_nodes.func_reorient_r,'format_string')
#	workflow.connect( nki_nodes.func_reorient , 'out_file' , datasink, '@rest_ro' )

	workflow.connect(nki_nodes.func_tstat, 'out_file',nki_nodes.func_tstat_o,'in_file')
	workflow.connect(renamer,'func_tstat_out_file',nki_nodes.func_tstat_o,'name')
	workflow.connect(nki_nodes.func_tstat, 'out_file',nki_nodes.func_tstat_r,'in_file')
	workflow.connect(nki_nodes.func_tstat_o,'out_file',nki_nodes.func_tstat_r,'format_string')
#	workflow.connect( nki_nodes.func_tstat , 'out_file' , datasink, '@rest_ro_mean' )

	workflow.connect(nki_nodes.func_volreg_1, 'oned_file',nki_nodes.func_volrego_o,'in_file')
	workflow.connect(renamer,'func_volreg_oned_file',nki_nodes.func_volrego_o,'name')
	workflow.connect(nki_nodes.func_volreg_1, 'oned_file',nki_nodes.func_volrego_r,'in_file')
	workflow.connect(nki_nodes.func_volrego_o,'out_file',nki_nodes.func_volrego_r,'format_string')
#	workflow.connect( nki_nodes.func_volreg , 'oned_file', datasink, '@rest_mc_1D')

	workflow.connect(nki_nodes.func_volreg_1, 'out_file',nki_nodes.func_volreg_o,'in_file')
	workflow.connect(renamer,'func_volreg_out_file',nki_nodes.func_volreg_o,'name')
	workflow.connect(nki_nodes.func_volreg_1, 'out_file',nki_nodes.func_volreg_r,'in_file')
	workflow.connect(nki_nodes.func_volreg_o,'out_file',nki_nodes.func_volreg_r,'format_string')
#	workflow.connect( nki_nodes.func_volreg , 'out_file', datasink, '@rest_mc')

	workflow.connect(nki_nodes.func_automask, 'out_file',nki_nodes.func_automask_o,'in_file')
	workflow.connect(renamer,'func_automask_out_file',nki_nodes.func_automask_o,'name')
	workflow.connect(nki_nodes.func_automask, 'out_file',nki_nodes.func_automask_r,'in_file')
	workflow.connect(nki_nodes.func_automask_o,'out_file',nki_nodes.func_automask_r,'format_string')
#	workflow.connect( nki_nodes.func_automask,'out_file', datasink,'@rest_mask')

	workflow.connect(nki_nodes.func_calcR, 'out_file',nki_nodes.func_calcR_o,'in_file')
	workflow.connect(renamer,'func_calcR_out_file',nki_nodes.func_calcR_o,'name')
	workflow.connect(nki_nodes.func_calcR, 'out_file',nki_nodes.func_calcR_r,'in_file')
	workflow.connect(nki_nodes.func_calcR_o,'out_file',nki_nodes.func_calcR_r,'format_string')
#	workflow.connect( nki_nodes.func_calcR,'out_file', datasink,'@rest_ss')

#	workflow.connect( nki_nodes.func_despike, 'out_file', datasink, '@rest_ds')
	workflow.connect(nki_nodes.func_mean, 'out_file',nki_nodes.func_mean_o,'in_file')
	workflow.connect(renamer,'func_mean_out_file',nki_nodes.func_mean_o,'name')
	workflow.connect(nki_nodes.func_mean, 'out_file',nki_nodes.func_mean_r,'in_file')
	workflow.connect(nki_nodes.func_mean_o,'out_file',nki_nodes.func_mean_r,'format_string')
#	workflow.connect( nki_nodes.func_smooth,'out_file' , datasink, '@rest_sm')

	workflow.connect(nki_nodes.func_scale, 'out_file',nki_nodes.func_scale_o,'in_file')
	workflow.connect(renamer,'func_scale_out_file',nki_nodes.func_scale_o,'name')
	workflow.connect(nki_nodes.func_scale, 'out_file',nki_nodes.func_scale_r,'in_file')
	workflow.connect(nki_nodes.func_scale_o,'out_file',nki_nodes.func_scale_r,'format_string')
#	workflow.connect( nki_nodes.func_scale, 'out_file', datasink, '@rest_gms')


	workflow.connect(nki_nodes.func_mask, 'out_file',nki_nodes.func_mask_o,'in_file')
	workflow.connect(renamer,'func_mask_out_file',nki_nodes.func_mask_o,'name')
	workflow.connect(nki_nodes.func_mask, 'out_file',nki_nodes.func_mask_r,'in_file')
	workflow.connect(nki_nodes.func_mask_o,'out_file',nki_nodes.func_mask_r,'format_string')
#	workflow.connect( nki_nodes.func_mask, 'out_file' , datasink, '@rest_pp_mask')

def rename_reg_outputs():

	workflow.connect(nki_nodes.reg_flirt, 'out_file',nki_nodes.reg_flirt_o,'in_file')
	workflow.connect(renamer,'reg_flirt_out_file',nki_nodes.reg_flirt_o,'name')
	workflow.connect(nki_nodes.reg_flirt, 'out_file',nki_nodes.reg_flirt_r,'in_file')
	workflow.connect(nki_nodes.reg_flirt_o,'out_file',nki_nodes.reg_flirt_r,'format_string')
	#workflow.connect(nki_nodes.reg_flirt, 'out_file',datasink,'reg.@example_func2highres')
	
	workflow.connect(nki_nodes.reg_flirt, 'out_matrix_file',nki_nodes.reg_flirto_o,'in_file')
	workflow.connect(renamer,'reg_flirt_out_matrix_file',nki_nodes.reg_flirto_o,'name')
	workflow.connect(nki_nodes.reg_flirt, 'out_matrix_file',nki_nodes.reg_flirto_r,'in_file')
	workflow.connect(nki_nodes.reg_flirto_o,'out_file',nki_nodes.reg_flirto_r,'format_string')
	#workflow.connect(nki_nodes.reg_flirt, 'out_matrix_file',datasink,'reg.@example_func2highresmat')
	
	workflow.connect(nki_nodes.reg_xfm1, 'out_file',nki_nodes.reg_xfm1_o,'in_file')
	workflow.connect(renamer,'reg_xfm1_out_file',nki_nodes.reg_xfm1_o,'name')
	workflow.connect(nki_nodes.reg_xfm1, 'out_file',nki_nodes.reg_xfm1_r,'in_file')
	workflow.connect(nki_nodes.reg_xfm1_o,'out_file',nki_nodes.reg_xfm1_r,'format_string')
	#workflow.connect(nki_nodes.reg_xfm1, 'out_file',datasink,'reg.@highres2example_funcmat')
	
	#workflow.connect(nki_nodes.reg_flirt1, 'out_file',datasinkAnat, 'reg.@highres2standard')
	#workflow.connect(nki_nodes.reg_flirt1, 'out_matrix_file',datasinkAnat, 'reg.@highres2standardmat')
	#workflow.connect(nki_nodes.reg_xfm2, 'out_file',datasinkAnat,'reg.@standard2highresmat')
	
	workflow.connect(nki_nodes.reg_xfm3, 'out_file',nki_nodes.reg_xfm3_o,'in_file')
	workflow.connect(renamer,'reg_xfm3_out_file',nki_nodes.reg_xfm3_o,'name')
	workflow.connect(nki_nodes.reg_xfm3, 'out_file',nki_nodes.reg_xfm3_r,'in_file')
	workflow.connect(nki_nodes.reg_xfm3_o,'out_file',nki_nodes.reg_xfm3_r,'format_string')
	#workflow.connect(nki_nodes.reg_xfm3, 'out_file',datasink,'reg.@example_func2standardmat')
	
	workflow.connect(nki_nodes.reg_flirt2, 'out_file',nki_nodes.reg_flirt2_o,'in_file')
	workflow.connect(renamer,'reg_flirt2_out_file',nki_nodes.reg_flirt2_o,'name')
	workflow.connect(nki_nodes.reg_flirt2, 'out_file',nki_nodes.reg_flirt2_r,'in_file')
	workflow.connect(nki_nodes.reg_flirt2_o,'out_file',nki_nodes.reg_flirt2_r,'format_string')
	#workflow.connect(nki_nodes.reg_flirt2, 'out_file',datasink,'reg.@example_func2standard')
	
	workflow.connect(nki_nodes.reg_xfm4, 'out_file',nki_nodes.reg_xfm4_o,'in_file')
	workflow.connect(renamer,'reg_xfm4_out_file',nki_nodes.reg_xfm4_o,'name')
	workflow.connect(nki_nodes.reg_xfm4, 'out_file',nki_nodes.reg_xfm4_r,'in_file')
	workflow.connect(nki_nodes.reg_xfm4_o,'out_file',nki_nodes.reg_xfm4_r,'format_string')
	#workflow.connect(nki_nodes.reg_xfm4, 'out_file',datasink,'reg.@standard2example_funcmat')
	
	#workflow.connect(nki_nodes.reg_fnt , 'warped_file',datasinkAnat,'reg.@highres2standard_NL')
	#workflow.connect(nki_nodes.reg_fnt, 'jacobian_file',datasinkAnat,'reg.@highres2standard_jac')
	#workflow.connect(nki_nodes.reg_fnt, 'fieldcoeff_file',datasinkAnat,'reg.@highres2standard_warp')
	workflow.connect(nki_nodes.reg_warp, 'out_file',nki_nodes.reg_warp_o,'in_file')
	workflow.connect(renamer,'reg_warp_out_file',nki_nodes.reg_warp_o,'name')
	workflow.connect(nki_nodes.reg_warp, 'out_file',nki_nodes.reg_warp_r,'in_file')
	workflow.connect(nki_nodes.reg_warp_o,'out_file',nki_nodes.reg_warp_r,'format_string')
	#workflow.connect(nki_nodes.reg_warp,'out_file',datasinkAnat,'reg.@example_func2standard_NL')

def rename_seg_outputs():

	global zeroSM
	#workflow.connect( nki_nodes.seg_segment, 'probability_maps',datasinkAnat,'segment.@seg_segment' )
	
	workflow.connect(nki_nodes.seg_flirt, 'out_file',nki_nodes.seg_flirt_o,'in_file')
	workflow.connect(renamer,'seg_flirt_out_file',nki_nodes.seg_flirt_o,'name')
	workflow.connect(nki_nodes.seg_flirt, 'out_file',nki_nodes.seg_flirt_r,'in_file')
	workflow.connect(nki_nodes.seg_flirt_o,'out_file',nki_nodes.seg_flirt_r,'format_string')
	#workflow.connect( nki_nodes.seg_flirt , 'out_file', datasink, 'segment.@seg_flirt' )

	workflow.connect(nki_nodes.seg_smooth, 'out_file',nki_nodes.seg_smooth_o,'in_file')
	workflow.connect(renamer,'seg_smooth_out_file',nki_nodes.seg_smooth_o,'name')
	workflow.connect(nki_nodes.seg_smooth, 'out_file',nki_nodes.seg_smooth_r,'in_file')
	workflow.connect(nki_nodes.seg_smooth_o,'out_file',nki_nodes.seg_smooth_r,'format_string')
#	workflow.connect( nki_nodes.seg_smooth, 'out_file', datasink, 'segment.@seg_smooth')

	workflow.connect(nki_nodes.seg_flirt1, 'out_file',nki_nodes.seg_flirt1_o,'in_file')
	workflow.connect(renamer,'seg_flirt1_out_file',nki_nodes.seg_flirt1_o,'name')
	workflow.connect(nki_nodes.seg_flirt1, 'out_file',nki_nodes.seg_flirt1_r,'in_file')
	workflow.connect(nki_nodes.seg_flirt1_o,'out_file',nki_nodes.seg_flirt1_r,'format_string')
#	workflow.connect( nki_nodes.seg_flirt1, 'out_file', datasink, 'segment.@seg_flirt1')

	workflow.connect(nki_nodes.seg_smooth1, 'out_file',nki_nodes.seg_smooth1_o,'in_file')
	workflow.connect(renamer,'seg_smooth1_out_file',nki_nodes.seg_smooth1_o,'name')
	workflow.connect(nki_nodes.seg_smooth1, 'out_file',nki_nodes.seg_smooth1_r,'in_file')
	workflow.connect(nki_nodes.seg_smooth1_o,'out_file',nki_nodes.seg_smooth1_r,'format_string')
#	workflow.connect( nki_nodes.seg_smooth1,'out_file', datasink,'segment.@seg_smooth1')

	workflow.connect(nki_nodes.seg_flirt2, 'out_file',nki_nodes.seg_flirt2_o,'in_file')
	workflow.connect(renamer,'seg_flirt2_out_file',nki_nodes.seg_flirt2_o,'name')
	workflow.connect(nki_nodes.seg_flirt2, 'out_file',nki_nodes.seg_flirt2_r,'in_file')
	workflow.connect(nki_nodes.seg_flirt2_o,'out_file',nki_nodes.seg_flirt2_r,'format_string')
#	workflow.connect( nki_nodes.seg_flirt2, 'out_file', datasink, 'segment.@seg_flirt2')

	workflow.connect(nki_nodes.seg_thresh, 'out_file',nki_nodes.seg_thresh_o,'in_file')
	workflow.connect(renamer,'seg_thresh_out_file',nki_nodes.seg_thresh_o,'name')
	workflow.connect(nki_nodes.seg_thresh, 'out_file',nki_nodes.seg_thresh_r,'in_file')
	workflow.connect(nki_nodes.seg_thresh_o,'out_file',nki_nodes.seg_thresh_r,'format_string')
#	workflow.connect( nki_nodes.seg_thresh,'out_file', datasink,'segment.@seg_thresh')

	workflow.connect(nki_nodes.seg_mask, 'out_file',nki_nodes.seg_mask_o,'in_file')
	workflow.connect(renamer,'seg_mask_out_file',nki_nodes.seg_mask_o,'name')
	workflow.connect(nki_nodes.seg_mask, 'out_file',nki_nodes.seg_mask_r,'in_file')
	workflow.connect(nki_nodes.seg_mask_o,'out_file',nki_nodes.seg_mask_r,'format_string')
#	workflow.connect( nki_nodes.seg_mask,'out_file', datasink,'segment.@seg_mask')

	workflow.connect(nki_nodes.seg_copy, 'out_file',nki_nodes.seg_copy_o,'in_file')
	workflow.connect(renamer,'seg_copy_out_file',nki_nodes.seg_copy_o,'name')
	workflow.connect(nki_nodes.seg_copy, 'out_file',nki_nodes.seg_copy_r,'in_file')
	workflow.connect(nki_nodes.seg_copy_o,'out_file',nki_nodes.seg_copy_r,'format_string')
#	workflow.connect( nki_nodes.seg_copy,'out_file',  datasink,'segment.@seg_copy')

	workflow.connect(nki_nodes.seg_flirt3, 'out_file',nki_nodes.seg_flirt3_o,'in_file')
	workflow.connect(renamer,'seg_flirt3_out_file',nki_nodes.seg_flirt3_o,'name')
	workflow.connect(nki_nodes.seg_flirt3, 'out_file',nki_nodes.seg_flirt3_r,'in_file')
	workflow.connect(nki_nodes.seg_flirt3_o,'out_file',nki_nodes.seg_flirt3_r,'format_string')
#	workflow.connect( nki_nodes.seg_flirt3, 'out_file', datasink, 'segment.@seg_flirt3')

	workflow.connect(nki_nodes.seg_smooth2, 'out_file',nki_nodes.seg_smooth2_o,'in_file')
	workflow.connect(renamer,'seg_smooth2_out_file',nki_nodes.seg_smooth2_o,'name')
	workflow.connect(nki_nodes.seg_smooth2, 'out_file',nki_nodes.seg_smooth2_r,'in_file')
	workflow.connect(nki_nodes.seg_smooth2_o,'out_file',nki_nodes.seg_smooth2_r,'format_string')
#	workflow.connect( nki_nodes.seg_smooth2,'out_file', datasink, 'segment.@seg_smooth2')

	workflow.connect(nki_nodes.seg_flirt4, 'out_file',nki_nodes.seg_flirt4_o,'in_file')
	workflow.connect(renamer,'seg_flirt4_out_file',nki_nodes.seg_flirt4_o,'name')
	workflow.connect(nki_nodes.seg_flirt4, 'out_file',nki_nodes.seg_flirt4_r,'in_file')
	workflow.connect(nki_nodes.seg_flirt4_o,'out_file',nki_nodes.seg_flirt4_r,'format_string')
#	workflow.connect( nki_nodes.seg_flirt4,'out_file', datasink, 'segment.@seg_flirt4')

	workflow.connect(nki_nodes.seg_prior1, 'out_file',nki_nodes.seg_prior1_o,'in_file')
	workflow.connect(renamer,'seg_prior1_out_file',nki_nodes.seg_prior1_o,'name')
	workflow.connect(nki_nodes.seg_prior1, 'out_file',nki_nodes.seg_prior1_r,'in_file')
	workflow.connect(nki_nodes.seg_prior1_o,'out_file',nki_nodes.seg_prior1_r,'format_string')
#	workflow.connect( nki_nodes.seg_prior1,'out_file', datasink, 'segment.@seg_prior1')

	workflow.connect(nki_nodes.seg_flirt5, 'out_file',nki_nodes.seg_flirt5_o,'in_file')
	workflow.connect(renamer,'seg_flirt5_out_file',nki_nodes.seg_flirt5_o,'name')
	workflow.connect(nki_nodes.seg_flirt5, 'out_file',nki_nodes.seg_flirt5_r,'in_file')
	workflow.connect(nki_nodes.seg_flirt5_o,'out_file',nki_nodes.seg_flirt5_r,'format_string')
#	workflow.connect( nki_nodes.seg_flirt5,'out_file', datasink,'segment.@seg_flirt5')

	workflow.connect(nki_nodes.seg_thresh1, 'out_file',nki_nodes.seg_thresh1_o,'in_file')
	workflow.connect(renamer,'seg_thresh1_out_file',nki_nodes.seg_thresh1_o,'name')
	workflow.connect(nki_nodes.seg_thresh1, 'out_file',nki_nodes.seg_thresh1_r,'in_file')
	workflow.connect(nki_nodes.seg_thresh1_o,'out_file',nki_nodes.seg_thresh1_r,'format_string')
#	workflow.connect( nki_nodes.seg_thresh1,'out_file', datasink,'segment.@seg_thresh1')
	
	workflow.connect(nki_nodes.seg_mask1, 'out_file',nki_nodes.seg_mask1_o,'in_file')
	workflow.connect(renamer,'seg_mask1_out_file',nki_nodes.seg_mask1_o,'name')
	workflow.connect(nki_nodes.seg_mask1, 'out_file',nki_nodes.seg_mask1_r,'in_file')
	workflow.connect(nki_nodes.seg_mask1_o,'out_file',nki_nodes.seg_mask1_r,'format_string')
#	workflow.connect( nki_nodes.seg_mask1,'out_file', datasink,'segment.@seg_mask1')

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

def rename_anat_outputs():

#	workflow.connect(nki_nodes.anat_refit, 'out_file',nki_nodes.anat_refit_o,'in_file')
#	workflow.connect(renamer,'anat_refit_out_file',nki_nodes.anat_refit_o,'name')
#	workflow.connect(nki_nodes.anat_refit, 'out_file',nki_nodes.anat_refit_r,'in_file')
#	workflow.connect(nki_nodes.anat_refit_o,'out_file',nki_nodes.anat_refit_r,'format_string')
	#workflow.connect(nki_nodes.anat_refit,'out_file',datasinkAnat,'@mprage')

#	workflow.connect(nki_nodes.anat_reorient, 'out_file',nki_nodes.anat_reorient_o,'in_file')
#	workflow.connect(renamer,'anat_reorient_out_file',nki_nodes.anat_reorient_o,'name')
#	workflow.connect(nki_nodes.anat_reorient, 'out_file',nki_nodes.anat_reorient_r,'in_file')
#	workflow.connect(nki_nodes.anat_reorient_o,'out_file',nki_nodes.anat_reorient_r,'format_string')
	#workflow.connect(nki_nodes.anat_reorient,'out_file',datasinkAnat,'@mprage_RPI')

	workflow.connect(nki_nodes.anat_skullstrip, 'out_file',nki_nodes.anat_skullstrip_o,'in_file')
	workflow.connect(renamer,'anat_skullstrip_out_file',nki_nodes.anat_skullstrip_o,'name')
	workflow.connect(nki_nodes.anat_skullstrip, 'out_file',nki_nodes.anat_skullstrip_r,'in_file')
	workflow.connect(nki_nodes.anat_skullstrip_o,'out_file',nki_nodes.anat_skullstrip_r,'format_string')
	#workflow.connect(nki_nodes.anat_skullstrip,'out_file', datasinkAnat,'@mprage_surf')
	
#	workflow.connect(nki_nodes.anat_calc, 'out_file',nki_nodes.anat_calc_o,'in_file')
#	workflow.connect(renamer,'anat_calc_out_file',nki_nodes.anat_calc_o,'name')
#	workflow.connect(nki_nodes.anat_calc, 'out_file',nki_nodes.anat_calc_r,'in_file')
#	workflow.connect(nki_nodes.anat_calc_o,'out_file',nki_nodes.anat_calc_r,'format_string')
	#workflow.connect(nki_nodes.anat_calc,'out_file',datasinkAnat,'@mprage_brain')

	workflow.connect(nki_nodes.reg_flirt1, 'out_matrix_file',nki_nodes.reg_flirt1_o,'in_file')
	workflow.connect(renamer,'reg_flirt1_out_matrix_file',nki_nodes.reg_flirt1_o,'name')
	workflow.connect(nki_nodes.reg_flirt1, 'out_matrix_file',nki_nodes.reg_flirt1_r,'in_file')
	workflow.connect(nki_nodes.reg_flirt1_o,'out_file',nki_nodes.reg_flirt1_r,'format_string')
	#workflow.connect(nki_nodes.reg_flirt1, 'out_matrix_file',datasinkAnat, 'reg.@highres2standardmat')
	
	workflow.connect(nki_nodes.reg_flirt1, 'out_file',nki_nodes.reg_flirt1o_o,'in_file')
	workflow.connect(renamer,'reg_flirt1_out_file',nki_nodes.reg_flirt1o_o,'name')
	workflow.connect(nki_nodes.reg_flirt1, 'out_file',nki_nodes.reg_flirt1o_r,'in_file')
	workflow.connect(nki_nodes.reg_flirt1o_o,'out_file',nki_nodes.reg_flirt1o_r,'format_string')


	workflow.connect(nki_nodes.reg_xfm2, 'out_file',nki_nodes.reg_xfm2_o,'in_file')
	workflow.connect(renamer,'reg_xfm2_out_file',nki_nodes.reg_xfm2_o,'name')
	workflow.connect(nki_nodes.reg_xfm2, 'out_file',nki_nodes.reg_xfm2_r,'in_file')
	workflow.connect(nki_nodes.reg_xfm2_o,'out_file',nki_nodes.reg_xfm2_r,'format_string')
	#workflow.connect(nki_nodes.reg_xfm2, 'out_file',datasinkAnat,'reg.@standard2highresmat')


	workflow.connect(nki_nodes.reg_fnt, 'warped_file',nki_nodes.reg_fnt_o,'in_file')
	workflow.connect(renamer,'reg_fnt_warped_file',nki_nodes.reg_fnt_o,'name')
	workflow.connect(nki_nodes.reg_fnt, 'warped_file',nki_nodes.reg_fnt_r,'in_file')
	workflow.connect(nki_nodes.reg_fnt_o,'out_file',nki_nodes.reg_fnt_r,'format_string')
	#workflow.connect(nki_nodes.reg_fnt , 'warped_file',datasinkAnat,'reg.@highres2standard_NL')


	workflow.connect(nki_nodes.reg_fnt, 'jacobian_file',nki_nodes.reg_fntj_o,'in_file')
	workflow.connect(renamer,'reg_fnt_jacobian_file',nki_nodes.reg_fntj_o,'name')
	workflow.connect(nki_nodes.reg_fnt, 'jacobian_file',nki_nodes.reg_fntj_r,'in_file')
	workflow.connect(nki_nodes.reg_fntj_o,'out_file',nki_nodes.reg_fntj_r,'format_string')
	#workflow.connect(nki_nodes.reg_fnt, 'jacobian_file',datasinkAnat,'reg.@highres2standard_jac')

	workflow.connect(nki_nodes.reg_fnt, 'fieldcoeff_file',nki_nodes.reg_fntf_o,'in_file')
	workflow.connect(renamer,'reg_fnt_fieldcoeff_file',nki_nodes.reg_fntf_o,'name')
	workflow.connect(nki_nodes.reg_fnt, 'fieldcoeff_file',nki_nodes.reg_fntf_r,'in_file')
	workflow.connect(nki_nodes.reg_fntf_o,'out_file',nki_nodes.reg_fntf_r,'format_string')
	#workflow.connect(nki_nodes.reg_fnt, 'fieldcoeff_file',datasinkAnat,'reg.@highres2standard_warp')

	#workflow.connect(nki_nodes.reg_warp,'out_file',datasinkAnat,'reg.@example_func2standard_NL')
	#workflow.connect( nki_nodes.seg_segment, 'probability_maps',datasinkAnat,'segment.@seg_segment' )



def renameOutputs():

	rename_anat_outputs()
	rename_func_outputs()
	rename_reg_outputs()
	rename_seg_outputs()
	rename_nuisance_outputs()

def getNodeOutputNames():

	names = ['reg_warp_out_file','reg_fnt_jacobian_file','reg_fnt_fieldcoeff_file','reg_fnt_warped_file','reg_flirt1_out_file','reg_flirt1_out_matrix_file','reg_xfm2_out_file','anat_refit_out_file','anat_reorient_out_file','anat_skullstrip_out_file','anat_calc_out_file','func_calc_out_file','func_refit_out_file','func_reorient_out_file','func_tstat_out_file','func_volreg_oned_file','func_volreg_out_file','func_automask_out_file','func_calcR_out_file','func_mean_out_file','func_calcI_out_file','func_despike_out_file','func_smooth_out_file','func_scale_out_file','func_filter_out_file','func_detrenda_out_file','func_detrendb_out_file','func_detrendc_out_file','func_mask_out_file','reg_flirt_out_file','reg_flirt_out_matrix_file','reg_xfm1_out_file','reg_xfm3_out_file','reg_flirt2_out_file','reg_xfm4_out_file','seg_flirt_out_file','seg_smooth_out_file','seg_flirt1_out_file','seg_smooth1_out_file','seg_flirt2_out_file','seg_thresh_out_file','seg_mask_out_file','seg_copy_out_file','seg_flirt3_out_file','seg_smooth2_out_file','seg_flirt4_out_file','seg_prior1_out_file','seg_flirt5_out_file','seg_thresh1_out_file','seg_mask1_out_file','nuisance_globalE_out_file','nuisance_csf_out_file','nuisance_wm_out_file','nuisance_featM_design_file','nuisance_fgls_residual4d','nuisance_stat_out_file','nuisance_calc_out_file','nuisance_warp_out_file','alff_cp_out_file','alff_mean_out_file','alff_pspec_out_file','alff_sum_out_file','alff_falff1_out_file','alff_Z_falff_out_file','alff_Z_alff_out_file','alff_warp_alff_out_file','alff_warp_falff_out_file','corrs','z_trans','register','nuisance_erosion_csf1_out_file','nuisance_erosion_wm1_out_file','nuisance_compcor_out_file','nuisance_MedianAngle_out_file','nuisance_warp_1_out_file']


	return names

def getRenamer():

	global renamer
	global rest_name
	global anat_name

	nodeOutputNames = getNodeOutputNames()
	renamer = pe.Node(interface=util.IdentityInterface(fields= nodeOutputNames), name="renamer")

	#renamer for anatomical preprocessing nodes
	renamer.inputs.anat_refit_out_file = anat_name + '.nii.gz'
	renamer.inputs.anat_reorient_out_file = anat_name + '_RPI.nii.gz'
	renamer.inputs.anat_skullstrip_out_file = anat_name + '_surf.nii.gz'
	renamer.inputs.anat_calc_out_file = anat_name + '_brain.nii.gz'
	renamer.inputs.reg_fnt_jacobian_file = 'highres2standard_jac.nii.gz'
	renamer.inputs.reg_fnt_fieldcoeff_file = 'highres2standard_warp.nii.gz'
	renamer.inputs.reg_fnt_warped_file = 'highres2standard_NL.nii.gz'
	renamer.inputs.reg_flirt1_out_matrix_file = 'highres2standard.mat'
	renamer.inputs.reg_flirt1_out_file = 'highres2standard.nii.gz'
	renamer.inputs.reg_xfm2_out_file = 'standard2highres.mat'
	#renamer for functional preprocessing nodes
	renamer.inputs.func_calc_out_file = rest_name +'_dr.nii.gz'
	renamer.inputs.func_refit_out_file = rest_name +'_dr.nii.gz'
	renamer.inputs.func_reorient_out_file = rest_name +'_ro.nii.gz'
	renamer.inputs.func_tstat_out_file = rest_name +'_ro_mean.nii.gz'
	renamer.inputs.func_volreg_oned_file = rest_name +'_mc.1D'
	renamer.inputs.func_volreg_out_file = rest_name +'_mc.nii.gz'
	renamer.inputs.func_automask_out_file = rest_name + '_mask.nii.gz'
	renamer.inputs.func_calcR_out_file = rest_name + '_ss.nii.gz'
	renamer.inputs.func_mean_out_file = 'example_func.nii.gz'
	renamer.inputs.func_scale_out_file = rest_name + '_pp.nii.gz'
	renamer.inputs.func_filter_out_file = rest_name + '_filt.nii.gz'
	renamer.inputs.func_mask_out_file = rest_name + '_pp_mask.nii.gz'

	#renamer for registration
	renamer.inputs.reg_flirt_out_file = 'example_func2highres.nii.gz'
	renamer.inputs.reg_flirt_out_matrix_file = 'example_func2highres.mat'
	renamer.inputs.reg_xfm1_out_file = 'highres2example_func.mat'
	renamer.inputs.reg_xfm3_out_file = 'example_func2standard.mat'
	renamer.inputs.reg_flirt2_out_file = 'example_func2standard.nii.gz'
	renamer.inputs.reg_xfm4_out_file = 'standard2example_func.mat'
	renamer.inputs.reg_warp_out_file = 'example_func2standard_NL.nii.gz'

	#renamer for segmentation
	renamer.inputs.seg_flirt_out_file = 'csf2func.nii.gz'
	renamer.inputs.seg_smooth_out_file = 'csf_sm.nii.gz'
	renamer.inputs.seg_flirt1_out_file = 'csf2standard.nii.gz'
	renamer.inputs.seg_smooth1_out_file = 'csf_masked.nii.gz'
	renamer.inputs.seg_flirt2_out_file = 'csf_native.nii.gz'
	renamer.inputs.seg_thresh_out_file = 'csf_bin.nii.gz'
	renamer.inputs.seg_mask_out_file = 'csf_mask.nii.gz'
	renamer.inputs.seg_copy_out_file = 'global_mask.nii.gz'
	renamer.inputs.seg_flirt3_out_file = 'wm2func.nii.gz'
	renamer.inputs.seg_smooth2_out_file = 'wm_sm.nii.gz'
	renamer.inputs.seg_flirt4_out_file = 'wm2standard.nii.gz'
	renamer.inputs.seg_prior1_out_file = 'wm_masked.nii.gz'
	renamer.inputs.seg_flirt5_out_file = 'wm_native.nii.gz'
	renamer.inputs.seg_thresh1_out_file = 'wm_bin.nii.gz'
	renamer.inputs.seg_mask1_out_file = 'wm_mask.nii.gz'

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

	#renamer for RSFC
	renamer.inputs.corrs = 'corrs' 
	renamer.inputs.z_trans = 'z_trans'
	renamer.inputs.register = 'register'

#def genPowerParams():

	

#def scrubbing():

#	genPowerParams()


def processS(sublist, analysisdirectory):

	from nipype.utils.config import config
	global infosource
	global datasource
	global workflow

	numCores = int(sys.argv[1])
	
	getInfoSource(sublist, analysisdirectory)
	getRenamer()
	gatherData(sublist, analysisdirectory)
	
	wfname =  'fcon1000'
	workflow = pe.Workflow(name=wfname)
	workflow.base_dir = working_dir
	anatpreproc()
	funcpreproc()
	registration()
	segment()
	nuisance()
	#scrubbing()
	renameOutputs()
	makeOutputConnections(analysisdirectory)

	#VMHC()
	workflow.run(plugin='MultiProc', plugin_args={'n_procs' : numCores})
	#workflow.run()
	#workflow.write_graph()

def checkProcessedSubject(dir,subject):


	checkList = ['mprage_brain.nii.gz']

	for file in checkList:
		if ('mprage' in file):
			if os.path.isdir(dir + '/' +subject + '/session_1/anat_1/'):

				fpath_1 = dir + '/' +subject + '/session_1/anat_1/' + file
				if os.path.isfile(fpath_1):
					return 1
				else:
					files = os.listdir(dir + '/' +subject + '/session_1/anat_1/')
					if(len(files) > 1):
						return 1
					else:
						return 0

			elif os.path.isdir(dir + '/' +subject + '/session_2/anat_1/'):	
				fpath_1 = dir + '/' +subject + '/session_2/anat_1/' + file
				if os.path.isfile(fpath_1):
					return 1
				else:
					files = os.listdir(dir + '/' +subject + '/session_2/anat_1/')
					if (len(files) > 1):
						return 1
					else:
						return 0
	return 1

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
				#if ( os.path.isdir(analysisdirectory + '/' + line + '/session_1/rest_1') or os.path.isdir(analysisdirectory + '/' + line + '/session_2/rest_1')):
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
#	l = line1[0]
#	print l
#	processSubjects((l.split(' '))[0],(l.split(' '))[1])

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
	parsermap = {}
	
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
	HP = parsermap['_hp_']
	LP = parsermap['_lp_']
	hp = parsermap['hp']
	lp = parsermap['lp']
	which_regression = parsermap['which_regression']
	zeroSM = parsermap['0smooth']

def main():

	if ( len(sys.argv) < 2 ):
		sys.stderr.write("./fcon_pipeline.py <Number_of_cores>")
	else:
	
		readDirSetup()
		readSubjects()

	

if __name__ == "__main__":

	sys.exit(main())
