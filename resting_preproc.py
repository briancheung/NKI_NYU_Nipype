#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python
import argparse
import e_afni
import sys
import os
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from base import (create_anat_preproc, create_func_preproc,
                    create_reg_preproc, create_seg_preproc,
                    create_alff_preproc, create_ifc_preproc,
                    create_vmhc_preproc, create_mnioutputs,
                    create_timeseries_preproc)

from utils import (create_anat_dataflow, create_func_dataflow,
                    create_alff_dataflow, create_ifc_dataflow,
                    create_vmhc_dataflow, create_parc_dataflow, 
                    create_mask_dataflow)


def getSubjectAndSeedLists(c):

    """
        Read Subject & Seed files to build
        corresponding lists.
    """
    subj_file = c.subj_file
    seed_file = c.seed_file

    reader_subj = open(subj_file, 'r')
    reader_seed = open(seed_file, 'r')

    subj_list = []
    seed_list = []

    for line in reader_subj.readlines():

        line = line.rstrip('\r\n')

        subj_list.append(line)

    for line in reader_seed.readlines():

        line = line.rstrip('\r\n')
        seed_list.append(line)

    return subj_list, seed_list

def getModelandSeedList(c):

    modelist=[]
    seedlist=[]

    try:
        fp = open(c.model_file,'r')
        models = fp.readlines()
        fp.close()

        for model in models:
            model = model.rstrip('\r\n')
            modelist.append(os.path.basename(model))
        
        print 'checking for models ', modelist
                
        f= open(c.seed_file, 'r')
        seeds=f.readlines()
        f.close()
        
        for seed in seeds:
            seed = os.path.basename(seed.rstrip('\r\n'))
            seed=os.path.splitext(os.path.splitext(seed)[0])[0]
            seedlist.append(seed)
            
        print 'checking for seeds ',seedlist
    except:
        raise
    
    return  modelist, seedlist

def get_seed_list(seed_file):

    f = open(seed_file, 'r')

    seed_list = []
    seeds = f.readlines()
    for seed in seeds:
        seed = seed.rstrip('\r\n')
        seed_list.append(seed)

    return seed_list



def get_workflow(wf_name, c):

    preproc = None
    """ 
        setup standard file paths
    """
    print 'inside get_wf ', '-->'+wf_name+'<--'
    prior_path = os.path.join(c.prior_dir, c.standard_res)
    PRIOR_CSF = os.path.join(prior_path, 'avg152T1_csf_bin.nii.gz')
    PRIOR_GRAY = os.path.join(prior_path, 'avg152T1_csf_bin.nii.gz')
    PRIOR_WHITE = os.path.join(prior_path, 'avg152T1_white_bin.nii.gz')
    standard_res_brain = os.path.join(c.FSLDIR, 'data/standard/MNI152_T1_%s_brain.nii.gz' % (c.standard_res))
    standard = os.path.join(c.FSLDIR, 'data/standard/MNI152_T1_%s.nii.gz' % (c.standard_res))
    standard_brain_mask_dil = os.path.join(c.FSLDIR, 'data/standard/MNI152_T1_%s_brain_mask_dil.nii.gz' % (c.standard_res))
    config_file = os.path.join(c.FSLDIR, 'etc/flirtsch/T1_2_MNI152_%s.cnf' % (c.standard_res))
    brain_symmetric = os.path.join(c.FSLDIR, 'data/standard/MNI152_T1_2mm_brain_symmetric.nii.gz')
    symm_standard = os.path.join(c.FSLDIR, 'data/standard/MNI152_T1_2mm_symmetric.nii.gz')
    twomm_brain_mask_dil = os.path.join(c.FSLDIR, 'data/standard/MNI152_T1_2mm_brain_mask_symmetric_dil.nii.gz')
    config_file_twomm = os.path.join(c.FSLDIR, 'etc/flirtsch/T1_2_MNI152_2mm.cnf')
    identity_matrix=os.path.join(c.FSLDIR,'etc/flirtsch/ident.mat')

    if wf_name.lower() == 'anat':
        preproc = create_anat_preproc()
        return preproc


    if wf_name.lower() == 'func':

        preproc = create_func_preproc()
        preproc.inputs.inputspec.start_idx = c.start_idx
        preproc.inputs.inputspec.stop_idx = c.stop_idx

        return preproc

    if wf_name.lower() == 'seg':

        preproc = create_seg_preproc()
        preproc.inputs.inputspec.PRIOR_CSF  = PRIOR_CSF
        preproc.inputs.inputspec.PRIOR_WHITE  = PRIOR_WHITE
        preproc.inputs.inputspec.PRIOR_GRAY  = PRIOR_GRAY
        preproc.inputs.inputspec.standard_res_brain = standard_res_brain

        return preproc

    if wf_name.lower() == 'reg':

        preproc = create_reg_preproc()
        preproc.inputs.inputspec.standard_res_brain = standard_res_brain
        preproc.inputs.inputspec.standard = standard
        preproc.inputs.inputspec.config_file = config_file
        preproc.inputs.inputspec.standard_brain_mask_dil = standard_brain_mask_dil

        return preproc

    if wf_name.lower() == 'alff':

        preproc = create_alff_preproc()
        preproc.inputs.inputspec.standard = standard
        preproc.inputs.hplp_input.hp = c.alff_hp
        preproc.inputs.hplp_input.lp = c.alff_lp
        preproc.inputs.fwhm_input.fwhm = c.fwhm
        preproc.get_node('hplp_input').iterables = ('hp', c.alff_HP, 'lp', c.alff_LP)
        preproc.get_node('fwhm_input').iterables = ('fwhm', c.fwhm)

        return preproc

    if wf_name.lower() == 'ifc':

        preproc = create_ifc_preproc()
        seeds = get_seed_list(c.seed_file)

        print seeds
        preproc.inputs.seed_list_input.seed_list = seeds
        preproc.get_node('seed_list_input').iterables = ('seed_list', seeds)
        preproc.inputs.fwhm_input.fwhm = c.fwhm
        preproc.get_node('fwhm_input').iterables = ('fwhm', c.fwhm)
        preproc.inputs.inputspec.standard = standard
        return preproc

    if wf_name.lower() == 'vmhc':

        preproc = create_vmhc_preproc()
        preproc.inputs.inputspec.brain_symmetric = brain_symmetric
        preproc.inputs.inputspec.symm_standard = symm_standard
        preproc.inputs.inputspec.twomm_brain_mask_dil = twomm_brain_mask_dil
        preproc.inputs.inputspec.config_file_twomm = config_file_twomm
        preproc.inputs.inputspec.standard = standard
        return preproc

    if wf_name.lower() =='mni':
        
        preproc= create_mnioutputs()
        preproc.inputs.inputspec.reference_file = standard
        return preproc
    
    if wf_name.lower() == 'ts' :
        
        preproc=create_timeseries_preproc(c.derivatives[4], c.derivatives[5], c.derivatives[6])
        preproc.inputs.inputspec.recon_subjects=c.reconSubjectsDirectory
        preproc.inputs.inputspec.standard = standard
        preproc.inputs.inputspec.identity_matrix=identity_matrix
        preproc.inputs.inputspec.unitTSOutputs=c.unitTSOutputs
        preproc.inputs.inputspec.voxelTSOutputs=c.voxelTSOutputs
        preproc.inputs.inputspec.verticesTSOutputs=c.verticesTSOutputs
        return preproc

def prep_workflow(c):

    wfname = 'resting_preproc'
    workflow = pe.Workflow(name=wfname)
    workflow.base_dir = c.working_dir
    workflow.crash_dir = c.crash_dir

    sublist, seed_list = getSubjectAndSeedLists(c)

    """
        BASIC and ALL preprocessing paths implemented below
    """
    if(c.analysis[0] or c.analysis[1]):

        """
            grab the subject data
        """
        flowAnat = create_anat_dataflow('anat_flow', sublist, c.subj_dir, c.anat_name, c.anat_template)
        flowFunc = create_func_dataflow('func_flow', sublist, c.subj_dir, c.rest_name, c.func_template)


        """
            Define workflows and set parameters
        """
        anatpreproc = get_workflow('anat', c)
        funcpreproc = get_workflow('func', c)
        regpreproc = get_workflow('reg', c)
        segpreproc = get_workflow('seg', c)

        """
            Make Connections
        """

        workflow.connect(flowAnat, 'anat', anatpreproc, 'inputspec.anat')
        workflow.connect(flowFunc, 'rest', funcpreproc, 'inputspec.rest')
        workflow.connect(funcpreproc, 'outputspec.example_func', regpreproc, 'inputspec.example_func')
        workflow.connect(anatpreproc, 'outputspec.brain', regpreproc, 'inputspec.brain')
        workflow.connect(anatpreproc, 'outputspec.reorient', regpreproc, 'inputspec.reorient')
        workflow.connect(funcpreproc, 'outputspec.preprocessed_mask', segpreproc, 'inputspec.preprocessed_mask')
        workflow.connect(anatpreproc, 'outputspec.brain', segpreproc, 'inputspec.brain')
        workflow.connect(funcpreproc, 'outputspec.example_func', segpreproc, 'inputspec.example_func')
        workflow.connect(regpreproc, 'outputspec.highres2example_func_mat', segpreproc, 'inputspec.highres2example_func_mat')
        workflow.connect(regpreproc, 'outputspec.stand2highres_warp', segpreproc, 'inputspec.stand2highres_warp')

        """
            Mni Outputs, SurfaceRegistration and timeseries
        """
        mnipreproc= get_workflow('mni', c)
        tspreproc=  get_workflow('ts',c)
        
        print mnipreproc
        print tspreproc
        """
            Connections
        """
        #getting mni outputs after segmentation
        merge = pe.Node(interface=util.Merge(3), name = 'outmni_concatnode')
        workflow.connect(segpreproc,'outputspec.partial_volume_files', merge, 'in1' )
        workflow.connect(segpreproc, 'outputspec.partial_volume_map', merge, 'in2')
        workflow.connect(segpreproc, 'outputspec.mixeltype', merge, 'in3')
        workflow.connect(merge, 'out',mnipreproc, 'inputspec.in_file' )
        workflow.connect(regpreproc,'outputspec.highres2standard_warp', mnipreproc, 'inputspec.warp_file')
        
        #timeseries preproc
        workflow.connect(flowFunc,'rest',tspreproc, 'inputspec.rest' )
        workflow.connect(anatpreproc,'outputspec.brain', tspreproc, 'inputspec.brain')
        workflow.connect(anatpreproc,'outputspec.reorient',tspreproc, 'inputspec.reorient')
        workflow.connect(funcpreproc, 'outputspec.motion_correct', tspreproc, 'inputspec.motion_correct')
        workflow.connect(regpreproc,'outputspec.highres2standard_warp',tspreproc, 'inputspec.warp_file')
        workflow.connect(regpreproc,'outputspec.example_func2highres',tspreproc, 'inputspec.premat')
        
        if c.derivatives[5]:
            pflow=create_parc_dataflow(c.unitDefinitionsDirectory)
            workflow.connect(pflow,'out_file',tspreproc, 'getparc.parcelations')
            
        if c.derivatives[6]:
            mflow=create_mask_dataflow(c.voxelMasksDirectory)
            workflow.connect(mflow, 'out_file', tspreproc,'getmask.masks')

        
        
    """
        ALFF Analysis
    """
    if((not c.analysis[0] and not c.analysis[1]) and c.analysis[4]):

        alffpreproc = get_workflow('alff', c)

        flowAlff, flowAlffWarp = create_alff_dataflow('alff_flow',
                                                      sublist,
                                                      c.subj_dir,
                                                      c.rest_name,
                                                      c.alff_template,
                                                      c.alff_warp_template)

        workflow.connect(flowAlff, 'rest_res', alffpreproc, 'inputspec.rest_res')
        workflow.connect(flowAlff, 'rest_mask', alffpreproc, 'inputspec.rest_mask')
        workflow.connect(flowAlff, 'rest_mask2standard', alffpreproc, 'inputspec.rest_mask2standard')
        workflow.connect(flowAlffWarp, 'premat', alffpreproc, 'inputspec.premat')
        workflow.connect(flowAlffWarp, 'fieldcoeff_file', alffpreproc, 'inputspec.fieldcoeff_file')

    """
        iFC Analysis
    """
    if((not c.analysis[0] and not c.analysis[1]) and c.analysis[5]):

        ifcpreproc = get_workflow('ifc', c)

        flowIfc, flowIfcWarp = create_ifc_dataflow('ifc_flow',
                                                      sublist,
                                                      c.subj_dir,
                                                      c.rest_name,
                                                      c.ifc_template,
                                                      c.ifc_warp_template)

        workflow.connect(flowIfc, 'ref', ifcpreproc, 'inputspec.ref')
        workflow.connect(flowIfcWarp, 'warp', ifcpreproc, 'inputspec.warp')
        workflow.connect(flowIfcWarp, 'postmat', ifcpreproc, 'inputspec.postmat')
        workflow.connect(flowIfcWarp, 'premat', ifcpreproc, 'inputspec.premat')
        workflow.connect(flowIfc, 'rest_res_filt', ifcpreproc, 'inputspec.rest_res_filt')
        workflow.connect(flowIfcWarp, 'fieldcoeff_file', ifcpreproc, 'inputspec.fieldcoeff_file')
#        workflow.connect(flowIfc, 'rest_res2standard', ifcpreproc, 'inputspec.rest_res2standard')
        workflow.connect(flowIfc, 'rest_mask2standard', ifcpreproc, 'inputspec.rest_mask2standard')


    """
        VMHC Analysis
    """
    if((not c.analysis[0] and not c.analysis[1]) and c.analysis[6]):

        vmhcpreproc = get_workflow('vmhc', c)

        flowVmhc_rest_res, flowVmhc_reorient, flowVmhc_example_func2highres_mat = create_vmhc_dataflow('vmhc_flow',
                                                      sublist,
                                                      c.subj_dir,
                                                      c.anat_name,
                                                      c.rest_name,
                                                      c.vmhc_rest_res_template,
                                                      c.vmhc_anat_reorient_template,
                                                      c.vmhc_example_func2highres_mat_template)

        workflow.connect(flowVmhc_rest_res, 'rest_res', vmhcpreproc, 'inputspec.rest_res')
        workflow.connect(flowVmhc_reorient, 'reorient', vmhcpreproc, 'inputspec.reorient')
        workflow.connect(flowVmhc_reorient, 'brain', vmhcpreproc, 'inputspec.brain')
        workflow.connect(flowVmhc_example_func2highres_mat, 'example_func2highres_mat', vmhcpreproc, 'inputspec.example_func2highres_mat')

    if(not c.run_on_grid):
        workflow.run(plugin='MultiProc', plugin_args={'n_procs': c.num_cores})
    else:
        workflow.run(plugin='SGE', plugin_args=dict(qsub_args=c.qsub_args))

def main():

    parser = argparse.ArgumentParser(description="example: \
                        run resting_preproc.py -c config.py")
    parser.add_argument('-c', '--config',
                        dest='config',
                        required=True,
                        help='location of config file'
                        )
    args = parser.parse_args()
    path, fname = os.path.split(os.path.realpath(args.config))
    sys.path.append(path)
    c = __import__(fname.split('.')[0])
    prep_workflow(c)

if __name__ == "__main__":

    sys.exit(main())
