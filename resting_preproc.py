import argparse
import e_afni
import sys
import os
import nipype.pipeline.engine as pe
from base import (create_anat_preproc, create_func_preproc,
                    create_reg_preproc, create_seg_preproc,
                    create_alff_preproc, create_sca_preproc,
                    create_vmhc_preproc, create_scrubbing_preproc,
                    mprage_in_mnioutputs, func_in_mnioutputs,
                    create_filter, create_timeseries_preproc,
                    create_group_analysis)

from utils import (create_anat_func_dataflow, create_seed_dataflow,
                    create_alff_dataflow, create_sca_dataflow,
                    create_vmhc_dataflow, selector_wf,
                    create_parc_dataflow, create_mask_dataflow,
                    create_gp_dataflow, create_datasink)

from base_nuisance import create_nuisance_preproc

from sink import (anat_sink, reg_sink, seg_sink, func_sink,
                  nuisance_sink, scrubbing_sink)

def getSubjectAndSeedLists(c):

    """
        Read Subject & Seed files to build
        corresponding lists.
    """
    
    def get_list(fname):
        flines = open(fname, 'r').readlines()
        return [fline.rstrip('\r\n') for fline in flines]
    
#    subj_file = c.subj_file
#    seed_file = c.seed_file
#    func_session_file = c.func_session_file
#    anat_session_file = c.anat_session_file
#    
#    reader_subj = open(subj_file, 'r')
#    reader_seed = open(seed_file, 'r')
#    reader_rest_session = open(func_session_file, 'r')
#
#    subj_list = []
#    seed_list = []
#    rest_session_list = []
#
#    for line in reader_subj.readlines():
#        line = line.rstrip('\r\n')
#        subj_list.append(line)
#
#    for line in reader_seed.readlines():
#        line = line.rstrip('\r\n')
#        seed_list.append(line)
#
#    for line in reader_rest_session.readlines():
#        line = line.rstrip('\r\n')
#        rest_session_list.append(line)
#
#
#    return subj_list, rest_session_list, seed_list

    return get_list(c.subj_file), get_list(c.func_session_file), get_list(c.anat_session_file), get_list(c.seed_file)

def get_seed_list(seed_file):

    f = open(seed_file, 'r')

    seed_list = []
    seeds = f.readlines()
    for seed in seeds:
        seed = seed.rstrip('\r\n')
        seed_list.append(seed)

    return seed_list


def getModelandSeedList(c):

    modelist = []
    seedlist = []

    try:
        fp = open(c.model_file, 'r')
        models = fp.readlines()
        fp.close()

        for model in models:
            model = model.rstrip('\r\n')
            modelist.append(os.path.basename(model))

        print 'checking for models ', modelist

        f = open(c.seed_file, 'r')
        seeds = f.readlines()
        f.close()

        for seed in seeds:
            seed = os.path.basename(seed.rstrip('\r\n'))
            seed = os.path.splitext(os.path.splitext(seed)[0])[0]
            seedlist.append(seed)

        print 'checking for seeds ', seedlist
    except:
        raise

    return  modelist, seedlist


def get_workflow(wf_name, c):

    preproc = None
    """ 
        setup standard file paths
    """
    print 'inside get_wf ', '-->'+wf_name+'<--'
    prior_path = os.path.join(c.prior_dir, c.standard_res)
    PRIOR_CSF = os.path.join(prior_path, 'avg152T1_csf_bin.nii.gz')
    PRIOR_GRAY = os.path.join(prior_path, 'avg152T1_gray_bin.nii.gz')
    PRIOR_WHITE = os.path.join(prior_path, 'avg152T1_white_bin.nii.gz')
    standard_res_brain = os.path.join(c.FSLDIR, 'data/standard/MNI152_T1_%s_brain.nii.gz' % (c.standard_res))
    standard = os.path.join(c.FSLDIR, 'data/standard/MNI152_T1_%s.nii.gz' % (c.standard_res))
    standard_brain_mask_dil = os.path.join(c.FSLDIR, 'data/standard/MNI152_T1_%s_brain_mask_dil.nii.gz' % (c.standard_res))
    config_file = os.path.join(c.FSLDIR, 'etc/flirtsch/T1_2_MNI152_%s.cnf' % (c.standard_res))
    brain_symmetric = os.path.join(c.FSLDIR, 'data/standard/MNI152_T1_2mm_brain_symmetric.nii.gz')
    symm_standard = os.path.join(c.FSLDIR, 'data/standard/MNI152_T1_2mm_symmetric.nii.gz')
    twomm_brain_mask_dil = os.path.join(c.FSLDIR, 'data/standard/MNI152_T1_2mm_brain_mask_symmetric_dil.nii.gz')
    config_file_twomm = os.path.join(c.FSLDIR, 'etc/flirtsch/T1_2_MNI152_2mm.cnf')
    identity_matrix = os.path.join(c.FSLDIR, 'etc/flirtsch/ident.mat')

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
        preproc.inputs.hplp_input.hp = c.highPassFreqALFF
        preproc.inputs.hplp_input.lp = c.lowPassFreqALFF
        preproc.inputs.fwhm_input.fwhm = c.fwhm
        preproc.get_node('hplp_input').iterables = ('hp', c.highPassFreqALFF, 'lp', c.lowPassFreqALFF)
        preproc.get_node('fwhm_input').iterables = ('fwhm', c.fwhm)

        return preproc

    if wf_name.lower() == 'sca':

        preproc = create_sca_preproc()
        #seeds = get_seed_list(c.seed_file)

        #print seeds
        #preproc.inputs.seed_list_input.seed_list = seeds
        #preproc.get_node('seed_list_input').iterables = ('seed_list', seeds)
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

    if wf_name.lower() == 'sc':

        preproc = create_scrubbing_preproc()
        return preproc

    if wf_name.lower() == 'select':

        preproc = selector_wf()
        preproc.inputs.run_scrubbing_input.run_scrubbing = c.scrubData
        preproc.get_node('run_scrubbing_input').iterables = ('run_scrubbing', c.scrubData)
        return preproc

    if wf_name.lower() == 'nuisance':

        preproc = create_nuisance_preproc()
        preproc.inputs.selector_input.selector = c.Corrections
        preproc.inputs.nc_input.nc = c.ncomponents
        preproc.inputs.target_angle_deg_input.target_angle_deg = c.target_angle_deg
        preproc.get_node('selector_input').iterables = ('selector', c.Corrections)
        preproc.get_node('nc_input').iterables = ('nc', c.ncomponents)
        preproc.get_node('target_angle_deg_input').iterables = ('target_angle_deg', c.target_angle_deg)

        return preproc

    if wf_name.lower() == 'mprage_in_mnioutputs':
        preproc = mprage_in_mnioutputs()
        preproc.inputs.inputspec.standard = standard

        return preproc


    if wf_name.lower() == 'func_in_mnioutputs':
        preproc = func_in_mnioutputs()
        preproc.inputs.inputspec.standard = standard

        return preproc

    if wf_name.lower() == 'freq_filter':


        preproc = create_filter(c.nuisanceHighPassFilter, c.nuisanceLowPassFilter)
        if c.nuisanceHighPassFilter:
            preproc.inputs.hp_input.hp = c.nuisanceHighPassLowCutOff
            preproc.get_node('hp_input').iterables = ('hp', c.nuisanceHighPassLowCutOff)

        if c.nuisanceLowPassFilter:
            preproc.inputs.lp_input.lp = c.nuisanceLowPassHighCutOff
            preproc.get_node('lp_input').iterables = ('lp', c.nuisanceLowPassHighCutOff)

        return preproc

    if wf_name.lower() == 'ts':

        preproc = create_timeseries_preproc(True, True, True)
        preproc.inputs.inputspec.recon_subjects = c.reconSubjectsDirectory
        preproc.inputs.inputspec.standard = standard
        preproc.inputs.inputspec.identity_matrix = identity_matrix
        preproc.inputs.inputspec.unitTSOutputs = [True, True]
        preproc.inputs.inputspec.voxelTSOutputs = [True, True]
        preproc.inputs.inputspec.verticesTSOutputs = [True, True]
        return preproc

    if wf_name.lower() == 'group_analysis':

        preproc = create_group_analysis(c.f_test)
        preproc.inputs.inputspec.z_threshold = c.z_threshold
        preproc.inputs.inputspec.p_threshold = c.p_threshold
        return preproc


def prep_workflow(c):

    wfname = 'resting_preproc'
    workflow = pe.Workflow(name=wfname)
    workflow.base_dir = c.working_dir
    workflow.crash_dir = c.crash_dir

    sublist, rest_session_list, anat_session_list, seed_list = getSubjectAndSeedLists(c)

    """
        BASIC and ALL preprocessing paths implemented below
    """

    """
        grab the subject data
    """
    flowAnatFunc = create_anat_func_dataflow( sublist,
                                              rest_session_list,
                                              anat_session_list,
                                              c.subj_dir,
                                              c.anat_template,
                                              c.func_template,
                                              c.anat_template_list,
                                              c.func_template_list)

    """
        grab the seeds data
    """
    seedFlow = create_seed_dataflow(c.seed_file)

    """
        grab parcellation data for time series extraction
    """
    pflow = create_parc_dataflow(c.unitDefinitionsDirectory)

    """
        grab mask data for time series extraction
    """
    mflow = create_mask_dataflow(c.voxelMasksDirectory)

    #modelist, seedlist = getModelandSeedList(c)
    #gp_flow = create_gp_dataflow(c.base_dir, modelist, seedlist, 'gpflow')
    """
        Define workflows and set parameters
    """
    anatpreproc = get_workflow('anat', c)
    funcpreproc = get_workflow('func', c)
    regpreproc = get_workflow('reg', c)
    segpreproc = get_workflow('seg', c)
    scpreproc = get_workflow('sc', c)
    select = get_workflow('select', c)
    nuisancepreproc = get_workflow('nuisance', c)
    mprage_mni = get_workflow('mprage_in_mnioutputs', c)
    func_in_mni = get_workflow('func_in_mnioutputs', c)
    freq_filter = get_workflow('freq_filter', c)
    alffpreproc = get_workflow('alff', c)
    scapreproc = get_workflow('sca', c)
    vmhcpreproc = get_workflow('vmhc', c)
    tspreproc = get_workflow('ts', c)
    gppreproc = get_workflow('group_analysis', c)
    """
        Make Connections
    """

    workflow.connect(flowAnatFunc, 'datasource.anat',
                     anatpreproc, 'inputspec.anat')
    workflow.connect(flowAnatFunc, 'datasource.rest',
                     funcpreproc, 'inputspec.rest')
    workflow.connect(funcpreproc, 'outputspec.example_func',
                     regpreproc, 'inputspec.example_func')
    workflow.connect(anatpreproc, 'outputspec.brain',
                     regpreproc, 'inputspec.brain')
    workflow.connect(anatpreproc, 'outputspec.reorient',
                     regpreproc, 'inputspec.reorient')
    workflow.connect(funcpreproc, 'outputspec.preprocessed_mask',
                     segpreproc, 'inputspec.preprocessed_mask')
    workflow.connect(anatpreproc, 'outputspec.brain',
                     segpreproc, 'inputspec.brain')
    workflow.connect(funcpreproc, 'outputspec.example_func',
                     segpreproc, 'inputspec.example_func')
    workflow.connect(regpreproc, 'outputspec.highres2example_func_mat',
                     segpreproc, 'inputspec.highres2example_func_mat')
    workflow.connect(regpreproc, 'outputspec.stand2highres_warp',
                     segpreproc, 'inputspec.stand2highres_warp')
    workflow.connect(flowAnatFunc, 'datasource.rest',
                     scpreproc, 'inputspec.rest')
    workflow.connect(funcpreproc, 'outputspec.movement_parameters',
                     scpreproc, 'inputspec.movement_parameters')
    workflow.connect(funcpreproc, 'outputspec.preprocessed',
                     select, 'inputspec.preprocessed')
    workflow.connect(scpreproc, 'outputspec.scrubbed_preprocessed',
                     select, 'inputspec.scrubbed_preprocessed')
    workflow.connect(funcpreproc, 'outputspec.movement_parameters',
                     select, 'inputspec.movement_parameters')
    workflow.connect(scpreproc, 'outputspec.scrubbed_movement_parameters',
                     select, 'inputspec.scrubbed_movement_parameters')
    workflow.connect(segpreproc, 'outputspec.wm_mask',
                     nuisancepreproc, 'inputspec.wm_mask')
    workflow.connect(segpreproc, 'outputspec.csf_mask',
                     nuisancepreproc, 'inputspec.csf_mask')
    workflow.connect(segpreproc, 'outputspec.gm_mask',
                     nuisancepreproc, 'inputspec.gm_mask')

    """
        Get T1 outputs in MNI
    """
    workflow.connect(regpreproc, 'outputspec.highres2standard_warp',
                     mprage_mni, 'inputspec.highres2standard_warp')
    workflow.connect(anatpreproc, 'outputspec.reorient',
                     mprage_mni, 'inputspec.reorient')
    workflow.connect(anatpreproc, 'outputspec.brain',
                     mprage_mni, 'inputspec.brain')
    workflow.connect(segpreproc, 'outputspec.probability_maps',
                     mprage_mni, 'inputspec.probability_maps')
    workflow.connect(segpreproc, 'outputspec.mixeltype',
                     mprage_mni, 'inputspec.mixeltype')
    workflow.connect(segpreproc, 'outputspec.partial_volume_map',
                     mprage_mni, 'inputspec.partial_volume_map')
    workflow.connect(segpreproc, 'outputspec.partial_volume_files',
                    mprage_mni, 'inputspec.partial_volume_files')

    """
         Nuisance
    """
    workflow.connect(select, 'outputspec.preprocessed_selector',
                     nuisancepreproc, 'inputspec.realigned_file')
    workflow.connect(select, 'outputspec.movement_parameters_selector',
                     nuisancepreproc, 'inputspec.motion_components')

    if(c.nuisanceLowPassFilter or c.nuisanceHighPassFilter):

        workflow.connect(nuisancepreproc, 'outputspec.residual_file',
                         freq_filter, 'inputspec.in_file')


    """
        Get Func outputs in MNI
    """
    workflow.connect(regpreproc, 'outputspec.highres2standard_warp',
                     func_in_mni, 'inputspec.highres2standard_warp')
    workflow.connect(regpreproc, 'outputspec.example_func2highres_mat',
                     func_in_mni, 'inputspec.premat')
    workflow.connect(funcpreproc, 'outputspec.preprocessed_mask',
                     func_in_mni, 'inputspec.preprocessed_mask')
    workflow.connect(nuisancepreproc, 'outputspec.residual_file',
                     func_in_mni, 'inputspec.residual_file')

    """
        ALFF/fALFF
    """
    workflow.connect(nuisancepreproc, 'outputspec.residual_file',
                     alffpreproc, 'inputspec.rest_res')
    workflow.connect(funcpreproc, 'outputspec.preprocessed_mask',
                     alffpreproc, 'inputspec.rest_mask')
    workflow.connect(func_in_mni, 'outputspec.preprocessed_mask_mni',
                     alffpreproc, 'inputspec.rest_mask2standard')
    workflow.connect(regpreproc, 'outputspec.example_func2highres_mat',
                     alffpreproc, 'inputspec.premat')
    workflow.connect(regpreproc, 'outputspec.highres2standard_warp',
                     alffpreproc, 'inputspec.fieldcoeff_file')

    """
        SCA (Seed Based Correlation Analysis)
    """
    workflow.connect(seedFlow, 'out_file',
                     scapreproc, 'seed_list_input.seed_list')
    workflow.connect(funcpreproc, 'outputspec.example_func',
                     scapreproc, 'inputspec.ref')
    workflow.connect(regpreproc, 'outputspec.stand2highres_warp',
                     scapreproc, 'inputspec.warp')
    workflow.connect(regpreproc, 'outputspec.highres2example_func_mat',
                     scapreproc, 'inputspec.postmat')
    workflow.connect(regpreproc, 'outputspec.example_func2highres_mat',
                     scapreproc, 'inputspec.premat')

    if(c.nuisanceLowPassFilter or c.nuisanceHighPassFilter):
        workflow.connect(freq_filter, 'outputspec.rest_res_filt',
                         scapreproc, 'inputspec.rest_res_filt')
    else:
       workflow.connect(nuisancepreproc, 'outputspec.residual_file',
                        scapreproc, 'inputspec.rest_res_filt')

    workflow.connect(regpreproc, 'outputspec.highres2standard_warp',
                     scapreproc, 'inputspec.fieldcoeff_file')
#        workflow.connect(func_in_mni, 'outputspec.residual_file_mni',
#                         scapreproc, 'inputspec.rest_res2standard')
    workflow.connect(func_in_mni, 'outputspec.preprocessed_mask_mni',
                     scapreproc, 'inputspec.rest_mask2standard')


    """
        VMHC (Voxel-Mirrored Homotopic Connectivity)
    """
    workflow.connect(nuisancepreproc, 'outputspec.residual_file',
                     vmhcpreproc, 'inputspec.rest_res')
    workflow.connect(anatpreproc, 'outputspec.reorient',
                     vmhcpreproc, 'inputspec.reorient')
    workflow.connect(anatpreproc, 'outputspec.brain',
                     vmhcpreproc, 'inputspec.brain')
    workflow.connect(regpreproc, 'outputspec.example_func2highres_mat',
                     vmhcpreproc, 'inputspec.example_func2highres_mat')

    """
        Time Series Extraction, Voxel Based and Vertices based
        Generates CSV and NUMPY files
    """
    #timeseries preproc
#    workflow.connect(anatpreproc, 'outputspec.brain',
#                     tspreproc, 'inputspec.brain')
#    workflow.connect(anatpreproc, 'outputspec.reorient',
#                     tspreproc, 'inputspec.reorient')
#    workflow.connect(funcpreproc, 'outputspec.motion_correct',
#                     tspreproc, 'inputspec.motion_correct')
#    workflow.connect(regpreproc, 'outputspec.highres2standard_warp',
#                     tspreproc, 'inputspec.warp_file')
#    workflow.connect(regpreproc, 'outputspec.example_func2highres_mat',
#                     tspreproc, 'inputspec.premat')
#
#    workflow.connect(pflow, 'out_file',
#                     tspreproc, 'getparc.parcelations')
#    workflow.connect(mflow, 'out_file',
#                     tspreproc, 'getmask.masks')


    """
        FSL Group Analysis
    """
#    workflow.connect(gp_flow, 'mat',
#                     gppreproc, 'inputspec.mat_file')
#    workflow.connect(gp_flow, 'con',
#                     gppreproc, 'inputspec.con_file')
#    workflow.connect(gp_flow, 'fts',
#                     gppreproc, 'inputspec.fts_file')
#    workflow.connect(gp_flow, 'grp',
#                     gppreproc, 'inputspec.grp_file')
#    workflow.connect(gp_flow, 'seedfiles',
#                     gppreproc, 'seed_files')

    """
        Calling datasink 
    """
    datasink = create_datasink(c.sink_dir)
    workflow.connect(flowAnatFunc, 'inputnode.subject_id', datasink, 'container')
    anat_sink(workflow, datasink, mprage_mni)
    func_sink(workflow, datasink, funcpreproc, func_in_mni)
    reg_sink(workflow, datasink, regpreproc)
    seg_sink(workflow, datasink, segpreproc, mprage_mni)
    nuisance_sink(workflow, datasink, nuisancepreproc, func_in_mni)
    scrubbing_sink(workflow, datasink, scpreproc)



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
