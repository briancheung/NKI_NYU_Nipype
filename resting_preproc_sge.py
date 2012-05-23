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
                    create_gp_dataflow, create_datasink,
                    symlink_creator)

from base_nuisance import create_nuisance_preproc


from sink import (anat_sink,
                  reg_sink,
                  seg_sink,
                  func_sink,
                  nuisance_sink,
                  scrubbing_sink,
                  sca_sink,
                  alff_sink,
                  vmhc_sink,
                  timeseries_sink,
                  group_analysis_sink)

from multiprocessing import Pool
from multiprocessing import Process

def getSubjectAndSeedLists(c):

    """
        Read Subject & Seed files to build
        corresponding lists.
    """

    def get_list(fname):
        flines = open(fname, 'r').readlines()
        return [fline.rstrip('\r\n') for fline in flines]

    def remove_subjects(sub_list, ex_sub_list):

        for sub in ex_sub_list:
            sub_list.remove(sub)

        return sub_list

    sublist = None
    func_session_list = None
    anat_session_list = None

    if c.subjectList is None:
        import commands
        commands.getoutput('ls -1 %s/ > CPAC_SUB_Auto.txt' % (c.subjectDirectory))

        sublist = get_list('CPAC_SUB_Auto.txt')
        commands.getoutput('rm -f CPAC_SUB_Auto.txt')
    else:
        sublist = get_list(c.subjectList)

    if not c.exclusionSubjectList is None:
            ex_sub_list = get_list(c.exclusionSubjectList)
            remove_subjects(sublist, ex_sub_list)

    func_session_list = get_list(c.funcSessionFile)
    anat_session_list = get_list(c.anatSessionFile)
    if c.derivatives[1]:
        return sublist,\
               func_session_list,\
               anat_session_list,\
               get_list(c.seedFile)
    else:
        return sublist,\
               func_session_list,\
               anat_session_list,\
               []


def get_seed_list(seed_file):

    f = open(seed_file, 'r')

    seed_list = []
    seeds = f.readlines()
    for seed in seeds:
        seed = seed.rstrip('\r\n')
        seed_list.append(seed)

    return seed_list



def getGroupAnalysisInputs(c):

    modelist = []
    dervlist = []
    labelist = []

    try:

        models = os.listdir(c.modelsDirectory)

        for model in models:
            model = model.rstrip('\r\n')
            modelist.append(model)

        print 'checking for models ', modelist

        f = open(c.labelFile, 'r')
        labels = f.readlines()
        f.close()
        for l in labels:
            labelist.append(l.split()[0])

        print 'checking for strategy list ', labelist

        dervlist = c.derivativeList
        print 'checking for derivatives ', dervlist

        template_dict = dict(mat=os.path.join(c.modelsDirectory, c.mat),
                             con=os.path.join(c.modelsDirectory, c.con),
                             fts=os.path.join(c.modelsDirectory, c.fts),
                             grp=os.path.join(c.modelsDirectory, c.grp),
                             derv=c.dervTemplate)

        print 'checking for template_dict', template_dict

        template_args = dict(mat=[c.matTemplateList],
                             con=[c.conTemplateList],
                             fts=[c.ftsTemplateList],
                             grp=[c.grpTemplateList],
                             derv=[c.dervTemplateList])

        print 'checking for template_args', template_args

    except:
        raise

    return  modelist, dervlist, labelist, template_dict, template_args



def get_workflow(wf_name, c):

    preproc = None
    """
        setup standard file paths
    """
    prior_path = os.path.join(c.priorDirectory, c.standardResolution)
    PRIOR_CSF = os.path.join(prior_path, 'avg152T1_csf_bin.nii.gz')
    PRIOR_GRAY = os.path.join(prior_path, 'avg152T1_gray_bin.nii.gz')
    PRIOR_WHITE = os.path.join(prior_path, 'avg152T1_white_bin.nii.gz')
    standard_res_brain = os.path.join(c.FSLDIR,
            'data/standard/MNI152_T1_%s_brain.nii.gz' % (c.standardResolution))
    standard = os.path.join(c.FSLDIR,
            'data/standard/MNI152_T1_%s.nii.gz' % (c.standardResolution))
    standard_brain_mask_dil = '/home2/ssikka/scripts1/templates/MNI152_T1_%s_brain_mask_dil.nii.gz' % (c.standardResolution)
    config_file = '/home2/ssikka/scripts1/templates/T1_2_MNI152_%s.cnf' % (c.standardResolution)
    brain_symmetric = '/home2/ssikka/scripts1/templates/MNI152_T1_2mm_brain_symmetric.nii.gz'
    symm_standard = '/home2/ssikka/scripts1/templates/MNI152_T1_2mm_symmetric.nii.gz'
    twomm_brain_mask_dil = '/home2/ssikka/scripts1/templates/MNI152_T1_2mm_brain_mask_symmetric_dil.nii.gz'
    config_file_twomm = '/home2/ssikka/scripts1/templates/T1_2_MNI152_2mm.cnf'
    identity_matrix = os.path.join(c.FSLDIR,
            'etc/flirtsch/ident.mat')


    if wf_name.lower() == 'anat':
        preproc = create_anat_preproc()
        return preproc

    if wf_name.lower() == 'func':

        preproc = create_func_preproc()
        preproc.inputs.inputspec.start_idx = c.startIdx
        preproc.inputs.inputspec.stop_idx = c.stopIdx

        return preproc

    if wf_name.lower() == 'seg':

        preproc = create_seg_preproc()
        preproc.inputs.inputspec.PRIOR_CSF = PRIOR_CSF
        preproc.inputs.inputspec.PRIOR_WHITE = PRIOR_WHITE
        preproc.inputs.inputspec.PRIOR_GRAY = PRIOR_GRAY
        preproc.inputs.inputspec.standard_res_brain = standard_res_brain
        preproc.inputs.csf_threshold.csf_threshold = \
                                        c.cerebralSpinalFluidThreshold
        preproc.inputs.wm_threshold.wm_threshold = \
                                        c.whiteMatterThreshold
        preproc.inputs.gm_threshold.gm_threshold = \
                                        c.grayMatterThreshold
        preproc.get_node('csf_threshold').iterables = ('csf_threshold',
                                        c.cerebralSpinalFluidThreshold)
        preproc.get_node('wm_threshold').iterables = ('wm_threshold',
                                        c.whiteMatterThreshold)
        preproc.get_node('gm_threshold').iterables = ('gm_threshold',
                                        c.grayMatterThreshold)

        return preproc

    if wf_name.lower() == 'reg':

        preproc = create_reg_preproc()
        preproc.inputs.inputspec.standard_res_brain = \
                                            standard_res_brain
        preproc.inputs.inputspec.standard = standard
        preproc.inputs.inputspec.config_file = config_file
        preproc.inputs.inputspec.standard_brain_mask_dil = \
                                            standard_brain_mask_dil

        return preproc

    if wf_name.lower() == 'alff':

        preproc = create_alff_preproc()
        preproc.inputs.inputspec.standard = standard
        preproc.inputs.hp_input.hp = c.highPassFreqALFF
        preproc.inputs.lp_input.lp = c.lowPassFreqALFF
        preproc.inputs.fwhm_input.fwhm = c.fwhm
        preproc.get_node('hp_input').iterables = ('hp',
                                                    c.highPassFreqALFF)
        preproc.get_node('lp_input').iterables = ('lp',
                                                    c.lowPassFreqALFF)
        preproc.get_node('fwhm_input').iterables = ('fwhm',
                                                    c.fwhm)

        return preproc

    if wf_name.lower() == 'sca':

        preproc = create_sca_preproc(c.correlationSpace)
        preproc.inputs.fwhm_input.fwhm = c.fwhm
        preproc.get_node('fwhm_input').iterables = ('fwhm',
                                                    c.fwhm)
        preproc.inputs.inputspec.standard = standard
        return preproc


    if wf_name.lower() == 'vmhc':

        preproc = create_vmhc_preproc()
        preproc.inputs.inputspec.brain_symmetric = \
                                        brain_symmetric
        preproc.inputs.inputspec.symm_standard = \
                                        symm_standard
        preproc.inputs.inputspec.twomm_brain_mask_dil = \
                                        twomm_brain_mask_dil
        preproc.inputs.inputspec.config_file_twomm = \
                                        config_file_twomm
        preproc.inputs.inputspec.standard = \
                                        standard
        return preproc

    if wf_name.lower() == 'sc':

        preproc = create_scrubbing_preproc()
        preproc.inputs.threshold_input.threshold = c.scrubbingThreshold
        preproc.get_node('threshold_input').iterables = ('threshold', c.scrubbingThreshold)
        return preproc

    if wf_name.lower() == 'select':

        preproc = selector_wf()
        preproc.inputs.run_scrubbing_input.run_scrubbing = c.scrubData
        preproc.get_node('run_scrubbing_input').iterables = \
                                    ('run_scrubbing', c.scrubData)
        return preproc

    if wf_name.lower() == 'nuisance':

        preproc = create_nuisance_preproc()
        preproc.inputs.selector_input.selector = c.Corrections
        preproc.inputs.nc_input.nc = c.nComponents
        preproc.inputs.target_angle_deg_input.target_angle_deg = \
                                            c.targetAngleDeg

        preproc.get_node('selector_input').iterables = \
                                            ('selector', c.Corrections)
        preproc.get_node('nc_input').iterables = \
                                            ('nc', c.nComponents)
        preproc.get_node('target_angle_deg_input').iterables = \
                                ('target_angle_deg', c.targetAngleDeg)

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

        from base_nuisance import bandpass_voxels
        import nipype.interfaces.utility as util
        preproc = pe.MapNode(util.Function(input_names=['realigned_file',
                                                        'sample_period',
                                                        'bandpass_freqs'],
                                           output_names=['bandpassed_file'],
                                           function=bandpass_voxels),
                                           name='bandpass_filter',
                                           iterfield=['realigned_file', 'bandpass_freqs'])
        preproc.inputs.bandpass_freqs = c.nuisanceBandpassFreq
        preproc.inputs.sample_period = c.TR
        return preproc


    if wf_name.lower() == 'ts':

        preproc = create_timeseries_preproc(True, True, True, c.runSurfaceRegistraion)
        preproc.inputs.inputspec.recon_subjects = c.reconSubjectsDirectory
        preproc.inputs.inputspec.standard = standard
        preproc.inputs.inputspec.identity_matrix = identity_matrix
        preproc.inputs.inputspec.unitTSOutputs = [True, True]
        preproc.inputs.inputspec.voxelTSOutputs = [True, True]
        preproc.inputs.inputspec.verticesTSOutputs = [True, True]
        return preproc

    if wf_name.lower() == 'group_analysis':

        preproc = create_group_analysis(c.fTest)
        preproc.inputs.inputspec.z_threshold = c.zThreshold
        preproc.inputs.inputspec.p_threshold = c.pThreshold
        preproc.inputs.inputspec.parameters = (c.FSLDIR,
                                               c.MNI)
        return preproc


def prep_workflow(sub, rest_session_list, anat_session_list, seed_list, c):



    wfname = 'resting_preproc_sge'
    workflow = pe.Workflow(name=wfname)
    workflow.base_dir = c.workingDirectory
    workflow.crash_dir = c.crashLogDirectory
    workflow.config['execution'] = {'hash_method': 'timestamp'}


    """
        BASIC and ALL preprocessing paths implemented below
    """

    """
        grab the subject data
    """
    flowAnatFunc = create_anat_func_dataflow(sub,
                                              rest_session_list,
                                              anat_session_list,
                                              c.subjectDirectory,
                                              c.anatTemplate,
                                              c.funcTemplate,
                                              c.anatTemplateList,
                                              c.funcTemplateList)


    if c.derivatives[6]:
        """
                grab the Group Analysis Data 
        """
        modelist, dervlist, labelist, template_dict, template_args = getGroupAnalysisInputs(c)

        gp_flow = create_gp_dataflow(c.sinkDirectory, modelist,
                                     dervlist, labelist, template_dict,
                                     template_args, 'gpflow')

    if c.derivatives[1]:
        """
            grab the seeds data
        """
        seedFlow = create_seed_dataflow(c.seedFile)

    if c.derivatives[4]:
        """
            grab parcellation data for time series extraction
        """
        pflow = create_parc_dataflow(c.unitDefinitionsDirectory)

    if c.derivatives[5]:
        pflow = create_parc_dataflow(c.unitDefinitionsDirectory)

    if c.derivatives[5]:
        """
            grab mask data for time series extraction
        """
        mflow = create_mask_dataflow(c.voxelMasksDirectory)

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

    if c.derivatives[0]:
        alffpreproc = get_workflow('alff', c)

    if c.derivatives[1]:
        scapreproc = get_workflow('sca', c)

    if c.derivatives[2]:
        vmhcpreproc = get_workflow('vmhc', c)

    if c.derivatives[4] or c.derivatives[5]:
        tspreproc = get_workflow('ts', c)

    if c.derivatives[6]:
        gppreproc = get_workflow('group_analysis', c)
    """
        Calling datasink
    """
    datasink = create_datasink(c.sinkDirectory)
    workflow.connect(flowAnatFunc, 'inputnode.subject_id',
    datasink, 'container')
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
    workflow.connect(funcpreproc, 'outputspec.preprocessed',
                     nuisancepreproc, 'inputspec.realigned_file')
    workflow.connect(funcpreproc, 'outputspec.movement_parameters',
                     nuisancepreproc, 'inputspec.motion_components')

    if(c.nuisanceBandpassFreq):

        workflow.connect(nuisancepreproc, 'outputspec.residual_file',
                         freq_filter, 'realigned_file')

        workflow.connect(freq_filter, 'bandpassed_file',
                         scpreproc, 'inputspec.preprocessed')
        workflow.connect(freq_filter, 'bandpassed_file',
                         select, 'inputspec.preprocessed')
    else:
        workflow.connect(nuisancepreproc, 'outputspec.residual_file',
                         scpreproc, 'inputspec.preprocessed')
        workflow.connect(nuisancepreproc, 'outputspec.residual_file',
                         select, 'inputspec.preprocessed')

    workflow.connect(flowAnatFunc, 'datasource.rest',
                     scpreproc, 'inputspec.rest')
    workflow.connect(funcpreproc, 'outputspec.movement_parameters',
                     scpreproc, 'inputspec.movement_parameters')

    workflow.connect(scpreproc, 'outputspec.scrubbed_preprocessed',
                     select, 'inputspec.scrubbed_preprocessed')
    workflow.connect(funcpreproc, 'outputspec.movement_parameters',
                     select, 'inputspec.movement_parameters')
    workflow.connect(scpreproc, 'outputspec.scrubbed_movement_parameters',
                     select, 'inputspec.scrubbed_movement_parameters')
    """
        Get Func outputs in MNI
    """
    workflow.connect(regpreproc, 'outputspec.highres2standard_warp',
                     func_in_mni, 'inputspec.highres2standard_warp')
    workflow.connect(regpreproc, 'outputspec.example_func2highres_mat',
                     func_in_mni, 'inputspec.premat')
    workflow.connect(funcpreproc, 'outputspec.preprocessed_mask',
                     func_in_mni, 'inputspec.preprocessed_mask')

    workflow.connect(select, 'outputspec.preprocessed_selector',
                     func_in_mni, 'inputspec.residual_file')


    anat_sink(workflow,
              datasink,
              mprage_mni)
    func_sink(workflow,
              datasink,
              funcpreproc,
              func_in_mni)
    reg_sink(workflow,
             datasink,
             regpreproc)
    seg_sink(workflow,
             datasink,
             segpreproc,
             mprage_mni)
    scrubbing_sink(workflow,
                   datasink,
                   scpreproc)
    nuisance_sink(workflow,
                  datasink,
                  nuisancepreproc)

    if c.derivatives[0]:
        """
                ALFF/fALFF
        """
        workflow.connect(select, 'outputspec.preprocessed_selector',
                         alffpreproc, 'inputspec.rest_res')
        workflow.connect(funcpreproc, 'outputspec.preprocessed_mask',
                         alffpreproc, 'inputspec.rest_mask')
        workflow.connect(func_in_mni, 'outputspec.preprocessed_mask_mni',
                         alffpreproc, 'inputspec.rest_mask2standard')
        workflow.connect(regpreproc, 'outputspec.example_func2highres_mat',
                         alffpreproc, 'inputspec.premat')
        workflow.connect(regpreproc, 'outputspec.highres2standard_warp',
                         alffpreproc, 'inputspec.fieldcoeff_file')
        alff_sink(workflow,
                  datasink,
                  alffpreproc)


    if c.derivatives[1]:
        """
            SCA (Seed Based Correlation Analysis)
        """
        workflow.connect(seedFlow, 'out_file',
                         scapreproc, 'seed_list_input.seed_list')
        workflow.connect(regpreproc, 'outputspec.example_func2highres_mat',
                         scapreproc, 'inputspec.premat')

        workflow.connect(select, 'outputspec.preprocessed_selector',
                         scapreproc, 'inputspec.residual_file')
        workflow.connect(select, 'outputspec.preprocessed_selector',
                         scapreproc, 'inputspec.rest_res_filt')

        workflow.connect(regpreproc, 'outputspec.highres2standard_warp',
                         scapreproc, 'inputspec.fieldcoeff_file')
        workflow.connect(func_in_mni, 'outputspec.preprocessed_mask_mni',
                         scapreproc, 'inputspec.rest_mask2standard')
        sca_sink(workflow,
                 datasink,
                 scapreproc)


    if c.derivatives[2]:
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

        vmhc_sink(workflow,
                  datasink,
                  vmhcpreproc)

    """
        Time Series Extraction, Voxel Based and Vertices based
        Generates CSV and NUMPY files
    """

    if c.derivatives[4] or c.derivatives[5]:
        workflow.connect(anatpreproc, 'outputspec.brain',
                         tspreproc, 'inputspec.brain')
        workflow.connect(anatpreproc, 'outputspec.reorient',
                         tspreproc, 'inputspec.reorient')
        workflow.connect(funcpreproc, 'outputspec.motion_correct',
                         tspreproc, 'inputspec.motion_correct')
        workflow.connect(regpreproc, 'outputspec.highres2standard_warp',
                         tspreproc, 'inputspec.warp_file')
        workflow.connect(regpreproc, 'outputspec.example_func2highres_mat',
                         tspreproc, 'inputspec.premat')

        workflow.connect(pflow, 'out_file',
                         tspreproc, 'getparc.parcelations')
        workflow.connect(mflow, 'out_file',
                         tspreproc, 'getmask.masks')

        timeseries_sink(workflow,
                        datasink,
                        tspreproc,
                        c.runSurfaceRegistraion)

    """
        FSL Group Analysis
    """
    if c.derivatives[6]:
        workflow.connect(gp_flow, 'gpflow.mat',
                     gppreproc, 'inputspec.mat_file')
        workflow.connect(gp_flow, 'gpflow.con',
                     gppreproc, 'inputspec.con_file')
        workflow.connect(gp_flow, 'gpflow.fts',
                     gppreproc, 'inputspec.fts_file')
        workflow.connect(gp_flow, 'gpflow.grp',
                     gppreproc, 'inputspec.grp_file')
        workflow.connect(gp_flow, 'gpflow.derv',
                     gppreproc, 'inputspec.zmap_files')


        gp_datasink = create_gp_datasink(c.sinkDirectory)
        workflow.connect(gp_flow, 'inputnode.label', gp_datasink, 'container')
        group_analysis_sink(workflow, gp_datasink, gppreproc)

    if(not c.runOnGrid):
        workflow.run(plugin='MultiProc',
                     plugin_args={'n_procs': c.numCoresPerSubject})
    else:
        template_str = """
#!/bin/bash
#$ -S /bin/bash
#$ -V
    """
        workflow.run(plugin='SGE',
                     plugin_args=dict(qsub_args=c.qsubArgs, template=template_str))

    workflow.write_graph()

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
    sublist, rest_session_list, anat_session_list, seed_list = \
                                            getSubjectAndSeedLists(c)

    prep_workflow(sublist, rest_session_list, anat_session_list, seed_list, c)

    symlink_creator(c.sinkDirectory, sublist)


if __name__ == "__main__":

    import os
    import commands
    cmd = 'bash /home2/ssikka/.bashrc'
    print cmd
    sys.stderr.write(commands.getoutput(cmd))

    sys.exit(main())
