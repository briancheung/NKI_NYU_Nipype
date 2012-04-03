import argparse
import e_afni
import sys
import os
import nipype.pipeline.engine as pe
from base import (create_anat_preproc, create_func_preproc,
                    create_reg_preproc, create_seg_preproc,
                    create_alff_preproc, create_RSFC_preproc)

from utils import (create_anat_dataflow, create_func_dataflow,
                    create_alff_dataflow, create_rsfc_dataflow)


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


def prep_workflow(c):

    wfname = 'resting_preproc'
    workflow = pe.Workflow(name=wfname)
    workflow.base_dir = c.working_dir
    workflow.crash_dir = c.crash_dir

    sublist, seed_list = getSubjectAndSeedLists(c)

    """ 
        setup standard file paths
    """

    prior_path = os.path.join(c.prior_dir, c.standard_res)
    PRIOR_CSF = os.path.join(prior_path, 'avg152T1_csf_bin.nii.gz')
    PRIOR_WHITE = os.path.join(prior_path, 'avg152T1_white_bin.nii.gz')
    standard_res_brain = os.path.join(c.FSLDIR, 'data/standard/MNI152_T1_%s_brain.nii.gz' % (c.standard_res))
    standard = os.path.join(c.FSLDIR, 'data/standard/MNI152_T1_%s.nii.gz' % (c.standard_res))
    standard_brain_mask_dil = os.path.join(c.FSLDIR, 'data/standard/MNI152_T1_%s_brain_mask_dil.nii.gz' % (c.standard_res))
    config_file = os.path.join(c.FSLDIR, 'etc/flirtsch/T1_2_MNI152_%s.cnf' % (c.standard_res))

    if(c.analysis[0] or c.analysis[1]):

        """
            grab the subject data
        """
        flowAnat = create_anat_dataflow('anat_flow', sublist, c.subj_dir, c.anat_name, c.anat_template)
        flowFunc = create_func_dataflow('func_flow', sublist, c.subj_dir, c.rest_name, c.func_template)


        """
            Define workflows and set parameters
        """
        anatpreproc = create_anat_preproc()

        funcpreproc = create_func_preproc()
        funcpreproc.inputs.inputspec.start_idx = c.start_idx
        funcpreproc.inputs.inputspec.stop_idx = c.stop_idx

        segpreproc = create_seg_preproc()
        segpreproc.inputs.inputspec.PRIOR_CSF  = PRIOR_CSF
        segpreproc.inputs.inputspec.PRIOR_WHITE  = PRIOR_WHITE
        segpreproc.inputs.inputspec.standard_res_brain = standard_res_brain

        regpreproc = create_reg_preproc()
        regpreproc.inputs.inputspec.standard_res_brain = standard_res_brain
        regpreproc.inputs.inputspec.standard = standard
        regpreproc.inputs.inputspec.config_file = config_file
        regpreproc.inputs.inputspec.standard_brain_mask_dil = standard_brain_mask_dil


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

    #if((not c.analysis[0] and not c.analysis[1]) and c.analysis[4]):
#
     #   alffpreproc = create_alff_preproc()
      #  flowAlff, flowAlffWarp = create_alff_dataflow('alff_flow', sublist, c.subj_dir, c.rest_name, c.func_template)
#
#		

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
