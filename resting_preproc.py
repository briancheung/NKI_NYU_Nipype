#!/frodo/shared/epd/bin/python
import e_afni
import sys
import os
import nipype.pipeline.engine as pe
from base import (create_anat_preproc, create_func_preproc,
                    create_reg_preproc, create_seg_preproc)




def processS():

    numCores=int(sys.argv[1])

    wfname='resting_preproc'
    workflow=pe.Workflow(name=wfname)

    anatpreproc=create_anat_preproc()
    anatpreproc.inputs.inputspec.anat=os.path.abspath('/home/sharad/Resources/0010052/session_1/anat_1/mprage.nii.gz')

    funcpreproc=create_func_preproc()
    funcpreproc.inputs.inputspec.rest=os.path.abspath('/home/sharad/Resources/0010052/session_1/rest_1/rest.nii.gz')
    funcpreproc.inputs.inputspec.start_idx=0
    funcpreproc.inputs.inputspec.stop_idx=175

    segpreproc=create_seg_preproc()
    segpreproc.inputs.inputspec.PRIOR_CSF = os.path.abspath('./tissuepriors/3mm/avg152T1_csf_bin.nii.gz')
    segpreproc.inputs.inputspec.PRIOR_WHITE = os.path.abspath('./tissuepriors/3mm/avg152T1_white_bin.nii.gz')
    segpreproc.inputs.inputspec.standard_res_brain=os.path.abspath('/frodo/shared/RH5_fsl/data/standard/MNI152_T1_3mm_brain.nii.gz')

    regpreproc=create_reg_preproc()
    regpreproc.inputs.inputspec.standard_res_brain=os.path.abspath('/frodo/shared/RH5_fsl/data/standard/MNI152_T1_3mm_brain.nii.gz')
    regpreproc.inputs.inputspec.standard=os.path.abspath('/frodo/shared/RH5_fsl/data/standard/MNI152_T1_3mm.nii.gz')
    regpreproc.inputs.inputspec.config_file=os.path.abspath('/frodo/shared/RH5_fsl/etc/flirtsch/T1_2_MNI152_3mm.cnf')
    regpreproc.inputs.inputspec.standard_brain_mask_dil=os.path.abspath('/frodo/shared/RH5_fsl/data/standard/MNI152_T1_3mm_brain_mask_dil.nii.gz')

    workflow.connect(funcpreproc, 'outputspec.example_func', regpreproc, 'inputspec.example_func')
    workflow.connect(anatpreproc, 'outputspec.brain', regpreproc, 'inputspec.brain')
    workflow.connect(anatpreproc, 'outputspec.reorient', regpreproc, 'inputspec.reorient')
    workflow.connect(funcpreproc, 'outputspec.preprocessed_mask', segpreproc, 'inputspec.preprocessed_mask')
    workflow.connect(anatpreproc, 'outputspec.brain', segpreproc, 'inputspec.brain')
    workflow.connect(funcpreproc, 'outputspec.example_func', segpreproc, 'inputspec.example_func')
    workflow.connect(regpreproc, 'outputspec.highres2example_func_mat', segpreproc, 'inputspec.highres2example_func_mat')
    workflow.connect(regpreproc, 'outputspec.stand2highres_warp', segpreproc, 'inputspec.stand2highres_warp')


#    workflow.add_nodes([anatpreproc, funcpreproc])
    workflow.base_dir='/home/sharad/nki_nyu_pipeline/working_dir'
    workflow.run(plugin='MultiProc', plugin_args={'n_procs': numCores})




def main():

    if (len(sys.argv) < 2 ):
        sys.stderr.write("./resting_preproc.py <Number_of_cores>")
    else:

        processS()


if __name__ == "__main__":

    sys.exit(main())
