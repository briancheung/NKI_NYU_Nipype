import sys
import e_afni
import os
import commands
import nipype.pipeline.engine as pe
import nipype.algorithms.rapidart as ra
import nipype.interfaces.afni as afni
import nipype.interfaces.fsl as fsl
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
from utils import *

def create_filter(high_pass_filter_bool, low_pass_filter_bool):

    preproc = pe.Workflow(name='frequency_filter')

    inputNode = pe.Node(util.IdentityInterface(fields=['in_file']),
                        name='inputspec')

    inputnode_hp = pe.Node(util.IdentityInterface(fields=['hp']),
                           name='hp_input')

    inputnode_lp = pe.Node(util.IdentityInterface(fields=['lp']),
                           name='lp_input')

    outputNode = pe.Node(util.IdentityInterface(fields=['rest_res_filt']),
                         name='outputspec')

    filter = pe.MapNode(interface=afni.Fourier(),
                        name='func_filter',
                        iterfield=['in_file'])
    filter.inputs.other = '-retrend'

    preproc.connect(inputNode, 'in_file',
                    filter, 'in_file')

    if(high_pass_filter_bool):

        preproc.connect(inputnode_hp, 'hp',
                        filter, 'highpass')

    if(low_pass_filter_bool):
        preproc.connect(inputnode_lp, 'lp',
                        filter, 'lowpass')

    preproc.connect(filter, 'out_file',
                    outputNode, 'rest_res_filt')

    return preproc


def create_fmri_mni():

    preproc = pe.Workflow(name='fmri_mnioutputs')

    inputNode = pe.Node(util.IdentityInterface(fields=['reference_file',
                                                       'warp_file',
                                                       'premat',
                                                       'in_file']),
                       name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['out_file']),
                       name='outputspec')

    apply_warp = pe.MapNode(interface=fsl.ApplyWarp(),
                            name='apply_warp',
                            iterfield=['in_file', 'premat'])

    preproc.connect(inputNode, 'in_file',
                    apply_warp, 'in_file')
    preproc.connect(inputNode, 'warp_file',
                    apply_warp, 'field_file')
    preproc.connect(inputNode, 'reference_file',
                    apply_warp, 'ref_file')
    preproc.connect(inputNode, 'premat',
                    apply_warp, 'premat')
    preproc.connect(apply_warp, 'out_file',
                    outputNode, 'out_file')

    return preproc


def create_Tone_mni():

    preproc = pe.Workflow(name='Tone_mnioutputs')

    inputNode = pe.Node(util.IdentityInterface(fields=['reference_file',
                                                       'warp_file',
                                                       'in_file']),
                       name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['out_file']),
                       name='outputspec')

    apply_warp = pe.MapNode(interface=fsl.ApplyWarp(),
                            name='apply_warp',
                            iterfield=['in_file'])

    preproc.connect(inputNode, 'in_file',
                    apply_warp, 'in_file')
    preproc.connect(inputNode, 'warp_file',
                    apply_warp, 'field_file')
    preproc.connect(inputNode, 'reference_file',
                    apply_warp, 'ref_file')
    preproc.connect(apply_warp, 'out_file',
                    outputNode, 'out_file')

    return preproc


def mprage_in_mnioutputs():

    preproc = pe.Workflow(name='mprage_in_mnioutputs')

    inputNode = pe.Node(util.IdentityInterface(fields=['reorient',
                                                'brain',
                                                'mixeltype',
                                                'probability_maps',
                                                'partial_volume_map',
                                                'partial_volume_files',
                                                'highres2standard_warp',
                                                'standard']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['reorient_mni',
                                                'brain_mni',
                                                'mixeltype_mni',
                                                'probability_maps_mni',
                                                'partial_volume_map_mni',
                                                'partial_volume_files_mni'
                                                 ]),
                       name='outputspec')

    r_mni = create_Tone_mni()
    b_mni = r_mni.clone('brain_mni')
    m_mni = r_mni.clone('mixeltype_mni')
    p_maps_mni = r_mni.clone('probability_maps_in_mni')
    pv_map_mni = r_mni.clone('partial_volume_map_mni')
    pv_files_mni = r_mni.clone('partial_volume_files_mni')

    preproc.connect(inputNode, 'reorient',
                    r_mni, 'inputspec.in_file')
    preproc.connect(inputNode, 'standard',
                    r_mni, 'inputspec.reference_file')
    preproc.connect(inputNode, 'highres2standard_warp',
                    r_mni, 'inputspec.warp_file')

    preproc.connect(inputNode, 'brain',
                    b_mni, 'inputspec.in_file')
    preproc.connect(inputNode, 'standard',
                    b_mni, 'inputspec.reference_file')
    preproc.connect(inputNode, 'highres2standard_warp',
                    b_mni, 'inputspec.warp_file')

    preproc.connect(inputNode, 'mixeltype',
                    m_mni, 'inputspec.in_file')
    preproc.connect(inputNode, 'standard',
                    m_mni, 'inputspec.reference_file')
    preproc.connect(inputNode, 'highres2standard_warp',
                    m_mni, 'inputspec.warp_file')

    preproc.connect(inputNode, 'probability_maps',
                    p_maps_mni, 'inputspec.in_file')
    preproc.connect(inputNode, 'standard',
                    p_maps_mni, 'inputspec.reference_file')
    preproc.connect(inputNode, 'highres2standard_warp',
                    p_maps_mni, 'inputspec.warp_file')

    preproc.connect(inputNode, 'partial_volume_map',
                    pv_map_mni, 'inputspec.in_file')
    preproc.connect(inputNode, 'standard',
                    pv_map_mni, 'inputspec.reference_file')
    preproc.connect(inputNode, 'highres2standard_warp',
                    pv_map_mni, 'inputspec.warp_file')

    preproc.connect(inputNode, 'partial_volume_files',
                    pv_files_mni, 'inputspec.in_file')
    preproc.connect(inputNode, 'standard',
                    pv_files_mni, 'inputspec.reference_file')
    preproc.connect(inputNode, 'highres2standard_warp',
                    pv_files_mni, 'inputspec.warp_file')

    preproc.connect(r_mni, 'outputspec.out_file',
                    outputNode, 'reorient_mni')
    preproc.connect(b_mni, 'outputspec.out_file',
                    outputNode, 'brain_mni')
    preproc.connect(m_mni, 'outputspec.out_file',
                    outputNode, 'mixeltype_mni')
    preproc.connect(p_maps_mni, 'outputspec.out_file',
                    outputNode, 'probability_maps_mni')
    preproc.connect(pv_map_mni, 'outputspec.out_file',
                    outputNode, 'partial_volume_map_mni')
    preproc.connect(pv_files_mni, 'outputspec.out_file',
                    outputNode, 'partial_volume_files_mni')

    return preproc


def func_in_mnioutputs():

    preproc = pe.Workflow(name='func_in_mnioutputs')

    inputNode = pe.Node(util.IdentityInterface(fields=['residual_file',
                                                       'highres2standard_warp',
                                                       'premat',
                                                       'standard',
                                                       'preprocessed_mask']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['residual_file_mni',
                                                    'preprocessed_mask_mni'
                                                       ]),
                         name='outputspec')

    residual_mni = create_fmri_mni()
    pm_mni = residual_mni.clone('pm_mni')

    preproc.connect(inputNode, 'residual_file',
                    residual_mni, 'inputspec.in_file')
    preproc.connect(inputNode, 'highres2standard_warp',
                    residual_mni, 'inputspec.warp_file')
    preproc.connect(inputNode, 'standard',
                    residual_mni, 'inputspec.reference_file')
    preproc.connect(inputNode, 'premat',
                    residual_mni, 'inputspec.premat')

    preproc.connect(inputNode, 'preprocessed_mask',
                    pm_mni, 'inputspec.in_file')
    preproc.connect(inputNode, 'highres2standard_warp',
                    pm_mni, 'inputspec.warp_file')
    preproc.connect(inputNode, 'standard',
                    pm_mni, 'inputspec.reference_file')
    preproc.connect(inputNode, 'premat',
                    pm_mni, 'inputspec.premat')

    preproc.connect(residual_mni, 'outputspec.out_file',
                    outputNode, 'residual_file_mni')

    preproc.connect(pm_mni, 'outputspec.out_file',
                    outputNode, 'preprocessed_mask_mni')

    return preproc


def create_anat_preproc():

    preproc = pe.Workflow(name='anatpreproc')

    #Node or Node ?
    inputNode = pe.Node(util.IdentityInterface(fields=['anat']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['refit',
                                                    'reorient',
                                                    'skullstrip',
                                                    'brain']),
                         name='outputspec')

    anat_refit = pe.Node(interface=e_afni.Threedrefit(),
                         name='anat_refit')
    anat_refit.inputs.deoblique = True

    anat_reorient = pe.Node(interface=afni.Resample(),
                            name='anat_reorient')
    anat_reorient.inputs.orientation = 'RPI'

    anat_skullstrip = pe.Node(interface=afni.SkullStrip(),
                              name='anat_skullstrip')
    anat_skullstrip.inputs.options = '-o_ply'

    anat_calc = pe.Node(interface=afni.Calc(),
                        name='anat_calc')
    anat_calc.inputs.expr = '\'a*step(b)\''

    preproc.connect(inputNode, 'anat',
                    anat_refit, 'in_file')
    preproc.connect(anat_refit, 'out_file',
                    anat_reorient, 'in_file')
    preproc.connect(anat_reorient, 'out_file',
                    anat_skullstrip, 'in_file')
    preproc.connect(anat_skullstrip, 'out_file',
                    anat_calc, 'in_file_b')
    preproc.connect(anat_reorient, 'out_file',
                    anat_calc, 'in_file_a')

    preproc.connect(anat_refit, 'out_file',
                    outputNode, 'refit')
    preproc.connect(anat_reorient, 'out_file',
                    outputNode, 'reorient')
    preproc.connect(anat_skullstrip, 'out_file',
                    outputNode, 'skullstrip')
    preproc.connect(anat_calc, 'out_file',
                    outputNode, 'brain')

    return preproc


def create_func_preproc():

    preproc = pe.Workflow(name='funcpreproc')
    inputNode = pe.Node(util.IdentityInterface(fields=['rest',
                                               'start_idx',
                                               'stop_idx']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['drop_tr',
                                                'refit',
                                                'reorient',
                                                'reorient_mean',
                                                'motion_correct',
                                                'motion_correct_ref',
                                                'movement_parameters',
                                                'max_displacement',
                                                'mask',
                                                'skullstrip',
                                                'example_func',
                                                'preprocessed',
                                                'preprocessed_mask']),

                          name='outputspec')

    func_calc = pe.MapNode(interface=afni.Calc(),
                           name='func_calc',
                           iterfield=['in_file_a'])
    func_calc.inputs.expr = '\'a\''

    func_refit = pe.MapNode(interface=e_afni.Threedrefit(),
                            name='func_refit',
                            iterfield=['in_file'])
    func_refit.inputs.deoblique = True

    func_reorient = pe.MapNode(interface=afni.Resample(),
                               name='func_reorient',
                               iterfield=['in_file'])

    func_reorient.inputs.orientation = 'RPI'

    func_tstat = pe.MapNode(interface=afni.TStat(),
                            name='func_tstat',
                            iterfield=['in_file'])
    func_tstat.inputs.options = '-mean'

    func_tstat_1 = func_tstat.clone('func_tstat_1')

    func_volreg = pe.MapNode(interface=e_afni.Threedvolreg(),
                             name='func_volreg',
                             iterfield=['in_file', 'basefile'])

    func_volreg.inputs.other = '-Fourier -twopass'
    func_volreg.inputs.zpad = '4'

    func_volreg_1 = func_volreg.clone('func_volreg_1')
    func_automask = pe.MapNode(interface=e_afni.ThreedAutomask(),
                               name='func_automask',
                               iterfield=['in_file'])
    func_automask.inputs.dilate = 1

    func_calcR = pe.MapNode(interface=afni.Calc(),
                            name='func_calcR',
                            iterfield=['in_file_a', 'in_file_b'])
    func_calcR.inputs.expr = '\'a*b\''

    func_mean = pe.MapNode(interface=afni.TStat(),
                           name='func_mean',
                           iterfield=['in_file'])
    func_mean.inputs.options = '-mean'

    func_scale = pe.MapNode(interface=fsl.ImageMaths(),
                            name='func_scale',
                            iterfield=['in_file'])
    func_scale.inputs.op_string = '-ing 10000'
    func_scale.inputs.out_data_type = 'float'

    func_mask = pe.MapNode(interface=fsl.ImageMaths(),
                           name='func_mask',
                           iterfield=['in_file'])
    func_mask.inputs.op_string = '-Tmin -bin'
    func_mask.inputs.out_data_type = 'char'

    preproc.connect(inputNode, 'rest',
                    func_calc, 'in_file_a')
    preproc.connect(inputNode, 'start_idx',
                    func_calc, 'start_idx')
    preproc.connect(inputNode, 'stop_idx',
                    func_calc, 'stop_idx')
    preproc.connect(func_calc, 'out_file',
                    func_refit, 'in_file')
    preproc.connect(func_refit, 'out_file',
                    func_reorient, 'in_file')
    preproc.connect(func_reorient, 'out_file',
                    func_tstat, 'in_file')
    preproc.connect(func_reorient, 'out_file',
                    func_volreg, 'in_file')
    preproc.connect(func_tstat, 'out_file',
                    func_volreg, 'basefile')
    preproc.connect(func_volreg, 'out_file',
                    func_tstat_1, 'in_file')
    preproc.connect(func_reorient, 'out_file',
                    func_volreg_1, 'in_file')
    preproc.connect(func_tstat_1, 'out_file',
                    func_volreg_1, 'basefile')
    preproc.connect(func_volreg_1, 'out_file',
                    func_automask, 'in_file')
    preproc.connect(func_volreg_1, 'out_file',
                    func_calcR, 'in_file_a')
    preproc.connect(func_automask, 'out_file',
                    func_calcR, 'in_file_b')
    preproc.connect(func_calcR, 'out_file',
                    func_mean, 'in_file')
    preproc.connect(func_calcR, 'out_file',
                    func_scale, 'in_file')
    preproc.connect(func_scale, 'out_file',
                    func_mask, 'in_file')

    preproc.connect(func_calc, 'out_file',
                    outputNode, 'drop_tr')
    preproc.connect(func_refit, 'out_file',
                    outputNode, 'refit')
    preproc.connect(func_reorient, 'out_file',
                    outputNode, 'reorient')
    preproc.connect(func_tstat_1, 'out_file',
                    outputNode, 'motion_correct_ref')
    preproc.connect(func_volreg_1, 'out_file',
                    outputNode, 'motion_correct')
    preproc.connect(func_volreg_1, 'md1d_file',
                    outputNode, 'max_displacement')
    preproc.connect(func_volreg_1, 'oned_file',
                    outputNode, 'movement_parameters')
    preproc.connect(func_automask, 'out_file',
                    outputNode, 'mask')
    preproc.connect(func_calcR, 'out_file',
                    outputNode, 'skullstrip')
    preproc.connect(func_mean, 'out_file',
                    outputNode, 'example_func')
    preproc.connect(func_scale, 'out_file',
                    outputNode, 'preprocessed')
    preproc.connect(func_mask, 'out_file',
                    outputNode, 'preprocessed_mask')

    return preproc


def create_reg_preproc():

    preproc = pe.Workflow(name='regpreproc')
    inputNode = pe.Node(util.IdentityInterface(fields=['example_func',
                                                'brain',
                                                'standard_res_brain',
                                                'reorient',
                                                'standard',
                                                'standard_brain_mask_dil',
                                                'config_file']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['example_func2highres',
                                                'example_func2highres_mat',
                                                'highres2example_func_mat',
                                                'highres2standard',
                                                'highres2standard_mat',
                                                'standard2highres_mat',
                                                'highres2standard_warp',
                                                'highres2standard_NL',
                                                'highres2standard_jac',
                                                'stand2highres_warp',
                                                'example_func2standard_NL']),
                        name='outputspec')

    reg_flirt = pe.MapNode(interface=fsl.FLIRT(),
                           name='reg_flirt',
                           iterfield=['in_file'])
    reg_flirt.inputs.cost = 'corratio'
    reg_flirt.inputs.dof = 6
    reg_flirt.inputs.interp = 'nearestneighbour'

    # Create mat file for conversion from subject's anatomical to functional
    reg_xfm1 = pe.MapNode(interface=fsl.ConvertXFM(),
                          name='reg_xfm1',
                          iterfield=['in_file'])
    reg_xfm1.inputs.invert_xfm = True

    ## T1->STANDARD
    ## NOTE THAT THIS IS Linear registration,
    ## you may want to use FNIRT (non-linear)
    reg_flirt1 = pe.Node(interface=fsl.FLIRT(),
                         name='reg_flirt1')
    reg_flirt1.inputs.cost = 'corratio'
    reg_flirt1.inputs.cost_func = 'corratio'
    reg_flirt1.inputs.dof = 12
    reg_flirt1.inputs.interp = 'nearestneighbour'

    ## Create mat file for conversion from standard to high res
    reg_xfm2 = pe.Node(interface=fsl.ConvertXFM(),
                       name='reg_xfm2')
    reg_xfm2.inputs.invert_xfm = True

    reg_inw = pe.Node(interface=e_afni.InvWarp(),
                      name='reg_inw')

    ## T1->STANDARD NONLINEAR
    # Perform nonlinear registration (higres to standard)
    reg_fnt = pe.Node(interface=fsl.FNIRT(),
                      name='reg_fnt')
    reg_fnt.inputs.fieldcoeff_file = True
    reg_fnt.inputs.jacobian_file = True
    reg_fnt.inputs.warp_resolution = (10, 10, 10)

    ## Apply nonlinear registration (func to standard)
    reg_warp = pe.MapNode(interface=fsl.ApplyWarp(),
                          name='reg_warp',
                          iterfield=['in_file',
                                     'premat'
                          ])

    preproc.connect(inputNode, 'example_func',
                    reg_flirt, 'in_file')
    preproc.connect(inputNode, 'brain',
                    reg_flirt, 'reference')
    preproc.connect(reg_flirt, 'out_matrix_file',
                    reg_xfm1, 'in_file')
    preproc.connect(inputNode, 'brain',
                    reg_flirt1, 'in_file')
    preproc.connect(inputNode, 'standard_res_brain',
                    reg_flirt1, 'reference')
    preproc.connect(reg_flirt1, 'out_matrix_file',
                    reg_xfm2, 'in_file')
    preproc.connect(inputNode, 'reorient',
                    reg_fnt, 'in_file')
    preproc.connect(reg_flirt1, 'out_matrix_file',
                    reg_fnt, 'affine_file')
    preproc.connect(inputNode, 'standard',
                    reg_fnt, 'ref_file')
    preproc.connect(inputNode, 'standard_brain_mask_dil',
                    reg_fnt, 'refmask_file')
    preproc.connect(inputNode, 'config_file',
                    reg_fnt, 'config_file')
    preproc.connect(reg_fnt, 'fieldcoeff_file',
                    reg_inw, 'in_file')
    preproc.connect(inputNode, 'brain',
                    reg_inw, 'ref_file')
    preproc.connect(inputNode, 'example_func',
                    reg_warp, 'in_file')
    preproc.connect(inputNode, 'standard',
                    reg_warp, 'ref_file')
    preproc.connect(reg_fnt, 'fieldcoeff_file',
                    reg_warp, 'field_file')
    preproc.connect(reg_flirt, 'out_matrix_file',
                    reg_warp, 'premat')

    preproc.connect(reg_flirt, 'out_file',
                    outputNode, 'example_func2highres')
    preproc.connect(reg_flirt, 'out_matrix_file',
                    outputNode, 'example_func2highres_mat')
    preproc.connect(reg_xfm1, 'out_file',
                    outputNode, 'highres2example_func_mat')
    preproc.connect(reg_warp, 'out_file',
                    outputNode, 'example_func2standard_NL')
    preproc.connect(reg_flirt1, 'out_file',
                    outputNode, 'highres2standard')
    preproc.connect(reg_flirt1, 'out_matrix_file',
                    outputNode, 'highres2standard_mat')
    preproc.connect(reg_inw, 'out_file',
                    outputNode, 'stand2highres_warp')
    preproc.connect(reg_fnt, 'jacobian_file',
                    outputNode, 'highres2standard_jac')
    preproc.connect(reg_fnt, 'fieldcoeff_file',
                    outputNode, 'highres2standard_warp')
    preproc.connect(reg_fnt, 'warped_file',
                    outputNode, 'highres2standard_NL')

    return preproc


def create_seg_preproc():

    preproc = pe.Workflow(name='segpreproc')
    inputNode = pe.Node(util.IdentityInterface(fields=['preprocessed_mask',
                                                'brain',
                                                'standard_res_brain',
                                                'example_func',
                                                'highres2example_func_mat',
                                                'stand2highres_warp',
                                                'PRIOR_CSF',
                                                'PRIOR_GRAY',
                                                'PRIOR_WHITE']),
                        name='inputspec')

    inputnode_csf_threshold = pe.Node(util.IdentityInterface(
                                    fields=['csf_threshold']),
                             name='csf_threshold')

    inputnode_wm_threshold = pe.Node(util.IdentityInterface(
                                    fields=['wm_threshold']),
                             name='wm_threshold')

    inputnode_gm_threshold = pe.Node(util.IdentityInterface(
                                    fields=['gm_threshold']),
                             name='gm_threshold')

    outputNode = pe.Node(util.IdentityInterface(fields=['csf_t12func',
                                                    'csf_mni2func',
                                                    'csf_combo',
                                                    'csf_bin',
                                                    'csf_mask',
                                                    'gm_t12func',
                                                    'gm_mni2func',
                                                    'gm_combo',
                                                    'gm_bin',
                                                    'gm_mask',
                                                    'global_mask',
                                                    'wm_t12func',
                                                    'wm_mni2func',
                                                    'wm_combo',
                                                    'wm_bin',
                                                    'probability_maps',
                                                    'mixeltype',
                                                    'partial_volume_map',
                                                    'partial_volume_files',
                                                    'wm_mask']),
                        name='outputspec')

    def form_threshold_string(threshold):

        return '-thr %f -bin ' % (threshold)

    seg_segment = pe.Node(interface=fsl.FAST(),
                          name='seg_segment')
    seg_segment.inputs.img_type = 1
    seg_segment.inputs.segments = True
    seg_segment.inputs.probability_maps = True
    seg_segment.inputs.out_basename = 'segment'

    seg_copy = pe.MapNode(interface=afni.Copy(),
                          name='seg_copy',
                          iterfield=['in_file'])

    seg_flirt = pe.MapNode(interface=fsl.FLIRT(),
                           name='seg_flirt',
                           iterfield=['reference',
                           'in_matrix_file'])
    seg_flirt.inputs.apply_xfm = True

    seg_warp = pe.MapNode(interface=fsl.ApplyWarp(),
                          name='seg_warp',
                          iterfield=['ref_file',
                          'postmat'])
    seg_warp.inputs.interp = 'nn'

    seg_warp1 = seg_warp.clone('seg_warp1')

    seg_smooth1 = pe.MapNode(interface=fsl.MultiImageMaths(),
                             name='seg_smooth1',
                             iterfield=['in_file',
                             'operand_files'])
    seg_str1 = '-mas %s '
    seg_smooth1.inputs.op_string = seg_str1

    seg_thresh = pe.MapNode(interface=fsl.ImageMaths(),
                            name='seg_thresh',
                            iterfield=['in_file'])

    seg_mask = pe.MapNode(interface=fsl.MultiImageMaths(),
                          name='seg_mask',
                          iterfield=['in_file',
                          'operand_files'])
    seg_str1 = '-mas %s '
    seg_mask.inputs.op_string = seg_str1

    seg_prior1 = pe.MapNode(interface=fsl.MultiImageMaths(),
                            name='seg_prior1',
                            iterfield=['in_file',
                            'operand_files'])
    seg_str1 = '-mas %s '
    seg_prior1.inputs.op_string = seg_str1

    seg_thresh1 = pe.MapNode(interface=fsl.ImageMaths(),
                             name='seg_thresh1',
                             iterfield=['in_file'])
    #seg_str1 = '-thr 0.66 -bin '
    #seg_thresh1.inputs.op_string = seg_str1

    seg_flirt3 = pe.MapNode(interface=fsl.FLIRT(),
                            name='seg_flirt3',
                            iterfield=['reference',
                            'in_matrix_file'])
    seg_flirt3.inputs.apply_xfm = True

    seg_mask1 = seg_mask.clone('seg_mask1')

    seg_flirt4 = seg_flirt3.clone('seg_flirt4')
    seg_thresh2 = seg_thresh.clone('seg_thresh2')
    #seg_thresh2.inputs.op_string = '-thr 0.2 -bin '
    seg_warp2 = seg_warp1.clone('seg_warp2')
    seg_prior2 = seg_prior1.clone('seg_prior2')
    seg_mask2 = seg_mask.clone('seg_mask2')

    preproc.connect(inputNode, 'brain',
                    seg_segment, 'in_files')
    preproc.connect(seg_segment, ('probability_maps', pick_wm_0),
                    seg_flirt, 'in_file')
    preproc.connect(inputNode, 'example_func',
                    seg_flirt, 'reference')
    preproc.connect(inputNode, 'highres2example_func_mat',
                    seg_flirt, 'in_matrix_file')
    preproc.connect(inputNode, 'example_func',
                    seg_warp, 'ref_file')
    preproc.connect(inputNode, 'stand2highres_warp',
                    seg_warp, 'field_file')
    preproc.connect(inputNode, 'PRIOR_CSF',
                    seg_warp, 'in_file')
    preproc.connect(inputNode, 'highres2example_func_mat',
                    seg_warp, 'postmat')
    preproc.connect(seg_flirt, 'out_file',
                    seg_smooth1, 'in_file')
    preproc.connect(seg_warp, 'out_file',
                    seg_smooth1, 'operand_files')
    preproc.connect(seg_smooth1, 'out_file',
                    seg_thresh, 'in_file')
    preproc.connect(inputnode_csf_threshold,
                    ('csf_threshold', form_threshold_string),
                    seg_thresh, 'op_string')
    preproc.connect(seg_thresh, 'out_file',
                    seg_mask, 'in_file')
    preproc.connect(inputNode, 'preprocessed_mask',
                    seg_copy, 'in_file')
    preproc.connect(seg_copy, 'out_file',
                    seg_mask, 'operand_files')
    preproc.connect(seg_segment,
                    ('probability_maps', pick_wm_2),
                    seg_flirt3, 'in_file')
    preproc.connect(inputNode, 'example_func',
                    seg_flirt3, 'reference')
    preproc.connect(inputNode, 'highres2example_func_mat',
                    seg_flirt3, 'in_matrix_file')
    preproc.connect(inputNode, 'example_func',
                    seg_warp1, 'ref_file')
    preproc.connect(inputNode, 'stand2highres_warp',
                    seg_warp1, 'field_file')
    preproc.connect(inputNode, 'PRIOR_WHITE',
                    seg_warp1, 'in_file')
    preproc.connect(inputNode, 'highres2example_func_mat',
                    seg_warp1, 'postmat')
    preproc.connect(seg_flirt3, 'out_file',
                    seg_prior1, 'in_file')
    preproc.connect(seg_warp1, 'out_file',
                    seg_prior1, 'operand_files')
    preproc.connect(seg_prior1, 'out_file',
                    seg_thresh1, 'in_file')
    preproc.connect(inputnode_wm_threshold,
                    ('wm_threshold', form_threshold_string),
                    seg_thresh1, 'op_string')
    preproc.connect(seg_thresh1, 'out_file',
                    seg_mask1, 'in_file')
    preproc.connect(seg_copy, 'out_file',
                    seg_mask1, 'operand_files')

    ##gray matter mask
    preproc.connect(seg_segment,
                    ('probability_maps', pick_wm_1),
                    seg_flirt4, 'in_file')
    preproc.connect(inputNode, 'example_func',
                    seg_flirt4, 'reference')
    preproc.connect(inputNode, 'highres2example_func_mat',
                    seg_flirt4, 'in_matrix_file')
    preproc.connect(inputNode, 'example_func',
                    seg_warp2, 'ref_file')
    preproc.connect(inputNode, 'stand2highres_warp',
                    seg_warp2, 'field_file')
    preproc.connect(inputNode, 'PRIOR_GRAY',
                    seg_warp2, 'in_file')
    preproc.connect(inputNode, 'highres2example_func_mat',
                    seg_warp2, 'postmat')
    preproc.connect(seg_flirt4, 'out_file',
                    seg_prior2, 'in_file')
    preproc.connect(seg_warp2, 'out_file',
                    seg_prior2, 'operand_files')
    preproc.connect(seg_prior2, 'out_file',
                    seg_thresh2, 'in_file')
    preproc.connect(inputnode_gm_threshold,
                    ('gm_threshold', form_threshold_string),
                    seg_thresh2, 'op_string')
    preproc.connect(seg_thresh2, 'out_file',
                    seg_mask2, 'in_file')
    preproc.connect(seg_copy, 'out_file',
                    seg_mask2, 'operand_files')

    preproc.connect(seg_segment, 'probability_maps',
                    outputNode, 'probability_maps')
    preproc.connect(seg_segment, 'mixeltype',
                    outputNode, 'mixeltype')
    preproc.connect(seg_segment, 'partial_volume_files',
                    outputNode, 'partial_volume_files')
    preproc.connect(seg_segment, 'partial_volume_map',
                    outputNode, 'partial_volume_map')
    preproc.connect(seg_flirt, 'out_file',
                    outputNode, 'csf_t12func')
    preproc.connect(seg_warp, 'out_file',
                    outputNode, 'csf_mni2func')
    preproc.connect(seg_smooth1, 'out_file',
                    outputNode, 'csf_combo')
    preproc.connect(seg_thresh, 'out_file',
                    outputNode, 'csf_bin')
    preproc.connect(seg_mask, 'out_file',
                    outputNode, 'csf_mask')
    preproc.connect(seg_copy, 'out_file',
                    outputNode, 'global_mask')
    preproc.connect(seg_flirt3, 'out_file',
                    outputNode, 'wm_t12func')
    preproc.connect(seg_warp1, 'out_file',
                    outputNode, 'wm_mni2func')
    preproc.connect(seg_prior1, 'out_file',
                    outputNode, 'wm_combo')
    preproc.connect(seg_thresh1, 'out_file',
                    outputNode, 'wm_bin')
    preproc.connect(seg_mask1, 'out_file',
                    outputNode, 'wm_mask')
    preproc.connect(seg_flirt4, 'out_file',
                    outputNode, 'gm_t12func')
    preproc.connect(seg_warp2, 'out_file',
                    outputNode, 'gm_mni2func')
    preproc.connect(seg_prior2, 'out_file',
                    outputNode, 'gm_combo')
    preproc.connect(seg_thresh2, 'out_file',
                    outputNode, 'gm_bin')
    preproc.connect(seg_mask2, 'out_file',
                    outputNode, 'gm_mask')

    return preproc


def create_scrubbing_preproc():


    sc = pe.Workflow(name='sc_preproc')
    inputNode = pe.Node(util.IdentityInterface(fields=['rest',
                                                    'movement_parameters',
                                                    'preprocessed'
                                                    ]),
                        name='inputspec')

    inputnode_threshold = pe.Node(util.IdentityInterface(fields=['threshold']),
                             name='threshold_input')

    outputNode = pe.Node(util.IdentityInterface(fields=['mean_deriv_sq_1D',
                                                            'mean_raw_sq_1D',
                                                            'scrubbed_preprocessed',
                                                            'temp_deriv_brik_file',
                                                            'temp_deriv_head_file',
                                                            'temp_deriv_sq_brik_file',
                                                            'temp_deriv_sq_head_file',
                                                            'raw_sq_brik_file',
                                                            'raw_sq_head_file',
                                                            'mask_brik_file',
                                                            'mask_head_file',
                                                            'FD_1D',
                                                            'sqrt_mean_deriv_sq_1D',
                                                            'sqrt_mean_raw_sq_1D',
                                                            'frames_ex_1D',
                                                            'frames_in_1D',
                                                            'pow_params',
                                                            'ftof_percent_change_1D',
                                                            'scrubbed_movement_parameters']),
                        name='outputspec')

    NVOLS = pe.Node(util.Function(input_names=['in_files'], 
                                  output_names=['nvols'],
                                  function=getImgNVols), 
                    name='NVOLS')

    sc_copy = pe.MapNode(util.Function(input_names=['in_file'], 
                                       output_names=['out_file'],
                                       function=scCopy), 
                         name='sc_copy', 
                         iterfield=["in_file"])

    sc_createSC = pe.MapNode(util.Function(input_names=['in_file'], 
                                           output_names=['out_file'],
                                           function=createSC), 
                             name='sc_createSC', 
                             iterfield=["in_file"])

    sc_MeanFD = pe.MapNode(util.Function(input_names=['infile_a', 'infile_b'], 
                                         output_names=['out_file'],
                                         function=setMeanFD), 
                           name='sc_MeanFD', 
                           iterfield=["infile_a", "infile_b"])

    sc_NumFD = pe.MapNode(util.Function(input_names=['infile_a', 'infile_b', 'threshold'], 
                                        output_names=['out_file'],
                                        function=setNumFD), 
                          name='sc_NumFD', 
                          iterfield=["infile_a", "infile_b"])

    sc_PercentFD = pe.MapNode(util.Function(input_names=['infile_a', 'infile_b', 'threshold'], 
                                            output_names=['out_file'],
                                            function=setPercentFD), 
                              name='sc_PercentFD', 
                              iterfield=["infile_a", "infile_b"])

    sc_FramesEx = pe.MapNode(util.Function(input_names=['in_file', 'threshold'], 
                                           output_names=['out_file'],
                                           function=setFramesEx), 
                             name='sc_FramesEx', 
                             iterfield=["in_file"])

    sc_FramesIN = pe.MapNode(util.Function(input_names=['in_file', 'threshold','exclude_list'], 
                                           output_names=['out_file'],
                                           function=setFramesIN), 
                             name='sc_FramesIN', 
                             iterfield=["in_file", "exclude_list"])

    sc_FramesInList = pe.MapNode(util.Function(input_names=['in_file'], 
                                               output_names=['out_file'],
                                               function=setFramesInList), 
                                 name='sc_FramesInList', 
                                 iterfield=["in_file"])

    sc_MeanDVARS = pe.MapNode(util.Function(input_names=['infile_a', 'infile_b'], 
                                            output_names=['out_file'],
                                            function=setMeanDVARS), 
                              name='sc_MeanDVARS', 
                              iterfield=["infile_a", "infile_b"])

    sc_NUM5 = pe.MapNode(util.Function(input_names=['infile_a', 'infile_b'], 
                                       output_names=['out_file'],
                                       function=setNUM5), 
                         name='sc_NUM5', 
                         iterfield=["infile_a", "infile_b"])

    sc_NUM10 = pe.MapNode(util.Function(input_names=['infile_a', 'infile_b'], 
                                        output_names=['out_file'],
                                        function=setNUM10), 
                          name='sc_NUM10', 
                          iterfield=["infile_a", "infile_b"])

    sc_NUMFD = pe.MapNode(util.Function(input_names=['in_file', 'threshold'], 
                                        output_names=['out_file'],
                                        function=setNUMFD), 
                          name='sc_NUMFD', 
                          iterfield=["in_file"])

    sc_ScrubbedMotion = pe.MapNode(util.Function(input_names=['infile_a', 'infile_b'], 
                                                 output_names=['out_file'],
                                                 function=setScrubbedMotion), 
                                   name='sc_ScrubbedMotion',
                                   iterfield=["infile_a", "infile_b"])

    sc_FtoFPercentChange = pe.MapNode(util.Function(input_names=['infile_a', 'infile_b'], 
                                                    output_names=['out_file'],
                                                    function=setFtoFPercentChange), 
                                      name='sc_FtoFPercentChange', 
                                      iterfield=["infile_a", "infile_b"])

    sc_SqrtMeanDeriv = pe.MapNode(util.Function(input_names=['in_file'], output_names=['out_file'],
                                                function=setSqrtMeanDeriv), 
                                  name='sc_SqrtMeanDeriv',
                                  iterfield=["in_file"])

    sc_SqrtMeanRaw = pe.MapNode(util.Function(input_names=['in_file'], 
                                              output_names=['out_file'],
                                              function=setSqrtMeanRaw), 
                                name='sc_SqrtMeanRaw', 
                                iterfield=["in_file"])

    sc_calc1 = pe.MapNode(interface=e_afni.Threedcalc(), 
                          name='sc_calc1',
                          iterfield=["infile_a", "stop_idx", "infile_b", "stop_idx2"])
    sc_calc1.inputs.start_idx = 4
    sc_calc1.inputs.start_idx2 = 3
    sc_calc1.inputs.expr = '\'(a-b)\''
    sc_calc1.inputs.out_file = 'temp_deriv'

    sc_calc2 = pe.MapNode(interface=e_afni.Threedcalc(), 
                          name='sc_calc2', 
                          iterfield=["infile_a"])
    sc_calc2.inputs.expr = '\'a*a\''
    sc_calc2.inputs.out_file = 'temp_deriv_sq'

    sc_calc3 = pe.MapNode(interface=e_afni.Threedcalc(), 
                          name='sc_calc3', 
                          iterfield=["infile_a", "stop_idx"])
    sc_calc3.inputs.start_idx = 3
    sc_calc3.inputs.expr = '\'a*a\''
    sc_calc3.inputs.out_file = 'raw_sq'

    sc_calc_scrub = pe.MapNode(interface=e_afni.Threedcalc(), 
                               name='sc_calc_scrub',
                               iterfield=["infile_a", "list_idx"] )
    sc_calc_scrub.inputs.expr = '\'a\''

    sc_automask = pe.MapNode(interface=e_afni.ThreedAutomask(), 
                             name='sc_automask', 
                             iterfield=["in_file"])
    sc_automask.inputs.dilate = 1
    sc_automask.inputs.genbrickhead = True
    sc_automask.inputs.out_file = './mask'


    sc_3dROIstats_1 = pe.MapNode(interface=e_afni.ThreedROIstats(), 
                                 name='sc_3dROIstats_1',
                                 iterfield=["in_file", "mask"])
    sc_3dROIstats_1.inputs.quiet = True

    sc_3dROIstats_2 = pe.MapNode(interface=e_afni.ThreedROIstats(), 
                                 name='sc_3dROIstats_2',
                                 iterfield=["in_file", "mask"])
    sc_3dROIstats_2.inputs.quiet = True

    sc.connect(inputNode, 'rest', NVOLS, 'in_files')
    sc.connect(inputNode, 'rest', sc_copy, 'in_file' )

    sc.connect(inputNode, 'rest', sc_calc1, 'infile_a')
    sc.connect(inputNode, 'rest', sc_calc1, 'infile_b')
    sc.connect(NVOLS, ('nvols', last_vol), sc_calc1, 'stop_idx')
    sc.connect(NVOLS, ('nvols', TRendminus1), sc_calc1, 'stop_idx2')

    sc.connect(sc_calc1, 'brik_file', sc_calc2, 'infile_a')

    sc.connect(inputNode, 'rest', sc_calc3, 'infile_a')
    sc.connect(NVOLS, ('nvols', TRendminus1), sc_calc3, 'stop_idx')

    sc.connect(inputNode, 'rest', sc_automask, 'in_file')

    sc.connect(sc_automask, 'brik_file', sc_3dROIstats_1, 'mask')
    sc.connect(sc_calc2, 'brik_file', sc_3dROIstats_1, 'in_file')

    sc.connect(sc_automask, 'brik_file', sc_3dROIstats_2, 'mask')
    sc.connect(sc_calc3, 'brik_file', sc_3dROIstats_2, 'in_file')

    sc.connect(sc_3dROIstats_1, 'stats', sc_SqrtMeanDeriv, 'in_file' )
    sc.connect(sc_3dROIstats_2, 'stats', sc_SqrtMeanRaw, 'in_file')

    sc.connect(sc_SqrtMeanDeriv, 'out_file', sc_FtoFPercentChange, 'infile_a' )
    sc.connect(sc_SqrtMeanRaw, 'out_file', sc_FtoFPercentChange, 'infile_b')

    ###Calculating mean Framewise Displacement
    sc.connect(inputNode, 'movement_parameters', sc_createSC, 'in_file' )
    sc.connect(sc_createSC, 'out_file', sc_MeanFD, 'infile_b' )
    sc.connect(sc_copy, 'out_file', sc_MeanFD, 'infile_a' )

    ##NUMBER OF FRAMES >0.5mm FD
    sc.connect(sc_MeanFD, 'out_file', sc_NumFD, 'infile_a')
    sc.connect(sc_createSC, 'out_file', sc_NumFD, 'infile_b')
    sc.connect(inputnode_threshold, 'threshold', sc_NumFD, 'threshold')

    ##NUMBER OF FRAMES >0.5mm FD as percentage of total num frames
    sc.connect(sc_NumFD, 'out_file', sc_PercentFD, 'infile_a')
    sc.connect(sc_createSC, 'out_file', sc_PercentFD, 'infile_b')
    sc.connect(inputnode_threshold, 'threshold', sc_PercentFD, 'threshold')
    
    ####Mean DVARS
    sc.connect(sc_PercentFD, 'out_file', sc_MeanDVARS, 'infile_a')
    sc.connect(sc_FtoFPercentChange, 'out_file', sc_MeanDVARS, 'infile_b')

    ###NUMBER OF relative FRAMES >5%
    sc.connect(sc_MeanDVARS, 'out_file', sc_NUM5, 'infile_a')
    sc.connect(sc_FtoFPercentChange, 'out_file', sc_NUM5, 'infile_b')

    ###NUMBER OF relative FRAMES >10%
    sc.connect(sc_NUM5, 'out_file', sc_NUM10, 'infile_a')
    sc.connect(sc_FtoFPercentChange, 'out_file', sc_NUM10, 'infile_b')

    ##WHAT FRAMES HAVE >0.5mm FD??
    ## FD timeseries starts at the second TR because it's a derivative
    ## but because images start at 0, it's ok to take the frame number directly from the FD file 
    ## (i.e., if the 5th number in the FD file indicates a bad frame, removing timepoint "5" will correctly remove the 6th frame).
    sc.connect(sc_createSC, 'out_file', sc_FramesEx, 'in_file')
    sc.connect(inputnode_threshold, 'threshold', sc_FramesEx, 'threshold')
    
    sc.connect(sc_FramesIN, 'out_file', sc_FramesInList, 'in_file')
    
    sc.connect(sc_createSC, 'out_file', sc_FramesIN, 'in_file')
    sc.connect(inputnode_threshold, 'threshold', sc_FramesIN, 'threshold')
    sc.connect(sc_FramesEx, 'out_file', sc_FramesIN, 'exclude_list')

    sc.connect(inputNode, 'preprocessed', sc_calc_scrub, 'infile_a')
    sc.connect(sc_FramesIN, ('out_file', getIndx), sc_calc_scrub, 'list_idx')

    sc.connect(inputNode, 'movement_parameters', sc_ScrubbedMotion, 'infile_b')
    sc.connect(sc_FramesInList, 'out_file', sc_ScrubbedMotion, 'infile_a' )

    sc.connect(sc_3dROIstats_1, 'stats', outputNode, 'mean_deriv_sq_1D')
    sc.connect(sc_3dROIstats_2, 'stats', outputNode, 'mean_raw_sq_1D')
    sc.connect(sc_calc_scrub, 'out_file', outputNode, 'scrubbed_preprocessed')
    sc.connect(sc_calc1, 'brik_file', outputNode, 'temp_deriv_brik_file')
    sc.connect(sc_calc1, 'head_file', outputNode, 'temp_deriv_head_file')
    sc.connect(sc_calc2, 'brik_file', outputNode, 'temp_deriv_sq_brik_file')
    sc.connect(sc_calc2, 'head_file', outputNode, 'temp_deriv_sq_head_file')
    sc.connect(sc_calc3, 'brik_file', outputNode, 'raw_sq_brik_file')
    sc.connect(sc_calc3, 'head_file', outputNode, 'raw_sq_head_file')
    sc.connect(sc_automask, 'brik_file', outputNode, 'mask_brik_file')
    sc.connect(sc_automask, 'head_file', outputNode, 'mask_head_file')
    sc.connect(sc_createSC, 'out_file', outputNode, 'FD_1D')
    sc.connect(sc_SqrtMeanDeriv, 'out_file', outputNode, 'sqrt_mean_deriv_sq_1D')
    sc.connect(sc_SqrtMeanRaw, 'out_file', outputNode, 'sqrt_mean_raw_sq_1D')
    sc.connect(sc_FramesEx, 'out_file', outputNode, 'frames_ex_1D')
    sc.connect(sc_FramesIN, 'out_file', outputNode, 'frames_in_1D')
    sc.connect(sc_NUM10, 'out_file', outputNode, 'pow_params')
    sc.connect(sc_FtoFPercentChange, 'out_file', outputNode, 'ftof_percent_change_1D')
    sc.connect(sc_ScrubbedMotion, 'out_file', outputNode, 'scrubbed_movement_parameters')

    return sc



def create_vmhc_preproc():

    vmhc = pe.Workflow(name='vmhc_preproc')
    inputNode = pe.Node(util.IdentityInterface(fields=['brain',
                                                'brain_symmetric',
                                                'rest_res',
                                                'reorient',
                                                'example_func2highres_mat',
                                                'symm_standard',
                                                'twomm_brain_mask_dil',
                                                'config_file_twomm',
                                                'standard']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['highres2symmstandard',
                                                'highres2symmstandard_mat',
                                                'highres2symmstandard_warp',
                                                'fnirt_highres2symmstandard',
                                                'highres2symmstandard_jac',
                                                'rest_res_2symmstandard',
                                                'VMHC_img',
                                                'VMHC_Z_img',
                                                'VMHC_Z_stat_img']),
                        name='outputspec')

    ## Linear registration of T1 --> symmetric standard
    flirt = pe.Node(interface=fsl.FLIRT(),
                    name='flirt')
    flirt.inputs.cost = 'corratio'
    flirt.inputs.cost_func = 'corratio'
    flirt.inputs.dof = 12
    flirt.inputs.interp = 'trilinear'

    ## Perform nonlinear registration
    ##(higres to standard) to symmetric standard brain
    fnt = pe.Node(interface=fsl.FNIRT(),
                  name='fnt')
    fnt.inputs.fieldcoeff_file = True
    fnt.inputs.jacobian_file = True
    fnt.inputs.warp_resolution = (10, 10, 10)

    ## Apply nonlinear registration (func to standard)
    warp = pe.MapNode(interface=fsl.ApplyWarp(),
                      name='warp',
                      iterfield=['in_file', 'premat'])

    ## copy and L/R swap file
    swap = pe.MapNode(interface=fsl.SwapDimensions(),
                      name='swap',
                      iterfield=['in_file'])
    swap.inputs.new_dims = ('-x', 'y', 'z')

    ## caculate vmhc
    corr = pe.MapNode(interface=afni.TCorrelate(),
                      name='corr',
                      iterfield=['xset', 'yset'])
    corr.inputs.pearson = True
    corr.inputs.polort = -1

    z_trans = pe.MapNode(interface=e_afni.Threedcalc(),
                         name='z_trans',
                         iterfield=['infile_a'])
    z_trans.inputs.expr = '\'log((1+a)/(1-a))/2\''
    z_trans.inputs.outputtype = 'NIFTI'
    z_stat = pe.MapNode(interface=e_afni.Threedcalc(),
                        name='z_stat',
                        iterfield=['infile_a',
                                   'expr'])
    z_stat.inputs.outputtype = 'NIFTI'

    NVOLS = pe.Node(util.Function(input_names=['in_files'],
                                  output_names=['nvols'],
                    function=getImgNVols),
                    name='NVOLS')

    generateEXP = pe.Node(util.Function(input_names=['nvols'],
                                        output_names=['expr'],
                          function=getEXP),
                          name='generateEXP')

    vmhc.connect(inputNode, 'brain',
                 flirt, 'in_file')
    vmhc.connect(inputNode, 'brain_symmetric',
                 flirt, 'reference')
    vmhc.connect(inputNode, 'reorient',
                 fnt, 'in_file')
    vmhc.connect(flirt, 'out_matrix_file',
                 fnt, 'affine_file')
    vmhc.connect(inputNode, 'symm_standard',
                 fnt, 'ref_file')
    vmhc.connect(inputNode, 'twomm_brain_mask_dil',
                 fnt, 'refmask_file')
    vmhc.connect(inputNode, 'config_file_twomm',
                 fnt, 'config_file')
    vmhc.connect(inputNode, 'rest_res',
                 warp, 'in_file')
    vmhc.connect(inputNode, 'symm_standard',
                 warp, 'ref_file')
    vmhc.connect(fnt, 'fieldcoeff_file',
                 warp, 'field_file')
    vmhc.connect(inputNode, 'example_func2highres_mat',
                 warp, 'premat')

    vmhc.connect(warp, 'out_file',
                 swap, 'in_file')
    vmhc.connect(warp, 'out_file',
                 corr, 'xset')
    vmhc.connect(swap, 'out_file',
                 corr, 'yset')
    vmhc.connect(corr, 'out_file',
                 z_trans, 'infile_a')
    vmhc.connect(swap, 'out_file',
                 NVOLS, 'in_files')
    vmhc.connect(NVOLS, 'nvols',
                 generateEXP, 'nvols')
    vmhc.connect(z_trans, 'out_file',
                 z_stat, 'infile_a')
    vmhc.connect(generateEXP, 'expr',
                 z_stat, 'expr')

    vmhc.connect(flirt, 'out_file',
                 outputNode, 'highres2symmstandard')
    vmhc.connect(flirt, 'out_matrix_file',
                 outputNode, 'highres2symmstandard_mat')
    vmhc.connect(fnt, 'jacobian_file',
                 outputNode, 'highres2symmstandard_jac')
    vmhc.connect(fnt, 'fieldcoeff_file',
                 outputNode, 'highres2symmstandard_warp')
    vmhc.connect(fnt, 'warped_file',
                 outputNode, 'fnirt_highres2symmstandard')
    vmhc.connect(warp, 'out_file',
                 outputNode, 'rest_res_2symmstandard')
    vmhc.connect(corr, 'out_file',
                 outputNode, 'VMHC_img')
    vmhc.connect(z_trans, 'out_file',
                 outputNode, 'VMHC_Z_img')
    vmhc.connect(z_stat, 'out_file',
                 outputNode, 'VMHC_Z_stat_img')

    return vmhc

def reho(in_file, mask_file):

    import nibabel as nb
    from utils import reho_statistic_filter

    lfo = np.load(in_file).get_data().as_type('float64')

    W = reho_statistic_filter(lfo, 0)

    out_file = None


    return out_file


def create_reho_preproc():

    preproc = pe.Workflow(name='reho_preproc')


    op_string = pe.MapNode(util.Function(input_names=['mean', 'std_dev'],
                                         output_names=['op_string'],
                           function=getOpString),
                           name='op_string',
                           iterfield=['mean', 'std_dev'])

    mat = pe.MapNode(util.Function(input_names=['in_file', 'mask_file'],
                                   output_names=['out_file'],
                     function=reho),
                     name='mat',
                     iterfield=['in_file',
                                'mask_file'])

    mean = pe.MapNode(interface=fsl.ImageStats(),
                      name='mean',
                      iterfield=['in_file',
                                 'mask_file'])
    mean.inputs.op_string = '-k %s -m'


    std = pe.MapNode(interface=fsl.ImageStats(),
                     name='std',
                     iterfield=['in_file',
                                'mask_file'])
    std.inputs.op_string = '-k %s -s'


    Z = pe.MapNode(interface=MultiImageMaths(),
                   name='Z',
                   iterfield=['in_file',
                              'operand_files',
                              'op_string'])


    #Registering Z-transformed ReHo to standard space using FLIRT
    Z_F = pe.MapNode(interface=fsl.FLIRT(),
                     name='Z_F',
                     iterfield=['in_file',
                                'in_matrix_file'])
    Z_F.inputs.apply_xfm = True

    #Registering Z-transformed ReHo to standard space using FNIRT
    Z_NF = pe.MapNode(interface=fsl.ApplyWarp(),
                      name='Z_NF',
                      iterfield=['in_file',
                                 'premat'])

    smooth = pe.MapNode(interface=MultiImageMaths(),
                        name='smooth',
                        iterfield=['in_file',
                                   'operand_files'])
    str = '-kernel gauss %f -fmean' % (sigma)
    smooth.inputs.op_string = str


    preproc.connect(nuisance_calc, 'out_file',
                    mat, 'in_file')
    preproc.connect(func_automask, 'out_file',
                    mat, 'mask_file')

    preproc.connect(mat, 'out_file',
                    mean, 'in_file')
    preproc.connect(func_automask, 'out_file',
                    mean, 'mask_file')
    preproc.connect(mat, 'out_file',
                    std, 'in_file')
    preproc.connect(func_automask, 'out_file',
                    std, 'mask_file')
    preproc.connect(mean, 'out_stat',
                    op_string, 'mean')
    preproc.connect(std, 'out_stat',
                    op_string, 'std_dev')
    preproc.connect(mat, 'out_file',
                    Z, 'in_file')
    preproc.connect(op_string, 'op_string',
                    Z, 'op_string')
    preproc.connect(func_automask, 'out_file',
                    Z, 'operand_files')

    preproc.connect(Z, 'out_file',
                    Z_F, 'in_file')
    preproc.connect(infosource, 'standard_res_brain',
                    Z_F, 'reference')
    preproc.connect(reg_xfm3, 'out_file',
                    Z_F, 'in_matrix_file')
    preproc.connect(Z_F, 'out_file',
                    smooth, 'in_file')


    preproc.connect(Z, 'out_file',
                    Z_NF, 'in_file')
    preproc.connect(infosource, 'standard',
                    Z_NF, 'ref_file')
    preproc.connect(reg_fnt, 'fieldcoeff_file',
                    Z_NF, 'field_file')
    preproc.connect(reg_flirt, 'out_matrix_file',
                    Z_NF, 'premat')
    preproc.connect(Z_NF, 'out_file',
                    smooth, 'in_file')

    return preproc


def create_sca_preproc_native():

    rsfc = pe.Workflow(name='sca_preproc_native')
    inputNode = pe.Node(util.IdentityInterface(fields=['ref',
                                                'warp',
                                                'postmat',
                                                'premat',
                                                'rest_res_filt',
                                                'fieldcoeff_file',
                                                'rest_mask2standard',
                                                'standard']),
                        name='inputspec')

    inputnode_fwhm = pe.Node(util.IdentityInterface(fields=['fwhm']),
                             name='fwhm_input')

    inputnode_seed_list = pe.Node(util.IdentityInterface(fields=['seed_list']),
                                  name='seed_list_input')

    outputNode = pe.Node(util.IdentityInterface(fields=['seed_mni2func',
                                                    'correlations',
                                                    'Z_trans_correlations',
                                                    'Z_2standard',
                                                    'Z_2standard_FWHM']),
                        name='outputspec')

    printToFile = pe.MapNode(util.Function(input_names=['time_series'],
                                           output_names=['ts_oneD'],
                             function=pToFile),
                             name='printToFile',
                             iterfield=['time_series'])

    ## 0. Register Seed template in to native space
    warp = pe.MapNode(interface=fsl.ApplyWarp(),
                      name='warp',
                      iterfield=['ref_file',
                      'postmat'])
    warp.inputs.interp = 'nn'
#    warp.iterables = ('in_file', inputnode_seed_list.inputs.seed_list)

    ## 1. Extract Timeseries
    time_series = pe.MapNode(interface=afni.ROIStats(),
                             name='time_series',
                             iterfield=['in_file',
                             'mask'])
    time_series.inputs.quiet = True
    time_series.inputs.mask_f2short = True

    ## 2. Compute voxel-wise correlation with Seed Timeseries
    corr = pe.MapNode(interface=afni.Fim(),
                      name='corr',
                      iterfield=['in_file',
                      'ideal_file'])
    corr.inputs.fim_thr = 0.0009
    corr.inputs.out = 'Correlation'

    ## 3. Z-transform correlations
    z_trans = pe.MapNode(interface=e_afni.Threedcalc(),
                         name='z_trans',
                         iterfield=['infile_a'])
    z_trans.inputs.expr = '\'log((1+a)/(1-a))/2\''

    ## 4. Register Z-transformed correlations to standard space
    register = pe.MapNode(interface=fsl.ApplyWarp(),
                          name='register',
                          iterfield=['premat',
                          'in_file'])

    smooth = pe.MapNode(interface=fsl.MultiImageMaths(),
                        name='smooth',
                        iterfield=['in_file',
                        'operand_files'])

    rsfc.connect(inputNode, 'ref',
                 warp, 'ref_file')
    rsfc.connect(inputNode, 'warp',
                 warp, 'field_file')
    rsfc.connect(inputNode, 'postmat',
                 warp, 'postmat')
    rsfc.connect(inputnode_seed_list, 'seed_list',
                 warp, 'in_file')

    rsfc.connect(inputNode, 'rest_res_filt',
                 time_series, 'in_file')
    rsfc.connect(warp, 'out_file',
                 time_series, 'mask')
    rsfc.connect(time_series, 'stats',
                 printToFile, 'time_series')
    rsfc.connect(printToFile, 'ts_oneD',
                 corr, 'ideal_file')
    rsfc.connect(inputNode, 'rest_res_filt',
                 corr, 'in_file')
    rsfc.connect(corr, 'out_file',
                 z_trans, 'infile_a')
    rsfc.connect(z_trans, 'out_file',
                 register, 'in_file')
    rsfc.connect(inputNode, 'standard',
                 register, 'ref_file')
    rsfc.connect(inputNode, 'fieldcoeff_file',
                 register, 'field_file')
    rsfc.connect(inputNode, 'premat',
                 register, 'premat')
    rsfc.connect(register, 'out_file',
                 smooth, 'in_file')
    rsfc.connect(inputnode_fwhm, ('fwhm', set_gauss),
                 smooth, 'op_string')
    rsfc.connect(inputNode, 'rest_mask2standard',
                 smooth, 'operand_files')

    rsfc.connect(warp, 'out_file',
                 outputNode, 'seed_mni2func')
    rsfc.connect(corr, 'out_file',
                 outputNode, 'correlations')
    rsfc.connect(z_trans, 'out_file',
                 outputNode, 'Z_trans_correlations')
    rsfc.connect(register, 'out_file',
                 outputNode, 'Z_2standard')
    rsfc.connect(smooth, 'out_file',
                 outputNode, 'Z_2standard_FWHM')

    return rsfc


def create_sca_preproc(corr_space):

    rsfc = pe.Workflow(name='sca_preproc')
    inputNode = pe.Node(util.IdentityInterface(fields=[
                                                'premat',
                                                'rest_res_filt',
                                                'fieldcoeff_file',
                                                'residual_file',
                                                'rest_mask2standard',
                                                'standard']),
                        name='inputspec')

    inputnode_fwhm = pe.Node(util.IdentityInterface(fields=['fwhm']),
                             name='fwhm_input')

    inputnode_seed_list = pe.Node(util.IdentityInterface(fields=['seed_list']),
                                  name='seed_list_input')

    outputNode = pe.Node(util.IdentityInterface(fields=[
                                                    'correlations',
                                                    'Z_trans_correlations',
                                                    'Z_2standard',
                                                    'Z_2standard_FWHM']),
                        name='outputspec')

    printToFile = pe.MapNode(util.Function(input_names=['time_series'],
                                           output_names=['ts_oneD'],
                             function=pToFile),
                             name='printToFile',
                             iterfield=['time_series'])

    warp = pe.MapNode(interface=fsl.ApplyWarp(),
                      name='warp',
                      iterfield=['in_file',
                                 'premat'])

    warp_filt = warp.clone('warp_filt')
    ## 1. Extract Timeseries

    time_series = pe.MapNode(interface=afni.ROIStats(),
                             name='time_series',
                             iterfield=['in_file'])
    time_series.inputs.quiet = True
    time_series.inputs.mask_f2short = True
    #time_series.iterables = ("mask",seed_list)


    ## 2. Compute voxel-wise correlation with Seed Timeseries
    corr = pe.MapNode(interface=afni.Fim(),
                      name='corr',
                      iterfield=['in_file',
                      'ideal_file'])
    corr.inputs.fim_thr = 0.0009
    corr.inputs.out = 'Correlation'

    ## 3. Z-transform correlations
    z_trans = pe.MapNode(interface=e_afni.Threedcalc(),
                         name='z_trans',
                         iterfield=['infile_a'])
    z_trans.inputs.expr = '\'log((1+a)/(1-a))/2\''

    ## 4. Register Z-transformed correlations to standard space
    register = pe.MapNode(interface=fsl.ApplyWarp(),
                          name='register',
                          iterfield=['premat',
                          'in_file'])

    smooth = pe.MapNode(interface=fsl.MultiImageMaths(),
                        name='smooth',
                        iterfield=['in_file',
                        'operand_files'])

    rsfc.connect(inputNode, 'rest_res_filt',
                 warp_filt, 'in_file')
    rsfc.connect(inputNode, 'standard',
                 warp_filt, 'ref_file')
    rsfc.connect(inputNode, 'fieldcoeff_file',
                 warp_filt, 'field_file')
    rsfc.connect(inputNode, 'premat',
                 warp_filt, 'premat')
    rsfc.connect(inputNode, 'residual_file',
                 warp, 'in_file')
    rsfc.connect(inputNode, 'standard',
                 warp, 'ref_file')
    rsfc.connect(inputNode, 'fieldcoeff_file',
                 warp, 'field_file')
    rsfc.connect(inputNode, 'premat',
                 warp, 'premat')
    rsfc.connect(warp, 'out_file',
                 time_series, 'in_file')
    rsfc.connect(time_series, 'stats',
                 printToFile, 'time_series')
    rsfc.connect(inputnode_seed_list, 'seed_list',
                time_series, 'mask')
    rsfc.connect(printToFile, 'ts_oneD',
                 corr, 'ideal_file')

    if corr_space == 'native':
        rsfc.connect(inputNode, 'rest_res_filt',
                     corr, 'in_file')
    else:
        rsfc.connect(warp_filt, 'out_file',
                     corr, 'in_file')
    rsfc.connect(corr, 'out_file',
                 z_trans, 'infile_a')
    rsfc.connect(z_trans, 'out_file',
                 register, 'in_file')
    rsfc.connect(inputNode, 'standard',
                 register, 'ref_file')
    rsfc.connect(inputNode, 'fieldcoeff_file',
                 register, 'field_file')
    rsfc.connect(inputNode, 'premat',
                 register, 'premat')
    rsfc.connect(register, 'out_file',
                 smooth, 'in_file')
    rsfc.connect(inputnode_fwhm, ('fwhm', set_gauss),
                 smooth, 'op_string')
    rsfc.connect(inputNode, 'rest_mask2standard',
                 smooth, 'operand_files')

    rsfc.connect(corr, 'out_file',
                 outputNode, 'correlations')
    rsfc.connect(z_trans, 'out_file',
                 outputNode, 'Z_trans_correlations')
    rsfc.connect(register, 'out_file',
                 outputNode, 'Z_2standard')
    rsfc.connect(smooth, 'out_file',
                 outputNode, 'Z_2standard_FWHM')

    return rsfc

def create_group_analysis(f_test):

    grp_analysis = pe.Workflow(name='group_analysis')

    inputnode = pe.Node(util.IdentityInterface(fields=['mat_file',
                                                        'con_file',
                                                        'fts_file',
                                                        'grp_file',
                                                        'zmap_files',
                                                        'z_threshold',
                                                        'p_threshold',
                                                        'parameters']),
                         name='inputspec')

    outputnode = pe.Node(util.IdentityInterface(fields=['merged',
                                                        'zstats',
                                                        'cluster_threshold',
                                                        'cluster_index',
                                                        'overlay_threshold',
                                                        'rendered_image']),
                         name='outputspec')


    inputnode_ftest = pe.Node(util.IdentityInterface(fields=['ftest']),
                             name='ftest_input')

    gp_fslmerge = pe.Node(interface=fsl.Merge(), name='gp_fslmerge')
    gp_fslmerge.inputs.dimension = 't'

    gp_fslmaths = pe.Node(interface=fsl.ImageMaths(), name='gp_fslmaths')
    gp_fslmaths.inputs.op_string = '-abs -Tmin -bin'

    gp_flameo = pe.Node(interface=fsl.FLAMEO(), name='gp_flameo')
    gp_flameo.inputs.run_mode = 'ols'


    gp_smooth_estimate = pe.MapNode(interface=fsl.SmoothEstimate(),
                                    name='gp_smooth_estimate',
                                    iterfield=['zstat_file'])

    gp_fslmultimaths = pe.MapNode(interface=fsl.MultiImageMaths(),
                                  name='gp_fslmultimaths',
                                  iterfield=['in_file'])
    gp_fslmultimaths.inputs.op_string = '-mas %s'

    gp_fslcpgeom = pe.MapNode(util.Function(input_names=['infile_a', 'infile_b'],
                                            output_names=['out_file'],
                                            function=copyGeom),
                              name='gp_fslcpgeom',
                              iterfield=['infile_a', 'infile_b'])

    gp_cluster = pe.MapNode(interface=fsl.Cluster(),
                            name='gp_cluster',
                            iterfield=['in_file', 'volume', 'dlh'])
    gp_cluster.inputs.out_index_file = True
    gp_cluster.inputs.out_threshold_file = True
    gp_cluster.inputs.out_localmax_txt_file = True
    gp_cluster.inputs.out_localmax_vol_file = True
    gp_cluster.inputs.out_size_file = True
    gp_cluster.inputs.out_max_file = True
    gp_cluster.inputs.out_mean_file = True
    gp_cluster.inputs.out_pval_file = True

    gp_fslstats = pe.MapNode(interface=fsl.ImageStats(),
                             name='gp_fslstats',
                             iterfield=['in_file'])
    gp_fslstats.inputs.op_string = '-R'

    gp_merge = pe.MapNode(util.Function(input_names=['infile_a', 'infile_b'],
                                        output_names=['out_file'],
                                        function=getTuple),
                          name='gp_merge', iterfield=['infile_b'])


    gp_overlay = pe.MapNode(interface=fsl.Overlay(),
                            name='gp_overlay',
                            iterfield=['stat_image', 'stat_thresh'])
    gp_overlay.inputs.transparency = True
    gp_overlay.inputs.auto_thresh_bg = True
    gp_overlay.inputs.out_type = 'float'



    gp_slicer = pe.MapNode(interface=fsl.Slicer(), name='gp_slicer',
                           iterfield=['in_file'])
    gp_slicer.inputs.image_width = 750
    gp_slicer.inputs.all_axial = True

    gp_fslmultimaths2 = pe.Node(interface=fsl.ImageMaths(), name='gp_fslmaths2')

    gp_nvols = pe.Node(util.Function(input_names=['in_file'],
                                     output_names=['out_file'], function=get_nvols),
                       name='gp_nvols')

    gp_getbackgroundimage = pe.MapNode(util.Function(input_names=['in_file',
                                                                  'file_parameters'],
                                                     output_names=['out_file'],
                            function=get_standard_background_img),
                            name='gp_getbackgroundimage', iterfield=['in_file'])

    gp_getbackgroundimage2 = pe.Node(util.Function(input_names=['in_file',
                                                                'file_parameters'],
                                                     output_names=['out_file'],
                            function=get_standard_background_img),
                            name='gp_getbackgroundimage2')

    grp_analysis.connect(inputnode, 'zmap_files', gp_fslmerge, 'in_files')
    grp_analysis.connect(gp_fslmerge, 'merged_file', gp_fslmaths, 'in_file')

    if f_test:
        grp_analysis.connect(gp_fslmerge, 'merged_file', gp_flameo, 'cope_file' )
        grp_analysis.connect(gp_fslmaths, 'out_file', gp_flameo, 'mask_file' )
        grp_analysis.connect(inputnode, 'mat_file', gp_flameo, 'design_file' )
        grp_analysis.connect(inputnode, 'con_file', gp_flameo, 't_con_file' )
        grp_analysis.connect(inputnode, 'grp_file', gp_flameo, 'cov_split_file' )
        grp_analysis.connect(inputnode, 'fts_file', gp_flameo, 'f_con_file' )
    else:
        grp_analysis.connect(gp_fslmerge, 'merged_file', gp_flameo, 'cope_file' )
        grp_analysis.connect(gp_fslmaths, 'out_file', gp_flameo, 'mask_file' )
        grp_analysis.connect(inputnode, 'mat_file', gp_flameo, 'design_file' )
        grp_analysis.connect(inputnode, 'con_file', gp_flameo, 't_con_file' )
        grp_analysis.connect(inputnode, 'grp_file', gp_flameo, 'cov_split_file' )


    grp_analysis.connect(gp_flameo, 'zstats', gp_smooth_estimate, 'zstat_file' )
    grp_analysis.connect(gp_fslmaths, 'out_file', gp_smooth_estimate, 'mask_file' )
    grp_analysis.connect(gp_flameo, 'zstats', gp_fslmultimaths, 'in_file')
    grp_analysis.connect(gp_fslmaths, 'out_file', gp_fslmultimaths, 'operand_files')
    grp_analysis.connect(gp_fslmultimaths, 'out_file', gp_getbackgroundimage, 'in_file' )
    grp_analysis.connect(inputnode, 'parameters', gp_getbackgroundimage, 'file_parameters')
    grp_analysis.connect(gp_getbackgroundimage, 'out_file', gp_fslcpgeom, 'infile_a' )
    grp_analysis.connect(gp_fslmultimaths, 'out_file', gp_fslcpgeom, 'infile_b')
    grp_analysis.connect(gp_fslcpgeom, 'out_file', gp_cluster, 'in_file')
    grp_analysis.connect(inputnode, 'z_threshold', gp_cluster, 'threshold')
    grp_analysis.connect(inputnode, 'p_threshold', gp_cluster, 'pthreshold')
    grp_analysis.connect(gp_smooth_estimate, 'volume', gp_cluster, 'volume')
    grp_analysis.connect(gp_smooth_estimate, 'dlh', gp_cluster, 'dlh')
    grp_analysis.connect(gp_cluster, 'threshold_file', gp_fslstats, 'in_file')
    grp_analysis.connect(gp_fslstats, 'out_stat', gp_merge, 'infile_b')
    grp_analysis.connect(inputnode, 'z_threshold', gp_merge, 'infile_a')

    grp_analysis.connect(gp_cluster, 'threshold_file', gp_overlay, 'stat_image')
    grp_analysis.connect(gp_merge, 'out_file', gp_overlay, 'stat_thresh')

    grp_analysis.connect(gp_fslmaths, 'out_file', gp_getbackgroundimage2, 'in_file' )
    grp_analysis.connect(inputnode, 'parameters', gp_getbackgroundimage2, 'file_parameters')
    grp_analysis.connect(gp_getbackgroundimage2, 'out_file', gp_overlay, 'background_image')
    grp_analysis.connect(gp_overlay, 'out_file', gp_slicer, 'in_file')

    grp_analysis.connect(gp_fslmerge, 'merged_file', gp_nvols, 'in_file')
    grp_analysis.connect(gp_fslmerge, 'merged_file', gp_fslmultimaths2, 'in_file' )
    grp_analysis.connect(gp_nvols, 'out_file', gp_fslmultimaths2, 'op_string')

    grp_analysis.connect(gp_fslmerge, 'merged_file', outputnode, 'merged')
    grp_analysis.connect(gp_flameo, 'zstats', outputnode, 'zstats')
    grp_analysis.connect(gp_cluster, 'threshold_file', outputnode, 'cluster_threshold')
    grp_analysis.connect(gp_cluster, 'index_file', outputnode, 'cluster_index')
    grp_analysis.connect(gp_overlay, 'out_file', outputnode, 'overlay_threshold')
    grp_analysis.connect(gp_slicer, 'out_file', outputnode, 'rendered_image')

    return grp_analysis



def create_alff_preproc():

    alff = pe.Workflow(name='alff_preproc')
    inputNode = pe.Node(util.IdentityInterface(fields=['rest_res',
                                                'rest_mask',
                                                'rest_mask2standard',
                                                'premat',
                                                'standard',
                                                'fieldcoeff_file',
                                                    ]),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['processed_alff',
                                            'processed_mean_alff',
                                            'power_spectrum_distribution',
                                            'alff_img',
                                            'falff_img',
                                            'alff_Z_img',
                                            'falff_Z_img',
                                            'alff_Z_2standard_img',
                                            'falff_Z_2standard_img',
                                            'alff_Z_2standard_fwhm_img',
                                            'falff_Z_2standard_fwhm_img']),
                          name='outputspec')
    inputnode_hp = pe.Node(util.IdentityInterface(fields=['hp']),
                             name='hp_input')

    inputnode_lp = pe.Node(util.IdentityInterface(fields=['lp']),
                             name='lp_input')

    inputnode_fwhm = pe.Node(util.IdentityInterface(fields=['fwhm']),
                             name='fwhm_input')

    TR = pe.Node(util.Function(input_names=['in_files'],
                               output_names=['TR'],
                 function=getImgTR), name='TR')
    NVOLS = pe.Node(util.Function(input_names=['in_files'],
                                  output_names=['nvols'],
                    function=getImgNVols),
                    name='NVOLS')

    cp = pe.MapNode(interface=fsl.ImageMaths(),
                    name='cp',
                    iterfield=['in_file'])

    mean = pe.MapNode(interface=fsl.ImageMaths(),
                      name='mean',
                      iterfield=['in_file'])
    mean.inputs.op_string = '-Tmean'

    roi = pe.MapNode(interface=fsl.ExtractROI(),
                     name='roi',
                     iterfield=['in_file',
                     't_size'])
    roi.inputs.t_min = 1

    concatnode = pe.MapNode(interface=util.Merge(2),
                            name='concatnode',
                            iterfield=['in1', 'in2'])

    selectnode = pe.MapNode(interface=util.Select(),
                            name='selectnode',
                            iterfield=['inlist', 'index'])

    pspec = pe.MapNode(interface=fsl.PowerSpectrum(),
                       name='pspec',
                       iterfield=['in_file'])

    ##compute sqrt of power spectrum
    sqrt = pe.MapNode(interface=fsl.ImageMaths(),
                      name='sqrt',
                      iterfield=['in_file'])
    sqrt.inputs.op_string = '-sqrt'

    calcN1 = pe.MapNode(util.Function(input_names=['nvols',
                                      'TR', 'HP'],
                                      output_names=['n1'],
                        function=getN1),
                        name='calcN1',
                        iterfield=['nvols', 'TR'])

    calcN2 = pe.MapNode(util.Function(input_names=['nvols',
                                      'TR', 'LP', 'HP'],
                                      output_names=['n2'],
                        function=getN2),
                        name='calcN2',
                        iterfield=['nvols', 'TR'])
    roi1 = pe.MapNode(interface=fsl.ExtractROI(),
                      name='roi1',
                      iterfield=['in_file',
                                 't_min', 't_size'])

    ## calculate ALFF as the _sum of the amplitudes
    ## in the low frequency band
    _sum = pe.MapNode(interface=fsl.ImageMaths(),
                      name='_sum',
                      iterfield=['in_file',
                      'op_string'])

    ## 4. Calculate fALFF
    falff = pe.MapNode(interface=fsl.ImageMaths(),
                       name='falff',
                       iterfield=['in_file',
                        'op_string'])

    falff1 = pe.MapNode(interface=fsl.MultiImageMaths(),
                        name='falff1',
                        iterfield=['in_file',
                        'operand_files'])
    falff1.inputs.op_string = '-div %s'

    ## 5. Z-normalisation across whole brain
    normM = pe.MapNode(interface=fsl.ImageStats(),
                       name='normM',
                       iterfield=['in_file',
                        'mask_file'])
    normM.inputs.op_string = '-k %s -m'

    normS = pe.MapNode(interface=fsl.ImageStats(),
                       name='normS',
                       iterfield=['in_file',
                       'mask_file'])
    normS.inputs.op_string = '-k %s -s'

    normM1 = pe.MapNode(interface=fsl.ImageStats(),
                        name='normM1',
                        iterfield=['in_file',
                        'mask_file'])
    normM1.inputs.op_string = '-k %s -m'

    normS1 = pe.MapNode(interface=fsl.ImageStats(),
                        name='normS1',
                        iterfield=['in_file',
                        'mask_file'])
    normS1.inputs.op_string = '-k %s -s'

    op_string = pe.MapNode(util.Function(input_names=['mean',
                                         'std_dev'],
                                         output_names=['op_string'],
                           function=getOpString),
                           name='alff_op_string',
                           iterfield=['mean',
                           'std_dev'])

    op_string1 = op_string.clone('op_string1')

    Z_alff = pe.MapNode(interface=fsl.MultiImageMaths(),
                        name='Z_alff',
                        iterfield=['in_file',
                        'operand_files',
                        'op_string'])

    Z_falff = pe.MapNode(interface=fsl.MultiImageMaths(),
                         name='Z_falff',
                         iterfield=['in_file',
                            'operand_files',
                            'op_string'])

    #Registering Z-transformed ALFF to standard space
    warp_alff = pe.MapNode(interface=fsl.ApplyWarp(),
                           name='warp_alff',
                           iterfield=['in_file',
                           'premat'])

    warp_falff = pe.MapNode(interface=fsl.ApplyWarp(),
                            name='warp_falff',
                            iterfield=['in_file',
                            'premat'])

    ## 3. Spatial Smoothing
    smooth = pe.MapNode(interface=fsl.MultiImageMaths(),
                        name='smooth',
                        iterfield=['in_file',
                        'operand_files'])

    fsmooth = pe.MapNode(interface=fsl.MultiImageMaths(),
                         name='fsmooth',
                         iterfield=['in_file',
                         'operand_files'])

    alff.connect(inputNode, 'rest_res',
                 TR, 'in_files')
    alff.connect(inputNode, 'rest_res',
                 NVOLS, 'in_files')
    alff.connect(inputNode, 'rest_res',
                 mean, 'in_file')
    alff.connect(inputNode, 'rest_res',
                 roi, 'in_file')
    alff.connect(NVOLS, 'nvols',
                 roi, 't_size')
    alff.connect(inputNode, 'rest_res',
                 cp, 'in_file')
    alff.connect(roi, 'roi_file',
                 concatnode, 'in1')
    alff.connect(cp, 'out_file',
                 concatnode, 'in2')
    alff.connect(concatnode, 'out',
                 selectnode, 'inlist')
    alff.connect(NVOLS, ('nvols', takemod),
                 selectnode, 'index')
    alff.connect(selectnode, 'out',
                 pspec, 'in_file')
    alff.connect(pspec, 'out_file',
                 sqrt, 'in_file')

    alff.connect(NVOLS, 'nvols',
                 calcN1, 'nvols')
    alff.connect(TR, 'TR',
                 calcN1, 'TR')
    alff.connect(inputnode_hp, 'hp',
                 calcN1, 'HP')

    alff.connect(NVOLS, 'nvols',
                 calcN2, 'nvols')
    alff.connect(TR, 'TR',
                 calcN2, 'TR')
    alff.connect(inputnode_lp, 'lp',
                 calcN2, 'LP')
    alff.connect(inputnode_hp, 'hp',
                 calcN2, 'HP')

    alff.connect(sqrt, 'out_file',
                 roi1, 'in_file')
    alff.connect(calcN1, 'n1',
                 roi1, 't_min')
    alff.connect(calcN2, 'n2',
                 roi1, 't_size')
    alff.connect(roi1, 'roi_file',
                 _sum, 'in_file')
    alff.connect(calcN2, ('n2', set_op_str),
                 _sum, 'op_string')

    alff.connect(sqrt, 'out_file',
                 falff, 'in_file')
    alff.connect(NVOLS, ('nvols', set_op1_str),
                 falff, 'op_string')
    alff.connect(_sum, 'out_file',
                 falff1, 'in_file')
    alff.connect(falff, 'out_file',
                 falff1, 'operand_files')

    alff.connect(_sum, 'out_file',
                 normM, 'in_file')
    alff.connect(inputNode, 'rest_mask',
                 normM, 'mask_file')
    alff.connect(_sum, 'out_file',
                 normS, 'in_file')
    alff.connect(inputNode, 'rest_mask',
                 normS, 'mask_file')
    alff.connect(falff1, 'out_file',
                 normM1, 'in_file')
    alff.connect(inputNode, 'rest_mask',
                 normM1, 'mask_file')
    alff.connect(falff1, 'out_file',
                 normS1, 'in_file')
    alff.connect(inputNode, 'rest_mask',
                 normS1, 'mask_file')

    alff.connect(normM, 'out_stat',
                 op_string, 'mean')
    alff.connect(normS, 'out_stat',
                 op_string, 'std_dev')
    alff.connect(op_string, 'op_string',
                 Z_alff, 'op_string')
    alff.connect(_sum, 'out_file',
                 Z_alff, 'in_file')
    alff.connect(inputNode, 'rest_mask',
                 Z_alff, 'operand_files')

    alff.connect(normM1, 'out_stat',
                 op_string1, 'mean')
    alff.connect(normS1, 'out_stat',
                 op_string1, 'std_dev')
    alff.connect(op_string1, 'op_string',
                 Z_falff, 'op_string')
    alff.connect(falff1, 'out_file',
                 Z_falff, 'in_file')
    alff.connect(inputNode, 'rest_mask',
                 Z_falff, 'operand_files')

    alff.connect(inputNode, 'standard',
                 warp_alff, 'ref_file')
    alff.connect(Z_alff, 'out_file',
                 warp_alff, 'in_file')
    alff.connect(inputNode, 'fieldcoeff_file',
                 warp_alff, 'field_file')
    alff.connect(inputNode, 'premat',
                 warp_alff, 'premat')

    alff.connect(inputNode, 'standard',
                 warp_falff, 'ref_file')
    alff.connect(Z_falff, 'out_file',
                 warp_falff, 'in_file')
    alff.connect(inputNode, 'fieldcoeff_file',
                 warp_falff, 'field_file')
    alff.connect(inputNode, 'premat',
                 warp_falff, 'premat')

    alff.connect(warp_alff, 'out_file',
                 smooth, 'in_file')
    alff.connect(inputnode_fwhm, ('fwhm', set_gauss),
                 smooth, 'op_string')
    alff.connect(inputNode, 'rest_mask2standard',
                 smooth, 'operand_files')

    alff.connect(warp_falff, 'out_file',
                 fsmooth, 'in_file')
    alff.connect(inputnode_fwhm, ('fwhm', set_gauss),
                 fsmooth, 'op_string')
    alff.connect(inputNode, 'rest_mask2standard',
                 fsmooth, 'operand_files')

    alff.connect(cp, 'out_file',
                 outputNode, 'processed_alff')
    alff.connect(mean, 'out_file',
                 outputNode, 'processed_mean_alff')
    alff.connect(pspec, 'out_file',
                 outputNode, 'power_spectrum_distribution')
    alff.connect(_sum, 'out_file',
                 outputNode, 'alff_img')
    alff.connect(falff1, 'out_file',
                 outputNode, 'falff_img')
    alff.connect(Z_alff, 'out_file',
                 outputNode, 'alff_Z_img')
    alff.connect(Z_falff, 'out_file',
                 outputNode, 'falff_Z_img')
    alff.connect(warp_alff, 'out_file',
                 outputNode, 'alff_Z_2standard_img')
    alff.connect(warp_falff, 'out_file',
                 outputNode, 'falff_Z_2standard_img')
    alff.connect(smooth, 'out_file',
                 outputNode, 'alff_Z_2standard_fwhm_img')
    alff.connect(fsmooth, 'out_file',
                 outputNode, 'falff_Z_2standard_fwhm_img')
    return alff

def create_timeseries_preproc(unit_time_series_extraction, voxel_time_series_extraction, vertices_time_series_extraction, run_surface_registraion):

    import nipype.interfaces.freesurfer as fs
    preproc = pe.Workflow(name='timeseries_preproc')

    inputNode = pe.Node(util.IdentityInterface(fields=['standard',
                                                       'recon_subjects',
                                                       'brain',
                                                       'reorient',
                                                       'motion_correct',
                                                       'warp_file',
                                                       'premat',
                                                       'identity_matrix',
                                                       'unitTSOutputs',
                                                       'voxelTSOutputs',
                                                       'verticesTSOutputs']),
                                                    name='inputspec')

    inputnode_getparc = pe.Node(util.IdentityInterface(fields=['parcelations']),
                             name='getparc')


    inputNode_getmask = pe.Node(util.IdentityInterface(fields=['masks']),
                              name='getmask')

    outputNode = pe.Node(util.IdentityInterface(fields=['mask_outputs',
                                                         'parc_outputs',
                                                         'surface_outputs'
                                                         ]),
                                                    name='outputspec')

    timeseries_reconall = pe.Node(interface=fs.ReconAll(), name="timeseries_reconall")
    timeseries_reconall.inputs.directive = 'all'

    timeseries_bbreg = pe.MapNode(interface=fs.BBRegister(init='fsl', contrast_type='t2', registered_file=True, out_fsl_file=True), name='timeseries_bbreg', iterfield=["source_file"] )

    timeseries_apply_warp = pe.MapNode(interface=fsl.ApplyWarp(), name='timeseries_apply_warp', iterfield=["in_file", "premat"])

    preproc.connect(inputNode, 'motion_correct', timeseries_apply_warp, 'in_file' )
    preproc.connect(inputNode, 'warp_file', timeseries_apply_warp, 'field_file')
    preproc.connect(inputNode, 'premat', timeseries_apply_warp, 'premat')
    preproc.connect(inputNode, 'standard', timeseries_apply_warp, 'ref_file' )


    timeseries_sampler_lh = pe.MapNode(interface=fs.SampleToSurface(hemi="lh"), name='timeseries_sampler_lh', iterfield=["source_file", "reg_file"])
    timeseries_sampler_lh.inputs.no_reshape = True
    timeseries_sampler_lh.inputs.interp_method = 'trilinear'
    timeseries_sampler_lh.inputs.sampling_method = "point"
    timeseries_sampler_lh.inputs.sampling_range = 0.5
    timeseries_sampler_lh.inputs.sampling_units = "frac"

    timeseries_sampler_rh = pe.MapNode(interface=fs.SampleToSurface(hemi="rh"), name='timeseries_sampler_rh', iterfield=["source_file", "reg_file"])
    timeseries_sampler_rh.inputs.no_reshape = True
    timeseries_sampler_rh.inputs.interp_method = 'trilinear'
    timeseries_sampler_rh.inputs.sampling_method = "point"
    timeseries_sampler_rh.inputs.sampling_range = 0.5
    timeseries_sampler_rh.inputs.sampling_units = "frac"

    timeseries_flirt = pe.MapNode(interface=fsl.FLIRT(), name='timeseries_flirt', iterfield=["in_file"])
    timeseries_flirt.inputs.interp = 'sinc'
    timeseries_flirt.inputs.apply_xfm = True

    timeseries_flirt1 = timeseries_flirt.clone('timeseries_flirt1')

    timeseries_gen_parc = pe.MapNode(util.Function(input_names=['data_file', 'template', 'unitTSOutputs'],
                                                  output_names=['out_file'],
                                                  function=gen_csv_for_parcelation),
                                                  name='timeseries_gen_parc',
                                                  iterfield=["data_file"])

    timeseries_gen_mask = pe.MapNode(util.Function(input_names=['data_file', 'template', 'voxelTSOutputs'],
                                                 output_names=['out_file'],
                                                 function=gen_csv_for_mask),
                                                 name='timeseries_gen_mask',
                                                 iterfield=["data_file"])

    timeseries_gen_surface = pe.MapNode(util.Function(input_names=['rh_surface_file', 'lh_surface_file', 'verticesTSOutputs'],
                                                    output_names=['out_file'],
                                                    function=gen_csv_for_surface),
                                                    name='timeseries_gen_surface',
                                                    iterfield=["rh_surface_file", "lh_surface_file"])


    """
        Surface Registration 
    """
    if run_surface_registraion:

        preproc.connect(inputNode, 'brain', timeseries_reconall, 'T1_files')
        preproc.connect(inputNode, ('reorient', extract_subjectID), timeseries_reconall, 'subject_id')
        preproc.connect(inputNode, 'recon_subjects', timeseries_reconall, 'subjects_dir')

        preproc.connect(timeseries_apply_warp, 'out_file', timeseries_bbreg, 'source_file' )
        preproc.connect(timeseries_reconall, 'subjects_dir', timeseries_bbreg, 'subjects_dir' )
        preproc.connect(timeseries_reconall, 'subject_id', timeseries_bbreg, 'subject_id' )

        preproc.connect(timeseries_bbreg, 'out_reg_file', timeseries_sampler_lh, 'reg_file' )
        preproc.connect(timeseries_apply_warp, 'out_file', timeseries_sampler_lh, 'source_file')

        preproc.connect(timeseries_bbreg, 'out_reg_file', timeseries_sampler_rh, 'reg_file' )
        preproc.connect(timeseries_apply_warp, 'out_file', timeseries_sampler_rh, 'source_file')

    """
        Time Series Extraction
    """

    if unit_time_series_extraction:

        preproc.connect(timeseries_apply_warp, 'out_file', timeseries_flirt, 'in_file')
        preproc.connect(inputNode, 'identity_matrix', timeseries_flirt, 'in_matrix_file')
        preproc.connect(inputnode_getparc, 'parcelations', timeseries_flirt, 'reference')

        preproc.connect(timeseries_flirt, 'out_file', timeseries_gen_parc, 'data_file')
        preproc.connect(inputNode, 'unitTSOutputs', timeseries_gen_parc, 'unitTSOutputs'  )
        preproc.connect(inputnode_getparc, 'parcelations', timeseries_gen_parc, 'template')

        preproc.connect(timeseries_gen_parc, 'out_file', outputNode, 'parc_outputs')

    if voxel_time_series_extraction:

        preproc.connect(timeseries_apply_warp, 'out_file', timeseries_flirt1, 'in_file')
        preproc.connect(inputNode, 'identity_matrix', timeseries_flirt1, 'in_matrix_file')
        preproc.connect(inputNode_getmask, 'masks', timeseries_flirt1, 'reference')

        preproc.connect(timeseries_flirt1, 'out_file', timeseries_gen_mask, 'data_file')
        preproc.connect(inputNode, 'voxelTSOutputs', timeseries_gen_mask, 'voxelTSOutputs')
        preproc.connect(inputNode_getmask, 'masks', timeseries_gen_mask, 'template')

        preproc.connect(timeseries_gen_mask, 'out_file', outputNode, 'mask_outputs')

    if vertices_time_series_extraction and run_surface_registraion:

        preproc.connect(timeseries_sampler_rh, 'out_file', timeseries_gen_surface, 'rh_surface_file' )
        preproc.connect(timeseries_sampler_lh, 'out_file', timeseries_gen_surface, 'lh_surface_file')
        preproc.connect(inputNode, 'verticesTSOutputs', timeseries_gen_surface, 'verticesTSOutputs' )

        preproc.connect(timeseries_gen_surface, 'out_file', outputNode, 'surface_outputs')

    return preproc
