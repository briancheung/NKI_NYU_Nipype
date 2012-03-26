#!/frodo/shared/epd/bin/python
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

def create_anat_preproc():

    preproc = pe.Workflow(name='anatpreproc')

    inputNode = pe.Node(util.IdentityInterface(fields=['anat']),
                                                    name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['refit',
                                                            'reorient',
                                                            'skullstrip',
                                                            'brain']),
                                                    name='outputspec')

    anat_refit = pe.MapNode(interface=afni.Refit(), name='anat_refit', iterfield=['in_file'])
    anat_refit.inputs.deoblique = True

    anat_reorient = pe.MapNode(interface=afni.Resample(), name='anat_reorient', iterfield=['in_file'])
    anat_reorient.inputs.orientation = 'RPI'

    anat_skullstrip = pe.MapNode(interface=afni.SkullStrip(), name='anat_skullstrip', iterfield=['in_file'])
    anat_skullstrip.inputs.options = '-o_ply'

    anat_calc = pe.MapNode(interface=afni.Calc(), name='anat_calc', iterfield=['infile_a', 'infile_b'])
    anat_calc.inputs.expr = '\'a*step(b)\''

    preproc.connect(inputNode, 'anat', anat_refit, 'in_file')
    preproc.connect(anat_refit, 'out_file', anat_reorient, 'in_file')
    preproc.connect(anat_reorient, 'out_file', anat_skullstrip, 'in_file')
    preproc.connect(anat_skullstrip, 'out_file', anat_calc, 'infile_b')
    preproc.connect(anat_reorient, 'out_file', anat_calc, 'infile_a')

    preproc.connect(anat_refit, 'out_file', outputNode, 'refit')
    preproc.connect(anat_reorient, 'out_file', outputNode, 'reorient')
    preproc.connect(anat_skullstrip, 'out_file', outputNode, 'skullstrip')
    preproc.connect(anat_calc, 'out_file', outputNode, 'brain')

    return preproc


def create_func_preproc():

    preproc = pe.Workflow(name='funcpreproc')
    inputNode = pe.Node(util.IdentityInterface(fields=['rest', 'start_idx', 'stop_idx']),
                                                name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['drop_tr',
                                                            'refit',
                                                            'reorient',
                                                            'reorient_mean',
                                                            'motion_correct',
                                                            'movement_parameters',
                                                            'max_displacement',
                                                            'mask',
                                                            'skullstrip',
                                                            'example_func',
                                                            'preprocessed',
                                                            'preprocessed_mask']),

                            name='outputspec')

    func_calc = pe.MapNode(interface=afni.Calc(), name='func_calc',
                                                    iterfield=['infile_a'])
    func_calc.inputs.expr = '\'a\''

    func_refit = pe.MapNode(interface=afni.Refit(), name='func_refit',
                                                    iterfield=['in_file'])
    func_refit.inputs.deoblique = True

    func_reorient = pe.MapNode(interface=afni.Resample(), name='func_reorient',
                                                        iterfield=['in_file'])

    func_reorient.inputs.orientation = 'RPI'

    func_tstat = pe.MapNode(interface=afni.TStat(), name='func_tstat',
                                                    iterfield=['in_file'])
    func_tstat.inputs.options = '-mean'

    func_tstat_1 = func_tstat.clone('func_tstat_1')

    func_volreg = pe.MapNode(interface=e_afni.Threedvolreg(), name='func_volreg',
                                                    iterfield=['in_file', 'basefile'])

    func_volreg.inputs.other = '-Fourier -twopass'
    func_volreg.inputs.zpad = '4'

    func_volreg_1 = func_volreg.clone('func_volreg_1')
    func_automask = pe.MapNode(interface=afni.Automask(), name='func_automask',
                                                        iterfield=['in_file'])
    func_automask.inputs.dilate = 1

    func_calcR = pe.MapNode(interface=afni.Calc(), name='func_calcR',
                                                    iterfield=['infile_a', 'infile_b'])
    func_calcR.inputs.expr = '\'a*b\''

    func_mean = pe.MapNode(interface=afni.TStat(), name='func_mean', iterfield=['in_file'])
    func_mean.inputs.options = '-mean'

#    func_calcI = pe.MapNode(interface=afni.Calc(), name='func_calcI', iterfield=['infile_a'])
#    func_calcI.inputs.single_idx = 4
#    func_calcI.inputs.expr = '\'a\''

    #func_despike = pe.MapNode(interface=afni.Despike(), name='func_despike',
     #                                                   iterfield=['in_file'])

#    func_smooth = pe.MapNode(interface=fsl.MultiImageMaths(), name='func_smooth',
#                                       iterfield=['in_file', 'operand_files'])
    #Note: Ask Satra about setting about building op_string iterfield
    #func_str1 = '-kernel gauss %f -fmean -mas' %(op_string)
    #func_smooth.inputs.op_string = func_str1 + ' %s'

    func_scale = pe.MapNode(interface=fsl.ImageMaths(), name='func_scale', iterfield=['in_file'])
    func_scale.inputs.op_string = '-ing 10000'
    func_scale.inputs.out_data_type = 'float'

    func_mask = pe.MapNode(interface=fsl.ImageMaths(), name='func_mask', iterfield=['in_file'])
    func_mask.inputs.op_string = '-Tmin -bin'
    func_mask.inputs.out_data_type = 'char'

    preproc.connect(inputNode, 'rest', func_calc, 'infile_a')
    preproc.connect(inputNode, 'start_idx', func_calc, 'start_idx')
    preproc.connect(inputNode, 'stop_idx', func_calc, 'stop_idx')
    preproc.connect(func_calc, 'out_file', func_refit, 'in_file')
    preproc.connect(func_refit, 'out_file', func_reorient, 'in_file')
    preproc.connect(func_reorient, 'out_file', func_tstat, 'in_file')
    preproc.connect(func_reorient, 'out_file', func_volreg, 'in_file')
    preproc.connect(func_tstat, 'out_file', func_volreg, 'basefile')
    preproc.connect(func_volreg, 'out_file', func_tstat_1, 'in_file')
    preproc.connect(func_reorient, 'out_file', func_volreg_1, 'in_file')
    preproc.connect(func_tstat_1, 'out_file', func_volreg_1, 'basefile')
    preproc.connect(func_volreg_1, 'out_file', func_automask, 'in_file')
    preproc.connect(func_volreg_1, 'out_file', func_calcR, 'infile_a')
    preproc.connect(func_automask, 'out_file', func_calcR, 'infile_b')
    preproc.connect(func_calcR, 'out_file', func_mean, 'in_file')
    preproc.connect(func_calcR, 'out_file', func_scale, 'in_file')
    preproc.connect(func_scale, 'out_file', func_mask, 'in_file')

    preproc.connect(func_calc, 'out_file', outputNode, 'drop_tr')
    preproc.connect(func_refit, 'out_file', outputNode, 'refit')
    preproc.connect(func_reorient, 'out_file', outputNode, 'reorient')
    preproc.connect(func_tstat_1, 'out_file', outputNode, 'reorient_mean')
    preproc.connect(func_volreg_1, 'out_file', outputNode, 'motion_correct')
    preproc.connect(func_volreg_1, 'md1d_file', outputNode, 'max_displacement')
    preproc.connect(func_volreg_1, 'oned_file', outputNode, 'movement_parameters')
    preproc.connect(func_automask, 'out_file', outputNode, 'mask')
    preproc.connect(func_calcR, 'out_file', outputNode, 'skullstrip')
    preproc.connect(func_mean, 'out_file', outputNode, 'example_func')
    preproc.connect(func_scale, 'out_file', outputNode, 'preprocessed')
    preproc.connect(func_mask, 'out_file', outputNode, 'preprocessed_mask')

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

    reg_flirt = pe.MapNode(interface=fsl.FLIRT(), name='reg_flirt', iterfield=['in_file', 'reference'])
    reg_flirt.inputs.cost = 'corratio'
    reg_flirt.inputs.dof = 6
    reg_flirt.inputs.interp = 'nearestneighbour'

    # Create mat file for conversion from subject's anatomical to functional
    reg_xfm1 = pe.MapNode(interface=fsl.ConvertXFM(), name='reg_xfm1', iterfield=['in_file'])
    reg_xfm1.inputs.invert_xfm = True

    ## T1->STANDARD
    ## NOTE THAT THIS IS Linear registration, you may want to use FNIRT (non-linear)
    reg_flirt1 = pe.MapNode(interface=fsl.FLIRT(), name='reg_flirt1', iterfield=['in_file'])
    reg_flirt1.inputs.cost = 'corratio'
    reg_flirt1.inputs.cost_func = 'corratio'
    reg_flirt1.inputs.dof = 12
    reg_flirt1.inputs.interp = 'nearestneighbour'

    ## Create mat file for conversion from standard to high res
    reg_xfm2 = pe.MapNode(interface=fsl.ConvertXFM(), name='reg_xfm2', iterfield=['in_file'])
    reg_xfm2.inputs.invert_xfm = True

    reg_inw = pe.MapNode(interface=e_afni.InvWarp(), name='reg_inw', iterfield=['in_file', 'ref_file'])

    ## T1->STANDARD NONLINEAR
    # Perform nonlinear registration (higres to standard)
    reg_fnt = pe.MapNode(interface=fsl.FNIRT(), name='reg_fnt', iterfield=['in_file',
                                                                        'affine_file'])
    reg_fnt.inputs.fieldcoeff_file = True
    reg_fnt.inputs.jacobian_file = True
    reg_fnt.inputs.warp_resolution = (10, 10, 10)

    ## Apply nonlinear registration (func to standard)
    reg_warp = pe.MapNode(interface=fsl.ApplyWarp(), name='reg_warp',
                                                    iterfield=['in_file',
                                                                'premat',
                                                                'field_file'])

    preproc.connect(inputNode, 'example_func', reg_flirt, 'in_file')
    preproc.connect(inputNode, 'brain', reg_flirt, 'reference')
    preproc.connect(reg_flirt, 'out_matrix_file', reg_xfm1, 'in_file')
    preproc.connect(inputNode, 'brain', reg_flirt1, 'in_file')
    preproc.connect(inputNode, 'standard_res_brain', reg_flirt1, 'reference')
    preproc.connect(reg_flirt1, 'out_matrix_file', reg_xfm2, 'in_file')
    preproc.connect(inputNode, 'reorient', reg_fnt, 'in_file')
    preproc.connect(reg_flirt1, 'out_matrix_file', reg_fnt, 'affine_file')
    preproc.connect(inputNode, 'standard', reg_fnt, 'ref_file')
    preproc.connect(inputNode, 'standard_brain_mask_dil', reg_fnt, 'refmask_file')
    preproc.connect(inputNode, 'config_file', reg_fnt, 'config_file')
    preproc.connect(reg_fnt, 'fieldcoeff_file', reg_inw, 'in_file')
    preproc.connect(inputNode, 'brain', reg_inw, 'ref_file')
    preproc.connect(inputNode, 'example_func', reg_warp, 'in_file')
    preproc.connect(inputNode, 'standard', reg_warp, 'ref_file')
    preproc.connect(reg_fnt, 'fieldcoeff_file', reg_warp, 'field_file')
    preproc.connect(reg_flirt, 'out_matrix_file', reg_warp, 'premat')

    preproc.connect(reg_flirt, 'out_matrix_file', outputNode, 'example_func2highres')
    preproc.connect(reg_flirt, 'out_file', outputNode, 'example_func2highres_mat')
    preproc.connect(reg_xfm1, 'out_file', outputNode, 'highres2example_func_mat')
    preproc.connect(reg_warp, 'out_file', outputNode, 'example_func2standard_NL')
    preproc.connect(reg_flirt1, 'out_file', outputNode, 'highres2standard')
    preproc.connect(reg_flirt1, 'out_matrix_file', outputNode, 'highres2standard_mat')
    preproc.connect(reg_inw, 'out_file', outputNode, 'stand2highres_warp')
    preproc.connect(reg_fnt, 'jacobian_file', outputNode, 'highres2standard_jac')
    preproc.connect(reg_fnt, 'fieldcoeff_file', outputNode, 'highres2standard_warp')
    preproc.connect(reg_fnt, 'warped_file', outputNode, 'highres2standard_NL')

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
                                                    'PRIOR_WHITE']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['csf_t12func',
                                                            'csf_mni2func',
                                                            'csf_combo',
                                                            'csf_bin',
                                                            'csf_mask',
                                                            'global_mask',
                                                            'wm_t12func',
                                                            'wm_mni2func',
                                                            'wm_combo',
                                                            'wm_bin',
                                                            'probability_maps',
                                                            'wm_mask']),
                        name='outputspec')

    seg_segment = pe.MapNode(interface=fsl.FAST(), name='seg_segment', iterfield=['in_files'])
    seg_segment.inputs.img_type = 1
    seg_segment.inputs.segments = True
    seg_segment.inputs.probability_maps = True
    seg_segment.inputs.out_basename = 'segment'

    seg_copy = pe.MapNode(interface=afni.Copy(), name='seg_copy', iterfield=['in_file'])

    seg_flirt = pe.MapNode(interface=fsl.FLIRT(), name='seg_flirt', iterfield=['reference', 'in_matrix_file'])
    seg_flirt.inputs.apply_xfm = True

    seg_warp = pe.MapNode(interface=fsl.ApplyWarp(), name='seg_warp', iterfield=['ref_file', 'postmat', 'field_file'])
    seg_warp.inputs.interp = 'nn'

    seg_warp1 = seg_warp.clone('seg_warp1')

    seg_smooth1 = pe.MapNode(interface=fsl.MultiImageMaths(), name='seg_smooth1', iterfield=['in_file'])
    seg_str1 = '-mas %s '
    seg_smooth1.inputs.op_string = seg_str1

    seg_thresh = pe.MapNode(interface=fsl.ImageMaths(), name='seg_thresh', iterfield=['in_file'])
    seg_str1 = '-thr 0.4 -bin '
    seg_thresh.inputs.op_string = seg_str1

    seg_mask = pe.MapNode(interface=fsl.MultiImageMaths(), name='seg_mask', iterfield=['in_file', 'operand_files'])
    seg_str1 = '-mas %s '
    seg_mask.inputs.op_string = seg_str1

    seg_prior1 = pe.MapNode(interface=fsl.MultiImageMaths(), name='seg_prior1', iterfield=['in_file'])
    seg_str1 = '-mas %s '
    seg_prior1.inputs.op_string = seg_str1

    seg_thresh1 = pe.MapNode(interface=fsl.ImageMaths(), name='seg_thresh1', iterfield=['in_file'])
    seg_str1 = '-thr 0.66 -bin '
    seg_thresh1.inputs.op_string = seg_str1

    seg_flirt3 = pe.MapNode(interface=fsl.FLIRT(), name='seg_flirt3', iterfield=['reference', 'in_matrix_file'])
    seg_flirt3.inputs.apply_xfm = True

    seg_mask1 = seg_mask.clone('seg_mask1')

    preproc.connect(inputNode, 'brain', seg_segment, 'in_files')
    preproc.connect(seg_segment, ('probability_maps', pick_wm_0), seg_flirt, 'in_file')
    preproc.connect(inputNode, 'example_func', seg_flirt, 'reference')
    preproc.connect(inputNode, 'highres2example_func_mat', seg_flirt, 'in_matrix_file')
    preproc.connect(inputNode, 'example_func', seg_warp, 'ref_file')
    preproc.connect(inputNode, 'stand2highres_warp', seg_warp, 'field_file')
    preproc.connect(inputNode, 'PRIOR_CSF', seg_warp, 'in_file')
    preproc.connect(inputNode, 'highres2example_func_mat', seg_warp, 'postmat')
    preproc.connect(seg_flirt, 'out_file', seg_smooth1, 'in_file')
    preproc.connect(seg_warp, 'out_file', seg_smooth1, 'operand_files')
    preproc.connect(seg_smooth1, 'out_file', seg_thresh, 'in_file')
    preproc.connect(seg_thresh, 'out_file', seg_mask, 'in_file')
    preproc.connect(inputNode, 'preprocessed_mask', seg_copy, 'in_file')
    preproc.connect(seg_copy, 'out_file', seg_mask, 'operand_files')
    preproc.connect(seg_segment, ('probability_maps', pick_wm_1), seg_flirt3, 'in_file')
    preproc.connect(inputNode, 'example_func', seg_flirt3, 'reference')
    preproc.connect(inputNode, 'highres2example_func_mat', seg_flirt3, 'in_matrix_file')
    preproc.connect(inputNode, 'example_func', seg_warp1, 'ref_file')
    preproc.connect(inputNode, 'stand2highres_warp', seg_warp1, 'field_file')
    preproc.connect(inputNode, 'PRIOR_WHITE', seg_warp1, 'in_file')
    preproc.connect(inputNode, 'highres2example_func_mat', seg_warp1, 'postmat')
    preproc.connect(seg_flirt3, 'out_file', seg_prior1, 'in_file')
    preproc.connect(seg_warp1, 'out_file', seg_prior1, 'operand_files')
    preproc.connect(seg_prior1, 'out_file', seg_thresh1, 'in_file')
    preproc.connect(seg_thresh1, 'out_file', seg_mask1, 'in_file')
    preproc.connect(seg_copy, 'out_file', seg_mask1, 'operand_files')

    preproc.connect(seg_segment, 'probability_maps', outputNode, 'probability_maps')
    preproc.connect(seg_flirt, 'out_file', outputNode, 'csf_t12func')
    preproc.connect(seg_warp, 'out_file', outputNode, 'csf_mni2func')
    preproc.connect(seg_smooth1, 'out_file', outputNode, 'csf_combo')
    preproc.connect(seg_thresh, 'out_file', outputNode, 'csf_bin')
    preproc.connect(seg_mask, 'out_file', outputNode, 'csf_mask')
    preproc.connect(seg_copy, 'out_file', outputNode, 'global_mask')
    preproc.connect(seg_flirt3, 'out_file', outputNode, 'wm_t12func')
    preproc.connect(seg_warp1, 'out_file', outputNode, 'wm_mni2func')
    preproc.connect(seg_prior1, 'out_file', outputNode, 'wm_combo')
    preproc.connect(seg_thresh1, 'out_file', outputNode, 'wm_bin')
    preproc.connect(seg_mask1, 'out_file', outputNode, 'wm_mask')

    return preproc
