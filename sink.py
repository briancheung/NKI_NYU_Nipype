import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
import nipype.pipeline.engine as pe
from utils import formatpath



def anat_sink(workflow, datasink, mprage_mni):        
        
        rename_anat= rename_anat_preproc()
        rename_anat.inputs.filenames.anat_calc= "mprage_brian.nii.gz"
        rename_anat.inputs.filenames.anat_reorient= "mprage_RPI.nii.gz"
        
        workflow.connect(mprage_mni, 'outputspec.brain_mni', rename_anat, 'inputspec.anat_calc' )
        workflow.connect(mprage_mni, 'outputspec.reorient_mni', rename_anat, 'inputspec.anat_reorient')
        
        workflow.connect(rename_anat, 'outputspec.anat_calc_out', datasink, 'anat')
        workflow.connect(rename_anat, 'outputspec.anat_reorient_out', datasink,'anat.@0')
        
        
def func_sink(workflow,datasink, funcpreproc, func_in_mni):
        
        rename= rename_func_preproc()
        
        rename.inputs.filenames.motion_correct = "rest_mc.nii.gz" 
        rename.inputs.filenames.motion_correct_ref = "rest_mc_init_mean.nii.gz"
        rename.inputs.filenames.movement_parameters = "rest_mc.1D"
        rename.inputs.filenames.max_displacement = "rest_maxdisp.1D"
        rename.inputs.filenames.preprocessed_mask = "rest_pp_mask.nii.gz"
        
        workflow.connect(funcpreproc, 'outputspec.motion_correct', rename, 'inputspec.motion_correct')
        workflow.connect(funcpreproc, 'outputspec.motion_correct_ref', rename, 'inputspec.motion_correct_ref')
        workflow.connect(funcpreproc, 'outputspec.movement_parameters', rename, 'inputspec.movement_parameters')
        workflow.connect(funcpreproc, 'outputspec.max_displacement', rename, 'inputspec.max_displacement')
        workflow.connect(func_in_mni, 'outputspec.preprocessed_mask_mni', rename, 'inputspec.preprocessed_mask')
         
         
        workflow.connect(rename, 'outputspec.motion_correct_out', datasink, 'func.@0')
        workflow.connect(rename, 'outputspec.motion_correct_ref_out', datasink, 'func.@1')
        workflow.connect(rename, 'outputspec.movement_parameters_out', datasink, 'func.@2')
        workflow.connect(rename, 'outputspec.max_displacement_out', datasink, 'func.@3')
        workflow.connect(rename, 'outputspec.preprocessed_mask_out', datasink, 'func.@4')


def reg_sink(workflow, datasink, regpreproc):
        
        rename_reg= rename_reg_preproc()
        
        rename_reg.inputs.filenames.highres2standard= "highres2standard.nii.gz" 
        rename_reg.inputs.filenames.example_func2highres_mat="example_func2highres.mat"
        rename_reg.inputs.filenames.example_func2highres="example_func2highres.nii.gz"
        rename_reg.inputs.filenames.example_func2standard_NL='example_func2standard_NL.nii.gz'
        rename_reg.inputs.filenames.highres2example_func_mat="highres2example_func.mat"
        
        workflow.connect(regpreproc, 'outputspec.highres2standard', rename_reg, 'inputspec.highres2standard')
        workflow.connect(regpreproc, 'outputspec.example_func2highres_mat', rename_reg, 'inputspec.example_func2highres_mat')
        workflow.connect(regpreproc, 'outputspec.example_func2highres', rename_reg, 'inputspec.example_func2highres')
        workflow.connect(regpreproc, 'outputspec.example_func2standard_NL', rename_reg, 'inputspec.example_func2standard_NL')
        workflow.connect(regpreproc, 'outputspec.highres2example_func_mat', rename_reg, 'inputspec.highres2example_func_mat')
        
        workflow.connect(rename_reg, 'outputspec.highres2standard_out', datasink, 'reg.@0')
        workflow.connect(rename_reg, 'outputspec.example_func2highres_mat_out', datasink, 'reg.@1')
        workflow.connect(rename_reg, 'outputspec.example_func2highres_out', datasink, 'reg.@2')
        workflow.connect(rename_reg, 'outputspec.example_func2standard_NL_out', datasink, 'reg.@3')
        workflow.connect(rename_reg, 'outputspec.highres2example_func_mat_out', datasink, 'reg.@4')


def seg_sink(workflow, datasink, segpreproc, mprage_mni):
           
        rename_seg=rename_seg_preproc()
        
        rename_seg.inputs.filenames.csf_t12func = "csf_t12func.nii.gz" 
        rename_seg.inputs.filenames.csf_mni2func = "csf_mni2func.nii.gz"
        rename_seg.inputs.filenames.csf_combo = "csf_combo.nii.gz"
        rename_seg.inputs.filenames.csf_bin ="csf_bin.nii.gz"
        rename_seg.inputs.filenames.csf_mask ="csf_mask.nii.gz"
        rename_seg.inputs.filenames.gm_t12func ="gm_t12func.nii.gz"
        rename_seg.inputs.filenames.gm_mni2func ="gm_mni2func.nii.gz"
        rename_seg.inputs.filenames.gm_combo = "gm_combo.nii.gz"
        rename_seg.inputs.filenames.gm_bin = "gm_bin.nii.gz"
        rename_seg.inputs.filenames.gm_mask = "gm_mask.nii.gz"
        rename_seg.inputs.filenames.global_mask = "global_mask.nii.gz"
        rename_seg.inputs.filenames.wm_t12func = "wm_t12func.nii.gz"
        rename_seg.inputs.filenames.wm_mni2func = "wm_mni2func.nii.gz"
        rename_seg.inputs.filenames.wm_combo = "wm_combo.nii.gz"
        rename_seg.inputs.filenames.wm_bin = "wm_bin.nii.gz"
        rename_seg.inputs.filenames.wm_mask = "wm_mask.nii.gz"
        
        workflow.connect(segpreproc, 'outputspec.csf_t12func', rename_seg, 'inputspec.csf_t12func')
        workflow.connect(segpreproc, 'outputspec.csf_mni2func',  rename_seg, 'inputspec.csf_mni2func')
        workflow.connect(segpreproc, 'outputspec.csf_combo',  rename_seg, 'inputspec.csf_combo')
        workflow.connect(segpreproc, 'outputspec.csf_bin',  rename_seg, 'inputspec.csf_bin')
        workflow.connect(segpreproc, 'outputspec.csf_mask',  rename_seg, 'inputspec.csf_mask')
        workflow.connect(segpreproc, 'outputspec.gm_t12func',  rename_seg, 'inputspec.gm_t12func')
        workflow.connect(segpreproc, 'outputspec.gm_mni2func',  rename_seg, 'inputspec.gm_mni2func')
        workflow.connect(segpreproc, 'outputspec.gm_combo',  rename_seg, 'inputspec.gm_combo')
        workflow.connect(segpreproc, 'outputspec.gm_bin',  rename_seg, 'inputspec.gm_bin')
        workflow.connect(segpreproc, 'outputspec.gm_mask',  rename_seg, 'inputspec.gm_mask')
        workflow.connect(segpreproc, 'outputspec.global_mask',  rename_seg, 'inputspec.global_mask')
        workflow.connect(segpreproc, 'outputspec.wm_t12func',  rename_seg, 'inputspec.wm_t12func')
        workflow.connect(segpreproc, 'outputspec.wm_mni2func',  rename_seg, 'inputspec.wm_mni2func')
        workflow.connect(segpreproc, 'outputspec.wm_combo',  rename_seg, 'inputspec.wm_combo')
        workflow.connect(segpreproc, 'outputspec.wm_bin',  rename_seg, 'inputspec.wm_bin')
        workflow.connect(segpreproc, 'outputspec.wm_mask',  rename_seg, 'inputspec.wm_mask')
        
        workflow.connect(rename_seg, 'outputspec.csf_t12func_out',  datasink, 'segment.@0')
        workflow.connect(rename_seg, 'outputspec.csf_mni2func_out',  datasink, 'segment.@1' )
        workflow.connect(rename_seg, 'outputspec.csf_combo_out',  datasink, 'segment.@2' )
        workflow.connect(rename_seg, 'outputspec.csf_bin_out',  datasink, 'segment.@3')
        workflow.connect(rename_seg, 'outputspec.csf_mask_out',  datasink, 'segment.@4')
        workflow.connect(rename_seg, 'outputspec.gm_t12func_out',  datasink, 'segment.@5')
        workflow.connect(rename_seg, 'outputspec.gm_mni2func_out',  datasink, 'segment.@6')
        workflow.connect(rename_seg, 'outputspec.gm_combo_out',  datasink, 'segment.@7')
        workflow.connect(rename_seg, 'outputspec.gm_bin_out',  datasink, 'segment.@8')
        workflow.connect(rename_seg, 'outputspec.gm_mask_out',  datasink, 'segment.@9')
        workflow.connect(rename_seg, 'outputspec.global_mask_out',  datasink, 'segment.@10')
        workflow.connect(rename_seg, 'outputspec.wm_t12func_out' ,  datasink, 'segment.@11')
        workflow.connect(rename_seg, 'outputspec.wm_mni2func_out',  datasink, 'segment.@12')
        workflow.connect(rename_seg, 'outputspec.wm_combo_out',  datasink, 'segment.@13')
        workflow.connect(rename_seg, 'outputspec.wm_bin_out',  datasink, 'segment.@14')
        workflow.connect(rename_seg, 'outputspec.wm_mask_out',  datasink, 'segment.@15')
        
        workflow.connect(segpreproc, 'outputspec.probability_maps', datasink, 'segment.@16')
        workflow.connect(segpreproc, 'outputspec.mixeltype', datasink, 'segment.@17')
        workflow.connect(segpreproc, 'outputspec.partial_volume_map', datasink, 'segment.@18')
        workflow.connect(segpreproc, 'outputspec.partial_volume_files', datasink, 'segment.@19')
        
        workflow.connect(mprage_mni, 'outputspec.mixeltype_mni', datasink, 'segment.@20')
        workflow.connect(mprage_mni, 'outputspec.probability_maps_mni', datasink, 'segment.@21')
        workflow.connect(mprage_mni, 'outputspec.partial_volume_map_mni',datasink, 'segment.@22')
        workflow.connect(mprage_mni, 'outputspec.partial_volume_files_mni', datasink, 'segment.@23')




def rename_connections(workflow, datasink, rename_list, sink_node):
    
    ncount = 0
    for rename in rename_list:
        din_file = sink_node + '.@' + str(ncount)
        if len(rename) == 4:
            rename_node = pe.MapNode(interface = util.Rename(), name = rename[2], 
                                     iterfield=['in_file', 'format_string'])
            rename_node.inputs.format_string = rename[3]
            
            workflow.connect(rename[0], rename[1], rename_node, 'in_file')
            workflow.connect(rename_node, 'out_file', datasink, din_file)
        else:
            workflow.connect(rename[0], rename[1], datasink, din_file)
        ncount+=1
        

def nuisance_sink(workflow,datasink, nuisancepreproc):
    rename_list = [(nuisancepreproc, 'outputspec.residual_file', 'rename', 'rest_residual.nii.gz'),
                   (nuisancepreproc, 'outputspec.median_angle_corrected_file')]
    
    rename_connections(workflow, datasink, rename_list, 'nuisance')


def scrubbing_sink(workflow, datasink, scpreproc):
    

    rename= rename_scrubbing_preproc()
        
    rename.inputs.filenames.mean_deriv_sq_1D = "mean_deriv_sq.1D" 
    rename.inputs.filenames.mean_raw_sq_1D = "mean_raw_sq.1D"
    rename.inputs.filenames.scrubbed_preprocessed = "rest_pp_scrubbed.nii.gz"
    
    workflow.connect(scpreproc, 'outputspec.mean_deriv_sq_1D', rename, 'inputspec.mean_deriv_sq_1D')
    workflow.connect(scpreproc, 'outputspec.mean_raw_sq_1D', rename, 'inputspec.mean_raw_sq_1D')
    workflow.connect(scpreproc, 'outputspec.scrubbed_preprocessed', rename, 'inputspec.scrubbed_preprocessed')
    
    workflow.connect(rename, 'outputspec.mean_deriv_sq_1D_out', datasink, 'scrubbing.@0')
    workflow.connect(rename, 'outputspec.mean_raw_sq_1D_out', datasink, 'scrubbing.@1')
    workflow.connect(rename, 'outputspec.scrubbed_preprocessed_out', datasink, 'scrubbing.@2')
    
    workflow.connect(scpreproc, 'outputspec.temp_deriv_brik_file', datasink, 'scrubbing.@3')
    workflow.connect(scpreproc, 'outputspec.temp_deriv_head_file', datasink, 'scrubbing.@4')
    workflow.connect(scpreproc, 'outputspec.temp_deriv_sq_brik_file', datasink, 'scrubbing.@5')
    workflow.connect(scpreproc, 'outputspec.temp_deriv_sq_head_file', datasink, 'scrubbing.@6')
    workflow.connect(scpreproc, 'outputspec.raw_sq_brik_file', datasink, 'scrubbing.@7')
    workflow.connect(scpreproc, 'outputspec.raw_sq_head_file', datasink, 'scrubbing.@8')
    workflow.connect(scpreproc, 'outputspec.mask_brik_file', datasink, 'scrubbing.@9')
    workflow.connect(scpreproc, 'outputspec.mask_head_file', datasink, 'scrubbing.@10')
    workflow.connect(scpreproc, 'outputspec.FD_1D', datasink, 'scrubbing.@11')
    workflow.connect(scpreproc, 'outputspec.sqrt_mean_deriv_sq_1D', datasink, 'scrubbing.@12')
    workflow.connect(scpreproc, 'outputspec.sqrt_mean_raw_sq_1D', datasink, 'scrubbing.@13')
    workflow.connect(scpreproc, 'outputspec.frames_ex_1D', datasink, 'scrubbing.@14')
    workflow.connect(scpreproc, 'outputspec.frames_in_1D', datasink, 'scrubbing.@15')
    workflow.connect(scpreproc, 'outputspec.pow_params', datasink, 'scrubbing.@16')
    workflow.connect(scpreproc, 'outputspec.ftof_percent_change_1D', datasink, 'scrubbing.@17')
    workflow.connect(scpreproc, 'outputspec.scrubbed_movement_parameters', datasink, 'scrubbing.@18')


def rename_anat_preproc():
    
    preproc = pe.Workflow(name='rename_anat_preproc')
    
    inputNode = pe.Node(util.IdentityInterface(fields=['anat_calc',
                                                       'anat_reorient']),
                        name='inputspec')
        
    outputNode = pe.Node(util.IdentityInterface(fields=['anat_calc_out',
                                                         'anat_reorient_out']),
                        name='outputspec')
    
    format_node = pe.Node(interface=util.Function(input_names=['in_file', 'filename'], 
                                                  output_names= ['format_string'], 
                                                  function= formatpath ), 
                          name='format_node')
    
    format_node1=format_node.clone('format_node1')

    rename_node = pe.MapNode(interface = util.Rename(), name = 'rename', 
                             iterfield=['in_file', 'format_string'])
    
    rename_node1=rename_node.clone('rename_node1')
     
    nodeOutputNames=['anat_calc', 'anat_reorient', 'anat_skullstrip']
    inputnode_names = pe.Node(interface=util.IdentityInterface(fields= nodeOutputNames), 
                              name="filenames")
    
    preproc.connect(inputnode_names, 'anat_calc', format_node, 'filename')
    preproc.connect(inputnode_names, 'anat_reorient', format_node1, 'filename')
    
    preproc.connect(inputNode, 'anat_calc', format_node, 'in_file')
    preproc.connect(inputNode, 'anat_reorient', format_node1, 'in_file')
    
    preproc.connect(format_node, 'format_string', rename_node, 'format_string')
    preproc.connect(format_node1, 'format_string', rename_node1, 'format_string')
    
    preproc.connect(inputNode, 'anat_calc', rename_node, 'in_file')
    preproc.connect(inputNode, 'anat_reorient', rename_node1, 'in_file')
    
    preproc.connect(rename_node, 'out_file', outputNode, 'anat_calc_out' )
    preproc.connect(rename_node1, 'out_file', outputNode, 'anat_reorient_out')
    
    return preproc


def rename_reg_preproc():
    
    preproc = pe.Workflow(name='rename_reg_preproc')
    
    inputNode = pe.Node(util.IdentityInterface(fields=['highres2standard',
                                                       'example_func2highres_mat',
                                                       'example_func2highres',
                                                       'example_func2standard_NL',
                                                       'highres2example_func_mat']),
                        name='inputspec')
        
    outputNode = pe.Node(util.IdentityInterface(fields=['highres2standard_out',
                                                        'example_func2highres_mat_out',
                                                        'example_func2highres_out',
                                                        'example_func2standard_NL_out',
                                                        'highres2example_func_mat_out']),
                        name='outputspec')
        
    nodeOutputNames=['highres2standard', 'example_func2highres_mat', 'example_func2highres', 
                     'example_func2standard_NL', 'highres2example_func_mat']
    
    inputnode_names = pe.Node(interface=util.IdentityInterface(fields= nodeOutputNames), 
                              name="filenames")
    
    
    format_node = pe.Node(interface=util.Function(input_names=['in_file', 'filename'], 
                                                  output_names= ['format_string'], 
                                                  function=formatpath ), 
                          name='format_node')
    
    format_node1=format_node.clone('format_node1')
    format_node2=format_node.clone('format_node2')
    format_node3=format_node.clone('format_node3')
    format_node4=format_node.clone('format_node4')
    
    rename_node = pe.MapNode(interface = util.Rename(), name = 'rename', 
                             iterfield=['in_file', 'format_string'])
    
    rename_node1=rename_node.clone('rename_node1')
    rename_node2=rename_node.clone('rename_node2')
    rename_node3=rename_node.clone('rename_node3')
    rename_node4=rename_node.clone('rename_node4')
    
    
    preproc.connect(inputnode_names, 'highres2standard', format_node, 'filename')
    preproc.connect(inputnode_names, 'example_func2highres_mat', format_node1, 'filename')
    preproc.connect(inputnode_names, 'example_func2highres', format_node2, 'filename')
    preproc.connect(inputnode_names, 'example_func2standard_NL', format_node3, 'filename')
    preproc.connect(inputnode_names, 'highres2example_func_mat', format_node4, 'filename')
    
    preproc.connect(inputNode, 'highres2standard', format_node, 'in_file')
    preproc.connect(inputNode, 'example_func2highres_mat', format_node1, 'in_file')
    preproc.connect(inputNode, 'example_func2highres', format_node2, 'in_file')
    preproc.connect(inputNode, 'example_func2standard_NL', format_node3, 'in_file')
    preproc.connect(inputNode, 'highres2example_func_mat', format_node4, 'in_file')
    
    preproc.connect(format_node, 'format_string', rename_node, 'format_string')
    preproc.connect(format_node1, 'format_string', rename_node1, 'format_string')
    preproc.connect(format_node2, 'format_string', rename_node2, 'format_string')
    preproc.connect(format_node3, 'format_string', rename_node3, 'format_string')
    preproc.connect(format_node4, 'format_string', rename_node4, 'format_string')
    
    preproc.connect(inputNode, 'highres2standard', rename_node, 'in_file')
    preproc.connect(inputNode, 'example_func2highres_mat', rename_node1, 'in_file')
    preproc.connect(inputNode, 'example_func2highres', rename_node2, 'in_file')
    preproc.connect(inputNode, 'example_func2standard_NL', rename_node3, 'in_file')
    preproc.connect(inputNode, 'highres2example_func_mat', rename_node4, 'in_file')
    
    preproc.connect(rename_node, 'out_file', outputNode,'highres2standard_out')
    preproc.connect(rename_node1, 'out_file', outputNode, 'example_func2highres_mat_out')
    preproc.connect(rename_node2, 'out_file', outputNode, 'example_func2highres_out')
    preproc.connect(rename_node3, 'out_file', outputNode, 'example_func2standard_NL_out')
    preproc.connect(rename_node4, 'out_file', outputNode, 'highres2example_func_mat_out')
    
    return preproc
    
    
def rename_func_preproc():
    
    preproc = pe.Workflow(name='rename_func_preproc')
    
    inputNode = pe.Node(util.IdentityInterface(fields=[ 'motion_correct',
                                                        'motion_correct_ref',
                                                        'movement_parameters',
                                                        'max_displacement',
                                                        'preprocessed_mask']),
                        name='inputspec')
        
    outputNode = pe.Node(util.IdentityInterface(fields=[ 'motion_correct_out',
                                                         'motion_correct_ref_out',
                                                         'movement_parameters_out',
                                                         'max_displacement_out',
                                                         'preprocessed_mask_out']),
                        name='outputspec')
        
    nodeOutputNames=['motion_correct','motion_correct_ref','movement_parameters',
                     'max_displacement','preprocessed_mask']
    inputnode_names = pe.Node(interface=util.IdentityInterface(fields= nodeOutputNames), 
                              name="filenames")
    
    rename_node = pe.MapNode(interface = util.Rename(), name = 'rename', 
                             iterfield=['in_file', 'format_string'])
    
    rename_node1=rename_node.clone('rename_node1')
    rename_node2=rename_node.clone('rename_node2')
    rename_node3=rename_node.clone('rename_node3')
    rename_node4=rename_node.clone('rename_node4')
    
    
    format_node = pe.Node(interface=util.Function(input_names=['in_file', 'filename'], 
                                                  output_names= ['format_string'], 
                                                  function=formatpath ), 
                          name='format_node')
    
    format_node1=format_node.clone('format_node1')
    format_node2=format_node.clone('format_node2')
    format_node3=format_node.clone('format_node3')
    format_node4=format_node.clone('format_node4')
    
    preproc.connect(inputnode_names, 'motion_correct', format_node, 'filename')
    preproc.connect(inputnode_names, 'motion_correct_ref', format_node1, 'filename')
    preproc.connect(inputnode_names, 'movement_parameters', format_node2, 'filename')
    preproc.connect(inputnode_names, 'max_displacement', format_node3, 'filename')
    preproc.connect(inputnode_names, 'preprocessed_mask', format_node4, 'filename')
    
    preproc.connect(inputNode, 'motion_correct', format_node, 'in_file')
    preproc.connect(inputNode, 'motion_correct_ref', format_node1, 'in_file')
    preproc.connect(inputNode, 'movement_parameters', format_node2, 'in_file')
    preproc.connect(inputNode, 'max_displacement', format_node3, 'in_file')
    preproc.connect(inputNode, 'preprocessed_mask', format_node4, 'in_file')
    
    preproc.connect(format_node, 'format_string', rename_node, 'format_string')
    preproc.connect(format_node1, 'format_string', rename_node1, 'format_string')
    preproc.connect(format_node2, 'format_string', rename_node2, 'format_string')
    preproc.connect(format_node3, 'format_string', rename_node3, 'format_string')
    preproc.connect(format_node4, 'format_string', rename_node4, 'format_string')
    
    preproc.connect(inputNode, 'motion_correct', rename_node, 'in_file')
    preproc.connect(inputNode, 'motion_correct_ref', rename_node1, 'in_file')
    preproc.connect(inputNode, 'movement_parameters', rename_node2, 'in_file')
    preproc.connect(inputNode, 'max_displacement', rename_node3, 'in_file')
    preproc.connect(inputNode, 'preprocessed_mask', rename_node4, 'in_file')
    
    preproc.connect(rename_node, 'out_file', outputNode, 'motion_correct_out')
    preproc.connect(rename_node1, 'out_file', outputNode, 'motion_correct_ref_out')
    preproc.connect(rename_node2, 'out_file', outputNode,'movement_parameters_out')
    preproc.connect(rename_node3, 'out_file', outputNode, 'max_displacement_out')
    preproc.connect(rename_node4, 'out_file', outputNode, 'preprocessed_mask_out')
    
    return preproc


def rename_seg_preproc():
    
    preproc = pe.Workflow(name='rename_seg_preproc')
    
    inputNode = pe.Node(util.IdentityInterface(fields=['csf_t12func',
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
                                                        'wm_mask'
                                                       ]),
                        name='inputspec')
        
    outputNode = pe.Node(util.IdentityInterface(fields=['csf_t12func_out',
                                                        'csf_mni2func_out',
                                                        'csf_combo_out',
                                                        'csf_bin_out',
                                                        'csf_mask_out',
                                                        'gm_t12func_out',
                                                        'gm_mni2func_out',
                                                        'gm_combo_out',
                                                        'gm_bin_out',
                                                        'gm_mask_out',
                                                        'global_mask_out',
                                                        'wm_t12func_out',
                                                        'wm_mni2func_out',
                                                        'wm_combo_out',
                                                        'wm_bin_out',
                                                        'wm_mask_out']),
                        name='outputspec')
    
    
    nodeOutputNames=['csf_t12func','csf_mni2func', 'csf_combo', 'csf_bin',
                     'csf_mask', 'gm_t12func', 'gm_mni2func', 'gm_combo',
                     'gm_bin', 'gm_mask', 'global_mask', 'wm_t12func',
                     'wm_mni2func', 'wm_combo', 'wm_bin', 'wm_mask']
    
    inputnode_names = pe.Node(interface=util.IdentityInterface(fields= nodeOutputNames), 
                              name="filenames")    
    
    
    rename_node = pe.MapNode(interface = util.Rename(), name = 'rename', 
                             iterfield=['in_file', 'format_string'])
    
    rename_node1=rename_node.clone('rename_node1')
    rename_node2=rename_node.clone('rename_node2')
    rename_node3=rename_node.clone('rename_node3')
    rename_node4=rename_node.clone('rename_node4')
    rename_node5=rename_node.clone('rename_node5')
    rename_node6=rename_node.clone('rename_node6')
    rename_node7=rename_node.clone('rename_node7')
    rename_node8=rename_node.clone('rename_node8')
    rename_node9=rename_node.clone('rename_node9')
    rename_node10=rename_node.clone('rename_node10')
    rename_node11=rename_node.clone('rename_node11')
    rename_node12=rename_node.clone('rename_node12')
    rename_node13=rename_node.clone('rename_node13')
    rename_node14=rename_node.clone('rename_node14')
    rename_node15=rename_node.clone('rename_node15')
    
    format_node = pe.Node(interface=util.Function(input_names=['in_file', 'filename'], 
                                                  output_names= ['format_string'], 
                                                  function=formatpath ), 
                          name='format_node')
    
    format_node1=format_node.clone('format_node1')
    format_node2=format_node.clone('format_node2')
    format_node3=format_node.clone('format_node3')
    format_node4=format_node.clone('format_node4')
    format_node5=format_node.clone('format_node5')
    format_node6=format_node.clone('format_node6')
    format_node7=format_node.clone('format_node7')
    format_node8=format_node.clone('format_node8')
    format_node9=format_node.clone('format_node9')
    format_node10=format_node.clone('format_node10')
    format_node11=format_node.clone('format_node11')
    format_node12=format_node.clone('format_node12')
    format_node13=format_node.clone('format_node13')
    format_node14=format_node.clone('format_node14')
    format_node15=format_node.clone('format_node15')
   

    preproc.connect(inputNode, 'csf_t12func', format_node, 'in_file')
    preproc.connect(inputNode, 'csf_mni2func', format_node1, 'in_file')
    preproc.connect(inputNode, 'csf_combo', format_node2, 'in_file')
    preproc.connect(inputNode, 'csf_bin', format_node3, 'in_file')
    preproc.connect(inputNode, 'csf_mask', format_node4, 'in_file')
    preproc.connect(inputNode, 'gm_t12func', format_node5, 'in_file')
    preproc.connect(inputNode, 'gm_mni2func', format_node6, 'in_file')
    preproc.connect(inputNode, 'gm_combo', format_node7, 'in_file')
    preproc.connect(inputNode, 'gm_bin', format_node8, 'in_file')
    preproc.connect(inputNode, 'gm_mask', format_node9, 'in_file')
    preproc.connect(inputNode, 'global_mask', format_node10, 'in_file')
    preproc.connect(inputNode, 'wm_t12func', format_node11, 'in_file')
    preproc.connect(inputNode, 'wm_mni2func', format_node12, 'in_file')
    preproc.connect(inputNode, 'wm_combo', format_node13, 'in_file')
    preproc.connect(inputNode, 'wm_bin', format_node14, 'in_file')
    preproc.connect(inputNode, 'wm_mask', format_node15, 'in_file')
   
   
    preproc.connect(inputnode_names, 'csf_t12func', format_node, 'filename')
    preproc.connect(inputnode_names, 'csf_mni2func', format_node1, 'filename')
    preproc.connect(inputnode_names, 'csf_combo', format_node2, 'filename')
    preproc.connect(inputnode_names, 'csf_bin', format_node3, 'filename')
    preproc.connect(inputnode_names, 'csf_mask', format_node4, 'filename')
    preproc.connect(inputnode_names, 'gm_t12func', format_node5, 'filename')
    preproc.connect(inputnode_names, 'gm_mni2func', format_node6, 'filename')
    preproc.connect(inputnode_names, 'gm_combo', format_node7, 'filename')
    preproc.connect(inputnode_names, 'gm_bin', format_node8, 'filename')
    preproc.connect(inputnode_names, 'gm_mask', format_node9, 'filename')
    preproc.connect(inputnode_names, 'global_mask', format_node10, 'filename')
    preproc.connect(inputnode_names, 'wm_t12func', format_node11, 'filename')
    preproc.connect(inputnode_names, 'wm_mni2func', format_node12, 'filename')
    preproc.connect(inputnode_names, 'wm_combo', format_node13, 'filename')
    preproc.connect(inputnode_names, 'wm_bin', format_node14, 'filename')
    preproc.connect(inputnode_names, 'wm_mask', format_node15, 'filename')
  
        
    preproc.connect(format_node, 'format_string', rename_node, 'format_string')
    preproc.connect(format_node1, 'format_string', rename_node1, 'format_string')
    preproc.connect(format_node2, 'format_string', rename_node2, 'format_string')
    preproc.connect(format_node3, 'format_string', rename_node3, 'format_string')
    preproc.connect(format_node4, 'format_string', rename_node4, 'format_string')
    preproc.connect(format_node5, 'format_string', rename_node5, 'format_string')
    preproc.connect(format_node6, 'format_string', rename_node6, 'format_string')
    preproc.connect(format_node7, 'format_string', rename_node7, 'format_string')
    preproc.connect(format_node8, 'format_string', rename_node8, 'format_string')
    preproc.connect(format_node9, 'format_string', rename_node9, 'format_string')
    preproc.connect(format_node10, 'format_string', rename_node10, 'format_string')
    preproc.connect(format_node11, 'format_string', rename_node11, 'format_string')
    preproc.connect(format_node12, 'format_string', rename_node12, 'format_string')
    preproc.connect(format_node13, 'format_string', rename_node13, 'format_string')
    preproc.connect(format_node14, 'format_string', rename_node14, 'format_string')
    preproc.connect(format_node15, 'format_string', rename_node15, 'format_string')
   
    
    preproc.connect(inputNode, 'csf_t12func', rename_node, 'in_file')
    preproc.connect(inputNode, 'csf_mni2func', rename_node1, 'in_file')
    preproc.connect(inputNode, 'csf_combo', rename_node2, 'in_file')
    preproc.connect(inputNode, 'csf_bin', rename_node3, 'in_file')
    preproc.connect(inputNode, 'csf_mask', rename_node4, 'in_file')
    preproc.connect(inputNode, 'gm_t12func', rename_node5, 'in_file')
    preproc.connect(inputNode, 'gm_mni2func', rename_node6, 'in_file')
    preproc.connect(inputNode, 'gm_combo', rename_node7, 'in_file')
    preproc.connect(inputNode, 'gm_bin', rename_node8, 'in_file')
    preproc.connect(inputNode, 'gm_mask', rename_node9, 'in_file')
    preproc.connect(inputNode, 'global_mask', rename_node10, 'in_file')
    preproc.connect(inputNode, 'wm_t12func', rename_node11, 'in_file')
    preproc.connect(inputNode, 'wm_mni2func', rename_node12, 'in_file')
    preproc.connect(inputNode, 'wm_combo', rename_node13, 'in_file')
    preproc.connect(inputNode, 'wm_bin', rename_node14, 'in_file')
    preproc.connect(inputNode,  'wm_mask', rename_node15, 'in_file')
 
        
      
    preproc.connect(rename_node, 'out_file', outputNode,'csf_t12func_out'  )
    preproc.connect(rename_node1, 'out_file', outputNode,'csf_mni2func_out' )
    preproc.connect(rename_node2, 'out_file', outputNode, 'csf_combo_out' )
    preproc.connect(rename_node3, 'out_file', outputNode, 'csf_bin_out')
    preproc.connect(rename_node4, 'out_file', outputNode, 'csf_mask_out')
    preproc.connect(rename_node5, 'out_file', outputNode, 'gm_t12func_out')
    preproc.connect(rename_node6, 'out_file', outputNode, 'gm_mni2func_out')
    preproc.connect(rename_node7, 'out_file', outputNode, 'gm_combo_out')
    preproc.connect(rename_node8, 'out_file', outputNode, 'gm_bin_out')
    preproc.connect(rename_node9, 'out_file', outputNode, 'gm_mask_out')
    preproc.connect(rename_node10, 'out_file', outputNode, 'global_mask_out')
    preproc.connect(rename_node11, 'out_file', outputNode, 'wm_t12func_out' )
    preproc.connect(rename_node12, 'out_file', outputNode, 'wm_mni2func_out')
    preproc.connect(rename_node13, 'out_file', outputNode, 'wm_combo_out')
    preproc.connect(rename_node14, 'out_file', outputNode, 'wm_bin_out')
    preproc.connect(rename_node15, 'out_file', outputNode, 'wm_mask_out')
    
    return preproc
  
    
def rename_nuisance_preproc():
    
    preproc = pe.Workflow(name='rename_nuisance_preproc')
    
    inputNode = pe.Node(util.IdentityInterface(fields=['residual_file']),
                        name='inputspec')
        
    outputNode = pe.Node(util.IdentityInterface(fields=['residual_file_out']),
                        name='outputspec')
    
    nodeOutputNames=['residual_file']
    inputnode_names = pe.Node(interface=util.IdentityInterface(fields= nodeOutputNames), 
                              name="filenames")    
    
    rename_node = pe.MapNode(interface = util.Rename(), name = 'rename', 
                             iterfield=['in_file', 'format_string'])
    
    format_node = pe.Node(interface=util.Function(input_names=['in_file', 'filename'], 
                                                  output_names= ['format_string'], 
                                                  function=formatpath ), 
                          name='format_node')
    
    preproc.connect(inputNode, 'residual_file', format_node, 'in_file')
    preproc.connect(inputnode_names, 'residual_file', format_node, 'filename')
    
    preproc.connect(format_node, 'format_string', rename_node, 'format_string')
    preproc.connect(inputNode, 'residual_file', rename_node, 'in_file')
    
    preproc.connect(rename_node, 'out_file', outputNode, 'residual_file_out')
    
    return preproc


def rename_scrubbing_preproc():
    
    preproc = pe.Workflow(name='rename_scrubbing_preproc')
    
    inputNode = pe.Node(util.IdentityInterface(fields=['mean_deriv_sq_1D', 
                                                       'mean_raw_sq_1D', 
                                                       'scrubbed_preprocessed']),
                        name='inputspec')
        
    outputNode = pe.Node(util.IdentityInterface(fields=['mean_deriv_sq_1D_out',
                                                        'mean_raw_sq_1D_out', 
                                                        'scrubbed_preprocessed_out']),
                        name='outputspec')
    
    nodeOutputNames=['mean_deriv_sq_1D', 'mean_raw_sq_1D', 'scrubbed_preprocessed']
    inputnode_names = pe.Node(interface=util.IdentityInterface(fields= nodeOutputNames), 
                              name="filenames")   
        
    
    rename_node = pe.MapNode(interface = util.Rename(), name = 'rename', 
                             iterfield=['in_file', 'format_string'])
    
    rename_node1=rename_node.clone('rename_node1')
    rename_node2=rename_node.clone('rename_node2')

    
    format_node = pe.Node(interface=util.Function(input_names=['in_file', 'filename'], 
                                                  output_names= ['format_string'], 
                                                  function=formatpath ), 
                          name='format_node')
    
    format_node1=format_node.clone('format_node1')
    format_node2=format_node.clone('format_node2')

    preproc.connect(inputNode, 'mean_deriv_sq_1D', format_node, 'in_file')
    preproc.connect(inputNode, 'mean_raw_sq_1D', format_node1, 'in_file')
    preproc.connect(inputNode, 'scrubbed_preprocessed', format_node2, 'in_file')
    
    preproc.connect(inputnode_names, 'mean_deriv_sq_1D', format_node, 'filename')
    preproc.connect(inputnode_names, 'mean_raw_sq_1D', format_node1, 'filename')
    preproc.connect(inputnode_names, 'scrubbed_preprocessed', format_node2, 'filename')
       
    preproc.connect(format_node, 'format_string', rename_node, 'format_string')
    preproc.connect(format_node1, 'format_string', rename_node1, 'format_string')
    preproc.connect(format_node2, 'format_string', rename_node2, 'format_string')
 
    preproc.connect(inputNode, 'mean_deriv_sq_1D', rename_node, 'in_file')
    preproc.connect(inputNode, 'mean_raw_sq_1D', rename_node1, 'in_file')
    preproc.connect(inputNode, 'scrubbed_preprocessed', rename_node2, 'in_file')
      
    preproc.connect(rename_node, 'out_file', outputNode, 'mean_deriv_sq_1D_out')
    preproc.connect(rename_node1, 'out_file', outputNode, 'mean_raw_sq_1D_out')
    preproc.connect(rename_node2, 'out_file', outputNode, 'scrubbed_preprocessed_out')


    return preproc


