import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
import nipype.pipeline.engine as pe


def rename_connections(workflow, datasink, rename_list, sink_node):

    ncount = 0
    rcount = 0
    for rename in rename_list:
        din_file = sink_node + '.@' + str(ncount)
        if len(rename) == 4:
            rename_node = pe.MapNode(interface=util.Rename(),
                                     name=rename[2] + str(rcount),
                                     iterfield=['in_file', 'format_string'])
            rename_node.inputs.format_string = rename[3]

            workflow.connect(rename[0], rename[1], rename_node, 'in_file')
            workflow.connect(rename_node, 'out_file', datasink, din_file)
            rcount += 1
        else:
            workflow.connect(rename[0], rename[1], datasink, din_file)
        ncount += 1


def anat_sink(workflow, datasink, mprage_mni):

    rename_list = [(mprage_mni, 'outputspec.brain_mni',
                    'anat_rename', 'mprage_brian.nii.gz'),
                   (mprage_mni, 'outputspec.reorient_mni',
                    'anat_rename', 'mprage_RPI.nii.gz')]
    rename_connections(workflow, datasink, rename_list, 'anat')


def func_sink(workflow, datasink, funcpreproc, func_in_mni):

    rename_list = [(funcpreproc, 'outputspec.motion_correct',
                    'func_rename', 'rest_mc.nii.gz'),
                   (funcpreproc, 'outputspec.motion_correct_ref',
                    'func_rename', 'rest_mc_init_mean.nii.gz'),
                   (funcpreproc, 'outputspec.movement_parameters',
                    'func_rename', 'rest_mc.1D'),
                   (funcpreproc, 'outputspec.max_displacement',
                    'func_rename', 'rest_maxdisp.1D'),
                   (funcpreproc, 'outputspec.preprocessed',
                    'func_rename', 'rest_pp.nii.gz'),
                   (func_in_mni, 'outputspec.preprocessed_mask_mni',
                    'func_rename', 'rest_pp_mask.nii.gz')]
    rename_connections(workflow, datasink, rename_list, 'func')


def reg_sink(workflow, datasink, regpreproc):

    rename_list = [(regpreproc, 'outputspec.highres2standard',
                    'reg_rename', 'highres2standard.nii.gz'),
                   (regpreproc, 'outputspec.example_func2highres_mat',
                    'reg_rename', 'example_func2highres.mat'),
                   (regpreproc, 'outputspec.example_func2highres',
                    'reg_rename', 'example_func2highres.nii.gz'),
                   (regpreproc, 'outputspec.example_func2standard_NL',
                    'reg_rename', 'example_func2standard_NL.nii.gz'),
                   (regpreproc, 'outputspec.highres2example_func_mat',
                    'reg_rename', 'highres2example_func.mat')]
    rename_connections(workflow, datasink, rename_list, 'reg')


def seg_sink(workflow, datasink, segpreproc, mprage_mni):

    rename_list = [(segpreproc, 'outputspec.csf_t12func',
                    'seg_rename', 'csf_t12func.nii.gz'),
                   (segpreproc, 'outputspec.csf_mni2func',
                    'seg_reanme', 'csf_mni2func.nii.gz'),
                   (segpreproc, 'outputspec.csf_combo',
                    'seg_rename', 'csf_combo.nii.gz'),
                   (segpreproc, 'outputspec.csf_bin',
                    'seg_rename', 'csf_bin.nii.gz'),
                   (segpreproc, 'outputspec.csf_mask',
                    'seg_rename', 'csf_mask.nii.gz'),
                   (segpreproc, 'outputspec.gm_t12func',
                    'seg_rename', 'gm_t12func.nii.gz'),
                   (segpreproc, 'outputspec.gm_mni2func',
                    'seg_rename', 'gm_mni2func.nii.gz'),
                   (segpreproc, 'outputspec.gm_combo',
                    'seg_rename', 'gm_combo.nii.gz'),
                   (segpreproc, 'outputspec.gm_bin',
                    'seg_rename', 'gm_bin.nii.gz'),
                   (segpreproc, 'outputspec.gm_mask',
                    'seg_rename', 'gm_mask.nii.gz'),
                   (segpreproc, 'outputspec.global_mask',
                    'seg_rename', 'global_mask.nii.gz'),
                   (segpreproc, 'outputspec.wm_t12func',
                    'seg_rename', 'wm_t12func.nii.gz'),
                   (segpreproc, 'outputspec.wm_mni2func',
                    'seg_rename', 'wm_mni2func.nii.gz'),
                   (segpreproc, 'outputspec.wm_combo',
                    'seg_rename', 'wm_combo.nii.gz'),
                   (segpreproc, 'outputspec.wm_bin',
                    'seg_rename', 'wm_bin.nii.gz'),
                   (segpreproc, 'outputspec.wm_mask',
                    'seg_rename', 'wm_mask.nii.gz'),
                   (segpreproc, 'outputspec.probability_maps'),
                   (segpreproc, 'outputspec.mixeltype'),
                   (segpreproc, 'outputspec.partial_volume_map'),
                   (segpreproc, 'outputspec.partial_volume_files'),
                   (mprage_mni, 'outputspec.mixeltype_mni'),
                   (mprage_mni, 'outputspec.probability_maps_mni'),
                   (mprage_mni, 'outputspec.partial_volume_map_mni'),
                   (mprage_mni, 'outputspec.partial_volume_files_mni')]
    rename_connections(workflow, datasink, rename_list, 'segment')


def nuisance_sink(workflow, datasink, nuisancepreproc):
    rename_list = [(nuisancepreproc, 'outputspec.residual_file',
                    'rename', 'rest_residual.nii.gz'),
                   (nuisancepreproc, 'outputspec.median_angle_corrected_file')]
    rename_connections(workflow, datasink, rename_list, 'nuisance')


def scrubbing_sink(workflow, datasink, scpreproc):

    rename_list = [(scpreproc, 'outputspec.mean_deriv_sq_1D',
                    'sc_rename', 'mean_deriv_sq.1D'),
                   (scpreproc, 'outputspec.mean_raw_sq_1D',
                    'sc_rename', 'mean_raw_sq.1D'),
                   (scpreproc, 'outputspec.scrubbed_preprocessed',
                    'sc_rename', 'rest_pp_scrubbed.nii.gz'),
                   (scpreproc, 'outputspec.temp_deriv_brik_file'),
                   (scpreproc, 'outputspec.temp_deriv_head_file'),
                   (scpreproc, 'outputspec.temp_deriv_sq_brik_file'),
                   (scpreproc, 'outputspec.temp_deriv_sq_head_file'),
                   (scpreproc, 'outputspec.raw_sq_brik_file'),
                   (scpreproc, 'outputspec.raw_sq_head_file'),
                   (scpreproc, 'outputspec.mask_brik_file'),
                   (scpreproc, 'outputspec.mask_head_file'),
                   (scpreproc, 'outputspec.FD_1D'),
                   (scpreproc, 'outputspec.sqrt_mean_deriv_sq_1D'),
                   (scpreproc, 'outputspec.sqrt_mean_raw_sq_1D'),
                   (scpreproc, 'outputspec.frames_ex_1D'),
                   (scpreproc, 'outputspec.frames_in_1D'),
                   (scpreproc, 'outputspec.pow_params'),
                   (scpreproc, 'outputspec.ftof_percent_change_1D'),
                   (scpreproc, 'outputspec.scrubbed_movement_parameters')]
    rename_connections(workflow, datasink, rename_list, 'scrubbing')


def sca_sink(workflow, datasink, scapreproc):
    rename_list = [
                   (scapreproc, 'outputspec.correlations',
                    'sca_rename', 'correlations.nii.gz'),
                   (scapreproc, 'outputspec.Z_trans_correlations',
                    'sca_rename', 'Z.nii.gz'),
                   (scapreproc, 'outputspec.Z_2standard',
                    'sca_rename', 'Z_2standard.nii.gz'),
                   (scapreproc, 'outputspec.Z_2standard_FWHM',
                    'sca_rename', 'Z_FWHM_2standard.nii.gz')]
    rename_connections(workflow, datasink, rename_list, 'sca')


def alff_sink(workflow, datasink, alffpreproc):
    rename_list = [(alffpreproc, 'outputspec.power_spectrum_distribution',
                    'alff_rename', 'power_spectrum_distribution.nii.gz'),
                   (alffpreproc, 'outputspec.alff_img',
                    'alff_rename', 'ALFF.nii.gz'),
                   (alffpreproc, 'outputspec.falff_img',
                    'alff_rename', 'fALFF.nii.gz'),
                   (alffpreproc, 'outputspec.alff_Z_img',
                    'alff_rename', 'ALFF_Z.nii.gz'),
                   (alffpreproc, 'outputspec.falff_Z_img',
                    'alff_rename', 'fALFF_Z.nii.gz'),
                   (alffpreproc, 'outputspec.alff_Z_2standard_img',
                    'alff_rename', 'ALFF_Z_2standard.nii.gz'),
                   (alffpreproc, 'outputspec.falff_Z_2standard_img',
                    'alff_rename', 'fALFF_Z_2standard.nii.gz'),
                   (alffpreproc, 'outputspec.alff_Z_2standard_fwhm_img',
                    'alff_rename', 'ALFF_Z_FWHM_2standard.nii.gz'),
                   (alffpreproc, 'outputspec.falff_Z_2standard_fwhm_img',
                    'alff_rename', 'fALFF_Z_FWHM_2standard.nii.gz')]
    rename_connections(workflow, datasink, rename_list, 'alff')


def vmhc_sink(workflow, datasink, vmhcpreproc):
    rename_list = [(vmhcpreproc, 'outputspec.VMHC_img', 'vmhc_rename', 'VMHC.nii.gz'),
                   (vmhcpreproc, 'outputspec.VMHC_Z_img', 'vmhc_rename', 'VMHC_Z.nii.gz'),
                   (vmhcpreproc, 'outputspec.VMHC_Z_stat_img', 'vmhc_rename', 'VMHC_Z_stat.nii.gz')]
    rename_connections(workflow, datasink, rename_list, 'vmhc')


def timeseries_sink(workflow, datasink, tspreproc, runSurfaceRegistration):

    rename_list = [(tspreproc, 'outputspec.mask_outputs'),
                  (tspreproc, 'outputspec.parc_outputs')]
    rename_connections(workflow, datasink, rename_list, 'timeseries')

    if runSurfaceRegistration:
        rename_connections(workflow, datasink,
                           [(tspreproc, 'outputspec.surface_outputs')],
                           'timeseries.@')


def group_analysis_sink(workflow, datasink, gppreproc):

    rename_connections(workflow, datasink,
                       [(gppreproc, 'outputspec.merged',
                        'gp_rename1', 'Merged.nii.gz')],
                       'Merged')
    rename_connections(workflow, datasink,
                       [(gppreproc, 'outputspec.zstats')],
                       'stats.unthreshold')
    rename_connections(workflow, datasink,
                       [(gppreproc, 'outputspec.cluster_threshold')],
                       'stats.threshold')
    rename_connections(workflow, datasink,
                       [(gppreproc, 'outputspec.cluster_index')],
                       'stats.clusterMap')
    rename_connections(workflow, datasink,
                       [(gppreproc, 'outputspec.overlay_threshold')],
                       'stats.Overlay')
    rename_connections(workflow, datasink,
                       [(gppreproc, 'outputspec.rendered_image')],
                       'rendered')
