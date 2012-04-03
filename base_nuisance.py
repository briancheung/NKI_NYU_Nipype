import nipype.pipeline.engine as pe 
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util

#Extraction commands are independent of nipype workflows, should be int heir own utils file.
def extract_compcor_components(num_components, realigned_file,
                             wm_mask, csf_mask):
    import os
    import nibabel as nb
    import scipy as sp
    import numpy as np
    from scipy.signal import detrend

    data = nb.load(realigned_file).get_data().astype('float64')
    wm_mask = nb.load(wm_mask).get_data().astype('float64')
    csf_mask = nb.load(csf_mask).get_data().astype('float64')
    print 'Data and masks loaded.'
    wmcsf_mask = (csf_mask+wm_mask).astype('bool')

    print 'Detrending and centering data'
    Y = detrend(data[wmcsf_mask], axis=1, type='linear').T
    Yc = Y - np.tile(Y.mean(0), (Y.shape[0], 1))

    print 'Calculating SVD decomposition of Y*Y\''
    U,S,Vh = np.linalg.svd(np.dot(Yc,Yc.T))
     
    components_file = os.path.join(os.getcwd(), 'noise_components.txt')
    print 'Saving components file:', components_file
    np.savetxt(components_file, U[:, :num_components])
    return components_file

def extract_global_component(realigned_file):
    import os
    import nibabel as nb
    import numpy as np
    from utils import mean_roi_signal

    data = nb.load(realigned_file).get_data().astype('float64')
    mask = (data != 0).sum(-1) != 0 #Global Mask 
    print 'Data loaded.'
#    Y = data[mask].T
#    Yc = Y - np.tile(Y.mean(0), (Y.shape[0], 1))
    glb_comp = mean_roi_signal(data, mask) 

    components_file = os.path.join(os.getcwd(), 'global_component.txt')
    print 'Saving components file:', components_file
    np.savetxt(components_file, glb_comp)

    return components_file 

def extract_wmcsf_components(realigned_file, wm_mask, csf_mask):
    import os
    import nibabel as nb
    import numpy as np
    from utils import mean_roi_signal

    data = nb.load(realigned_file).get_data().astype('float64')
    wm_mask = nb.load(wm_mask).get_data().astype('float64')
    csf_mask = nb.load(csf_mask).get_data().astype('float64')
    print 'Data and masks loaded.'
    wm_comp = mean_roi_signal(data, wm_mask.astype('bool'))
    csf_comp = mean_roi_signal(data, csf_mask.astype('bool')) 

    wmcsf_comp = np.vstack((wm_comp,csf_comp)).T

    components_file = os.path.join(os.getcwd(), 'wmcsf_components.txt') 
    print 'Saving components file:', components_file
    np.savetxt(components_file, wmcsf_comp)

    return components_file

def extract_firstprinc_component(realigned_file):
    import os
    import nibabel as nb
    import numpy as np

    data = nb.load(realigned_file).get_data().astype('float64')
    mask = (data != 0).sum(-1) != 0 #Global Mask 
    print 'Data loaded.'
    Y = data[mask].T
    Yc = Y - np.tile(Y.mean(0), (Y.shape[0], 1))

    print 'Calculating SVD decomposition of Y'
    U,S,Vh = np.linalg.svd(Yc, full_matrices=False)
    
    components_file = os.path.join(os.getcwd(), 'firstprinc_component.txt')
    print 'Saving components file:', components_file
    np.savetxt(components_file, U[:,0])

    return components_file 
 
#Nuisance selection structure based on https://github.com/satra/BrainImagingPipelines/tree/master/fmri
def create_filter_matrix(global_component, compcor_components, 
                         wmcsf_components, firstprinc_component,
                         selector):
    import numpy as np
    import os

    def try_import(fname):
        try:
            a = np.genfromtxt(fname)
            return a
        except:
            return np.array([])
    
    options = np.array([global_component, compcor_components, wmcsf_components, firstprinc_component])
    fieldnames = ['global', 'compcor', 'wmcsf', 'firstprinc']

    selector = np.array(selector) #Use selector as an index mask
    #Grab component filenames of according to selector
    filenames = [fieldnames[i] for i, val in enumerate(selector) if val]
    filter_file = os.path.abspath('filter_%s.txt' % '_'.join(filenames))

    z = None
    for i,opt in enumerate(options[selector]):
        a = try_import(opt)
        if len(a.shape) == 1:
            a = np.array([a]).T
        print a.shape
        if i == 0:
            z = a
        else:
            z = np.hstack((z,a))


    print 'Writing filter design matrix of size', z.shape, 'to file', filter_file
    np.savetxt(filter_file, z)
    return filter_file 

def create_compcor_extraction(name='compcor_extract'):
    inputspec = pe.Node(util.IdentityInterface(fields=['num_components',
                                                       'realigned_file',
                                                       'wm_mask',
                                                       'csf_mask']),
                        name='inputspec')
    outputspec = pe.Node(util.IdentityInterface(fields=['noise_components']),
                         name='outputspec')

    compcor = pe.MapNode(util.Function(input_names=['num_components',
                                                    'realigned_file',
                                                    'wm_mask',
                                                    'csf_mask'],
                                       output_names=['noise_components'],
                                       function=extract_compcor_components),
                                       name='compcor_components',
                                       iterfield=['realigned_file',
                                                  'wm_mask',
                                                  'csf_mask'])
    comproc = pe.Workflow(name=name) 
    comproc.connect(inputspec, 'num_components', compcor, 'num_components')
    comproc.connect(inputspec, 'realigned_file', compcor, 'realigned_file')
    comproc.connect(inputspec, 'wm_mask', compcor, 'wm_mask')
    comproc.connect(inputspec, 'csf_mask', compcor, 'csf_mask')
    comproc.connect(compcor, 'noise_components', output_spec, 'noise_components')
 
    return comproc 

def create_nuisance_preproc(name='nuisance_preproc'):
    inputspec = pe.Node(util.IdentityInterface(fields=['num_components',
                                                       'realigned_file',
                                                       'wm_mask',
                                                       'csf_mask',
                                                       'selector']),
                        name='inputspec')
    outputspec = pe.Node(util.IdentityInterface(fields=['residual_file']),
                         name='outputspec')

    nuisance_preproc = pe.Workflow(name=name)

    compcor = pe.MapNode(util.Function(input_names=['num_components',
                                                    'realigned_file',
                                                    'wm_mask',
                                                    'csf_mask'],
                                       output_names=['noise_components'],
                                       function=extract_compcor_components),
                                       name='compcor',
                                       iterfield=['realigned_file',
                                                  'wm_mask',
                                                  'csf_mask'])
    glb_sig = pe.MapNode(util.Function(input_names=['realigned_file'],
                                       output_names=['global_component'],
                                       function=extract_global_component),
                                       name='glb_sig',
                                       iterfield=['realigned_file'])
    wmcsf_sig = pe.MapNode(util.Function(input_names=['realigned_file',
                                                      'wm_mask',
                                                      'csf_mask'],
                                         output_names=['wmcsf_components'],
                                         function=extract_wmcsf_components),
                                         name='wmcsf_sig',
                                         iterfield=['realigned_file',
                                                    'wm_mask',
                                                    'csf_mask'])
    fp1_sig = pe.MapNode(util.Function(input_names=['realigned_file'],
                                       output_names=['firstprinc_component'],
                                       function=extract_firstprinc_component),
                                       name='fp1_sig',
                                       iterfield=['realigned_file'])

    addoutliers = pe.MapNode(util.Function(input_names=['global_component', 
                                                        'compcor_components', 
                                                        'wmcsf_components',
                                                        'firstprinc_component',
                                                        'selector'],
                                           output_names=['filter_file'],
                                           function=create_filter_matrix),
                                           name='create_nuisance_filter',
                                           iterfield=['global_component', 
                                                      'compcor_components',
                                                      'wmcsf_components',
                                                      'firstprinc_component'])

    remove_noise = pe.MapNode(fsl.FilterRegressor(filter_all=True),
                              name='regress_nuisance',
                              iterfield=['design_file','in_file'])

    nuisance_preproc.connect(inputspec, 'realigned_file',
                             compcor, 'realigned_file')
    nuisance_preproc.connect(inputspec, 'num_components',
                             compcor, 'num_components')
    nuisance_preproc.connect(inputspec, 'wm_mask',
                             compcor, 'wm_mask')
    nuisance_preproc.connect(inputspec, 'csf_mask',
                             compcor, 'csf_mask')
    nuisance_preproc.connect(inputspec, 'realigned_file',
                             glb_sig, 'realigned_file') 
    nuisance_preproc.connect(inputspec, 'realigned_file',
                             wmcsf_sig, 'realigned_file')
    nuisance_preproc.connect(inputspec, 'wm_mask',
                             wmcsf_sig, 'wm_mask')
    nuisance_preproc.connect(inputspec, 'csf_mask',
                             wmcsf_sig, 'csf_mask')
    nuisance_preproc.connect(inputspec, 'realigned_file',
                             fp1_sig, 'realigned_file')

    nuisance_preproc.connect(glb_sig, 'global_component',
                             addoutliers, 'global_component')
    nuisance_preproc.connect(compcor, 'noise_components',
                             addoutliers, 'compcor_components')
    nuisance_preproc.connect(wmcsf_sig, 'wmcsf_components',
                             addoutliers, 'wmcsf_components')
    nuisance_preproc.connect(fp1_sig, 'firstprinc_component',
                             addoutliers, 'firstprinc_component')
    nuisance_preproc.connect(inputspec, 'selector',
                             addoutliers, 'selector')
    nuisance_preproc.connect(addoutliers, 'filter_file',
                             remove_noise, 'design_file')
    nuisance_preproc.connect(inputspec, 'realigned_file',
                             remove_noise, 'in_file')
    nuisance_preproc.connect(remove_noise, 'out_file',
                             outputspec, 'residual_file')

    return nuisance_preproc
