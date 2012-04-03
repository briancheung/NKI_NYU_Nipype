def create_compcor_extraction():
    inputspec = pe.Node(util.IdentityInterface(fields=['num_components',
                                                       'realigned_file',
                                                       'wm_mask',
                                                       'csf_mask']),
                        name='inputspec')
    outputspec = pe.Node(util.IdentityInterface(fields=['noise_components',
                                                        'residual_file']),
                         name='outputspec')


def create_global_extraction():

def create_wm_extraction():

def create_grey_extraction():

def create_csf_extraction():

def create_pc1_extraction():

def create_motion_extraction():

 

def create_nuisance_matrix(motion_time, composit_norm, compcorr_components):

def create_nuisance_preproc():
    
