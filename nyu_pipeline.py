import os
import nipype.interfaces.io as nio
import nipype.pipeline.engine as pe 
import re
import glob

from nipype.interfaces.utility import IdentityInterface
from nipype.interfaces import fsl

subject_paths = glob.glob('/home/data/Originals/NYU_TRT/NYU_TRT_session1/sub*')
subject_ids = [re.search('(sub\w+)',sp).group(1) for sp in subject_paths]
session_ids = ['NYU_TRT_session1', 'NYU_TRT_session2']
data_path = os.path.abspath('/home/data/Originals/NYU_TRT/')
output_path = os.path.abspath('/home/data/Projects/nuisance_reliability_paper/')

inputnode = pe.Node(interface=IdentityInterface(fields=['session_id', 'subject_id'],
                                                mandatory_inputs=True),
                    name='inputnode')
inputnode.iterables = [('session_id', session_ids),
                       ('subject_id', subject_ids)]

datasource = pe.Node(interface=nio.DataGrabber(infields=['session_id','subject_id'],
                                               outfields=['func', 'struct']),
                                               name='datasource')

datasource.inputs.base_directory = data_path
datasource.inputs.template = '*' 
datasource.inputs.field_template = dict(func='%s/%s/func/lfo.nii.gz',
                                        struct='%s/%s/anat/mprage_anonymized.nii.gz')

datasource.inputs.template_args = dict(func=[['session_id',  'subject_id']],
                                       struct=[['session_id', 'subject_id']])

st = pe.Node(interface=fsl.SliceTimer(), name='st')
st.iterables = [('time_repetition', [2.0, 3.0])]

datasink = pe.Node(nio.DataSink(infields=['container.@subject']), name='sinker')
datasink.inputs.base_directory = output_path
#datasink.inputs.substitutions = [('session\d_subject_id_', '')]

workflow = pe.Workflow(name='wf')
workflow.base_dir = output_path

workflow.connect(inputnode, 'session_id', datasource, 'session_id')
workflow.connect(inputnode, 'subject_id', datasource, 'subject_id')
workflow.connect(datasource, 'func', st, 'in_file')
workflow.connect(inputnode, 'session_id', datasink, 'container')
workflow.connect(st, 'slice_time_corrected_file', datasink, 'slicecorrected')
