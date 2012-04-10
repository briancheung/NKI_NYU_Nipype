# Utility Functions ---------------------------------------------------------
import e_afni
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util

def pToFile(time_series):

    import os
    import re
    import commands
    import sys

    dir  = os.path.dirname(time_series)

    dir1 = re.sub(r'(.)+_subject_id_(.)+/_mask', '', dir)

    dir1 = dir1.split('..')
    dir1 = dir1[len(dir1) -1]
    dir1 = dir1.split('.nii.gz')
    dir1 = dir1[0]

    ts_oneD = os.path.join(os.getcwd(), dir1 + '.1D')
    cmd = "cp %s %s" % (time_series, ts_oneD)
    print cmd

    sys.stderr.write('\n'+ commands.getoutput(cmd))
    return os.path.abspath(ts_oneD)


def mean_roi_signal(data_volume, roi_mask):
    import numpy as np
    Y = data_volume[roi_mask].T
    Yc = Y - np.tile(Y.mean(0), (Y.shape[0], 1))
    return Yc.mean(1)

def get_standard_background_img(in_file):
    import os

    from nibabel import load
    img = load(in_file)
    hdr = img.get_header()
    group_mm = int(hdr.get_zooms()[2])
    print "gorup_mm -> ", group_mm
    path = '/usr/local/fsl' + '/data/standard/MNI152_T1_%smm_brain.nii.gz' % (group_mm)
    print "path ->", path
    return os.path.abspath(path)

def get_nvols(in_file):

    from nibabel import load
    img = load(in_file)
    hdr = img.get_header()
    n_vol = int(hdr.get_data_shape()[3])
    op_string = '-abs -bin -Tmean -mul %d' % (n_vol)
    return op_string

def copyGeom(infile_a, infile_b):
    import subprocess as sb
    out_file = infile_b
    cmd = sb.Popen(['fslcpgeom', infile_a, out_file], stdin=sb.PIPE, stdout=sb.PIPE,)
    stdout_value, stderr_value = cmd.communicate()
    return out_file

def getTuple(infile_a, infile_b):

    print "inisde getTuple"
    print "infile_a -> ", infile_a
    print "infile_b -> ", infile_b[1]

    return (infile_a, infile_b[1])

def getOpString(mean, std_dev):

    str1 = "-sub %f -div %f" % (float(mean), float(std_dev))

    op_string = str1 + " -mas %s"

    return op_string

def pick_wm_0(probability_maps):

    import sys
    import os


    if(isinstance(probability_maps, list)):

        if(len(probability_maps) == 1):
            probability_maps = probability_maps[0]
        for file in probability_maps:
            print file
            if file.endswith("prob_0.nii.gz"):

                return file
    return None


def pick_wm_1(probability_maps):

    import sys
    import os


    if(isinstance(probability_maps, list)):

        if(len(probability_maps) == 1):
            probability_maps = probability_maps[0]
        for file in probability_maps:
            print file
            if file.endswith("prob_1.nii.gz"):

                return file
    return None


def pick_wm_2(probability_maps):

    import sys
    import os
    if(isinstance(probability_maps, list)):

        if(len(probability_maps) == 1):
            probability_maps = probability_maps[0]
        for file in probability_maps:
            print file
            if file.endswith("prob_2.nii.gz"):

                return file
    return None

def getImgNVols(in_files):

    out = []
    from nibabel import load
    if(isinstance(in_files, list)):
        for in_file in in_files:
            img = load(in_file)
            hdr = img.get_header()
            out.append(int(hdr.get_data_shape()[3]))
        return out

    else:
        img = load(in_files)
        hdr = img.get_header()
        n_vol = int(hdr.get_data_shape()[3])
        return [n_vol]

def getImgTR(in_files):

    out = []
    from nibabel import load
    if(isinstance(in_files,list)):
        for in_file in in_files:
            img = load(in_file)
            hdr = img.get_header()
            tr = float(hdr.get_zooms()[3])
            if tr > 10:
                tr = float(float(tr)/1000.0)
            out.append(tr)
        return out
    else:
        img = load(in_files)
        hdr = img.get_header()
        tr = float(hdr.get_zooms()[3])
        if tr > 10:
            tr = float(float(tr)/1000.0)
        return [tr]


def getN1(TR, nvols, HP):

    n_lp = float(HP) * float(int(nvols)) * float (TR)
    n1 = int("%1.0f" % (float(n_lp - 1.0)))

    print '>>>n_lp ', n_lp
    print '>>>n1 ', n1

    return n1

def getN2(TR, nvols, LP, HP):

    n_lp = float(HP) * float(int(nvols)) * float (TR)
    n_hp = float(LP) * float(int(nvols)) * float (TR)
    n2 = int("%1.0f" % (float(n_hp - n_lp + 1.0)))

    return n2

def takemod(nvols):

    decisions = []
    for vol in nvols:
        mod = int (int(vol) % 2)

        if mod == 1:
            decisions.append([0])
            #return [0]
        else:
            decisions.append([1])
            #return [1]

    return decisions


def set_op_str(n2):

    strs = []
    for n in n2:
        str = "-Tmean -mul %f" % (n)
        strs.append(str)
    return strs


def set_op1_str(nvols):

    strs = []
    for vol in nvols:
        str = '-Tmean -mul %d -div 2' % (int(vol))
        strs.append(str)

    return strs

def set_gauss(fwhm):


    op_string = ""

    fwhm = float(fwhm)

    sigma = float(fwhm/2.3548)

    op = "-kernel gauss %f -fmean -mas " % (sigma) + "%s"
    op_string = op

    return op_string

def getEXP(nvols):

    expr = []
    for vol in nvols:
        vol = int(vol)
        expr.append("'a*sqrt('%d'-3)'" % vol)

    return expr

def compcor():


    preproc = pe.Workflow(name='compcor')
    inputNode = pe.Node(util.IdentityInterface(fields=['csf_mask',
                                                    'wm_mask',
                                                    'preprocessed',
                                                    'ncomponents']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['csf_mask_erosion',
                                                            'wm_mask_erosion',
                                                            'preprocessed_compcor']),
                        name='outputspec')



    erosion_csf = pe.MapNode(interface=e_afni.Threedcalc(), name='erosion_csf', iterfield=['infile_a'])
    erosion_csf.inputs.infile_b_prime = 'a+i -c a-i -d a+j -e a-j -f a+k -g a-k'
    erosion_csf.inputs.expr = "'a*(1-amongst(0,b,c,d,e,f,g))'"


    erosion_csf1 = erosion_csf.clone('erosion_csf1')

    erosion_wm = pe.MapNode(interface=e_afni.Threedcalc(), name='erosion_wm', iterfield=['infile_a'])
    erosion_wm.inputs.infile_b_prime = 'a+i -c a-i -d a+j -e a-j -f a+k -g a-k'
    erosion_wm.inputs.expr = "'a*(1-amongst(0,b,c,d,e,f,g))'"

    erosion_wm1 = erosion_wm.clone('erosion_wm1')

    compcor = pe.MapNode(interface=e_afni.compcor(), name='compcor', iterfield=['in_file', 'wmmask', 'csfmask'])


    preproc.connect(inputNode, 'csf_mask', erosion_csf, 'infile_a')
    preproc.connect(erosion_csf, 'out_file', erosion_csf1, 'infile_a')

    preproc.connect(inputNode, 'wm_mask', erosion_wm, 'infile_a')
    preproc.connect(erosion_wm, 'out_file', erosion_wm1, 'infile_a')

    preproc.connect(inputNode, 'preprocessed', compcor, 'in_file')
    preproc.connect(erosion_csf1, 'out_file', compcor, 'csfmask')
    preproc.connect(erosion_wm1, 'out_file', compcor, 'wmmask')
    preproc.connect(inputNode, 'ncomponents', compcor, 'ncomponents')

    preproc.connect(erosion_csf1, 'out_file', outputNode, 'csf_mask_erosion')
    preproc.connect(erosion_wm1, 'out_file', outputNode, 'wm_mask_erosion')
    preproc.connect(compcor, 'out_file', outputNode, 'preprocessed_compcor')

    return preproc

def median_angle_correction():


    preproc = pe.Workflow(name='median_angle_correction')
    inputNode = pe.Node(util.IdentityInterface(fields=['angle',
                                                    'preprocessed']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['preprocessed_median_angle']),
                        name='outputspec')

    MedianAngle = pe.MapNode(interface=e_afni.MedianAngle(), name='MedianAngle', iterfield=['in_file'])

    preproc.connect(inputNode, 'preprocessed', MedianAngle, 'in_file')
    preproc.connect(inputNode, 'angle', MedianAngle, 'angle')

    preproc.connect(MedianAngle, 'out_file', outputNode, 'preprocessed_median_angle')

    return preproc

def createMC(in_file, pp):

    import os
    import sys
    import commands
    EV_Lists = []
    idx = 0

    if not(isinstance(in_file, list)):
        in_file = [in_file]

    for file in in_file:
        print file
        parent_path = os.path.dirname(pp[idx])

        EV_List = []
        cmd = "awk '{print $1}' %s > %s/mc1.1D" %(file, parent_path)
        print cmd
        sys.stderr.write(commands.getoutput(cmd))

        cmd = "awk '{print $2}' %s > %s/mc2.1D" %(file, parent_path)
        print cmd
        sys.stderr.write(commands.getoutput(cmd))

        cmd = "awk '{print $3}' %s > %s/mc3.1D" %(file, parent_path)
        print cmd
        sys.stderr.write(commands.getoutput(cmd))

        cmd = "awk '{print $4}' %s > %s/mc4.1D" %(file, parent_path)
        print cmd
        sys.stderr.write(commands.getoutput(cmd))

        cmd = "awk '{print $5}' %s > %s/mc5.1D" %(file, parent_path)
        print cmd
        sys.stderr.write(commands.getoutput(cmd))

        cmd = "awk '{print $6}' %s > %s/mc6.1D" %(file, parent_path)
        print cmd
        sys.stderr.write(commands.getoutput(cmd))

        EV_List.append(os.path.abspath(parent_path+'/mc1.1D'))
        EV_List.append(os.path.abspath(parent_path+'/mc2.1D'))
        EV_List.append(os.path.abspath(parent_path+'/mc3.1D'))
        EV_List.append(os.path.abspath(parent_path+'/mc4.1D'))
        EV_List.append(os.path.abspath(parent_path+'/mc5.1D'))
        EV_List.append(os.path.abspath(parent_path+'/mc6.1D'))

        EV_Lists.append(EV_List)
        idx += 1

    return EV_Lists

def createFSF(nuisance_template, rest_pp, TR, n_vols):

    import os
    import sys
    import commands
    nuisance_files = []
    idx = 0
    for rest_p in rest_pp:
        parent_path = os.path.dirname(rest_p)

        cmd = "sed -e s:nuisance_dir:\"%s\":g <%s >%s/temp1"\
            %(parent_path, nuisance_template, parent_path)
        print cmd
        sys.stderr.write(commands.getoutput(cmd))


        cmd = "sed -e s:nuisance_model_outputdir:\"%s/residuals.feat\":g <%s/temp1 >%s/temp2"\
            % (parent_path, parent_path, parent_path)
        print cmd
        sys.stderr.write(commands.getoutput(cmd))

        cmd = "sed -e s:nuisance_model_TR:\"%s\":g <%s/temp2 >%s/temp3" %(TR[idx],parent_path,parent_path)
        print cmd
        sys.stderr.write(commands.getoutput(cmd))

        cmd = "sed -e s:nuisance_model_numTRs:\"%s\":g <%s/temp3 >%s/temp4" %(n_vols[idx],parent_path,parent_path)
        print cmd
        sys.stderr.write(commands.getoutput(cmd))

        cmd = "sed -e s:nuisance_model_input_data:\"%s\":g <%s/temp4 >%s/nuisance.fsf" %(rest_p,parent_path,parent_path)
        print cmd
        sys.stderr.write(commands.getoutput(cmd))

    #   cmd = "rm %s/temp*" %(parent_path)
    #   print cmd
    #   sys.stderr.write(commands.getoutput(cmd))
        nuisance_files.append(os.path.abspath("%s/nuisance.fsf"%(parent_path)))

        idx += 1


    return nuisance_files


def copyStuff(EV_Lists, nuisance_files, global1Ds, csf1Ds, wm1Ds, regressors):

    idx = 0
    import os
    import sys
    import commands

    global1ds = []
    csf1ds = []
    wm1ds = []
    EV_final_lists_set = []
    for nuisance_file in nuisance_files:

        parent_path = os.path.dirname(nuisance_file)
        if (isinstance(global1Ds, list)):
            cmd = "cp %s %s/global.1D" % (global1Ds[idx], parent_path)
            print cmd
            sys.stderr.write(commands.getoutput(cmd))
            EV_Lists[idx].append("%s/global.1D" % (parent_path))

            cmd = "cp %s %s/csf.1D" % (csf1Ds[idx], parent_path)
            print cmd
            sys.stderr.write(commands.getoutput(cmd))
            EV_Lists[idx].append("%s/csf.1D" % (parent_path))

            cmd = "cp %s %s/wm.1D" % (wm1Ds[idx], parent_path)
            print cmd
            sys.stderr.write(commands.getoutput(cmd))

            EV_Lists[idx].append("%s/wm.1D" % (parent_path))

            path = os.path.join(parent_path, os.path.basename((regressors[idx])[0]))
            cmd = "cp %s %s" % (os.path.abspath((regressors[idx])[0]), path)
            print cmd
            sys.stderr.write(commands.getoutput(cmd))

            EV_Lists[idx].append("%s" % (path))

            path = os.path.join(parent_path, os.path.basename((regressors[idx])[1]))
            cmd = "cp %s %s" % (os.path.abspath((regressors[idx])[1]), path)
            print cmd
            sys.stderr.write(commands.getoutput(cmd))

            EV_Lists[idx].append("%s" % (path))

        elif (isinstance(csf1Ds, list)):
            cmd = "cp %s %s/csf.1D" % (csf1Ds[idx], parent_path)
            print cmd
            sys.stderr.write(commands.getoutput(cmd))
            EV_Lists[idx].append("%s/csf.1D" % (parent_path))

            cmd = "cp %s %s/wm.1D" % (wm1Ds[idx], parent_path)
            print cmd
            sys.stderr.write(commands.getoutput(cmd))

            EV_Lists[idx].append("%s/wm.1D" % (parent_path))

            path = os.path.join(parent_path, os.path.basename((regressors[idx])[0]))
            cmd = "cp %s %s" % (os.path.abspath((regressors[idx])[0]), path)
            print cmd
            sys.stderr.write(commands.getoutput(cmd))

            EV_Lists[idx].append("%s" % (path))

            path = os.path.join(parent_path, os.path.basename((regressors[idx])[1]))
            cmd = "cp %s %s" % (os.path.abspath((regressors[idx])[1]), path)
            print cmd
            sys.stderr.write(commands.getoutput(cmd))

            EV_Lists[idx].append("%s" % (path))

        EV_final_lists_set.append(EV_Lists[idx])
        idx += 1

    return EV_final_lists_set

def getStatsDir(in_files):

    import os

    stats_dir = []

    for in_file in in_files:

        parent_path = os.path.dirname(in_file)

        stats_dir.append(parent_path + '/stats')

    return stats_dir

def create_anat_dataflow(name, sublist, analysisdirectory, anat_name, at):

    import nipype.pipeline.engine as pe
    import nipype.interfaces.io as nio

    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'], outfields=['anat']), name=name)
    datasource.inputs.base_directory = analysisdirectory
    #datasource.inputs.template = '%s/*/%s.nii.gz'
    datasource.inputs.template = at
    datasource.inputs.template_args = dict(anat=[['subject_id', anat_name]])
    datasource.iterables = ('subject_id', sublist)

    return datasource

def create_func_dataflow(name, sublist, analysisdirectory, rest_name, rt):

    import nipype.pipeline.engine as pe
    import nipype.interfaces.io as nio

    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'], outfields=['rest']), name=name)
    datasource.inputs.base_directory = analysisdirectory
    #datasource.inputs.template = '%s/*/%s.nii.gz'
    datasource.inputs.template = rt
    datasource.inputs.template_args['rest'] = [['subject_id', rest_name]]
    datasource.iterables = ('subject_id', sublist)

    return datasource

def create_alff_dataflow(name, sublist, analysisdirectory, rest_name, at, atw):

    import nipype.pipeline.engine as pe
    import nipype.interfaces.io as nio

    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                                   outfields=['rest_res', 'rest_mask', 'rest_mask2standard' ]),
                         name=name)
    datasource.inputs.base_directory = analysisdirectory
    #datasource.inputs.template = '%s/*/%s.nii.gz'
    datasource.inputs.template = at
    datasource.inputs.template_args = dict(rest_res=[['subject_id', rest_name +'_res' ]],
                                           rest_mask=[['subject_id', rest_name + '_mask']],
                                           rest_mask2standard=[['subject_id', rest_name + '_mask2standard']])
    datasource.iterables = ('subject_id', sublist)

    datasource_warp = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                                        outfields=['premat', 'fieldcoeff_file' ]),
                              name=name + '_warp')

    datasource_warp.inputs.base_directory = analysisdirectory
    #datasource_warp.inputs.template = '%s/*/%s.nii.gz'
    datasource_warp.inputs.template = atw
    datasource_warp.inputs.template_args = dict(premat=[['subject_id', 'example_func2highres.mat' ]],
                                                fieldcoeff_file=[['subject_id', 'highres2standard_warp.nii.gz']])
    datasource_warp.iterables = ('subject_id', sublist)

    return datasource, datasource_warp

def create_ifc_dataflow(name, sublist, analysisdirectory, rt, rtw):

    import nipype.pipeline.engine as pe
    import nipype.interfaces.io as nio

    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                                   outfields=['rest_res2standard',
                                                              'rest_res_filt',
                                                              'rest_mask2standard',
                                                              'ref' ]),
                         name=name)
    datasource.inputs.base_directory = analysisdirectory
    #datasource.inputs.template = '%s/*/%s.nii.gz'
    datasource.inputs.template = rt
    datasource.inputs.template_args = dict(rest_res2standard=[['subject_id', rest_name + '_res2standard']],
                                           rest_res_filt=[['subject_id', rest_name + '_res_filt']],
                                           rest_mask2standard=[['subject_id', rest_name + '_mask2standard']],
                                           ref=[['subject_id', 'example_func' ]])
    datasource.iterables = ('subject_id', sublist)


    datasource_warp = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                                        outfields=['premat',
                                                                   'fieldcoeff_file',
                                                                   'warp',
                                                                   'postmat']),
                              name=name+'_warp')
    datasource_warp.inputs.base_directory = analysisdirectory
    #datasource_warp.inputs.template = '%s/*/%s.nii.gz'
    datasource_warp.inputs.template = rtw
    datasource_warp.inputs.template_args = dict(premat=[['subject_id', 'example_func2highres.mat' ]],
                                                fieldcoeff_file=[['subject_id', 'highres2standard_warp.nii.gz']],
                                                warp=[['subject_id', 'stand2highres_warp.nii.gz' ]],
                                                postmat=[['subject_id', 'highres2example_func.mat' ]])
    datasource_warp.iterables = ('subject_id', sublist)

    return datasource, datasource_warp



def create_vmhc_dataflow(name, sublist, analysisdirectory, anat_name, rest_name, vt, vta, vtw):

    import nipype.pipeline.engine as pe
    import nipype.interfaces.io as nio

    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                                   outfields=['rest_res']),
                         name=name+'_res')
    datasource.inputs.base_directory = analysisdirectory
    #datasource.inputs.template = '%s/*/%s.nii.gz'
    datasource.inputs.template = vt
    datasource.inputs.template_args = dict(rest_res=[['subject_id', rest_name + '_res']])
    datasource.iterables = ('subject_id', sublist)


    da = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                           outfields=['reorient', 'brain']),
                 name=name+'_RPI')
    da.inputs.base_directory = analysisdirectory
    da.inputs.template = vta
    da.inputs.template_args = dict(reorient=[['subject_id', anat_name + '_RPI']], brain=[['subject_id', anat_name + '_brain']])
    da.iterables = ('subject_id', sublist)


    datasource_warp = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                                        outfields=['example_func2highres_mat']),
                              name=name+'_example_func')
    datasource_warp.inputs.base_directory = analysisdirectory
    #datasource_warp.inputs.template = '%s/*/%s.nii.gz'
    datasource_warp.inputs.template = vtw
    datasource_warp.inputs.template_args = dict(example_func2highres_mat=[['subject_id', 'example_func2highres.mat' ]])
    datasource_warp.iterables = ('subject_id', sublist)

    return datasource, da, datasource_warp

def create_gp_dataflow(base_dir, modelist, seedlist, wfname= "datasource"):
    import nipype.pipeline.engine as pe
    import nipype.interfaces.io as nio
    
    #define your subject directory structure inside the base directory folder
    subject_dir='*/*/%s.nii.gz'
    
    datasource = pe.Node( interface = nio.DataGrabber(infields=['model_name','seed'], outfields = ['mat','con','fts','grp','seedfiles']) , name= 'datasource')
    datasource.inputs.base_directory = base_dir
    datasource.inputs.template= '*'
    datasource.inputs.field_template = dict(mat='*/%s/%s.mat',
                                        con='*/%s/%s.con',
                                        fts='*/%s/%s.fts',
                                        grp='*/%s/%s.grp',
                                        seedfiles=subject_dir)
    
    datasource.inputs.template_args = dict(mat=[['model_name', 'model_name']],
                                           con=[['model_name', 'model_name']],
                                           fts=[['model_name', 'model_name']],
                                           grp=[['model_name', 'model_name']],
                                           seedfiles=[['seed']])
    
    datasource.iterables = [('model_name',modelist),('seed',seedlist)]
    
    return datasource


def extract_subjectID(out_file):
    
    import sys
    
    outs = (out_file.split('_subject_id_'))[1]
    out_file = (outs.split('/'))[0]
    print "extract_subjectID out_file", out_file
    sys.stderr.write('\n>>>>>D<<<<< ' + out_file)

    return out_file


def create_parc_dataflow(unitDefinitionsDirectory):
    
    import nipype.interfaces.io as nio
    import os
    
    unitlist=[os.path.splitext(os.path.splitext(f)[0])[0] for f in os.listdir(unitDefinitionsDirectory)]
    print "unitList ->", unitlist
    datasource = pe.Node(interface=nio.DataGrabber(infields=['parcelation'],
                                                   outfields=['out_file']),
                         name="datasource_parc")
    datasource.inputs.base_directory = unitDefinitionsDirectory
    datasource.inputs.template='*'
    datasource.inputs.field_template=dict(out_file='%s.nii.gz')
    datasource.inputs.template_args=dict(out_file=[['parcelation']])
    datasource.iterables= ('parcelation', unitlist)


    return datasource

def create_mask_dataflow(voxelMasksDirectory):
    
    import nipype.interfaces.io as nio
    import os
    
    masklist=[os.path.splitext(os.path.splitext(f)[0])[0] for f in os.listdir(voxelMasksDirectory)]
    print masklist
    datasource = pe.Node(interface=nio.DataGrabber(infields=['mask'],
                                                   outfields=['out_file']),
                         name="datasource_mask")
    datasource.inputs.base_directory = voxelMasksDirectory
    datasource.inputs.template='*'
    datasource.inputs.field_template=dict(out_file='%s.nii.gz')
    datasource.inputs.template_args=dict(out_file=[['mask']])
    datasource.iterables= ('mask', masklist)

    return datasource

def gen_csv_for_parcelation(data_file, template, csv_file_name, unitTSOutputs):
    
    import nibabel as nib
    import csv
    import numpy as np
    import os

    unit=nib.load(template)
    unit_data=unit.get_data()
    datafile= nib.load(data_file)
    img_data=datafile.get_data()
    vol=img_data.shape[3]
    
    nodes=np.unique(unit_data).tolist()
    sorted_list=[]
    node_dict={}
    out_list=[]
    
    for n in nodes:
        if n>0:
            node_array=img_data[unit_data==n]
            node_array=node_array.T
            time_points, no_of_voxels=node_array.shape
            list1=[n]
            node_str='node %s' %(n)
            node_dict[node_str]=node_array
            for t in range(0,time_points):
                avg=node_array[t].mean()
                avg=np.round(avg,3)
                list1.append(avg)
            sorted_list.append(list1)
    
    sub_file = os.path.splitext(os.path.basename(csv_file_name))[0]
    sub_file= os.path.splitext(sub_file)[0]
    tmp_file=os.path.splitext(os.path.basename(template))[0]
    tmp_file=os.path.splitext(tmp_file)[0]
    csv_file= os.path.abspath(sub_file+'_'+tmp_file+'.csv')
    numpy_file=os.path.abspath(sub_file+'_'+tmp_file+'.npz')
    
    if unitTSOutputs[0]==True:
        f= open(csv_file, 'wt')
        writer = csv.writer(f, delimiter=',', quoting=csv.QUOTE_MINIMAL) 
        headers=['node/volume']
        for i in range (0, vol):
            headers.append(i)   
        writer.writerow(headers)
        writer.writerows(sorted_list)
        f.close()
        out_list.append(csv_file)
   
    if unitTSOutputs[1]==True: 
        np.savez(numpy_file, **dict(node_dict))
        out_list.append(numpy_file)
        
    return out_list

def gen_csv_for_mask(data_file, template, csv_file_name,voxelTSOutputs):
    
    import nibabel as nib
    import numpy as np
    import csv
    import os

    
    unit=nib.load(template)
    unit_data=unit.get_data()
    datafile= nib.load(data_file)
    img_data=datafile.get_data()
    header_data=datafile.get_header()
    qform=header_data.get_qform()
    sorted_list=[]
    vol_dict={}
    out_list=[]
        
    node_array=img_data[unit_data!=0]
    node_array=node_array.T
    time_points=node_array.shape[0]
    for t in range(0,time_points):
        str= 'vol %s' %(t)
        vol_dict[str]=node_array[t]
        val=node_array[t].tolist()
        val.insert(0,t)
        sorted_list.append(val)
        
    sub_file = os.path.splitext(os.path.basename(csv_file_name))[0]
    sub_file= os.path.splitext(sub_file)[0]
    tmp_file=os.path.splitext(os.path.basename(template))[0]
    tmp_file=os.path.splitext(tmp_file)[0]
    csv_file= os.path.abspath(sub_file+'_'+tmp_file+'.csv')
    numpy_file=os.path.abspath(sub_file+'_'+tmp_file+'.npz')
    
    if voxelTSOutputs[0]==True:
        f= open(csv_file, 'wt')
        writer = csv.writer(f, delimiter=',', quoting=csv.QUOTE_MINIMAL)  
        one=np.array([1])
        headers=['volume/xyz']
        cordinates=np.argwhere(unit_data!=0)
        for val in range(np.alen(cordinates)):
            ijk_mat=np.concatenate([cordinates[val],one])
            ijk_mat=ijk_mat.T
            product=np.dot(qform,ijk_mat)
            val=tuple(product.tolist()[0:3])
            headers.append(val)
        writer.writerow(headers)
        writer.writerows(sorted_list)
        f.close()
        out_list.append(csv_file)
    
    if voxelTSOutputs[1]==True:
        np.savez(numpy_file, **dict(vol_dict))
        out_list.append(numpy_file)
    
    
    
    return out_list

def gen_csv_for_surface(rh_surface_file, lh_surface_file, verticesTSOutputs):
    
    import gradunwarp
    import numpy as np
    import os
    out_list=[]
    
    if verticesTSOutputs[0]==True:
        rh_file= os.path.splitext(os.path.basename(rh_surface_file))[0] +'_rh.csv'
        mghobj1 = gradunwarp.mgh.MGH()
        
        mghobj1.load(rh_surface_file)
        vol=mghobj1.vol
        (x,y)=vol.shape
        print "rh shape" ,x,y
        
        np.savetxt(rh_file, vol, delimiter=',')
        out_list.append(rh_file)
        
        lh_file= os.path.splitext(os.path.basename(lh_surface_file))[0] +'_lh.csv'
        mghobj2 = gradunwarp.mgh.MGH()
        
        mghobj2.load(lh_surface_file)
        vol=mghobj2.vol
        (x,y)=vol.shape
        print "lh shape" ,x,y
        
        np.savetxt(lh_file, vol, delimiter=',')
        out_list.append(lh_file)
    
    return out_list
