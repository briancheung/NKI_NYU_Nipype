#!/frodo/shared/epd/bin/python
import e_afni
import sys
import os
import commands
import nipype.pipeline.engine as pe
import nipype.algorithms.rapidart as ra
import nipype.interfaces.afni as afni
import nipype.interfaces.fsl as fsl
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
from nipype.interfaces.utility import Rename
from nipype.interfaces.afni.base import Info, AFNITraitedSpec, AFNICommand
from nipype.interfaces.base import Bunch, TraitedSpec, File, Directory, InputMultiPath
from nipype.utils.filemanip import fname_presuffix, list_to_filename, split_filename
from nipype.utils.docparse import get_doc
from nipype.interfaces.fsl import MultiImageMaths
from nipype.interfaces.fsl import ExtractROI
from nipype.interfaces.fsl.maths import DilateImage
from multiprocessing import Process
from multiprocessing import Pool
from ConfigParser import SafeConfigParser


FWHM = 6.0
sigma = 2.5479870901
#n_lp = float(HP) * float(int(n_vols)) * float (TR)
#n1 = float("%1.0f" %(float(n_lp - 1.0)))
#n_hp = float(LP) * float(int(n_vols)) * float (TR)
#n2 = float("%1.0f" %(float(n_hp - n_lp + 1.0)))
seed_list = []


def setParameters(_fwhm):

    global FWHM
    global sigma
    FWHM = float(_fwhm)
    sigma = FWHM / 2.3548


def setSeedList(seed_file):

    global seed_list

    f = open(seed_file, 'r')
    lines = f.readlines()
    f.close()

    for line in lines:
        line = line.rstrip('\r\n')
        seed_list.append(line)

    print str(seed_list)


def getN1(TR, nvols, HP):

    n_lp = float(HP) * float(int(nvols)) * float(TR)
    n1 = int("%1.0f" % (float(n_lp - 1.0)))

    return n1


def getN2(TR, nvols, LP, HP):

    n_lp = float(HP) * float(int(nvols)) * float(TR)
    n_hp = float(LP) * float(int(nvols)) * float(TR)
    n2 = int("%1.0f" % (float(n_hp - n_lp + 1.0)))

    return n2


def getImgTR(in_files):

    out = []
    from nibabel import load
    if(isinstance(in_files, list)):
        for in_file in in_files:
            img = load(in_file)
            hdr = img.get_header()
            tr = float(hdr.get_zooms()[3])
            if tr > 10:
                tr = float(float(tr) / 1000.0)
            out.append(tr)
        return out
    else:
        img = load(in_files)
        hdr = img.get_header()
        tr = float(hdr.get_zooms()[3])
        if tr > 10:
            tr = float(float(tr) / 1000.0)
        return [tr]


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
        cmd = "awk '{print $1}' %s > %s/mc1.1D" % (file, parent_path)
        print cmd
        sys.stderr.write(commands.getoutput(cmd))

        cmd = "awk '{print $2}' %s > %s/mc2.1D" % (file, parent_path)
        print cmd
        sys.stderr.write(commands.getoutput(cmd))

        cmd = "awk '{print $3}' %s > %s/mc3.1D" % (file, parent_path)
        print cmd
        sys.stderr.write(commands.getoutput(cmd))

        cmd = "awk '{print $4}' %s > %s/mc4.1D" % (file, parent_path)
        print cmd
        sys.stderr.write(commands.getoutput(cmd))

        cmd = "awk '{print $5}' %s > %s/mc5.1D" % (file, parent_path)
        print cmd
        sys.stderr.write(commands.getoutput(cmd))

        cmd = "awk '{print $6}' %s > %s/mc6.1D" % (file, parent_path)
        print cmd
        sys.stderr.write(commands.getoutput(cmd))

        cmd = "awk '{print $5}' %s > %s/mc5.1D" % (file, parent_path)
        print cmd
        sys.stderr.write(commands.getoutput(cmd))

        cmd = "awk '{print $6}' %s > %s/mc6.1D" % (file, parent_path)
        print cmd
        sys.stderr.write(commands.getoutput(cmd))

        EV_List.append(os.path.abspath(parent_path + '/mc1.1D'))
        EV_List.append(os.path.abspath(parent_path + '/mc2.1D'))
        EV_List.append(os.path.abspath(parent_path + '/mc3.1D'))
        EV_List.append(os.path.abspath(parent_path + '/mc4.1D'))
        EV_List.append(os.path.abspath(parent_path + '/mc5.1D'))
        EV_List.append(os.path.abspath(parent_path + '/mc6.1D'))

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

        cmd = "sed -e s:nuisance_dir:\"%s\":g <%s >%s/temp1" % (parent_path, nuisance_template, parent_path)
        print cmd
        sys.stderr.write(commands.getoutput(cmd))

        cmd = "sed -e s:nuisance_model_outputdir:\"%s/residuals.feat\":g <%s/temp1 >%s/temp2" % (parent_path, parent_path, parent_path)
        print cmd
        sys.stderr.write(commands.getoutput(cmd))

        cmd = "sed -e s:nuisance_model_TR:\"%s\":g <%s/temp2 >%s/temp3" % (TR[idx], parent_path, parent_path)
        print cmd
        sys.stderr.write(commands.getoutput(cmd))

        cmd = "sed -e s:nuisance_model_numTRs:\"%s\":g <%s/temp3 >%s/temp4" % (n_vols[idx], parent_path, parent_path)
        print cmd
        sys.stderr.write(commands.getoutput(cmd))

        cmd = "sed -e s:nuisance_model_input_data:\"%s\":g <%s/temp4 >%s/nuisance.fsf" % (rest_p, parent_path, parent_path)
        print cmd
        sys.stderr.write(commands.getoutput(cmd))

    #    cmd = "rm %s/temp*" %(parent_path)
    #    print cmd
    #    sys.stderr.write(commands.getoutput(cmd))
        nuisance_files.append(os.path.abspath("%s/nuisance.fsf" % (parent_path)))

        idx += 1

    return nuisance_files


def get_EV_LIST(EV_Lists, global1Ds, csf1Ds, wm1Ds):

    idx = 0
    EV_final_lists_set = []
    for EV_List in EV_Lists:
        print '>>>>><<<<', EV_List
        EV_List.append(global1Ds[idx])
        EV_List.append(csf1Ds[idx])
        EV_List.append(wm1Ds[idx])
        EV_final_lists_set.append(EV_List)

        idx += 1

    return EV_final_lists_set


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


def getOpString(mean, std_dev):

    str1 = "-sub %f -div %f" % (float(mean), float(std_dev))

    op_string = str1 + " -mas %s"

    return op_string


def pToFile(time_series):

    import os
    import re
    import commands
    import sys

#    ts_oneD =   'abcd.1D'

    dir = os.path.dirname(time_series)

    dir1 = re.sub(r'(.)+_subject_id_(.)+/_mask', '', dir)

    dir1 = dir1.split('..')
    dir1 = dir1[len(dir1) - 1]
    dir1 = dir1.split('.nii.gz')
    dir1 = dir1[0]

    ts_oneD = dir + '/' + dir1 + '.1D'
    cmd = "mv %s %s" % (time_series, ts_oneD)
    print cmd

    sys.stderr.write('\n' + commands.getoutput(cmd))
    return os.path.abspath(ts_oneD)


RSFC_printToFile = pe.MapNode(util.Function(input_names=['time_series'], output_names=['ts_oneD'], function=pToFile), name='RSFC_printToFile', iterfield=["time_series"])


def mPath(rest_path, strategy, sublist):

    import re

    def modifyPath(path, strategy, sublist):

        rpath = ''
        parts = path.split('-')

        idx = 0

        if strategy in path:
            return path

        for part in parts:

            if not part in sublist:
                rpath += '-' + part
            else:
                rpath += '-' + strategy + '-' + part
            idx += 1

        return rpath

    result_path = []

    if(isinstance(rest_path, list)):
        for r_p in rest_path:
            path = re.sub(r'\/', '-', r_p)
            path = modifyPath(path, strategy, sublist)
            result_path.append(path)

        print str(result_path)
    else:
        path = re.sub(r'\/', '-', rest_path)
        path = modifyPath(path, strategy, sublist)
        result_path.append(path)
        print str(result_path)
    return result_path


def mAnat(anat, strategy, sublist):

    import os
    parts = anat.split('/')

    result_path = ''
    idx = 0

    if strategy in anat:
        result_path = os.path.dirname(anat)
        return result_path

    for part in parts:

        if not part.endswith('.nii.gz'):
            if not part in sublist:
                result_path += '/' + part
            else:
                result_path += '/' + strategy + '/' + part
        idx += 1

    return result_path

#TR and n_vols nodes

adjustPath = pe.Node(util.Function(input_names=['rest_path', 'strategy', 'sublist'], output_names=['result_path'], function=mPath), name='adjustPath')

adjustPath1 = adjustPath.clone('adjustPath1')

adjustAnat = pe.Node(util.Function(input_names=['anat', 'strategy', 'sublist'], output_names=['result_path'], function=mAnat), name='adjustAnat')


TR = pe.Node(util.Function(input_names=['in_files'], output_names=['TR'], function=getImgTR), name='TR')
NVOLS = pe.Node(util.Function(input_names=['in_files'], output_names=['nvols'], function=getImgNVols), name='NVOLS')

calcN1 = pe.MapNode(util.Function(input_names=['nvols', 'TR', 'HP'], output_names=['n1'], function=getN1), name='calcN1', iterfield=["nvols", "TR"])

calcN2 = pe.MapNode(util.Function(input_names=['nvols', 'TR', 'LP', 'HP'], output_names=['n2'], function=getN2), name='calcN2', iterfield=["nvols", "TR"])

alff_op_string = pe.MapNode(util.Function(input_names=['mean', 'std_dev'], output_names=['op_string'], function=getOpString), name='alff_op_string', iterfield=["mean", "std_dev"])

alff_op_string1 = pe.MapNode(util.Function(input_names=['mean', 'std_dev'], output_names=['op_string'], function=getOpString), name='alff_op_string1', iterfield=["mean", "std_dev"])


#create Motion Correction ID files for all six motion parameters
MC = pe.Node(util.Function(input_names=['in_file', 'pp'], output_names=['EV_Lists'], function=createMC), name='MC')


FSF = pe.Node(util.Function(input_names=['nuisance_template', 'rest_pp', 'TR', 'n_vols'], output_names=['nuisance_files'], function=createFSF), name='FSF')

copyS = pe.Node(util.Function(input_names=['EV_Lists', 'nuisance_files', 'global1Ds', 'csf1Ds', 'wm1Ds', 'regressors'], output_names=['EV_final_lists_set'], function=copyStuff), name='copyS')

EV_List = pe.Node(util.Function(input_names=['EV_Lists', 'global1Ds', 'csf1Ds', 'wm1Ds'], output_names=['EV_final_lists_set'], function=get_EV_LIST), name='EV_List')

getStats = pe.Node(util.Function(input_names=['in_files'], output_names=['stats_dir'], function=getStatsDir), name='getStats')


lifeSaver = pe.MapNode(interface=util.Rename(), name='lifeSaver', iterfield=["in_file", "format_string"])

lifeSaverALFF = pe.MapNode(interface=util.Rename(), name='lifeSaverALFF', iterfield=["in_file", "format_string"])

lifeSaverNuisancepp = pe.MapNode(interface=util.Rename(), name='lifeSaverNuisancepp', iterfield=["in_file", "format_string"])

lifeSaverNuisanceMask = pe.MapNode(interface=util.Rename(), name='lifeSaverNuisanceMask', iterfield=["in_file", "format_string"])

lifeSaver_reg_flirt = pe.MapNode(interface=util.Rename(), name='lifeSaver_reg_flirt', iterfield=["in_file", "format_string"])

lifeSaver_func2highresmat = pe.MapNode(interface=util.Rename(), name='lifeSaver_func2highresmat', iterfield=["in_file", "format_string"])

lifeSaver_highres2example_funcmat = pe.MapNode(interface=util.Rename(), name='lifeSaver_highres2example_funcmat', iterfield=["in_file", "format_string"])

lifeSaver_func2standardmat = pe.MapNode(interface=util.Rename(), name='lifeSaver_func2standardmat', iterfield=["in_file", "format_string"])

lifeSaver_func2standard = pe.MapNode(interface=util.Rename(), name='lifeSaver_func2standard', iterfield=["in_file", "format_string"])

lifeSaver_standard2example_funcmat = pe.MapNode(interface=util.Rename(), name='lifeSaver_standard2example_funcmat', iterfield=["in_file", "format_string"])

lifeSaver_example_func2standard_NL = pe.MapNode(interface=util.Rename(), name='lifeSaver_example_func2standard_NL', iterfield=["in_file", "format_string"])

lifeSaver_seg_flirt = pe.MapNode(interface=util.Rename(), name='lifeSaver_seg_flirt', iterfield=["in_file", "format_string"])

lifeSaver_seg_warp = pe.MapNode(interface=util.Rename(), name='lifeSaver_seg_warp', iterfield=["in_file", "format_string"])


lifeSaver_seg_smooth1 = pe.MapNode(interface=util.Rename(), name='lifeSaver_seg_smooth1', iterfield=["in_file", "format_string"])

lifeSaver_seg_thresh = pe.MapNode(interface=util.Rename(), name='lifeSaver_seg_thresh', iterfield=["in_file", "format_string"])

lifeSaver_seg_mask = pe.MapNode(interface=util.Rename(), name='lifeSaver_seg_mask', iterfield=["in_file", "format_string"])

lifeSaver_seg_copy = pe.MapNode(interface=util.Rename(), name='lifeSaver_seg_copy', iterfield=["in_file", "format_string"])

lifeSaver_seg_flirt3 = pe.MapNode(interface=util.Rename(), name='lifeSaver_seg_flirt3', iterfield=["in_file", "format_string"])

lifeSaver_seg_warp1 = pe.MapNode(interface=util.Rename(), name='lifeSaver_seg_warp1', iterfield=["in_file", "format_string"])

lifeSaver_seg_prior1 = pe.MapNode(interface=util.Rename(), name='lifeSaver_seg_prior1', iterfield=["in_file", "format_string"])

lifeSaver_seg_thresh1 = pe.MapNode(interface=util.Rename(), name='lifeSaver_seg_thresh1', iterfield=["in_file", "format_string"])

lifeSaver_seg_mask1 = pe.MapNode(interface=util.Rename(), name='lifeSaver_seg_mask1', iterfield=["in_file", "format_string"])


lifeSaver_alff_cp = pe.MapNode(interface=util.Rename(), name='lifeSaver_alff_cp', iterfield=["in_file", "format_string"])

lifeSaver_alff_mean = pe.MapNode(interface=util.Rename(), name='lifeSaver_alff_mean', iterfield=["in_file", "format_string"])

lifeSaver_alff_pspec = pe.MapNode(interface=util.Rename(), name='lifeSaver_alff_pspec', iterfield=["in_file", "format_string"])

lifeSaver_alff_sum = pe.MapNode(interface=util.Rename(), name='lifeSaver_alff_sum', iterfield=["in_file", "format_string"])

lifeSaver_alff_falff1 = pe.MapNode(interface=util.Rename(), name='lifeSaver_alff_falff1', iterfield=["in_file", "format_string"])

lifeSaver_alff_Z_falff = pe.MapNode(interface=util.Rename(), name='lifeSaver_alff_Z_falff', iterfield=["in_file", "format_string"])

lifeSaver_alff_Z_alff = pe.MapNode(interface=util.Rename(), name='lifeSaver_alff_Z_alff', iterfield=["in_file", "format_string"])

lifeSaver_alff_warp_alff = pe.MapNode(interface=util.Rename(), name='lifeSaver_alff_warp_alff', iterfield=["in_file", "format_string"])

lifeSaver_alff_warp_falff = pe.MapNode(interface=util.Rename(), name='lifeSaver_alff_warp_falff', iterfield=["in_file", "format_string"])

lifeSaver_alff_smooth = pe.MapNode(interface=util.Rename(), name='lifeSaver_alff_smooth', iterfield=["in_file", "format_string"])

lifeSaver_falff_smooth = pe.MapNode(interface=util.Rename(), name='lifeSaver_falff_smooth', iterfield=["in_file", "format_string"])


def saveSeg(in_file, pp):

    import os
    out_file = ""

    fname = os.path.basename(in_file)

    dir = (in_file.split(fname))[0]

    split_path = pp.split('-')

    for index in range(0, len(split_path) - 1):

        if index == 0:

            out_file += dir + '-'
        else:
            out_file += split_path[index] + '-'

    out_file += 'segment' + '-' + fname

    return out_file


def saveNuisance(in_file, pp):

    import os
    out_file = ""

    fname = os.path.basename(in_file)

    dir = (in_file.split(fname))[0]

    split_path = pp.split('-')

    for index in range(0, len(split_path) - 1):

        if index == 0:

            out_file += dir + '-'
        else:
            out_file += split_path[index] + '-'

    out_file += 'nuisance' + '-' + fname

    return out_file


def saveRSFC_seed_ts(in_file, pp):

    import os
    out_file = ""

    fname = os.path.basename(in_file)

    dir = (in_file.split(fname))[0]

    split_path = pp.split('-')

    for index in range(0, len(split_path) - 1):

        if index == 0:

            out_file += dir + '-'
        else:
            out_file += split_path[index] + '-'

    out_file += 'seed_ts' + '-' + fname

    return out_file


def saveRSFC(in_file, pp):

    import os
    out_file = ""

    fname = os.path.basename(in_file)

    dir = (in_file.split(fname))[0]

    split_path = pp.split('-')

    for index in range(0, len(split_path) - 1):

        if index == 0:

            out_file += dir + '-'
        else:
            out_file += split_path[index] + '-'

    out_file += 'RSFC' + '-' + fname

    return out_file


def saveNuisance_fgls(in_file, pp):

    import os
    out_file = ""

    fname = os.path.basename(in_file)

    dir = (in_file.split(fname))[0]

    split_path = pp.split('-')

    for index in range(0, len(split_path) - 1):

        if index == 0:

            out_file += dir + '-'
        else:
            out_file += split_path[index] + '-'

    out_file += 'nuisance' + '-stats-' + fname

    return out_file


def saveNuisance_calc(in_file, pp):

    import os
    out_file = ""

    fname = os.path.basename(in_file)

    dir = (in_file.split(fname))[0]

    split_path = pp.split('-')

    for index in range(0, len(split_path) - 1):

        if index == 0:

            out_file += dir + '-'
        else:
            out_file += split_path[index] + '-'

    out_file += fname

    return out_file

Saver_nuisance_featM = pe.MapNode(util.Function(input_names=['in_file', 'pp'], output_names=['out_file'], function=saveNuisance), name='Saver_nuisance_featM', iterfield=["in_file", "pp"])

Saver_nuisance_erosion_csf1 = pe.MapNode(util.Function(input_names=['in_file', 'pp'], output_names=['out_file'], function=saveNuisance), name='Saver_nuisance_erosion_csf1', iterfield=["in_file", "pp"])

Saver_nuisance_erosion_wm1 = pe.MapNode(util.Function(input_names=['in_file', 'pp'], output_names=['out_file'], function=saveNuisance), name='Saver_nuisance_erosion_wm1', iterfield=["in_file", "pp"])

Saver_nuisance_fgls = pe.MapNode(util.Function(input_names=['in_file', 'pp'], output_names=['out_file'], function=saveNuisance_fgls), name='Saver_nuisance_fgls', iterfield=["in_file", "pp"])

Saver_nuisance_fglsd = pe.MapNode(util.Function(input_names=['in_file', 'pp'], output_names=['out_file'], function=saveNuisance_fgls), name='Saver_nuisance_fglsd', iterfield=["in_file", "pp"])

Saver_nuisance_stat = pe.MapNode(util.Function(input_names=['in_file', 'pp'], output_names=['out_file'], function=saveNuisance_fgls), name='Saver_nuisance_stat', iterfield=["in_file", "pp"])

Saver_nuisance_calc = pe.MapNode(util.Function(input_names=['in_file', 'pp'], output_names=['out_file'], function=saveNuisance_calc), name='Saver_nuisance_calc', iterfield=["in_file", "pp"])

Saver_func_filter = pe.MapNode(util.Function(input_names=['in_file', 'pp'], output_names=['out_file'], function=saveNuisance_calc), name='Saver_func_filter', iterfield=["in_file", "pp"])

Saver_nuisance_warp = pe.MapNode(util.Function(input_names=['in_file', 'pp'], output_names=['out_file'], function=saveNuisance_calc), name='Saver_nuisance_warp', iterfield=["in_file", "pp"])

Saver_nuisance_warp_1 = pe.MapNode(util.Function(input_names=['in_file', 'pp'], output_names=['out_file'], function=saveNuisance_calc), name='Saver_nuisance_warp_1', iterfield=["in_file", "pp"])

Saver_seg_flirt = pe.MapNode(util.Function(input_names=['in_file', 'pp'], output_names=['out_file'], function=saveSeg), name='Saver_seg_flirt', iterfield=["in_file", "pp"])

Saver_seg_warp = pe.MapNode(util.Function(input_names=['in_file', 'pp'], output_names=['out_file'], function=saveSeg), name='Saver_seg_warp', iterfield=["in_file", "pp"])

Saver_seg_smooth1 = pe.MapNode(util.Function(input_names=['in_file', 'pp'], output_names=['out_file'], function=saveSeg), name='Saver_seg_smooth1', iterfield=["in_file", "pp"])

Saver_seg_thresh = pe.MapNode(util.Function(input_names=['in_file', 'pp'], output_names=['out_file'], function=saveSeg), name='Saver_seg_thresh', iterfield=["in_file", "pp"])

Saver_seg_mask = pe.MapNode(util.Function(input_names=['in_file', 'pp'], output_names=['out_file'], function=saveSeg), name='Saver_seg_mask', iterfield=["in_file", "pp"])

Saver_seg_flirt3 = pe.MapNode(util.Function(input_names=['in_file', 'pp'], output_names=['out_file'], function=saveSeg), name='Saver_seg_flirt3', iterfield=["in_file", "pp"])

Saver_seg_warp1 = pe.MapNode(util.Function(input_names=['in_file', 'pp'], output_names=['out_file'], function=saveSeg), name='Saver_seg_warp1', iterfield=["in_file", "pp"])

Saver_seg_prior1 = pe.MapNode(util.Function(input_names=['in_file', 'pp'], output_names=['out_file'], function=saveSeg), name='Saver_seg_prior1', iterfield=["in_file", "pp"])

Saver_seg_thresh1 = pe.MapNode(util.Function(input_names=['in_file', 'pp'], output_names=['out_file'], function=saveSeg), name='Saver_seg_thresh1', iterfield=["in_file", "pp"])

Saver_seg_mask1 = pe.MapNode(util.Function(input_names=['in_file', 'pp'], output_names=['out_file'], function=saveSeg), name='Saver_seg_mask1', iterfield=["in_file", "pp"])

Saver_RSFC_printToFile = pe.MapNode(util.Function(input_names=['in_file', 'pp'], output_names=['out_file'], function=saveRSFC_seed_ts), name='Saver_RSFC_printToFile', iterfield=["in_file", "pp"])

Saver_RSFC_corr = pe.MapNode(util.Function(input_names=['in_file', 'pp'], output_names=['out_file'], function=saveRSFC), name='Saver_RSFC_corr', iterfield=["in_file", "pp"])

Saver_RSFC_z_trans = pe.MapNode(util.Function(input_names=['in_file', 'pp'], output_names=['out_file'], function=saveRSFC), name='Saver_RSFC_z_trans', iterfield=["in_file", "pp"])

Saver_RSFC_register = pe.MapNode(util.Function(input_names=['in_file', 'pp'], output_names=['out_file'], function=saveRSFC), name='Saver_RSFC_register', iterfield=["in_file", "pp"])

Saver_RSFC_smooth = pe.MapNode(util.Function(input_names=['in_file', 'pp'], output_names=['out_file'], function=saveRSFC), name='Saver_RSFC_smooth', iterfield=["in_file", "pp"])
#anatomical nodes


anat_refit = pe.Node(interface=e_afni.Threedrefit(), name='anat_refit')
anat_refit.inputs.deoblique = True

anat_reorient = pe.Node(interface=e_afni.Threedresample(), name='anat_reorient')
anat_reorient.inputs.orientation = 'RPI'

anat_skullstrip = pe.Node(interface=e_afni.ThreedSkullStrip(), name='anat_skullstrip')
anat_skullstrip.inputs.options = '-o_ply'

anat_calc = pe.Node(interface=e_afni.Threedcalc(), name='anat_calc')
anat_calc.inputs.expr = '\'a*step(b)\''
anat_calc.inputs.out_file = 'mprage_brain.nii.gz'


def rename(in_file, name):

    out_file = ""
    print 'IN_FILE>>>>> ' + in_file
    print 'name>>>>> ' + name
    file_path = in_file.split('-')

    for index in range(0, len(file_path) - 1):

        out_file += file_path[index] + '-'

    out_file += name

    return out_file


def renameRSFC(in_file, name, whichNode):
    import sys

    out_file = ""
    sys.stderr.write('IN_FILE>>>>> ' + str(in_file))
    sys.stderr.write('name>>>>> ' + str(name))
    file_path = in_file.split('/')
    nameL = name.split('/')
    name = nameL[len(nameL) - 1]
    if(whichNode == 'corrs'):
        name = (name.split('.'))[0] + '_corr.nii.gz'

    elif (whichNode == 'z_trans'):
        name = (name.split('.'))[0] + '_Z.nii.gz'

    elif (whichNode == 'register'):
        name = (name.split('.'))[0] + '_Z_2standard.nii.gz'

    elif('Z_2standard_FWHM' in whichNode):

        name = (name.split('.'))[0] + whichNode + '.nii.gz'

    sys.stderr.write('name>>>>> ' + str(name))

    for index in range(0, len(file_path) - 1):

        out_file += file_path[index] + '/'

    out_file += name

    return out_file


def renameAnat(in_file, name):

    out_file = ""
    print 'IN_FILE>>>>> ' + in_file
    print 'name>>>>> ' + name
    file_path = in_file.split('/')

    for index in range(0, len(file_path) - 1):

        out_file += file_path[index] + '/'

    out_file += name

    return out_file


#anatpreproc Nodes

anat_refit_r = pe.Node(interface=util.Rename(), name='anat_refit_r')
anat_refit_o = pe.Node(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=renameAnat), name='anat_refit_o')

anat_reorient_r = pe.Node(interface=util.Rename(), name='anat_reorient_r')
anat_reorient_o = pe.Node(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=renameAnat), name='anat_reorient_o')

anat_reorient_r = pe.Node(interface=util.Rename(), name='anat_reorient_r')
anat_reorient_o = pe.Node(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=renameAnat), name='anat_reorient_o')

anat_skullstrip_r = pe.Node(interface=util.Rename(), name='anat_skullstrip_r')
anat_skullstrip_o = pe.Node(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=renameAnat), name='anat_skullstrip_o')

anat_calc_r = pe.Node(interface=util.Rename(), name='anat_calc_r')
anat_calc_o = pe.Node(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=renameAnat), name='anat_calc_o')

reg_flirt1_r = pe.Node(interface=util.Rename(), name='reg_flirt1_r')
reg_flirt1_o = pe.Node(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=renameAnat), name='reg_flirt1_o')

reg_flirt1o_r = pe.Node(interface=util.Rename(), name='reg_flirt1o_r')
reg_flirt1o_o = pe.Node(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=renameAnat), name='reg_flirt1o_o')

reg_xfm2_r = pe.Node(interface=util.Rename(), name='reg_xfm2_r')
reg_xfm2_o = pe.Node(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=renameAnat), name='reg_xfm2_o')

reg_fnt_r = pe.Node(interface=util.Rename(), name='reg_fnt_r')
reg_fnt_o = pe.Node(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=renameAnat), name='reg_fnt_o')
reg_inw_r = reg_fnt_r.clone('reg_inw_r')
reg_inw_o = reg_fnt_o.clone('reg_inw_o')

reg_fntj_r = pe.Node(interface=util.Rename(), name='reg_fntj_r')
reg_fntj_o = pe.Node(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=renameAnat), name='reg_fntj_o')

reg_fntf_r = pe.Node(interface=util.Rename(), name='reg_fntf_r')
reg_fntf_o = pe.Node(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=renameAnat), name='reg_fntf_o')


#funcpreproc Nodes
func_calc_r = pe.MapNode(interface=util.Rename(), name='func_calc_r', iterfield=["in_file", "format_string"])
func_calc_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='func_calc_o', iterfield=["in_file"])

func_refit_r = pe.MapNode(interface=util.Rename(), name='func_refit_r', iterfield=["in_file", "format_string"])
func_refit_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='func_refit_o', iterfield=["in_file"])

func_reorient_r = pe.MapNode(interface=util.Rename(), name='func_reorient_r', iterfield=["in_file", "format_string"])
func_reorient_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='func_reorient_o', iterfield=["in_file"])

func_tstat_r = pe.MapNode(interface=util.Rename(), name='func_tstat_r', iterfield=["in_file", "format_string"])
func_tstat_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='func_tstat_o', iterfield=["in_file"])

func_volrego_r = pe.MapNode(interface=util.Rename(), name='func_volrego_r', iterfield=["in_file", "format_string"])
func_volrego_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='func_volrego_o', iterfield=["in_file"])

func_volreg_r = pe.MapNode(interface=util.Rename(), name='func_volreg_r', iterfield=["in_file", "format_string"])
func_volreg_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='func_volreg_o', iterfield=["in_file"])

func_automask_r = pe.MapNode(interface=util.Rename(), name='func_automask_r', iterfield=["in_file", "format_string"])
func_automask_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='func_automask_o', iterfield=["in_file"])

func_calcR_r = pe.MapNode(interface=util.Rename(), name='func_calcR_r', iterfield=["in_file", "format_string"])
func_calcR_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='func_calcR_o', iterfield=["in_file"])

func_calcI_r = pe.MapNode(interface=util.Rename(), name='func_calcI_r', iterfield=["in_file", "format_string"])
func_calcI_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='func_calcI_o', iterfield=["in_file"])

func_despike_r = pe.MapNode(interface=util.Rename(), name='func_despike_r', iterfield=["in_file", "format_string"])
func_despike_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='func_despike_o', iterfield=["in_file"])

func_mean_r = pe.MapNode(interface=util.Rename(), name='func_mean_r', iterfield=["in_file", "format_string"])
func_mean_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='func_mean_o', iterfield=["in_file"])


func_scale_r = pe.MapNode(interface=util.Rename(), name='func_scale_r', iterfield=["in_file", "format_string"])
func_scale_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='func_scale_o', iterfield=["in_file"])

func_filter_r = pe.MapNode(interface=util.Rename(), name='func_filter_r', iterfield=["in_file", "format_string"])
func_filter_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='func_filter_o', iterfield=["in_file"])

func_detrenda_r = pe.MapNode(interface=util.Rename(), name='func_detrenda_r', iterfield=["in_file", "format_string"])
func_detrenda_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='func_detrenda_o', iterfield=["in_file"])

func_detrendb_r = pe.MapNode(interface=util.Rename(), name='func_detrendb_r', iterfield=["in_file", "format_string"])
func_detrendb_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='func_detrendb_o', iterfield=["in_file"])

func_detrendc_r = pe.MapNode(interface=util.Rename(), name='func_detrendc_r', iterfield=["in_file", "format_string"])
func_detrendc_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='func_detrendc_o', iterfield=["in_file"])

func_mask_r = pe.MapNode(interface=util.Rename(), name='func_mask_r', iterfield=["in_file", "format_string"])
func_mask_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='func_mask_o', iterfield=["in_file"])

reg_flirt_r = pe.MapNode(interface=util.Rename(), name='reg_flirt_r', iterfield=["in_file", "format_string"])
reg_flirt_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='reg_flirt_o', iterfield=["in_file"])

reg_flirto_r = pe.MapNode(interface=util.Rename(), name='reg_flirto_r', iterfield=["in_file", "format_string"])
reg_flirto_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='reg_flirto_o', iterfield=["in_file"])

reg_xfm1_r = pe.MapNode(interface=util.Rename(), name='reg_xfm1_r', iterfield=["in_file", "format_string"])
reg_xfm1_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='reg_xfm1_o', iterfield=["in_file"])

reg_xfm3_r = pe.MapNode(interface=util.Rename(), name='reg_xfm3_r', iterfield=["in_file", "format_string"])
reg_xfm3_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='reg_xfm3_o', iterfield=["in_file"])

reg_flirt2_r = pe.MapNode(interface=util.Rename(), name='reg_flirt2_r', iterfield=["in_file", "format_string"])
reg_flirt2_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='reg_flirt2_o', iterfield=["in_file"])

reg_xfm4_r = pe.MapNode(interface=util.Rename(), name='reg_xfm4_r', iterfield=["in_file", "format_string"])
reg_xfm4_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='reg_xfm4_o', iterfield=["in_file"])

reg_warp_r = pe.MapNode(interface=util.Rename(), name='reg_warp_r', iterfield=["in_file", "format_string"])
reg_warp_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='reg_warp_o', iterfield=["in_file"])

seg_flirt_r = pe.MapNode(interface=util.Rename(), name='seg_flirt_r', iterfield=["in_file", "format_string"])
seg_flirt_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='seg_flirt_o', iterfield=["in_file"])

seg_warp_r = pe.MapNode(interface=util.Rename(), name='seg_warp_r', iterfield=["in_file", "format_string"])
seg_warp_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='seg_warp_o', iterfield=["in_file"])

seg_smooth1_r = pe.MapNode(interface=util.Rename(), name='seg_smooth1_r', iterfield=["in_file", "format_string"])
seg_smooth1_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='seg_smooth1_o', iterfield=["in_file"])

seg_thresh_r = pe.MapNode(interface=util.Rename(), name='seg_thresh_r', iterfield=["in_file", "format_string"])
seg_thresh_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='seg_thresh_o', iterfield=["in_file"])

seg_mask_r = pe.MapNode(interface=util.Rename(), name='seg_mask_r', iterfield=["in_file", "format_string"])
seg_mask_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='seg_mask_o', iterfield=["in_file"])

seg_copy_r = pe.MapNode(interface=util.Rename(), name='seg_copy_r', iterfield=["in_file", "format_string"])
seg_copy_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='seg_copy_o', iterfield=["in_file"])

seg_flirt3_r = pe.MapNode(interface=util.Rename(), name='seg_flirt3_r', iterfield=["in_file", "format_string"])
seg_flirt3_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='seg_flirt3_o', iterfield=["in_file"])

seg_warp1_r = pe.MapNode(interface=util.Rename(), name='seg_warp1_r', iterfield=["in_file", "format_string"])
seg_warp1_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='seg_warp1_o', iterfield=["in_file"])

seg_prior1_r = pe.MapNode(interface=util.Rename(), name='seg_prior1_r', iterfield=["in_file", "format_string"])
seg_prior1_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='seg_prior1_o', iterfield=["in_file"])

seg_thresh1_r = pe.MapNode(interface=util.Rename(), name='seg_thresh1_r', iterfield=["in_file", "format_string"])
seg_thresh1_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='seg_thresh1_o', iterfield=["in_file"])

seg_mask1_r = pe.MapNode(interface=util.Rename(), name='seg_mask1_r', iterfield=["in_file", "format_string"])
seg_mask1_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='seg_mask1_o', iterfield=["in_file"])

nuisance_globalE_r = pe.MapNode(interface=util.Rename(), name='nuisance_globalE_r', iterfield=["in_file", "format_string"])
nuisance_globalE_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='nuisance_globalE_o', iterfield=["in_file"])

nuisance_erosion_csf1_r = pe.MapNode(interface=util.Rename(), name='nuisance_erosion_csf1_r', iterfield=["in_file", "format_string"])
nuisance_erosion_csf1_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='nuisance_erosion_csf1_o', iterfield=["in_file"])

nuisance_erosion_wm1_r = pe.MapNode(interface=util.Rename(), name='nuisance_erosion_wm1_r', iterfield=["in_file", "format_string"])
nuisance_erosion_wm1_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='nuisance_erosion_wm1_o', iterfield=["in_file"])

nuisance_compcor_r = pe.MapNode(interface=util.Rename(), name='nuisance_compcor_r', iterfield=["in_file", "format_string"])
nuisance_compcor_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='nuisance_compcor_o', iterfield=["in_file"])

nuisance_MedianAngle_r = pe.MapNode(interface=util.Rename(), name='nuisance_MedianAngle_r', iterfield=["in_file", "format_string"])
nuisance_MedianAngle_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='nuisance_MedianAngle_o', iterfield=["in_file"])

nuisance_csf_r = pe.MapNode(interface=util.Rename(), name='nuisance_csf_r', iterfield=["in_file", "format_string"])
nuisance_csf_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='nuisance_csf_o', iterfield=["in_file"])

nuisance_wm_r = pe.MapNode(interface=util.Rename(), name='nuisance_wm_r', iterfield=["in_file", "format_string"])
nuisance_wm_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='nuisance_wm_o', iterfield=["in_file"])

nuisance_featM_r = pe.MapNode(interface=util.Rename(), name='nuisance_featM_r', iterfield=["in_file", "format_string"])
nuisance_featM_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='nuisance_featM_o', iterfield=["in_file"])

nuisance_fgls_r = pe.MapNode(interface=util.Rename(), name='nuisance_fgls_r', iterfield=["in_file", "format_string"])
nuisance_fgls_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='nuisance_fgls_o', iterfield=["in_file"])

nuisance_stat_r = pe.MapNode(interface=util.Rename(), name='nuisance_stat_r', iterfield=["in_file", "format_string"])
nuisance_stat_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='nuisance_stat_o', iterfield=["in_file"])

nuisance_calc_r = pe.MapNode(interface=util.Rename(), name='nuisance_calc_r', iterfield=["in_file", "format_string"])
nuisance_calc_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='nuisance_calc_o', iterfield=["in_file"])

nuisance_warp_r = pe.MapNode(interface=util.Rename(), name='nuisance_warp_r', iterfield=["in_file", "format_string"])
nuisance_warp_1_r = pe.MapNode(interface=util.Rename(), name='nuisance_warp_1_r', iterfield=["in_file", "format_string"])
nuisance_warp_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='nuisance_warp_o', iterfield=["in_file"])

nuisance_warp_1_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='nuisance_warp_1_o', iterfield=["in_file"])

alff_cp_r = pe.MapNode(interface=util.Rename(), name='alff_cp_r', iterfield=["in_file", "format_string"])
alff_cp_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='alff_cp_o', iterfield=["in_file"])

alff_mean_r = pe.MapNode(interface=util.Rename(), name='alff_mean_r', iterfield=["in_file", "format_string"])
alff_mean_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='alff_mean_o', iterfield=["in_file"])

alff_pspec_r = pe.MapNode(interface=util.Rename(), name='alff_pspec_r', iterfield=["in_file", "format_string"])
alff_pspec_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='alff_pspec_o', iterfield=["in_file"])

alff_sum_r = pe.MapNode(interface=util.Rename(), name='alff_sum_r', iterfield=["in_file", "format_string"])
alff_sum_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='alff_sum_o', iterfield=["in_file"])

alff_falff1_r = pe.MapNode(interface=util.Rename(), name='alff_falff1_r', iterfield=["in_file", "format_string"])
alff_falff1_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='alff_falff1_o', iterfield=["in_file"])

alff_Z_falff_r = pe.MapNode(interface=util.Rename(), name='alff_Z_falff_r', iterfield=["in_file", "format_string"])
alff_Z_falff_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='alff_Z_falff_o', iterfield=["in_file"])

alff_Z_alff_r = pe.MapNode(interface=util.Rename(), name='alff_Z_alff_r', iterfield=["in_file", "format_string"])
alff_Z_alff_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='alff_Z_alff_o', iterfield=["in_file"])

alff_warp_alff_r = pe.MapNode(interface=util.Rename(), name='alff_warp_alff_r', iterfield=["in_file", "format_string"])
alff_warp_alff_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='alff_warp_alff_o', iterfield=["in_file"])

alff_warp_falff_r = pe.MapNode(interface=util.Rename(), name='alff_warp_falff_r', iterfield=["in_file", "format_string"])
alff_warp_falff_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='alff_warp_falff_o', iterfield=["in_file"])

alff_smooth_r = pe.MapNode(interface=util.Rename(), name='alff_smooth_r', iterfield=["in_file", "format_string"])
alff_smooth_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='alff_smooth_o', iterfield=["in_file"])

falff_smooth_r = pe.MapNode(interface=util.Rename(), name='falff_smooth_r', iterfield=["in_file", "format_string"])
falff_smooth_o = pe.MapNode(util.Function(input_names=['in_file', 'name'], output_names=['out_file'], function=rename), name='falff_smooth_o', iterfield=["in_file"])

RSFC_corr_r = pe.MapNode(interface=util.Rename(), name='RSFC_corr_r', iterfield=["in_file", "format_string"])
RSFC_corr_o = pe.MapNode(util.Function(input_names=['in_file', 'name', 'whichNode'], output_names=['out_file'], function=renameRSFC), name='RSFC_corr_o', iterfield=["in_file", "name"])

RSFC_z_trans_r = pe.MapNode(interface=util.Rename(), name='RSFC_z_trans_r', iterfield=["in_file", "format_string"])
RSFC_z_trans_o = pe.MapNode(util.Function(input_names=['in_file', 'name', 'whichNode'], output_names=['out_file'], function=renameRSFC), name='RSFC_z_trans_o', iterfield=["in_file", "name"])

RSFC_register_r = pe.MapNode(interface=util.Rename(), name='RSFC_register_r', iterfield=["in_file", "format_string"])
RSFC_register_o = pe.MapNode(util.Function(input_names=['in_file', 'name', 'whichNode'], output_names=['out_file'], function=renameRSFC), name='RSFC_register_o', iterfield=["in_file", "name"])


RSFC_smooth_r = pe.MapNode(interface=util.Rename(), name='RSFC_smooth_r', iterfield=["in_file", "format_string"])
RSFC_smooth_o = pe.MapNode(util.Function(input_names=['in_file', 'name', 'whichNode'], output_names=['out_file'], function=renameRSFC), name='RSFC_smooth_o', iterfield=["in_file", "name"])


func_calc = pe.MapNode(interface=e_afni.Threedcalc(), name='func_calc', iterfield=["infile_a", "stop_idx"])
func_calc.inputs.start_idx = 0
#calc.inputs.infile_a = analysisdirectory + '/' + subject + '/' + 'rest_1/rest.nii.gz[%d..%d]' %(int(first_vol),int(last_vol))
func_calc.inputs.expr = '\'a\''
#func_calc.inputs.out_file = 'rest_dr.nii.gz'

func_refit = pe.MapNode(interface=e_afni.Threedrefit(), name='func_refit', iterfield=["in_file"])
func_refit.inputs.deoblique = True

func_reorient = pe.MapNode(interface=e_afni.Threedresample(), name='func_reorient', iterfield=["in_file"])
func_reorient.inputs.orientation = 'RPI'
#reorient.inputs.infile = anat_name+'.nii.gz'
#reorient.inputs.out_file = os.path.abspath( analysisdirectory + '/' + subject + '/' + 'rest_1/rest_ro.nii.gz')

func_tstat = pe.MapNode(interface=e_afni.ThreedTstat(), name='func_tstat', iterfield=["in_file"])
func_tstat.inputs.options = "-mean"
#func_tstat.inputs.out_file = 'rest_ro_mean.nii.gz'

func_volreg = pe.MapNode(interface=e_afni.Threedvolreg(), name='func_volreg', iterfield=["in_file", "basefile"])
func_volreg.inputs.other = '-Fourier -twopass'
func_volreg.inputs.zpad = '4'
#func_volreg.inputs.oned_file = 'rest_mc.1D'
#func_volreg.inputs.out_file = 'rest_mc.nii.gz'


func_tstat_1 = pe.MapNode(interface=e_afni.ThreedTstat(), name='func_tstat_1', iterfield=["in_file"])
func_tstat_1.inputs.options = "-mean"
#func_tstat.inputs.out_file = 'rest_ro_mean.nii.gz'

func_volreg_1 = pe.MapNode(interface=e_afni.Threedvolreg(), name='func_volreg_1', iterfield=["in_file", "basefile"])
func_volreg_1.inputs.other = '-Fourier -twopass'
func_volreg_1.inputs.zpad = '4'
#func_volreg_1.inputs.oned_file = 'rest_mc.1D'
#func_volreg_1.inputs.out_file = 'rest_mc.nii.gz'

func_automask = pe.MapNode(interface=e_afni.ThreedAutomask(), name='func_automask', iterfield=["in_file"])
func_automask.inputs.dilate = 1
#automask.inputs.out_file = os.path.abspath(analysisdirectory + '/' + subject + '/' + 'rest_1/rest_mask.nii.gz')


func_calcR = pe.MapNode(interface=e_afni.Threedcalc(), name='func_calcR', iterfield=["infile_a", "infile_b"])
#func_calcR.inputs.infile_a = analysisdirectory + '/' + subject + '/' + 'rest_1/rest_mc.nii.gz'
#func_calcR.inputs.infile_b = analysisdirectory + '/' + subject + '/' + 'rest_1/rest_mask.nii.gz'
func_calcR.inputs.expr = '\'a*b\''
#func_calcR.inputs.out_file =  'rest_ss.nii.gz'

func_mean = pe.MapNode(interface=e_afni.ThreedTstat(), name='func_mean', iterfield=['in_file'])
func_mean.inputs.options = "-mean"

func_calcI = pe.MapNode(interface=e_afni.Threedcalc(), name='func_calcI', iterfield=["infile_a"])
#func_calcI.inputs.infile_a = analysisdirectory + '/' + subject + '/' + 'rest_1/rest_ss.nii.gz[7]'
func_calcI.inputs.single_idx = 4
func_calcI.inputs.expr = '\'a\''
#func_calcI.inputs.out_file = 'example_func.nii.gz'

func_despike = pe.MapNode(interface=e_afni.ThreedDespike(), name='func_despike', iterfield=["in_file"])
#func_despike.inputs.out_file = 'rest_ds.nii.gz'

func_smooth = pe.MapNode(interface=MultiImageMaths(), name='func_smooth', iterfield=["in_file", "operand_files"])
func_str1 = "-kernel gauss %f -fmean -mas" % (sigma)
func_smooth.inputs.op_string = func_str1 + " %s"
#func_smooth.inputs.out_file = os.path.abspath(analysisdirectory + '/' + subject + '/' + 'rest_1/rest_sm.nii.gz')

func_scale = pe.MapNode(interface=fsl.ImageMaths(), name='func_scale', iterfield=["in_file"])
func_scale.inputs.op_string = '-ing 10000'
func_scale.inputs.out_data_type = 'float'
#func_scale.inputs.out_file = os.path.abspath(analysisdirectory + '/' + subject + '/' + 'rest_1/rest_gms.nii.gz')

func_filter = pe.MapNode(interface=e_afni.ThreedFourier(), name='func_filter', iterfield=["in_file"])
#func_filter.inputs.highpass = str(hp)
#func_filter.inputs.lowpass = str(lp)
func_filter.inputs.other = '-retrend'
#func_filter.inputs.out_file = 'rest_filt.nii.gz'

func_detrenda = pe.MapNode(interface=e_afni.ThreedTstat(), name='func_detrenda', iterfield=["in_file"])
func_detrenda.inputs.options = '-mean'
#func_detrenda.inputs.out_file = 'rest_filt_mean.nii.gz'

func_detrendb = pe.MapNode(interface=e_afni.ThreedDetrend(), name='func_detrendb', iterfield=["in_file"])
func_detrendb.inputs.options = '-polort 2'
#func_detrendb.inputs.out_file = 'rest_dt.nii.gz'


func_detrendc = pe.MapNode(interface=e_afni.Threedcalc(), name='func_detrendc', iterfield=["infile_a", "infile_b"])
func_detrendc.inputs.expr = '\'a+b\''
#func_detrendc.inputs.out_file = 'rest_pp.nii.gz'

func_mask = pe.MapNode(interface=fsl.ImageMaths(), name='func_mask', iterfield=["in_file"])
func_mask.inputs.op_string = '-Tmin -bin'
func_mask.inputs.out_data_type = 'char'
#func_mask.inputs.out_file = os.path.abspath(analysisdirectory + '/' + subject + '/' + 'rest_1/rest_pp_mask.nii.gz')


##registration
reg_flirt = pe.MapNode(interface=fsl.FLIRT(), name='reg_flirt', iterfield=["in_file"])
#reg_flirt.inputs.in_file  = os.path.abspath(func_reg_dir+'/example_func.nii.gz')
#reg_flirt.inputs.reference = os.path.abspath(anat_reg_dir+'/highres.nii.gz')
#reg_flirt.inputs.out_file = 'example_func2highres.nii.gz'
#reg_flirt.inputs.out_matrix_file = os.path.abspath('example_func2highres.mat')
reg_flirt.inputs.cost = 'corratio'
reg_flirt.inputs.dof = 6
reg_flirt.inputs.interp = 'nearestneighbour'

# Create mat file for conversion from subject's anatomical to functional
reg_xfm1 = pe.MapNode(interface=fsl.ConvertXFM(), name='reg_xfm1', iterfield=["in_file"])
#reg_xfm1.inputs.out_file = 'highres2example_func.mat'
reg_xfm1.inputs.invert_xfm = True

## 4. T1->STANDARD
## NOTE THAT THIS IS Linear registration, you may want to use FNIRT (non-linear)
reg_flirt1 = pe.Node(interface=fsl.FLIRT(), name='reg_flirt1')
#reg_flirt1.inputs.in_file  = os.path.abspath(anat_reg_dir+'/highres.nii.gz')
#reg_flirt1.inputs.reference = os.path.abspath(anat_reg_dir+'/standard.nii.gz')
#reg_flirt1.inputs.out_file = 'highres2standard.nii.gz'
#reg_flirt1.inputs.out_matrix_file ='highres2standard.mat'
reg_flirt1.inputs.cost = 'corratio'
reg_flirt1.inputs.cost_func = 'corratio'
reg_flirt1.inputs.dof = 12
reg_flirt1.inputs.interp = 'nearestneighbour'

## Create mat file for conversion from standard to high res
reg_xfm2 = pe.Node(interface=fsl.ConvertXFM(), name='reg_xfm2')
#reg_xfm2.inputs.in_file = os.path.abspath(reg_dir+'/highres2standard.mat')
#reg_xfm2.inputs.out_file = 'standard2highres.mat'
reg_xfm2.inputs.invert_xfm = True

## 5. FUNC->STANDARD
## Create mat file for registration of functional to standard
#reg_xfm3.inputs.in_file = os.path.abspath(reg_dir+'/example_func2highres.mat')
#reg_xfm3.inputs.in_file2 = os.path.abspath(reg_dir+'/highres2standard.mat')
#reg_xfm3.inputs.out_file = 'example_func2standard.mat'
#reg_flirt2.inputs.in_file  = os.path.abspath(func_reg_dir+'/example_func.nii.gz')
#reg_flirt2.inputs.reference = os.path.abspath(anat_reg_dir+'/standard.nii.gz')
#reg_flirt2.inputs.out_file = 'example_func2standard.nii.gz'
#reg_flirt2.inputs.in_matrix_file = os.path.abspath(reg_dir+'/example_func2standard.mat')
#reg_xfm4.inputs.in_file = os.path.abspath(reg_dir+'/example_func2standard.mat')
#reg_xfm4.inputs.out_file = 'standard2example_func.mat'

reg_inw = pe.Node(interface=e_afni.InvWarp(), name='reg_inw')

## 6. T1->STANDARD NONLINEAR
# Perform nonlinear registration (higres to standard)
reg_fnt = pe.Node(interface=fsl.FNIRT(), name='reg_fnt')
#reg_fnt.inputs.in_file = os.path.abspath(anat_reg_dir+'/highres_head.nii.gz')
#reg_fnt.inputs.affine_file = os.path.abspath(anat_reg_dir+'/highres2standard.mat')
reg_fnt.inputs.fieldcoeff_file = True
#reg_fnt.inputs.warped_file = os.path.abspath(anat_reg_dir+'/highres2standard_NL.nii.gz')
reg_fnt.inputs.jacobian_file = True
reg_fnt.inputs.fieldcoeff_file = True
#reg_fnt.inputs.config_file = os.path.abspath(FSLDIR+'/etc/flirtsch/T1_2_MNI152_%s.cnf'%(standard_res) )
#reg_fnt.inputs.ref_file = os.path.abspath(FSLDIR+'/data/standard/MNI152_T1_%s.nii.gz' %(standard_res))
#reg_fnt.inputs.refmask_file = os.path.abspath(FSLDIR+'/data/standard/MNI152_T1_%s_brain_mask_dil.nii.gz' %(standard_res))
reg_fnt.inputs.warp_resolution = (10, 10, 10)

reg_warp = pe.MapNode(interface=fsl.ApplyWarp(), name='reg_warp', iterfield=["in_file", "premat"])
#reg_warp.inputs.in_file = os.path.abspath(func_reg_dir+'/example_func.nii.gz')
#reg_warp.inputs.ref_file = os.path.abspath(FSLDIR + '/data/standard/MNI152_T1_%s.nii.gz' %(standard_res))
#reg_warp.inputs.field_file = 'highres2standard_warp.nii.gz'
#reg_warp.inputs.premat = os.path.abspath(func_reg_dir + '/example_func2highres.mat')
#reg_warp.inputs.out_file = os.path.abspath(func_reg_dir + '/example_func2standard_NL.nii.gz')


#segmentation

seg_segment = pe.Node(interface=fsl.FAST(), name='seg_segment')
#seg_segment.inputs.in_files = os.path.abspath(anat_dir+'/mprage_brain.nii.gz')
seg_segment.inputs.img_type = 1
seg_segment.inputs.segments = True
seg_segment.inputs.probability_maps = True
seg_segment.inputs.out_basename = 'segment'

## 4. Copy functional mask from FSLpreproc step 5 - this is the global signal mask
seg_copy = pe.MapNode(interface=e_afni.Threedcopy(), name='seg_copy', iterfield=["in_file"])
#seg_copy.inputs.in_file = os.path.abspath(func_dir+'/rest_pp_mask.nii.gz')
#seg_copy.inputs.out_file = 'global_mask.nii.gz'

## CSF
## 5. Register csf to native space
seg_flirt = pe.MapNode(interface=fsl.FLIRT(), name='seg_flirt', iterfield=["reference", "in_matrix_file"])
#seg_flirt.inputs.in_file  = os.path.abspath(anat_dir+'/segment_prob_0.nii.gz')
#seg_flirt.inputs.reference = os.path.abspath(func_reg_dir+'/example_func.nii.gz')
#seg_flirt.inputs.out_file = os.path.abspath(segment_dir+'/csf2func.nii.gz')
#seg_flirt.inputs.in_matrix_file = os.path.abspath(func_reg_dir+'/highres2example_func.mat')
seg_flirt.inputs.apply_xfm = True

## 6 register CSF template (in MNI space) to native space
seg_warp = pe.MapNode(interface=fsl.ApplyWarp(), name='seg_warp', iterfield=["ref_file", "postmat"])
seg_warp.inputs.interp = 'nn'

seg_warp1 = seg_warp.clone('seg_warp1')
## 7. register to standard

## 8. find overlap with prior
seg_smooth1 = pe.MapNode(interface=MultiImageMaths(), name='seg_smooth1', iterfield=["in_file"])
seg_str1 = "-mas %s "
#seg_smooth1.inputs.in_file = os.path.abspath(segment_dir+'/csf2standard.nii.gz')
seg_smooth1.inputs.op_string = seg_str1
#seg_smooth1.inputs.operand_files = os.path.abspath(PRIOR_CSF)
#seg_smooth1.inputs.out_file = os.path.abspath(segment_dir+'/csf_masked.nii.gz')

## 9. revert back to functional space

## 10. Threshold and binarize probability map of csf
seg_thresh = pe.MapNode(interface=fsl.ImageMaths(), name='seg_thresh', iterfield=["in_file"])
seg_str1 = "-thr 0.4 -bin "
#seg_thresh.inputs.in_file = os.path.abspath(segment_dir+'/csf_native.nii.gz')
seg_thresh.inputs.op_string = seg_str1
#seg_thresh.inputs.out_file = os.path.abspath(segment_dir+'/csf_bin.nii.gz')

## 11. Mask again by the subject's functional
seg_mask = pe.MapNode(interface=fsl.MultiImageMaths(), name='seg_mask', iterfield=["in_file", "operand_files"])
seg_str1 = "-mas %s "
#seg_mask.inputs.in_file = os.path.abspath(segment_dir+'/csf_bin.nii.gz')
seg_mask.inputs.op_string = seg_str1
#seg_mask.inputs.operand_files = os.path.abspath(segment_dir+'global_mask.nii.gz')
#seg_mask.inputs.out_file = os.path.abspath(segment_dir+'/csf_mask.nii.gz')

## WM
## 12. Register wm to native space
seg_flirt3 = pe.MapNode(interface=fsl.FLIRT(), name='seg_flirt3', iterfield=["reference", "in_matrix_file"])
#seg_flirt3.inputs.in_file  = os.path.abspath(anat_dir+'/segment_prob_2.nii.gz')
#seg_flirt3.inputs.reference = os.path.abspath(func_reg_dir+'/example_func.nii.gz')
#seg_flirt3.inputs.out_file = os.path.abspath(segment_dir+'/wm2func.nii.gz')
#seg_flirt3.inputs.in_matrix_file = os.path.abspath(func_reg_dir+'/highres2example_func.mat')
seg_flirt3.inputs.apply_xfm = True

## 13. Smooth image to match smoothing on functional
seg_smooth2 = pe.MapNode(interface=fsl.ImageMaths(), name='seg_smooth2', iterfield=["in_file"])
seg_str1 = "-kernel gauss %f -fmean " % (sigma)
#seg_smooth2.inputs.in_file = os.path.abspath(segment_dir+'/wm2func.nii.gz')
seg_smooth2.inputs.op_string = seg_str1
#seg_smooth2.inputs.out_file = os.path.abspath(segment_dir+'/wm_sm.nii.gz')

## 14. register to standard

## 15. find overlap with prior
seg_prior1 = pe.MapNode(interface=fsl.MultiImageMaths(), name='seg_prior1', iterfield=["in_file"])
seg_str1 = "-mas %s "
#seg_prior1.inputs.in_file = os.path.abspath(segment_dir+'/wm2standard.nii.gz')
seg_prior1.inputs.op_string = seg_str1
#seg_prior1.inputs.operand_files = os.path.abspath(PRIOR_WHITE)
#seg_prior1.inputs.out_file = os.path.abspath(segment_dir+'/wm_masked.nii.gz')


## 17. Threshold and binarize probability map of wm
seg_thresh1 = pe.MapNode(interface=fsl.ImageMaths(), name='seg_thresh1', iterfield=["in_file"])
seg_str1 = "-thr 0.66 -bin "
#seg_thresh1.inputs.in_file = os.path.abspath(segment_dir+'/wm_native.nii.gz')
seg_thresh1.inputs.op_string = seg_str1
#seg_thresh1.inputs.out_file = os.path.abspath(segment_dir+'/wm_bin.nii.gz')

## 18. Mask again by the subject's functional
seg_mask1 = pe.MapNode(interface=fsl.MultiImageMaths(), name='seg_mask1', iterfield=["in_file", "operand_files"])
seg_str1 = "-mas %s "
#seg_mask1.inputs.in_file = os.path.abspath(segment_dir+'/wm_bin.nii.gz')
seg_mask1.inputs.op_string = seg_str1
#seg_mask1.inputs.operand_files = os.path.abspath(segment_dir+'global_mask.nii.gz')
#seg_mask1.inputs.out_file = os.path.abspath(segment_dir+'/wm_mask.nii.gz')


lifeSaver_nuisance_compcor = pe.MapNode(interface=util.Rename(), name='lifeSaver_nuisance_compcor', iterfield=["in_file", "format_string"])

lifeSaver_nuisance_MedianAngle = pe.MapNode(interface=util.Rename(), name='lifeSaver_nuisance_MedianAngle', iterfield=["in_file", "format_string"])

lifeSaver_nuisance_erosion_csf1 = pe.MapNode(interface=util.Rename(), name='lifeSaver_nuisance_erosion_csf1', iterfield=["in_file", "format_string"])

lifeSaver_nuisance_erosion_wm1 = pe.MapNode(interface=util.Rename(), name='lifeSaver_nuisance_erosion_wm1', iterfield=["in_file", "format_string"])

lifeSaver_nuisance_globalE = pe.MapNode(interface=util.Rename(), name='lifeSaver_nuisance_globalE', iterfield=["in_file", "format_string"])

lifeSaver_nuisance_csf = pe.MapNode(interface=util.Rename(), name='lifeSaver_nuisance_csf', iterfield=["in_file", "format_string"])

lifeSaver_nuisance_wm = pe.MapNode(interface=util.Rename(), name='lifeSaver_nuisance_wm', iterfield=["in_file", "format_string"])

lifeSaver_nuisance_featM = pe.MapNode(interface=util.Rename(), name='lifeSaver_nuisance_featM', iterfield=["in_file", "format_string"])

lifeSaver_nuisance_fgls = pe.MapNode(interface=util.Rename(), name='lifeSaver_nuisance_fgls', iterfield=["in_file", "format_string"])

lifeSaver_nuisance_fglsd = pe.MapNode(interface=util.Rename(), name='lifeSaver_nuisance_fglsd', iterfield=["in_file", "format_string"])

lifeSaver_nuisance_stat = pe.MapNode(interface=util.Rename(), name='lifeSaver_nuisance_stat', iterfield=["in_file", "format_string"])

lifeSaver_nuisance_calc = pe.MapNode(interface=util.Rename(), name='lifeSaver_nuisance_calc', iterfield=["in_file", "format_string"])

lifeSaver_func_filter = pe.MapNode(interface=util.Rename(), name='lifeSaver_func_filter', iterfield=["in_file", "format_string"])

lifeSaver_nuisance_warp = pe.MapNode(interface=util.Rename(), name='lifeSaver_nuisance_warp', iterfield=["in_file", "format_string"])

lifeSaver_nuisance_warp_1 = pe.MapNode(interface=util.Rename(), name='lifeSaver_nuisance_warp_1', iterfield=["in_file", "format_string"])

lifeSaver_RSFC_printToFile = pe.MapNode(interface=util.Rename(), name='lifeSaver_RSFC_printToFile', iterfield=["in_file", "format_string"])

lifeSaverRSFC = pe.MapNode(interface=util.Rename(), name='lifeSaver_RSFC', iterfield=["in_file", "format_string"])

lifeSaverRSFC1 = pe.MapNode(interface=util.Rename(), name='lifeSaver_RSFC1', iterfield=["in_file", "format_string"])

lifeSaver_RSFC_corr = pe.MapNode(interface=util.Rename(), name='lifeSaver_RSFC_corr', iterfield=["in_file", "format_string"])

lifeSaver_RSFC_z_trans = pe.MapNode(interface=util.Rename(), name='lifeSaver_RSFC_z_trans', iterfield=["in_file", "format_string"])

lifeSaver_RSFC_register = pe.MapNode(interface=util.Rename(), name='lifeSaver_RSFC_register', iterfield=["in_file", "format_string"])

lifeSaver_RSFC_smooth = pe.MapNode(interface=util.Rename(), name='lifeSaver_RSFC_smooth', iterfield=["in_file", "format_string"])


def getPolyRegressors(nvols, nPoly):

    import sys
    import os
    import commands

    cmd = 'poly_trends.py %d %d' % (int(nvols), int(nPoly))

    sys.stderr.write(commands.getoutput(cmd))

    if os.path.exists(os.path.join(os.getcwd(), 'poly_detrend_1.txt')) and os.path.exists(os.path.join(os.getcwd(), 'poly_detrend_2.txt')):
        regressors = [os.path.join(os.getcwd(), 'poly_detrend_1.txt'), os.path.join(os.getcwd(), 'poly_detrend_2.txt')]
        return regressors
    else:
        sys.stderr.write('\nError: output files from poly_trends.py not found. check! \n')
        regressors = []
        return regressors

nuisance_poly = pe.MapNode(util.Function(input_names=['nvols', 'nPoly'], output_names=['regressors'], function=getPolyRegressors), name='nuisance_poly', iterfield=['nvols'])
nuisance_poly.inputs.nPoly = 2

nuisance_erosion_csf = pe.MapNode(interface=e_afni.Threedcalc(), name='nuisance_erosion_csf', iterfield=["infile_a"])
nuisance_erosion_csf.inputs.infile_b_prime = 'a+i -c a-i -d a+j -e a-j -f a+k -g a-k'
nuisance_erosion_csf.inputs.expr = "'a*(1-amongst(0,b,c,d,e,f,g))'"

nuisance_erosion_csf1 = pe.MapNode(interface=e_afni.Threedcalc(), name='nuisance_erosion_csf1', iterfield=["infile_a"])
nuisance_erosion_csf1.inputs.infile_b_prime = 'a+i -c a-i -d a+j -e a-j -f a+k -g a-k'
nuisance_erosion_csf1.inputs.expr = "'a*(1-amongst(0,b,c,d,e,f,g))'"

nuisance_erosion_wm = pe.MapNode(interface=e_afni.Threedcalc(), name='nuisance_erosion_wm', iterfield=["infile_a"])
nuisance_erosion_wm.inputs.infile_b_prime = 'a+i -c a-i -d a+j -e a-j -f a+k -g a-k'
nuisance_erosion_wm.inputs.expr = "'a*(1-amongst(0,b,c,d,e,f,g))'"

nuisance_erosion_wm1 = pe.MapNode(interface=e_afni.Threedcalc(), name='nuisance_erosion_wm1', iterfield=["infile_a"])
nuisance_erosion_wm1.inputs.infile_b_prime = 'a+i -c a-i -d a+j -e a-j -f a+k -g a-k'
nuisance_erosion_wm1.inputs.expr = "'a*(1-amongst(0,b,c,d,e,f,g))'"

nuisance_globalE = pe.MapNode(interface=e_afni.ThreedMaskave(), name='nuisance_globalE', iterfield=["in_file", "mask"])
#nuisance_globalE.inputs.in_file = os.path.abspath(func_dir + '/' + 'rest_pp.nii.gz')
#nuisance_globalE.inputs.mask = os.path.abspath(segment_dir + '/' + 'global_mask.nii.gz')
nuisance_globalE.inputs.quiet = True
#nuisance_globalE.inputs.out_file = 'global.1D'

## 4. csf
nuisance_csf = pe.MapNode(interface=e_afni.ThreedMaskave(), name='nuisance_csf', iterfield=["in_file", "mask"])
#nuisance_csf.inputs.in_file = os.path.abspath(func_dir + '/' + 'rest_pp.nii.gz')
#nuisance_csf.inputs.mask = os.path.abspath(segment_dir + '/' + 'csf_mask.nii.gz')
nuisance_csf.inputs.quiet = True
#nuisance_csf.inputs.out_file =  'csf.1D'

## 5. wm
nuisance_wm = pe.MapNode(interface=e_afni.ThreedMaskave(), name='nuisance_wm', iterfield=["in_file", "mask"])
#nuisance_wm.inputs.in_file = os.path.abspath(func_dir + '/' + 'rest_pp.nii.gz')
#nuisance_wm.inputs.mask = os.path.abspath(segment_dir + '/' + 'wm_mask.nii.gz')
nuisance_wm.inputs.quiet = True
#nuisance_wm.inputs.out_file = 'wm.1D'

nuisance_compcor = pe.MapNode(interface=e_afni.compcor(), name='nuisance_compcor', iterfield=["in_file", "wmmask", "csfmask"])
nuisance_compcor.inputs.ncomponents = 5


nuisance_MedianAngle = pe.MapNode(interface=e_afni.MedianAngle(), name='nuisance_MedianAngle', iterfield=["in_file"])
nuisance_MedianAngle.inputs.angle = 90.0

nuisance_concatnode = pe.MapNode(interface=util.Merge(3), name='nuisance_concatnode', iterfield=["in1", "in2", "in3"])

nuisance_selectnode = pe.MapNode(interface=util.Select(), name='nuisance_selectnode', iterfield=["inlist", "index"])


def makeRegressionDecision(pp, which_regression):

    decisions = []

    if which_regression.lower() == "default":
        print '\n\n INSIDE 0'
        decisions.append(0)
    elif which_regression.lower() == "compcor":
        print '\n\n INSIDE 1'
        decisions.append(1)
    elif which_regression.lower() == "median_angle":
        print '\n\n INSIDE 2'
        decisions.append(2)
    else:
        print '\n\n INSIDE 3'
        decisions.append(0)

    print '\n\n REGRESSION  -->%s DECISION: ' % (which_regression.lower()) + str(decisions)
    return decisions


def makeTemplateDecision(which_regression):

    decisions = 0
    if which_regression.lower() == "default":

        decisions = 0
    else:
        decisions = 1

    return decisions


nuisance_decision = pe.MapNode(util.Function(input_names=['pp', 'which_regression'], output_names=['decisions'], function=makeRegressionDecision), name='nuisance_decision', iterfield=["pp"])


nuisance_template_concatnode = pe.Node(interface=util.Merge(2), name='nuisance_template_concatnode', iterfield=["in1", "in2"])

nuisance_template_selectnode = pe.Node(interface=util.Select(), name='nuisance_template_selectnode', )


nuisance_template_decision = pe.Node(util.Function(input_names=['which_regression'], output_names=['decisions'], function=makeTemplateDecision), name='nuisance_template_decision')

## feat model
nuisance_featM = pe.MapNode(interface=fsl.FEATModel(), name='nuisance_featM', iterfield=["fsf_file", "ev_files"])
#nuisance_featM.inputs.fsf_file = os.path.abspath(nuisance_dir + '/nuisance.fsf')
#nuisance_featM.inputs.ev_files = EV_List
#nuisance_featM.inputs.design_file = os.path.abspath(nuisance_dir + '/nuisance.mat')

nuisance_brick = pe.MapNode(interface=e_afni.ThreedBrickStat(), name='nuisance_brick', iterfield=["in_file", "mask"])
#nuisance_brick.inputs.in_file = os.path.abspath(func_dir + '/rest_pp.nii.gz')
#nuisance_brick.inputs.mask = os.path.abspath(func_dir + '/rest_pp_mask.nii.gz')
nuisance_brick.inputs.min = True

## 7. Get residuals
nuisance_fgls = pe.MapNode(interface=fsl.FILMGLS(), name='nuisance_fgls', iterfield=["in_file", "results_dir", "design_file", "threshold"])

#nuisance_fgls.inputs.in_file = os.path.abspath(func_dir + '/rest_pp.nii.gz')
#nuisance_fgls.inputs.results_dir = nuisance_dir + '/stats'
nuisance_fgls.inputs.mask_size = 5
nuisance_fgls.inputs.smooth_autocorr = True
nuisance_fgls.inputs.autocorr_noestimate = True
#nuisance_fgls.inputs.design_file = os.path.abspath(nuisance_dir + '/nuisance.mat')

## 8. Demeaning residuals and ADDING 100
nuisance_stat = pe.MapNode(interface=e_afni.ThreedTstat(), name='nuisance_stat', iterfield=["in_file"])
#nuisance_stat.inputs.in_file = os.path.abspath(nuisance_dir + '/stats/res4d.nii.gz')
nuisance_stat.inputs.options = ' -mean '
#nuisance_stat.inputs.out_file = 'res4d_mean.nii.gz'

nuisance_calc = pe.MapNode(interface=e_afni.Threedcalc(), name='nuisance_calc', iterfield=["infile_a", "infile_b"])
#nuisance_calc.inputs.infile_a = os.path.abspath(nuisance_dir + '/stats/res4d.nii.gz')
#nuisance_calc.inputs.infile_b = os.path.abspath(nuisance_dir + '/stats/res4d_mean.nii.gz')
nuisance_calc.inputs.expr = '\'(a-b)+100\''
#nuisance_calc.inputs.out_file = 'rest_res.nii.gz'

## 9. Resampling residuals to MNI space
nuisance_warp = pe.MapNode(interface=fsl.ApplyWarp(), name='nuisance_warp', iterfield=["in_file", "premat"])
#nuisance_warp.inputs.in_file = os.path.abspath(func_dir+'/rest_res.nii.gz')
#nuisance_warp.inputs.ref_file = os.path.abspath( FSLDIR + '/data/standard/MNI152_T1_%s.nii.gz' %(standard_res) )
#nuisance_warp.inputs.field_file = os.path.abspath(anat_reg_dir + '/highres2standard_warp.nii.gz')
#nuisance_warp.inputs.premat = os.path.abspath(func_reg_dir + '/example_func2highres.mat')
#nuisance_warp.inputs.out_file = os.path.abspath(func_dir + '/rest_res2standard.nii.gz')

nuisance_warp_1 = pe.MapNode(interface=fsl.ApplyWarp(), name='nuisance_warp_1', iterfield=["in_file", "premat"])

alff_detrend = pe.MapNode(interface=e_afni.ThreedTcat(), name='alff_detrend', iterfield=["in_file"])
#alff_detrend.inputs.in_file = os.path.abspath(func_dir + '/rest_ds.nii.gz')
#alff_detrend.inputs.out_file = 'prefiltered_func_data_rlt.nii.gz'
alff_detrend.inputs.rlt = '+'

## 3. Spatial Smoothing
alff_smooth = pe.MapNode(interface=MultiImageMaths(), name='alff_smooth', iterfield=["in_file", "operand_files"])
alff_str1 = "-kernel gauss %f -fmean -mas " % (sigma)
alff_smooth.inputs.op_string = alff_str1 + " %s"
#alff_smooth.inputs.in_file = os.path.abspath(alff_dir+'/prefiltered_func_data_rlt.nii.gz')
#alff_smooth.inputs.operand_files = os.path.abspath(func_dir + '/rest_mask.nii.gz')
#alff_smooth.inputs.out_file = os.path.abspath(alff_dir+'/prefiltered_func_data_smooth.nii.gz')

falff_smooth = pe.MapNode(interface=MultiImageMaths(), name='falff_smooth', iterfield=["in_file", "operand_files"])
falff_str1 = "-kernel gauss %f -fmean -mas" % (sigma)
falff_smooth.inputs.op_string = falff_str1 + " %s"
#func_smooth.inputs.out_file = os.path.abspath(analysisdirectory + '/' + subject + '/' + 'rest_1/rest_sm.nii.gz')
## 4. Grand mean Scaling
alff_scale = pe.MapNode(interface=fsl.ImageMaths(), name='alff_scale', iterfield=["in_file"])
alff_scale.inputs.op_string = '-ing 10000'
#alff_scale.inputs.in_file = os.path.abspath(alff_dir+'/prefiltered_func_data_smooth.nii.gz')
#alff_scale.inputs.out_file = os.path.abspath(alff_dir+'/prefiltered_func_data_inn_volsorm.nii.gz')

alff_cp = pe.MapNode(interface=fsl.ImageMaths(), name='alff_cp', iterfield=["in_file"])
#alff_cp.inputs.in_file = os.path.abspath(alff_dir+'/prefiltered_func_data_inn_volsorm.nii.gz')
#alff_cp.inputs.out_file = os.path.abspath(alff_dir+'/rest_alff_pp.nii.gz')

alff_mean = pe.MapNode(interface=fsl.ImageMaths(), name='alff_mean', iterfield=["in_file"])
#alff_mean.inputs.in_file = os.path.abspath(alff_dir+'/rest_alff_pp.nii.gz')
alff_mean.inputs.op_string = '-Tmean'
#alff_mean.inputs.out_file = os.path.abspath(alff_dir+'/mean_rest_alff_pp.nii.gz')

## 4. Calculate fALFF
alff_falff = pe.MapNode(interface=fsl.ImageMaths(), name='alff_falff', iterfield=["in_file", "op_string"])
#alff_falff.inputs.op_string = '-Tmean -mul %s -div 2' %(n_vols)
#alff_falff.inputs.in_file = os.path.abspath(alff_dir+'/prealff_func_data_ps_sqrt.nii.gz')
#alff_falff.inputs.out_file = os.path.abspath(alff_dir+'/prealff_func_data_ps_sum.nii.gz')

alff_falff1 = pe.MapNode(interface=MultiImageMaths(), name='falff1', iterfield=["in_file", "operand_files"])
alff_falff1.inputs.op_string = "-div %s"
#alff_falff1.inputs.in_file = os.path.abspath(alff_dir+'/prealff_func_ps_alff4slow.nii.gz')
#alff_falff1.inputs.operand_files = os.path.abspath(alff_dir+'/prealff_func_data_ps_sum.nii.gz')
#alff_falff1.inputs.out_file = os.path.abspath(alff_dir+'/fALFF.nii.gz')

## 5. Z-normalisation across whole brain
## 4. Calculate fALFF
alff_falff = pe.MapNode(interface=fsl.ImageMaths(), name='alff_falff', iterfield=["in_file", "op_string"])
#alff_falff.inputs.op_string = '-Tmean -mul %s -div 2' %(n_vols)
#alff_falff.inputs.in_file = os.path.abspath(alff_dir+'/prealff_func_data_ps_sqrt.nii.gz')
#alff_falff.inputs.out_file = os.path.abspath(alff_dir+'/prealff_func_data_ps_sum.nii.gz')

alff_falff1 = pe.MapNode(interface=MultiImageMaths(), name='falff1', iterfield=["in_file", "operand_files"])
alff_falff1.inputs.op_string = "-div %s"
#alff_falff1.inputs.in_file = os.path.abspath(alff_dir+'/prealff_func_ps_alff4slow.nii.gz')
#alff_falff1.inputs.operand_files = os.path.abspath(alff_dir+'/prealff_func_data_ps_sum.nii.gz')
#alff_falff1.inputs.out_file = os.path.abspath(alff_dir+'/fALFF.nii.gz')

## 5. Z-normalisation across whole brain
alff_normM = pe.MapNode(interface=fsl.ImageStats(), name='alff_normM', iterfield=["in_file", "mask_file"])
#alff_normM.inputs.in_file = os.path.abspath(alff_dir+'/ALFF.nii.gz')
#alff_normM.inputs.mask_file = os.path.abspath(func_dir+'/rest_mask.nii.gz')
alff_normM.inputs.op_string = "-k %s -m"

alff_normS = pe.MapNode(interface=fsl.ImageStats(), name='alff_normS', iterfield=["in_file", "mask_file"])
#alff_normS.inputs.in_file = os.path.abspath(alff_dir+'/ALFF.nii.gz')
#alff_normS.inputs.mask_file = os.path.abspath(func_dir+'/rest_mask.nii.gz')
alff_normS.inputs.op_string = "-k %s -s"

alff_normM1 = pe.MapNode(interface=fsl.ImageStats(), name='alff_normM1', iterfield=["in_file", "mask_file"])
#alff_normM1.inputs.in_file = os.path.abspath(alff_dir+'/fALFF.nii.gz')
#alff_normM1.inputs.mask_file = os.path.abspath(func_dir+'/rest_mask.nii.gz')
alff_normM1.inputs.op_string = "-k %s -m"

alff_normS1 = pe.MapNode(interface=fsl.ImageStats(), name='alff_normS1', iterfield=["in_file", "mask_file"])
#alff_normS1.inputs.in_file = os.path.abspath(alff_dir+'/fALFF.nii.gz')
#alff_normS1.inputs.mask_file = os.path.abspath(func_dir+'/rest_mask.nii.gz')
alff_normS1.inputs.op_string = "-k %s -s"

alff_Z_alff = pe.MapNode(interface=MultiImageMaths(), name='alff_Z_alff', iterfield=["in_file", "operand_files", "op_string"])
#alff_str = "-sub %f -div %f" %(mean,std_dev)
#alff_str1 = str + " -mas %s "
#alff_Z_alff.inputs.in_file = os.path.abspath(alff_dir+'/ALFF.nii.gz')
#alff_Z_alff.inputs.op_string = alff_str1
#alff_Z_alff.inputs.operand_files = os.path.abspath(func_dir+'/rest_mask.nii.gz')
#alff_Z_alff.inputs.out_file = os.path.abspath(alff_dir+'/ALFF_Z.nii.gz')

alff_Z_falff = pe.MapNode(interface=MultiImageMaths(), name='alff_Z_falff', iterfield=["in_file", "operand_files", "op_string"])
#alff_str = "-sub %f -div %f" %(mean1,std_dev1)
#alff_str1 = str + " -mas %s "
#alff_Z_falff.inputs.in_file = os.path.abspath(alff_dir+'/fALFF.nii.gz')
#alff_Z_falff.inputs.op_string = alff_str1

#alff_Z_falff.inputs.operand_files = os.path.abspath(func_dir+'/rest_mask.nii.gz')
#alff_Z_falff.inputs.out_file = os.path.abspath(alff_dir+'/fALFF_Z.nii.gz')


#Registering Z-transformed ALFF to standard space
alff_warp_alff = pe.MapNode(interface=fsl.ApplyWarp(), name='alff_warp_alff', iterfield=["in_file", "premat"])
#alff_warp_alff.inputs.in_file = os.path.abspath(alff_dir+'/ALFF_Z.nii.gz')
#alff_warp_alff.inputs.ref_file = os.path.abspath( FSLDIR + '/data/standard/MNI152_T1_%s.nii.gz' %(standard_res))
#alff_warp_alff.inputs.field_file = os.path.abspath(anat_reg_dir + '/highres2standard_warp.nii.gz')
#alff_warp_alff.inputs.premat = os.path.abspath(func_reg_dir + '/example_func2highres.mat')
#alff_warp_alff.inputs.out_file = os.path.abspath(alff_dir + '/ALFF_Z_2standard.nii.gz')

alff_warp_falff = pe.MapNode(interface=fsl.ApplyWarp(), name='alff_warp_falff', iterfield=["in_file", "premat"])
#alff_warp_falff.inputs.in_file = os.path.abspath(alff_dir+'/fALFF_Z.nii.gz')
#alff_warp_falff.inputs.ref_file = os.path.abspath( FSLDIR + '/data/standard/MNI152_T1_%s.nii.gz' %(standard_res))
#alff_warp_falff.inputs.field_file = os.path.abspath(anat_reg_dir + '/highres2standard_warp.nii.gz')
#alff_warp_falff.inputs.premat = os.path.abspath(func_reg_dir + '/example_func2highres.mat')
#alff_warp_falff.inputs.out_file = os.path.abspath(alff_dir + '/fALFF_Z_2standard.nii.gz')

alff_roi = pe.MapNode(interface=fsl.ExtractROI(), name="alff_roi", iterfield=["in_file", "t_size"])
#alff_roi.inputs.in_file = os.path.abspath(alff_dir+'/rest_alff_pp.nii.gz')
#alff_roi.inputs.roi_file = os.path.abspath(alff_dir+'/prealff_func_data.nii.gz')
alff_roi.inputs.t_min = 1
#alff_roi.inputs.t_size = int(n_vols)
#Registering Z-transformed fALFF to standard space
alff_cp1 = pe.MapNode(interface=fsl.ImageMaths(), name='alff_cp1', iterfield=["in_file"])
#alff_cp1.inputs.in_file = os.path.abspath(alff_dir+'/rest_alff_pp.nii.gz')
#alff_cp1.inputs.out_file = os.path.abspath(alff_dir+'/prealff_func_data.nii.gz')

alff_concatnode = pe.MapNode(interface=util.Merge(2), name='alff_concatnode', iterfield=["in1", "in2"])

alff_selectnode = pe.MapNode(interface=util.Select(), name='alff_selectnode', iterfield=["inlist", "index"])

alff_pspec = pe.MapNode(interface=fsl.PowerSpectrum(), name='alff_pspec', iterfield=["in_file"])
#alff_pspec.inputs.in_file = os.path.abspath(alff_dir+'/prealff_func_data.nii.gz')
#alff_pspec.inputs.out_file = os.path.abspath(alff_dir+'/prealff_func_data_ps.nii.gz')

##compute sqrt of power spectrum
alff_sqrt = pe.MapNode(interface=fsl.ImageMaths(), name='alff_sqrt', iterfield=["in_file"])
alff_sqrt.inputs.op_string = '-sqrt'
#alff_sqrt.inputs.in_file = os.path.abspath(alff_dir+'/prealff_func_data_ps.nii.gz')
#alff_sqrt.inputs.out_file = os.path.abspath(alff_dir+'/prealff_func_data_ps_sqrt.nii.gz')

alff_roi1 = pe.MapNode(interface=fsl.ExtractROI(), name="alff_roi1", iterfield=["in_file", "t_min", "t_size"])
#alff_roi1.inputs.in_file = os.path.abspath(alff_dir+'/prealff_func_data_ps_sqrt.nii.gz')
#alff_roi1.inputs.roi_file = os.path.abspath(alff_dir+'/prealff_func_ps_slow.nii.gz')
#alff_roi1.inputs.t_min = int(n1)
#alff_roi1.inputs.t_size = int(n2)

## calculate ALFF as the sum of the amplitudes in the low frequency band
alff_sum = pe.MapNode(interface=fsl.ImageMaths(), name='alff_sum', iterfield=["in_file", "op_string"])
#alff_sum.inputs.op_string = '-Tmean -mul %f' %(n2)
#alff_sum.inputs.in_file = os.path.abspath(alff_dir+'/prealff_func_ps_slow.nii.gz')
#alff_sum.inputs.out_file = os.path.abspath(alff_dir+'/prealff_func_ps_alff4slow.nii.gz')


## 0. Register Seed template in to native space
RSFC_warp = pe.MapNode(interface=fsl.ApplyWarp(), name='RSFC_warp', iterfield=["ref_file", "postmat"])
RSFC_warp.inputs.interp = 'nn'
RSFC_warp.iterables = ("in_file", seed_list)
## 1. Extract Timeseries
RSFC_time_series = pe.MapNode(interface=e_afni.ThreedROIstats(), name='RSFC_time_series', iterfield=["in_file", "mask"])
#RSFC_time_series.inputs.in_file = os.path.abspath(func_dir+'/rest_res2standard.nii.gz')
#RSFC_time_series.inputs.mask = os.path.abspath(seed)
RSFC_time_series.inputs.quiet = True
RSFC_time_series.inputs.mask_f2short = True
#RSFC_time_series.iterables = ("mask",seed_list)
#RSFC_(time_series.run()).outputs.stats

## 2. Compute voxel-wise correlation with Seed Timeseries
RSFC_corr = pe.MapNode(interface=e_afni.Threedfim(), name='RSFC_corr', iterfield=["in_file", "ideal_file"])
#RSFC_corr.inputs.in_file = os.path.abspath(func_dir + '/rest_res.nii.gz')
#RSFC_corr.inputs.ideal_file = os.path.abspath(seed_ts_dir + '/%s.1D'%(seed_name))
RSFC_corr.inputs.fim_thr = 0.0009
RSFC_corr.inputs.out = 'Correlation'
RSFC_corr.inputs.out_file = 's_corr.nii.gz'

## 3. Z-transform correlations
RSFC_z_trans = pe.MapNode(interface=e_afni.Threedcalc(), name='RSFC_z_trans', iterfield=["infile_a"])
#RSFC_z_trans.inputs.infile_a = os.path.abspath(RSFC_dir+'/%s_corr.nii.gz'%(seed_name))
RSFC_z_trans.inputs.expr = '\'log((1+a)/(1-a))/2\''
RSFC_z_trans.inputs.out_file = 's_Z.nii.gz'

## 4. Register Z-transformed correlations to standard space
RSFC_register = pe.MapNode(interface=fsl.ApplyWarp(), name='RSFC_register', iterfield=["premat", "in_file"])
#RSFC_register.inputs.in_file = os.path.abspath(RSFC_dir+'/%s_Z.nii.gz'%(seed_name))
#RSFC_register.inputs.ref_file = os.path.abspath(FSLDIR + '/data/standard/MNI152_T1_%s.nii.gz' %(standard_res))
#RSFC_register.inputs.field_file = os.path.abspath(anat_reg_dir + '/highres2standard_warp.nii.gz')
#RSFC_register.inputs.premat = os.path.abspath(func_reg_dir + '/example_func2highres.mat')
#RSFC_register.inputs.out_file = os.path.abspath(RSFC_dir + '/%s_Z_2standard.nii.gz' %(seed_name))

RSFC_smooth = pe.MapNode(interface=MultiImageMaths(), name='RSFC_smooth', iterfield=["in_file", "operand_files"])
RSFC_str1 = "-kernel gauss %f -fmean -mas " % (sigma)
RSFC_smooth.inputs.op_string = RSFC_str1 + " %s"

## Linear registration of T1 --> symmetric standard
VMHC_flirt = pe.Node(interface=fsl.FLIRT(), name='VMHC_flirt')
#VMHC_flirt.inputs.in_file  = os.path.abspath(${anat_dir}/brain.nii.gz)
#VMHC_flirt.inputs.reference = os.path.abspath(${symm_standard_brain})
#VMHC_flirt.inputs.out_file = '${anat_reg_dir}/highres2symmstandard.nii.gz'
#VMHC_flirt.inputs.out_matrix_file = os.path.abspath('${anat_reg_dir}/highres2symmstandard.mat')
VMHC_flirt.inputs.cost = 'corratio'
VMHC_flirt.inputs.cost_func = 'corratio'
VMHC_flirt.inputs.dof = 12
VMHC_flirt.inputs.interp = 'trilinear'

## Perform nonlinear registration (higres to standard) to symmetric standard brain
VMHC_fnt = pe.Node(interface=fsl.FNIRT(), name='VMHC_fnt')
#VMHC_fnt.inputs.in_file = os.path.abspath('${anat_dir}/head.nii.gz')
#VMHC_fnt.inputs.affine_file = os.path.abspath('${anat_reg_dir}/highres2symmstandard.mat')
VMHC_fnt.inputs.fieldcoeff_file = True
VMHC_fnt.inputs.jacobian_file = True
#VMHC_fnt.inputs.warped_file = os.path.abspath(${anat_reg_dir}/highres2symmstandard_warp.nii.gz)
#VMHC_fnt.inputs.jacobian_file = os.path.abspath(${anat_reg_dir}/highres2symmstandard_jac)
#VMHC_fnt.inputs.config_file = os.path.abspath(${T1_2_MNI152_2mm_symmetric} )
#VMHC_fnt.inputs.ref_file = os.path.abspath(${symm_standard})
#VMHC_fnt.inputs.refmask_file = os.path.abspath(${FSLDIR}/data/standard/MNI152_T1_2mm_brain_mask_symmetric_dil.nii.gz)
VMHC_fnt.inputs.warp_resolution = (10, 10, 10)

## Apply nonlinear registration (func to standard)
VMHC_warp = pe.MapNode(interface=fsl.ApplyWarp(), name='VMHC_warp', iterfield=["in_file", "premat"])
#VMHC_warp.inputs.in_file = os.path.abspath(${func_dir}/${image2use}.nii.gz)
#VMHC_warp.inputs.ref_file = os.path.abspath(${symm_standard})
#VMHC_warp.inputs.field_file = '${anat_reg_dir}/highres2symmstandard_warp.nii.gz'
#VMHC_warp.inputs.premat = os.path.abspath(${reg_dir}/example_func2highres.mat)
#VMHC_warp.inputs.out_file = os.path.abspath('${vmhc_dir}/${image2use}2symmstandard.nii.gz')


## copy and L/R swap file
VMHC_swap = pe.MapNode(interface=fsl.SwapDimensions(), name='VMHC_swap', iterfield=["in_file"])
#VMHC_swap.inputs.in_file = '${vmhc_dir}/${image2use}.nii.gz'
VMHC_swap.inputs.new_dims = ('-x', 'y', 'z')
VMHC_swap.inputs.out_file = 'tmp_LRflipped.nii.gz'

## caculate vmhc
VMHC_corr = pe.MapNode(interface=e_afni.ThreedTcorrelate(), name='VMHC_corr', iterfield=["xset"])
VMHC_corr.inputs.pearson = True
VMHC_corr.inputs.polort = -1
VMHC_corr.inputs.out_file = 'VMHC.nii.gz'
#VMHC_corr.inputs.xset = '${vmhc_dir}/${image2use}.nii.gz'
#VMHC_corr.inputs.yset = 'tmp_LRflipped.nii.gz'

VMHC_z_trans = pe.MapNode(interface=e_afni.Threedcalc(), name='VMHC_z_trans', iterfield=["infile_a"])
#VMHC_z_trans.inputs.infile_a = os.path.abspath('VMHC.nii.gz')
VMHC_z_trans.inputs.expr = '\'log((a+1)/(a-1))/2\''
VMHC_z_trans.inputs.out_file = 'VMHC_Z.nii.gz'

VMHC_z_stat = pe.MapNode(interface=e_afni.Threedcalc(), name='VMHC_z_stat', iterfield=["infile_a", "expr"])
VMHC_z_stat.inputs.out_file = 'VMHC_z_stat.nii.gz'
