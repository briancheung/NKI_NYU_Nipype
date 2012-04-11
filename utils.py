# Utility Functions ---------------------------------------------------------
import e_afni
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util

def getStartIdx(in_file):
       print "getStart Idx in_file -----------> ", in_file

       start_indx = []

       if(isinstance(in_file, list)):
               for file in in_file:
                       f = open(file, 'r')
                       line = f.readline()
                       start_indx.append(int(line.split(",")[0]))
                       f.close()
               print "start_indx ", start_indx
               return start_indx
       else:
               f = open(in_file, 'r')
               myList = []
               line = f.readline()
               myList.extend(line.split(","))
               f.close()
               return int (myList[0])

def getStopIdx(in_file):

       print "getStop Idx in_file ------------> ", in_file

       stop_indx = []
       if (isinstance(in_file, list)):
               for file in in_file:
                       f = open(file, 'r')
                       line = f.readline()
                       myList = []
                       myList.extend(line.split(","))
                       length = len(myList)
                       f.close()

                       if str(myList[length-1]) == "":
                               stop_indx.append(int(myList[length-2]))
                       else:
                               stop_indx.append(int(myList[length-1]))
               print "stop_indx ", stop_indx
               return stop_indx
       else:
               f = open(in_file, 'r')
               line = f.readline()
               myList = []
               myList.extend(line.split(","))
               length = len(myList)
               f.close()
               if str(myList[length-1]) == "":
                       return int(myList[length-2])
               else:
                       return int(myList[length-1])

def last_vol(vols):

    v = []
    for vol in vols:
        v.append(int(vol) - 1)

    return v

def TRendminus1(vols):
    v = []
    for vol in vols:
        v.append(int(vol) - 2)

    return v


def scCopy(in_file):

    import os

    cwd = os.getcwd()
    print "****************************** l -->", cwd

    out_file = 'pow_params.txt'
    path = os.path.join(cwd, out_file)
    print "##########path --->", path
    out_file = path

    dir_name = os.path.dirname(in_file)
    rest_name = os.path.basename(dir_name)
    subject_id = os.path.basename(os.path.dirname(dir_name))

    f = open(out_file, 'a')

    f.write("%s," % subject_id)
    f.write("%s," % rest_name)

    f.close()

    return out_file

def createSC(in_file):
    import subprocess as sb
    import os

    out_file = os.path.join(os.getcwd(), 'FD.1D')


    print "outfile ", out_file

    cmd1 = sb.Popen(['awk', '{x=$4} {y=$5} {z=$6} {a=$1} {b=$2} {c=$3} '+
                    '{print 2*3.142*50*(a/360),2*3.142*50*(b/360), 2*3.142*50*(c/360), x, y, z}',
                    in_file], stdin=sb.PIPE, stdout=sb.PIPE,)
    cmd2 = sb.Popen(['awk', '{a=$1} {b=$2} {c=$3} {x=$4} {y=$5} {z=$6} '+
                    'NR>1{print a-d, b-e, c-f, x-u, y-v, z-w}{d=a} {e=b} {f=c} {u=x} {v=y} {w=z}'],
                    stdin=cmd1.stdout, stdout=sb.PIPE,)
    cmd3 = sb.Popen(['awk', '{ for (i=1; i<=NF; i=i+1) {if ($i < 0) $i = -$i} print}'],
                    stdin=cmd2.stdout, stdout=sb.PIPE,)
    cmd4 = sb.Popen(['awk', '{a=$1+$2+$3+$4+$5+$6} {print a}'],
                    stdin=cmd3.stdout, stdout=sb.PIPE,)


    stdout_value, stderr_value = cmd4.communicate()

    output = stdout_value.split()
    f = open(out_file, 'w')

    for out in output:
      print >>f, float(out)

    f.close()

    print "stdout in createSC --> ", stdout_value
    print "stderr in createSC ---> ", stderr_value

    return out_file


def setMeanFD(infile_a, infile_b):

    import subprocess as sb
    import os

    out_file = os.path.join(os.getcwd(), 'pow_params.txt')
    copycmd = sb.Popen(['cp', infile_a, out_file], stdin=sb.PIPE, stdout=sb.PIPE )
    out_val, error_val = copycmd.communicate()
    print "outval ---> ", out_val
    print "error_val -->", error_val


    cmd = sb.Popen(['awk', '{a += $1} END {print a/NR}', infile_b],
                   stdin=sb.PIPE, stdout=sb.PIPE,)
    stdout_value, stderr_value = cmd.communicate()

    output = stdout_value.split()
    f = open(out_file, 'a')

    for out in output:
        f.write('%.4f,' % float(out))

    f.close()

    print "stdout in setMeanFD --> ", stdout_value
    print "stderr in setMeanFD---> ", stderr_value

    return out_file

def setNumFD(infile_a, infile_b):

    import subprocess as sb
    import os

    out_file = os.path.join(os.getcwd(), 'pow_params.txt')

    copycmd = sb.Popen(['cp', infile_a, out_file], stdin=sb.PIPE, stdout=sb.PIPE )
    out_val, error_val = copycmd.communicate()
    print "outval ---> ", out_val
    print "error_val -->", error_val

    cmd = sb.Popen(['awk', '{ if($1>=0.5) {a += 1}} END {print a}',
                    infile_b], stdin=sb.PIPE, stdout=sb.PIPE,)

    stdout_value, stderr_value = cmd.communicate()
    output = stdout_value.split()

    f = open(out_file, 'a')

    for out in output:
        f.write('%.4f,' % float(out))

    f.close()

    print "stdout in setNumFD --> ", stdout_value
    print "stderr in setNumFD---> ", stderr_value

    return out_file

def setPercentFD(infile_a, infile_b):

    import subprocess as sb
    import os

    out_file = os.path.join(os.getcwd(), 'pow_params.txt')

    copycmd = sb.Popen(['cp', infile_a, out_file], stdin=sb.PIPE, stdout=sb.PIPE )
    out_val, error_val = copycmd.communicate()
    print "outval ---> ", out_val
    print "error_val -->", error_val

    cmd = sb.Popen(['awk', '{ if($1>=0.5) {a += 1}} END {print (a/(NR+1)*100)}',
                   infile_b], stdin=sb.PIPE, stdout=sb.PIPE,)
    stdout_value, stderr_value = cmd.communicate()

    output = stdout_value.split()
    f = open(out_file, 'a')

    for out in output:
       f.write('%.4f,' % float(out))

    f.close()

    print "stdout in setPercentFD --> ", stdout_value
    print "stderr in setPercentFD---> ", stderr_value

    return out_file

def setFramesEx(in_file):

    import subprocess as sb
    import os

    out_file = os.path.join(os.getcwd(), 'frames_ex.1D')

    cmd = sb.Popen(['awk', '{ if($1>=0.5) {print NR}}', in_file],
                   stdin=sb.PIPE, stdout=sb.PIPE,)

    stdout_value, stderr_value = cmd.communicate()

    output = stdout_value.split()
    f = open(out_file, 'a')

    for out in output:
       f.write('%s,' % int(out))

    f.close()

    print "stdout in setFramesEx --> ", stdout_value
    print "stderr in setFramesEx---> ", stderr_value

    return out_file

def setFramesIN(in_file):

    import subprocess as sb
    import os

    out_file = os.path.join(os.getcwd(), 'frames_in.1D')

    cmd = sb.Popen(['awk', '{ if($1<0.5) {print NR}}', in_file],
                   stdin=sb.PIPE, stdout=sb.PIPE,)

    stdout_value, stderr_value = cmd.communicate()

    output = stdout_value.split()
    f = open(out_file, 'a')

    for out in output:
       f.write('%s,' % int(out))

    f.close()

    print "stdout in setFramesIn --> ", stdout_value
    print "stderr in setFramesIn---> ", stderr_value

    return out_file


def setFramesInList(in_file):

    import subprocess as sb
    import os

    out_file = os.path.join(os.getcwd(), 'frames_in.1D')

    cmd = sb.Popen(['awk', '{ if($1<0.5) {print NR}}', in_file],
                   stdin=sb.PIPE, stdout=sb.PIPE,)

    stdout_value, stderr_value = cmd.communicate()

    output = stdout_value.split()
    f = open(out_file, 'a')

    for out in output:
       print >>f, int(out)

    f.close()

    print "stdout in setFramesInList --> ", stdout_value
    print "stderr in setFramesInList---> ", stderr_value

    return out_file

def setMeanDVARS(infile_a, infile_b):

    import subprocess as sb
    import os

    out_file = os.path.join(os.getcwd(), 'pow_params.txt')
    copycmd = sb.Popen(['cp', infile_a, out_file], stdin=sb.PIPE, stdout=sb.PIPE )
    out_val, error_val = copycmd.communicate()
    print "outval ---> ", out_val
    print "error_val -->", error_val

    cmd = sb.Popen(['awk', '{a += $1} END {print a/NR}',
                    infile_b], stdin=sb.PIPE, stdout=sb.PIPE,)

    stdout_value, stderr_value = cmd.communicate()

    output = stdout_value.split()
    f = open(out_file, 'a')

    for out in output:
      f.write('%.4f,' % float(out))

    f.close()

    print "stdout in setMeanDVARS --> ", stdout_value
    print "stderr in setMeanDVARS---> ", stderr_value

    return out_file

def setNUM5(infile_a, infile_b):

    import subprocess as sb
    import os

    out_file = os.path.join(os.getcwd(), 'pow_params.txt')
    copycmd = sb.Popen(['cp', infile_a, out_file], stdin=sb.PIPE, stdout=sb.PIPE )
    out_val, error_val = copycmd.communicate()
    print "outval ---> ", out_val
    print "error_val -->", error_val

    cmd = sb.Popen(['awk', '{ if($1>=5) {a += 1}} END {print a}',
                    infile_b], stdin=sb.PIPE, stdout=sb.PIPE,)

    stdout_value, stderr_value = cmd.communicate()

    output = stdout_value.split()
    f = open(out_file, 'a')

    for out in output:
      f.write('%.4f,' % float(out))

    f.close()

    print "stdout in setNUM5 --> ", stdout_value
    print "stderr in setNUM5---> ", stderr_value

    return out_file

def setNUM10(infile_a, infile_b):

    import subprocess as sb
    import os

    out_file = os.path.join(os.getcwd(), 'pow_params.txt')
    copycmd = sb.Popen(['cp', infile_a, out_file], stdin=sb.PIPE, stdout=sb.PIPE )
    out_val, error_val = copycmd.communicate()
    print "outval ---> ", out_val
    print "error_val -->", error_val

    cmd = sb.Popen(['awk', '{ if($1>=10) {a += 1}} END {print a}',
                    infile_b], stdin=sb.PIPE, stdout=sb.PIPE,)

    stdout_value, stderr_value = cmd.communicate()

    output = stdout_value.split()
    f = open(out_file, 'a')

    for out in output:
      f.write('%.4f,' % float(out))

    f.close()

    print "stdout in setNUM5 --> ", stdout_value
    print "stderr in setNUM5---> ", stderr_value

    return out_file


def setSqrtMeanDeriv (in_file):
    import subprocess as sb
    import os

    out_file = os.path.join(os.getcwd(), 'sqrt_mean_deriv_sq.1D')

    cmd = sb.Popen(['awk', '{x=$1} {printf( "%.6f ", sqrt(x))}',
                    in_file], stdin=sb.PIPE, stdout=sb.PIPE,)

    stdout_value, stderr_value = cmd.communicate()

    output = stdout_value.split()
    f = open(out_file, 'a')

    for out in output:
      print >>f, float(out)

    f.close()

    print "stdout in setSqrtMeanDeriv --> ", stdout_value
    print "stderr in setSqrtMeanDeriv---> ", stderr_value

    return out_file

def setSqrtMeanRaw(in_file):

    import subprocess as sb
    import os

    out_file = os.path.join(os.getcwd(), 'sqrt_mean_raw_sq.1D')

    cmd = sb.Popen(['awk', '{x=$1} {printf( "%.6f ", sqrt(x))}',
                    in_file], stdin=sb.PIPE, stdout=sb.PIPE,)

    stdout_value, stderr_value = cmd.communicate()

    output = stdout_value.split()
    f = open(out_file, 'a')

    for out in output:
      print >>f, float(out)

    f.close()

    print "stdout in setSqrtMeanRaw --> ", stdout_value
    print "stderr in setSqrtMeanRaw---> ", stderr_value

    return out_file


def setFtoFPercentChange(infile_a, infile_b):

    import subprocess as sb
    import os

    out_file = os.path.join(os.getcwd(), 'ftof_percent_change.1D')

    cmd1 = sb.Popen(['awk', 'NR==FNR{a[NR]=$1; next} {print a[FNR], $1}',
                     infile_a, infile_b], stdin=sb.PIPE, stdout=sb.PIPE,)

    cmd2 = sb.Popen(['awk', '{x=$1} {y=$2} {printf( "%.6f ", ((x/y)*100))}'],
                    stdin=cmd1.stdout, stdout=sb.PIPE,)

    stdout_value, stderr_value = cmd2.communicate()

    output = stdout_value.split()
    f = open(out_file, 'a')

    for out in output:
      print >>f, float(out)

    f.close()

    print "stdout in ftofPercentChange --> ", stdout_value
    print "stderr in ftofPercentChange---> ", stderr_value

    return out_file



def setNUMFD(in_file):

    import subprocess as sb
    import os

    out_file = os.path.join(os.getcwd(), 'numFD')

    cmd = sb.Popen(['awk', '{ if($1>=0.5) {a += 1}} END {print a}',
                   in_file], stdin=sb.PIPE, stdout=sb.PIPE,)

    stdout_value, stderr_value = cmd.communicate()

    output = stdout_value.split()
    f = open(out_file, 'a')

    for out in output:
      f.write('%.4f,' % float(out))

    f.close()

    print "stdout in setNUMFD --> ", stdout_value
    print "stderr in setNUMFD---> ", stderr_value

    return out_file

def setScrubbedMotion(infile_a, infile_b):

    import subprocess as sb
    import os

    out_file = os.path.join(os.getcwd(), 'rest_mc_scrubbed.1D')

    cmd = sb.Popen(['awk', 'FNR==NR{a[$1];next}(FNR in a){print}',
                   infile_a, infile_b], stdin=sb.PIPE, stdout=sb.PIPE,)

    stdout_value, stderr_value = cmd.communicate()

    f = open(out_file, 'a')

    f.write(stdout_value)

    f.close()

    print "stdout in setScrubbedMotion --> ", stdout_value
    print "stderr in setScrubbedMotion---> ", stderr_value

    return out_file


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


def select(pp, run_scrubbing, sc_pp):

    pp_s = None
    if (run_scrubbing):

        pp_s = sc_pp

    else:

        pp_s = pp

    return pp_s


def selector_wf():

    swf = pe.Workflow(name='scrubbing_selector')

    inputNode = pe.Node(util.IdentityInterface(fields=['preprocessed',
                                                       'scrubbed_preprocessed'
                                                      ]),
                        name='inputspec')

    iRun = pe.Node(util.IdentityInterface(fields=['run_scrubbing']),
                             name='run_scrubbing_input')


    outputNode = pe.Node(util.IdentityInterface(fields=['preprocessed_selector'
                                                       ]),
                        name='outputspec')


    selectP = pe.MapNode(util.Function(input_names=['pp',
                                                  'run_scrubbing',
                                                  'sc_pp'],
                                      output_names=['pp_s'], function=select),
                        name='selectP', iterfield=['pp', 'sc_pp'])

    swf.connect(inputNode, 'preprocessed', selectP, 'pp')
    swf.connect(inputNode, 'scrubbed_preprocessed', selectP, 'sc_pp')
    swf.connect(iRun, 'run_scrubbing', selectP, 'run_scrubbing')

    swf.connect(selectP, 'pp_s', outputNode, 'preprocessed_selector')


    return swf




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
