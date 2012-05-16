import os
import sys
import re
import commands


def get_vectors(strat):

    paths = []
    def dfs(val_list, path):

        if val_list == []:
            paths.append(list(path.split(',')[1:]))

        else:
            vals = val_list.pop()
            for val in vals:
                dfs(list(val_list), path + ',' + val)




    key_list = []
    val_list = []
    for k, v in sorted(strat.items()):
        key_list.append(k)
        val_list.append(v)

    dfs(val_list, '')

    return paths

def get_strategies(K, iterables):

    p_s_m = {}
    sm = {}
    strat = {}
    for key in K:

        vals = key.split('/')

        vals.remove('')
        for val in vals:

            if not val in p_s_m:
                p_s_m[val] = 1

    k = sorted(p_s_m.keys())
    for val in k:

        match = 0
        for variable in iterables:

            if variable in val and val.split(variable)[0] == '':

                match = 1

                it = val.split(variable)[1]
                if variable in strat:
                    item = strat[variable]
                    item.append(val)
                else:
                    strat[variable] = [val]

        if match == 0:
            strat[val] = [val]


    strat_vect = get_vectors(strat)
    return strat_vect

def get_match(sink_dir, p, strategies):

    p1 = p.split(sink_dir)[1:]
    p1 = p1[0]
    p_vals = p1.split('/')[2:]

    print p_vals
    flag = 0
    strats = []
    for strat in strategies:
        flag = 0
        for val in p_vals:

            if not val in strat:
                flag = 1
                break
        if flag == 0:
            strats.append(str(strat))

    if strats == []:
        print "WARNING : No strategy found for ", p
        return None

    return strats

def smash(sink_dir, maps, strategies):

    new_maps = {}

    for p in sorted(maps.keys()):

        strats = get_match(sink_dir, p, strategies)


        for strat in strats:
            if strat in new_maps:
                list1 = list(new_maps[strat])
                new_maps[strat] = sorted(list(set(maps[p] + list1)))
            else:
                new_maps[strat] = sorted(list(maps[p]))

    return dict(new_maps)

def make_links(new_maps, sink_dir, sub):

    sym_path = os.path.join(sink_dir, 'sym_links')
    cmd = 'mkdir -p %s' % (sym_path)
    print cmd
    commands.getoutput(cmd)

    idx = 0
    f = open('%s/label_linkage.txt' % (sym_path), 'w')

    labels = {}
    strats = sorted(new_maps.keys())
    orig_sub_dir = os.path.join(sink_dir, sub)
    wfs = sorted(os.listdir(orig_sub_dir))
    sessions = os.listdir(os.path.join(orig_sub_dir, wfs[1]))

    for strat in strats:
        idx += 1
        print >>f, 'label_'+str(idx), ' ', strat
        labels[strat] = 'label_'+str(idx)
        l_path = os.path.join(sym_path, 'label_'+str(idx))

        cmd = 'mkdir -p %s' % (l_path)
        print cmd
        commands.getoutput(cmd)

        sub_path = os.path.join(l_path, sub)
        cmd = 'mkdir -p %s' % (sub_path)
        print cmd
        commands.getoutput(cmd)


        for wf in wfs:
            wf_path = os.path.join(sub_path, wf)
            cmd = 'mkdir -p %s' % (wf_path)
            print cmd
            commands.getoutput(cmd)

            for file in new_maps[strat]:

                if '/' + wf + '/' in file:
                    cmd = 'ln -s %s %s' %(file, os.path.join(wf_path, os.path.basename(file)))
                    print cmd
                    commands.getoutput(cmd)

    f.close()


def make_sym_links(sink_dir, strategies, subj_list):


    for subject in subj_list:
        maps = {}

        dfs_string = commands.getoutput("find %s/ -type f | perl -lne 'print tr:/::, \" $_\"' | sort -n | cut -d' ' -f2" % (os.path.join(sink_dir, subject)))

        dfs_files = dfs_string.split('\n')

        for file in dfs_files:
            pdir = os.path.dirname(file)

            val = None
            if pdir in maps:
                val = maps[pdir]
                val.append(file)
            else:
                val = [file]
            maps[os.path.dirname(file)] = val

        new_maps = smash(sink_dir, maps, strategies)

        make_links(new_maps, sink_dir, subject)

def main():

    sink_dir = '/home/sharad/nki_nyu_pipeline/results/'
    subj_list = ['s1001', 's1004']
    subjectDir = '/home/sharad/nki_nyu_pipeline/results/%s' % (subj_list[0])

    iterables = ['threshold', 'csf_threshold', 'fwhm', 'gm_threshold', 'hp', 'lp', 'nc', 'run_scrubbing', 'seeds', 'selector', 'session_id', 'target_angle', 'wm_threshold']
    wfs = os.listdir(subjectDir)

    labels = {}
    label_id = 0

    dfs_string = commands.getoutput("find %s/ -type f | perl -lne 'print tr:/::, \" $_\"' | sort -n | cut -d' ' -f2" % (subjectDir))

    dfs_files = dfs_string.split('\n')

    maps = {}

    m1 = {}
    for file in dfs_files:
        f = str(file)
        spath = os.path.join(sink_dir, subj_list[0])
        vmhc_path = os.path.join(spath, 'vmhc')
        alff_path = os.path.join(spath, 'alff')
        nuisance_path = os.path.join(spath, 'nuisance')
        sca_path = os.path.join(spath, 'sca')
        scrubbing_path = os.path.join(spath, 'scrubbing')
        anat_path = os.path.join(spath, 'anat')
        func_path = os.path.join(spath, 'func')
        seg_path = os.path.join(spath, 'segment')
        reg_path = os.path.join(spath, 'reg')

        file = re.sub(r'%s' % (vmhc_path), '', file)
        file = re.sub(r'%s' % (alff_path), '', file)
        file = re.sub(r'%s' % (nuisance_path), '', file)
        file = re.sub(r'%s' % (sca_path), '', file)
        file = re.sub(r'%s' % (scrubbing_path), '', file)
        file = re.sub(r'%s' % (anat_path), '', file)
        file = re.sub(r'%s' % (func_path), '', file)
        file = re.sub(r'%s' % (reg_path), '', file)
        file = re.sub(r'%s' % (seg_path), '', file)
        #print file
        pdir = os.path.dirname(file)

        val = None
        v = None
        if pdir in maps:
            val = maps[pdir]
            v = m1[pdir]
            val.append(file)
            v.append(f)
        else:
            val = [file]
            v = [f]
        maps[os.path.dirname(file)] = val
        m1[os.path.dirname(file)] = v

    strategies = list(get_strategies(maps.keys(), iterables))

    print strategies

    make_sym_links(sink_dir, strategies, subj_list)


if __name__ == "__main__":

    sys.exit(main())
