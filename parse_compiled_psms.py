from itertools import islice

#infile = "Compiled\M_fortuitum_specific_rawfile_cluster_psm_rawinfo_mapped.txt"

def get_header_idx(infile):
    with open(infile) as file:
        for i in islice(file, 0, 1):
            split_i = i.rstrip().split('\t')
            try:
                peptide = split_i.index('Peptide')
                rawfile = split_i.index('Raw File')
                scan = split_i.index('Scan')
                mz = split_i.index('m/z')
                z = split_i.index('Charge(z)')
                rt = split_i.index('RT (min)')

            except:
                peptide = split_i.index('"Peptide"')
                rawfile = split_i.index('"Raw File"')
                scan = split_i.index('"Scan"')
                mz = split_i.index('"m/z"')
                z = split_i.index('"Charge(z)"')
                rt = split_i.index('"RT (min)"')
                
            return peptide, rawfile, scan, z, rt, split_i

def psm_rawinfo(infile):
    a = get_header_idx(infile)
    dicts = {}
    with open(infile) as file:
        for i in islice(file, 1, None):
            split_i = i.rstrip().split('\t')
            pep = split_i[a[0]]
            rawfile = split_i[a[1]]
            scan = split_i[a[2]]
            z = split_i[a[3]]
            rt = split_i[a[4]]
            
            if pep != 'NA':
                if pep not in dicts:
                    dicts[pep] = [rawfile + '@' + scan + '@' + z + '@' + rt]
                else:
                    dicts[pep].append(rawfile + '@' + scan + '@' + z + '@' + rt)

    return dicts

def rawinfo(infile):
    a = get_header_idx(infile)
    dicts = {}
    with open(infile) as file:
        for i in islice(file, 1, None):
            split_i = i.rstrip().split('\t')
            rawfile = split_i[a[1]]
            scan = split_i[a[2]]
            
            if rawfile + '@' + scan not in dicts:
                dicts[rawfile + '@' + scan] = [split_i]
            else:
                dicts[rawfile + '@' + scan].append(split_i)

    return dicts
