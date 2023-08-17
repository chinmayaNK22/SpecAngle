import numpy as np
from itertools import groupby

#infile = "M_fortuitum_SSPPs_SpecLib.msp"

def sort_peptides(infile):
    with open(infile) as file:
        
            
        #split_i = i.rstrip()
            
        peptide_pred_spec = (pred_spec[1] for pred_spec in groupby(file, lambda x: x.startswith('Name:')))

        for spec in peptide_pred_spec:
            peptide = next(spec)
            info = [s.rstrip() for s in next(peptide_pred_spec)]

            yield peptide, info


def extract_msp(infile):
    msp_file = ""
    if infile.split('.')[-1] == 'msp':
        msp_file = infile
    else:
        raise Exception("Input file is not in MSP format")
    
    predicted_peptide_info = sort_peptides(msp_file)
    for peptide in predicted_peptide_info:
        pep = peptide[0]
        info = peptide[1]

        sequence = pep.split(':')[1].strip().split('/')[0]
        z = pep.split(':')[1].strip().split('/')[1]

        attributes = {"m/z":[],"MODS":[],"Mod_Peptide":[],"iRT":[],"proteotypicity":[]}
        mz_values = np.array([], dtype=float)
        mz_ions = np.array([], dtype='object')
        mz_intense = np.array([], dtype='f')
        for i in info:
            if i.rstrip().startswith('Comment:'):
                    
                attrs = i.rstrip().split(':')[1].strip().split(' ')

                try:
                    attributes["m/z"].append(float(attrs[0].split('=')[1]))
                    attributes["MODS"].append(attrs[2].split('=')[1])
                    attributes["Mod_Peptide"].append(attrs[3].split('=')[1].split('/')[0])
                    attributes["iRT"].append(float(attrs[4].split('=')[1]))
                    attributes["proteotypicity"].append(float(attrs[5].split('=')[1]))
                    
                except:
                    print (f'Found a modified peptide precursor {pep.rstrip()}')
                    for attribute in attrs:
                        if 'iRT' in attribute:
                            attributes["iRT"].append(float(attribute.split('=')[1]))
                        if 'proteotypicity' in attribute:
                            attributes["proteotypicity"].append(float(attribute.split('=')[1]))
                
            if i.rstrip()[0].isdigit():
            
                spectra = i.rstrip().split('\t')

                ions = spectra[2].strip('"').split('/')[0]
            
                mz_values = np.append(mz_values, float(spectra[0]))
                mz_ions = np.append(mz_ions, ions)
                mz_intense = np.append(mz_intense, float(spectra[1]))

        msms_spec = {}
        msms_ions = {}
        for idx, mz in np.ndenumerate(mz_values):
            #print (sequence, idx, mz, mz_ions[idx])
            if mz not in msms_spec:
                msms_ions[mz] = [mz_ions[idx]]
                msms_spec[mz] = [mz_intense[idx]]
            else:
                msms_ions[mz].append(mz_ions[idx])
                msms_spec[mz].append(mz_intense[idx])

        sorted_mz_values = np.array([], dtype=float)
        sorted_mz_ions = np.array([], dtype='object')
        sorted_mz_intense = np.array([], dtype='f')
        for mz in sorted(msms_ions):
            #print (mz, msms_ions[mz], msms_spec[mz])
            sorted_mz_values = np.append(mz_values, mz)
            sorted_mz_ions = np.append(mz_ions, msms_ions[mz][0])
            sorted_mz_intense = np.append(mz_intense, msms_spec[mz][0])
                    

        yield sequence, z, attributes, sorted_mz_values, sorted_mz_ions, sorted_mz_intense

    
##if __name__== "__main__":
##    if infile.split('.')[-1] == 'msp':
##        dissociation_msp_info(infile)
##    else:
##        raise Exception("Input file is not in MSP format")
    
