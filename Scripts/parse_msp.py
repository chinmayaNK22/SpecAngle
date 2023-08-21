import numpy as np
from itertools import groupby

def sort_peptides(infile):
    with open(infile) as file:
            
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

        mod_strings = 0
        
        attributes = {"m/z":[],"MODS":[],"Mod_Peptide":[],"iRT":[],"proteotypicity":[]}
        mz_values = np.array([], dtype=float)
        mz_ions = np.array([], dtype='object')
        mz_intense = np.array([], dtype='d')
        
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
                    mod_strings += 1
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

        print (f'INFO: Found a {mod_strings} modified peptide precursors in the library.')
        
        msms_spec = {}
        msms_ions = {}
        for idx, mz in np.ndenumerate(mz_values):
            
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
            
            sorted_mz_values = np.append(mz_values, mz)
            sorted_mz_ions = np.append(mz_ions, msms_ions[mz][0])
            sorted_mz_intense = np.append(mz_intense, msms_spec[mz][0])
                    

        yield sequence, z, attributes, sorted_mz_values, sorted_mz_ions, sorted_mz_intense

    
