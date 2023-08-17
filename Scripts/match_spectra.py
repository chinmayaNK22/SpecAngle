import numpy as np


def mz_range(tolerance_val, mz_val, unit):
    mz_range = 0
    if unit.lower() == 'da':
        mz_range += tolerance_val
    elif unit.lower() == 'ppm':
        mz_range += (tolerance_val/10**6)*mz_val

    return mz_range
        

def peak_match(ex_frag_mz, theo_mz, ex_frag_intensity, theo_frag_intensity, ion_type, tolerance, unit):

    matching_values = np.array([], dtype=float)
    matching_ions = np.array([], dtype='object')
    matching_exp_intense = np.array([], dtype='f')
    matching_pred_intense = np.array([], dtype='f')
    
    # Loop through the values in thoeretical mz array and identify the matching experimental mz values within given fragment tolerance
    for index, value in np.ndenumerate(theo_mz):
        if np.any(np.abs(ex_frag_mz - value) <= float(mz_range(tolerance, value, unit))):
            exp_mz = ex_frag_mz[np.where(np.abs(ex_frag_mz - value) <= float(mz_range(tolerance, value, unit)))]

            #print (exp_mz)
            
            if len(exp_mz) > 1:
                close_exp_idx = np.abs(exp_mz - value).argmin()
                new_exp_mz = exp_mz[close_exp_idx]
                if len(np.where(matching_values == new_exp_mz)[0]) == 0: ### Exclude the mz values which already present the the final numpy array
                    matching_values = np.append(matching_values, new_exp_mz)
                    matching_ions = np.append(matching_ions, ion_type[index])
                    matching_exp_intense = np.append(matching_exp_intense, ex_frag_intensity[close_exp_idx])
                    matching_pred_intense = np.append(matching_pred_intense, theo_frag_intensity[index])
                    #if close_exp_idx == len(ex_frag_intensity)-1:
                    #    matching_intense = np.append(matching_intense, ex_frag_intensity[-1])
                    #
                    #else:
                    #    matching_intense = np.append(matching_intense, ex_frag_intensity[close_exp_idx])
                        
            else:
                exp_mz_idx = np.where(np.abs(ex_frag_mz - value) <= float(mz_range(tolerance, value, unit)))[0]
                if len(np.where(matching_values == exp_mz)[0]) == 0: ### Exclude the mz values which already present the the final numpy array
                    matching_values = np.append(matching_values, exp_mz)
                    matching_ions = np.append(matching_ions, ion_type[index])
                    #print (exp_mz_idx, exp_mz, value, ion_type[index], len(ex_frag_intensity), len(ex_frag_mz), ex_frag_mz, ex_frag_intensity[exp_mz_idx], ex_frag_intensity)

                    #if exp_mz_idx == len(ex_frag_intensity)-1 or exp_mz_idx == len(ex_frag_intensity):
                    #    matching_intense = np.append(matching_intense, ex_frag_intensity[-1])
                    #else:
                    
                    matching_exp_intense = np.append(matching_exp_intense, ex_frag_intensity[exp_mz_idx])
                    matching_pred_intense = np.append(matching_pred_intense, theo_frag_intensity[index])
                

    return matching_values, matching_exp_intense, matching_pred_intense, matching_ions
