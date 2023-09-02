import parse_msp
import parse_compiled_psms
import numpy as np
import statsmodels.api as sm
import pylab
import seaborn as sns
import logging
logging.warning('Watch out!')  # will print a message to the console
logging.info('I told you so')

logger = logging.getLogger(__name__)

print (logger)

msp_file = "..\SSPP_Prosit_MSPs\M_vaccae_Peptides_Prosit_HCD.zip"

infile = "M_vaccae_PSMs_formatted.txt"

predicted_spectra = parse_msp.extract_msp(msp_file)

psms = parse_compiled_psms.psm_rawinfo(infile)

pep_raw_scan_mapped = {}
for predicted_pep in predicted_spectra:
    peptide = predicted_pep[0]
    z = predicted_pep[1]
    attrbs = predicted_pep[2]
    pred_mzs = predicted_pep[3]
    pred_intensities = predicted_pep[5]
    pred_ions = predicted_pep[4]

    #print (peptide, z, attrbs)

    irt = attrbs['iRT'][0]

    if peptide in psms:
        nr_raw_scans = {}
        for rawfile_scan in psms[peptide]:
            exp_rawfile = rawfile_scan.split('@')[0]
            exp_scan = rawfile_scan.split('@')[1]
            exp_z = rawfile_scan.split('@')[2]
            exp_rt = rawfile_scan.split('@')[3]
            nr_raw_scans[exp_rawfile + '@' + exp_scan] = [peptide + '@' + exp_z + '@' + exp_rt + '@' + str(irt)]

    for f, s in nr_raw_scans.items():
        rawf = f.split('@')[0]
        raw_scan = f.split('@')[1]
        for psm in s:
            pep_sequnece = psm.split('@')[0]
            pep_z = psm.split('@')[1]
            pep_rt = psm.split('@')[2]
            pep_irt = psm.split('@')[3] 
            if rawf not in pep_raw_scan_mapped:
                pep_raw_scan_mapped[rawf] = [[raw_scan + '@' + pep_sequnece + '@' + pep_z + '@' + pep_rt + '@' + pep_irt]]
            else:
                pep_raw_scan_mapped[rawf].append([raw_scan + '@' + pep_sequnece + '@' + pep_z + '@' + pep_rt + '@' + pep_irt])


def plot_scatter(file_name, x, y, loess_x, loess_y):

    sns.set_style("darkgrid")
    pylab.rc("figure", figsize=(16, 8))
    pylab.rc("font", size=14)

    #Plot the fit line
    fig, ax = pylab.subplots()

    #ax.scatter(x, y)
    ax.scatter(x, y, color='blue', s=2, linewidths=None, label='Predicted RT (min)')
    #ax.scatter(x, x,color= 'red', label='Observed RT (min)')
    ax.plot(loess_x, loess_y, c="k")
    
    pylab.autoscale(enable=True, axis="x", tight=True)
    pylab.title("Experimental vs Predicted Retention Time Plot", fontsize=25)
    pylab.xlabel('Experimental Retention Time (min)', fontsize=18)
    pylab.ylabel('Predicted Retention Time (min)', fontsize=18)
    pylab.savefig(file_name + 'RT_vs_iRT.pdf', dpi = 300, bbox_inches="tight", transparent=True)
    pylab.savefig(file_name + 'RT_vs_iRT.png', dpi = 300, bbox_inches="tight", transparent=True)
    pylab.show()

def loess_reg(peps, x, y):
   
    lowess = sm.nonparametric.lowess
    #z = lowess(y, x, frac= 1./3, it=0)
    lowess_smoothed = lowess(y, x, frac=1./4)
    
    #print (lowess_smoothed[:,0], lowess_smoothed[:,1])

    return lowess_smoothed[:,0], lowess_smoothed[:,1]

smoothed_irts = np.array([], dtype=float)
smoothed_rts = np.array([], dtype=float)

all_pep_RTs = np.array([], dtype=float)
all_pep_iRTs = np.array([], dtype=float)
all_pep_sequences = np.array([], dtype='object')

store_irt = {}
for rawfile, pep_scan_info in pep_raw_scan_mapped.items():
    peptide_RTs = np.array([], dtype=float)
    peptide_sequence = np.array([], dtype='object')
    peptide_iRTs = np.array([], dtype=float)
    for pep_scan in pep_scan_info:

        '''Storing the experimental retention time (min) and predicted iRT (min) to numpy array'''
        scan_num = pep_scan[0].split('@')[0]
        RT = float(pep_scan[0].split('@')[3])
        iRT = float(pep_scan[0].split('@')[4])       
        sequence = pep_scan[0].split('@')[1]

        peptide_RTs = np.append(peptide_RTs, RT)
        peptide_iRTs = np.append(peptide_iRTs, iRT)
        peptide_sequence = np.append(peptide_sequence, sequence)

        RT_diff = RT-iRT
        if rawfile + '@' + scan_num not in store_irt:
            store_irt[rawfile + '@' + scan_num] = [str(iRT) + '@' +  str(RT_diff) + '@' + str(abs(RT_diff))]
        else:
            store_irt[rawfile + '@' + scan_num].append(str(iRT) + '@' +  str(RT_diff) + '@' + str(abs(RT_diff)))

    ''' Perform the loess smoothing for experimental RT and predicted iRT values. Both the values will be provided in the numpy array format.'''
    smoothed_irt, smoothed_rt = loess_reg(peptide_sequence, peptide_RTs, peptide_iRTs)

    print (f'Loess smoothing of experimental and predicted iRT values were performed for {len(peptide_sequence)} PSMs of Raw File {rawfile}.')

    smoothed_irts = np.append(smoothed_irts, smoothed_irt)
    smoothed_rts = np.append(smoothed_rts, smoothed_rt)

    all_pep_RTs = np.append(all_pep_RTs, peptide_RTs)
    all_pep_iRTs = np.append(all_pep_iRTs, peptide_iRTs)
    all_pep_sequences = np.append(all_pep_sequences, peptide_sequence)

print (f'Completed the Loess smoothing of experimental retention time (min) and predicted iRT (min) of {len(all_pep_sequences)} PSMs.')

#print (len(all_pep_sequences), len(all_pep_RTs), len(all_pep_iRTs), len(smoothed_rts), len(smoothed_irts))

plot_scatter(infile, all_pep_RTs, all_pep_iRTs, smoothed_rts, smoothed_irts)


''' Map the experimental and predicted iRT differences of peptide spectrum matches to the input file'''
all_rawscans = parse_compiled_psms.rawinfo(infile)
output = []
for raw_scan, scan_info in all_rawscans.items():
    if raw_scan in store_irt:
        RT_info = store_irt[raw_scan][0].split('@')
        for info in scan_info:
            output.append(info + RT_info)
    else:
        for info in scan_info:
            output.append(info + ['No iRT', 'NA','NA'])

outfile = "{0}_RTinfo.txt".format(infile.rstrip('txt').rstrip('txt'))
header = parse_compiled_psms.get_header_idx(infile)[-1] + ["Predicted iRT (min)", "RT Difference (min)", "Absolute RT Difference (min)"]

with open(outfile, 'w') as outf:
    outf.write('\t'.join(header) + '\n')
    outf.writelines('\t'.join(i) + '\n' for i in output)
