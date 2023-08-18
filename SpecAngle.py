import numpy as np
from pyteomics import mzml
import math
from Scripts import parse_msp
from Scripts import parse_psm_scan
import os
from Scripts import match_spectra
import argparse

parser = argparse.ArgumentParser(description='''Calculates the Normalized Spectral Angle Score for Pepetide Spectrum Matches by comparing to the Prosit predicted library''')

parser.add_argument('infiles', metavar='-i', type=str, nargs='+', help='Input file consisting peptide spectrum matches with raw file and scan details in tab delimitted format')

parser.add_argument('msp_file', metavar='-msp', type=str, nargs='+', help='Prosit predicted peptide spectral library in MSP format')

parser.add_argument('raw_path', metavar='-raw', type=str, nargs='+', help='Path to a location where DDA derived MS/MS data files in mzML format')

args = parser.parse_args()

#sspp_msp_file = "../SSPP_Prosit_MSPs/M_abscessus_SSPPs_SpecLib.msp"

#predicted_pep_spec = parse_msp.extract_msp(sspp_msp_file)

#rawfile_path = "."

def store_rawfiles(raw_path):
    rawfiles = {}
    for files in os.listdir(raw_path):
        if os.path.isfile(os.path.join(raw_path, files)):
            if files.split('.')[-1] == 'mzML':
                if files not in rawfiles:
                    rawfiles[files] = [os.path.join(raw_path, files)]

    return rawfiles

#rawscan_psm = ["../PSM_RawScan_Mapped/M_abscessus_In_house_specific_rawfile_cluster_psm_rawinfo_mapped.txt", "../PSM_RawScan_Mapped/M_abscessus_specific_rawfile_cluster_psm_rawinfo_mapped.txt"]

#psms = parse_psm_scan.psm_rawinfo(rawscan_psm)

def sqrt_norm(ms2_intensities):

    # Square root transformation of intensities
    sqrt_intensity = np.sqrt(ms2_intensities)

    norm_intensity = np.linalg.norm(sqrt_intensity)

    # Normalize the vector to unit length
    normalized_vector = sqrt_intensity / norm_intensity
    
    return normalized_vector

def dot_product(vector1, vector2):
    # Ensure both vectors have the same length
    if len(vector1) != len(vector2):
        raise ValueError("Both vectors must have the same length.")

    # Compute the dot product using numpy's dot() function

    dot_product_result = np.dot(vector1, vector2)

    return dot_product_result

def spec_contrast_angle(dotp):
    spec_angle = 0
    if abs(dotp) <= 1:
        spec_angle += 1 - (2 * math.acos(dotp)/math.pi)
    else:
        spec_angle = 0

    return spec_angle

def process_spectrum(sequence, scan_num, exp_spectra, pred_mzs, pred_intensities, pred_ions):
    #print (scan_num, info['m/z array'], info['intensity array'])
    spec_mz = exp_spectra['m/z array']
    spec_intensity = exp_spectra['intensity array']
    
    ''' Normalize the intensity values of all MS/MS peaks from experimental spectra'''
    norm_exp_intensity = sqrt_norm(spec_intensity)
    #norm_pred_intensity = sqrt_norm(pred_intensities)

    tolerance = 20.0 ## Fragment mass tolerance for matching m/z values from experimental and Prosit predicted spectra
    tolerance_unit = 'ppm' ## Fragment mass tolerance unit either in 'parts per million' (ppm) or in Dalton (Da)

    ''' Perform m/z peak match of experimental raw spectra (MS/MS) of a peptide sequence with its respective Prosit predicted fragment spectra.
    The peak match will be performed with the fragmenta match tolerance either in 'ppm' or 'Da'.'''

    matched_values, matched_exp_intense, matched_pred_intense, matched_ions = match_spectra.peak_match(spec_mz, pred_mzs, norm_exp_intensity, pred_intensities, pred_ions, tolerance, tolerance_unit)

    #exp_spec_dot = dot_product(spec_mz, norm_intensity)
    #norm_exp_intensity = sqrt_norm(matched_exp_intense)
    #prosit_spec_dot = dot_product(pred_mzs, pred_intensities)
    #print (matched_values, matched_exp_intense, matched_pred_intense, matched_ions, norm_exp_intensity)
    
    dotp = dot_product(matched_pred_intense, matched_exp_intense)
    #print (dotp, matched_pred_intense, matched_ions, norm_exp_intensity) 
    
    norm_spec_angle = spec_contrast_angle(dotp)
    
    if norm_spec_angle == 0:
        print (f"Normalized spectral contrast could not be calculated for the raw spectra {scan_num} of {rawfile} file for {sequence}. Returning Normalized spectral contrast angle score as 0.")
        return scan_num, dotp, norm_spec_angle
    else:
        return scan_num, dotp, norm_spec_angle

def list_infiles(infiles):
    files = []
    split_i = infiles.split(' ')

    print (split_i)
    for i in split_i:
        if len(i) != 0:
            print (i)

            files.append(i)

    return files

def SpecAngle_Calc(infiles, msp_file, raw_path):

    psms = parse_psm_scan.psm_rawinfo(list_infiles(infiles))

    predicted_pep_spec = parse_msp.extract_msp(msp_file)

    pep_raw_scan_mapped = {}
    for predicted_pep in predicted_pep_spec:
        peptide = predicted_pep[0]
        z = predicted_pep[1]
        attrbs = predicted_pep[2]
        pred_mzs = predicted_pep[3]
        pred_intensities = predicted_pep[5]
        pred_ions = predicted_pep[4]

        if peptide in psms:
            nr_raw_scans = {}
            for rawfile_scan in psms[peptide]:
                if z == rawfile_scan.split('@')[2]:
                    exp_rawfile = rawfile_scan.split('@')[0] + '.mzML'
                    exp_scan = rawfile_scan.split('@')[1]
                    nr_raw_scans[exp_rawfile + '@' + exp_scan] = [peptide + '@' + z]

            for f, s in nr_raw_scans.items():
                rawf = f.split('@')[0]
                raw_scan = f.split('@')[1]
                if rawf not in pep_raw_scan_mapped:
                    pep_raw_scan_mapped[rawf] = [[raw_scan + '@' + peptide + '@' + z] + [pred_mzs] + [pred_intensities] + [pred_ions]]
                else:
                    pep_raw_scan_mapped[rawf].append([raw_scan + '@' + peptide + '@' + z] + [pred_mzs] + [pred_intensities] + [pred_ions])

    rawfiles = store_rawfiles(raw_path)

    pep_spectrum_scores = {}
    for rawfile, info in pep_raw_scan_mapped.items():  
    
        try:
            r_file = rawfiles[rawfile][0]

        except:
            raise ValueError(f'No raw file information found for {peptide} with charge {z}')
    
        rawinfo = mzml.read(r_file, read_schema=True)
        scans = {info['id'].split(' ')[-1].strip('scan='):info for info in rawinfo}
    
        scan_count = 0
        for scan_pep in info:
            exp_scan_num = scan_pep[0].split('@')[0]
            sequence = scan_pep[0].split('@')[1]
            charge = scan_pep[0].split('@')[2]

            #print (rawfile, exp_scan_num, scan_pep)
            if exp_scan_num in scans:
            
                ''' Performs pair-wise spectral comparison of experimentally derived MS/MS spectra of a peptide sequence with its Prosit predicted MS/MS spectra
                to get the Dot Product and Normalized spectral contrast angle scores.'''

                scan_number, dotp, spec_angle = process_spectrum(sequence, exp_scan_num, scans[exp_scan_num], scan_pep[1], scan_pep[2], scan_pep[3])

                #print (f'The Dot product and Normalized spectral contrast angle score for peptide {sequence} with respect to its Prosit predicted spectra is {dotp} and {spec_angle}')
        
                pep_spectrum_scores[rawfile + '@' + exp_scan_num + '@' + sequence + '@' + charge] = [str(dotp) + '@' + str(spec_angle)]

            else:
                print (sequence, exp_scan_num)
        
            scan_count += 1
        
        print (f"Processed {scan_count} SSPP spectrum matches from {rawfile}.")

    outfile = "{0}_Dotp_NSCA_scores.txt".format(sspp_msp_file.rstrip('txt').rstrip('.'))
    with open(outfile, 'w') as outf:
        outf.write('Raw File\tScan Num\tPeptide Sequence\tCharge\tDot Product Score\tNorm. Spectral Contract Angle\n')
        for pep_spec, scores in pep_spectrum_scores.items():
            outf.write('\t'.join(pep_spec.split('@') + scores[0].split('@'))+ '\n')

if __name__== "__main__":
    print (args.infiles[0], args.msp_file[0], args.raw_path[0])
    
    SpecAngle_Calc(args.infiles[0], args.msp_file[0], args.raw_path[0])

