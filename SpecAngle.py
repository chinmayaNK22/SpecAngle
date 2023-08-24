import numpy as np
from pyteomics import mzml
import math
from Scripts import parse_msp
from Scripts import parse_psm_scan
import os
from Scripts import match_spectra
import argparse

parser = argparse.ArgumentParser(description='''Calculates the Normalized Spectral Angle Score for Pepetide Spectrum Matches by comparing to the Prosit predicted library''')

parser.add_argument('infiles', metavar='-i', type=str, nargs='*', help='Input file consisting peptide spectrum matches with raw file and scan details in tab delimitted format')

parser.add_argument('msp_file', metavar='-msp', type=str, nargs='+', help='Prosit predicted peptide spectral library in MSP format')

parser.add_argument('raw_path', metavar='-raw', type=str, nargs='+', help='Path to a location where DDA derived MS/MS data files in mzML format')

parser.add_argument('tol', metavar='-t', type=str, nargs='+', help="Set the fragment tolerance in 'ppm' or 'Da' for matching experimental and theoretical MS/MS peaks")

args = parser.parse_args()

def store_rawfiles(raw_path):
    rawfiles = {}
    for files in os.listdir(raw_path):
        if os.path.isfile(os.path.join(raw_path, files)):
            if files.split('.')[-1] == 'mzML':
                if files not in rawfiles:
                    rawfiles[files] = [os.path.join(raw_path, files)]

    return rawfiles


def L2_norm(intensity_array):
    '''Performs L2 normalization of an intensity array.
    Refered the code from https://github.com/wilhelm-lab/spectrum_fundamentals/blob/development/spectrum_fundamentals/metrics/similarity.py'''

    #L2 normalization
    norm_intensity = np.linalg.norm(intensity_array, 2)

    # Normalize the vector to unit length
    normalized_vector = intensity_array / norm_intensity

    return normalized_vector


def sqrt_norm(ms2_intensities):
    '''Perform square root intensity transformation and normalized intensity vectors
    to unit length as defined in the 2014 MCP paper by Toprak et al. (2014)'''
    
    # Square root transformation of intensities
    sqrt_intensity = np.sqrt(ms2_intensities)

    norm_intensity = np.linalg.norm(sqrt_intensity, 2)

    # Normalize the vector to unit length
    normalized_vector = sqrt_intensity / norm_intensity
    
    return normalized_vector

def calc_dot_product(vector1, vector2):
    ''' This function was generated by refering the dot-product code Statistics::Angle of Statistics::NormalizedContrastAngle
    provided in github repo of pwiz/pwiz_tools/Skyline/Util/Statistics.cs'''

    # Ensure both vectors have the same length
    if len(vector1) != len(vector2):
        raise ValueError("The length of intensity array of matched fragment m/z values from experimental and predicted spectra are not the same")
    
    else:
        # Compute the dot product using numpy's dot() function
        
        '''Dot-product'''
        dot_product = np.dot(vector1, vector2)

        '''Cosine Similary Angle'''
        Cos_Angle = dot_product/math.sqrt(sum(np.square(vector1)) * sum(np.square(vector2)))

        return dot_product, Cos_Angle

def NormSpectralContrastAngle(dotp):
    '''Toprak UH, Gillet LC, Maiolica A, Navarro P, Leitner A, Aebersold R.
    Conserved peptide fragmentation as a benchmarking tool for mass spectrometers
    and a discriminating feature for targeted proteomics.
    Molecular & Cellular Proteomics. 2014 Aug 1;13(8):2056-71.'''
    
    spec_angle = 0
    if abs(dotp) <= 1:
        spec_angle += 1 - (2 * math.acos(dotp)/math.pi)
    else:
        spec_angle = 0

    return spec_angle

def process_spectrum(sequence, scan_num, rawfile, exp_spectra, pred_mzs, pred_intensities, pred_ions, tolerance):
    #print (scan_num, info['m/z array'], info['intensity array'])
    spec_mz = exp_spectra['m/z array']
    spec_intensity = exp_spectra['intensity array']
    
    ''' Normalize the intensity values of all MS/MS peaks from experimental spectra'''
    #norm_exp_intensity = sqrt_norm(spec_intensity)
    norm_pred_intensity = L2_norm(pred_intensities)
    norm_exp_intensity = L2_norm(spec_intensity)

    tolerance_val = float(''.join(_ for _ in tolerance if _.isdigit())) ## Fragment mass tolerance for matching m/z values from experimental and Prosit predicted spectra
    tolerance_unit = ''.join(_ for _ in tolerance if not _.isdigit()) ## Fragment mass tolerance unit either in 'parts per million' (ppm) or in Dalton (Da)

    ''' Perform m/z peak match of experimental raw spectra (MS/MS) of a peptide sequence with its respective Prosit predicted fragment spectra.
    The peak match will be performed with the fragmenta match tolerance either in 'ppm' or 'Da'.'''

    matched_values, matched_exp_intense, matched_pred_intense, matched_ions = match_spectra.peak_match(spec_mz, pred_mzs, norm_exp_intensity, norm_pred_intensity, pred_ions, tolerance_val, tolerance_unit)
    
    '''Calculate the Dot-Product and Cosine Similarity Angle as defined in the 2011 MCP paper by Yen et al. (2011)'''
    dotp, cos_angle = calc_dot_product(matched_pred_intense, matched_exp_intense)
    
    '''The Normalized spectral contrast angle will be calculated as defined in the 2014 MCP paper by Toprak et al. (2014)'''
    norm_spec_angle = NormSpectralContrastAngle(dotp)

    if norm_spec_angle == 0:
        print (f"Normalized spectral contrast angle could not be calculated for the raw spectra {scan_num} of {rawfile} file for {sequence}. Returning Normalized spectral contrast angle score as 0.")
        return scan_num, dotp, cos_angle, norm_spec_angle
    else:
        return scan_num, dotp, cos_angle, norm_spec_angle

def SpecAngle_Calc(infiles, msp_file, raw_path, tolerance):

    psms = parse_psm_scan.psm_rawinfo(infiles)

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
        scans = {rinfo['id'].split(' ')[-1].strip('scan='):rinfo for rinfo in rawinfo}
    
        scan_count = 0
        for scan_pep in info:
            exp_scan_num = scan_pep[0].split('@')[0]
            sequence = scan_pep[0].split('@')[1]
            charge = scan_pep[0].split('@')[2]

            #print (rawfile, exp_scan_num, scan_pep)
            if exp_scan_num in scans:
            
                ''' Performs pair-wise spectral comparison of experimentally derived MS/MS spectra of a peptide sequence with its Prosit predicted MS/MS spectra
                to get the Dot Product and Normalized spectral contrast angle scores.'''

                scan_number, dotp, cos_angle, Norm_spec_angle = process_spectrum(sequence, exp_scan_num, r_file, scans[exp_scan_num], scan_pep[1], scan_pep[2], scan_pep[3], tolerance)

                #print (f'The Dot product and Normalized spectral contrast angle score for peptide {sequence} with respect to its Prosit predicted spectra is {dotp} and {spec_angle}')
        
                pep_spectrum_scores[rawfile + '@' + exp_scan_num + '@' + sequence + '@' + charge] = [str(dotp) + '@' + str(cos_angle) +'@' + str(Norm_spec_angle)]

            else:
                print (sequence, exp_scan_num)
        
            scan_count += 1
        
        print (f"Processed {scan_count} SSPP spectrum matches from {rawfile}.")

    outfile = "{0}_Dotp_NSCA_scores.txt".format(msp_file.rstrip('msp').rstrip('.'))
    with open(outfile, 'w') as outf:
        outf.write('Raw File\tScan Num\tPeptide Sequence\tCharge\tDot Product\tCosine Angle\tNorm. Spectral Contract Angle\n')
        for pep_spec, scores in pep_spectrum_scores.items():
            outf.write('\t'.join(pep_spec.split('@') + scores[0].split('@'))+ '\n')

if __name__== "__main__":
    
    if len(args.infiles) > 1:
        print (f'INFO: There are {len(args.infiles)} input files provided')
    
    SpecAngle_Calc(args.infiles, args.msp_file[0], args.raw_path[0], args.tol[0])


