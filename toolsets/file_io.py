import pandas as pd 
# import toolsets.spectra_operations as so
# from toolsets.spectra_operations import clean_spectrum, break_spectra, pack_spectra,convert_arr_to_string
from toolsets.search import num_search, string_search
import re
import os
import toolsets.spectra_operations as so
from toolsets.msp_file import read_one_spectrum
from toolsets.filename import smart_io
import re
import shutil
import numpy as np
from tqdm import tqdm
from fuzzywuzzy import fuzz
import pymzml
import json
def save_df(df, save_path):
    data = df.copy()
    cols = []
    for c in df.columns:
        if isinstance(df.iloc[0][c], np.ndarray):
            if np.shape(df.iloc[0][c])[1]==2:
                cols.append(c)
    print(cols)
    if save_path.endswith('.csv') == False:
        save_path = save_path+'.csv'
    for col in cols:
        specs = []
        for index, row in tqdm(data.iterrows(), total = len(data)):
            specs.append(so.convert_arr_to_string(row[col]))
        data[col]=specs
    data.to_csv(save_path, index = False)
def read_df(path):
    df = pd.read_csv(path)
    print('done read in df')
    cols = []
    for c in df.columns:
        if check_pattern(df.iloc[0][c]):
            cols.append(c)
    for col in cols:
        df[col] = [so.convert_string_to_arr(y[col]) for x,y in df.iterrows()]
    return(df)
def check_pattern(input_string):
    # Regular expression to match pairs of floats in standard or scientific notation separated by a tab,
    # with each pair on a new line
    if isinstance(input_string, str):
        if '\t' in input_string:
            return True
    return False
def split_pos_neg(all_folder, tail = '.mzML'):
    print(all_folder)
    pos_folder = os.path.join(all_folder, 'pos')
    neg_folder = os.path.join(all_folder, 'neg')
    for folder in [pos_folder, neg_folder]:
        if os.path.exists(folder)==False:
            os.makedirs(folder)
    for root, dirs, files in os.walk(all_folder):
        for file in files:
            # print(file)
            if file.endswith(tail):
                # print('i am in if')
                if len(file.split('.'))==2:
                    if file[-(len(tail)+1)]=='P':
                        # print('i am in other if')
                        shutil.move(os.path.join(all_folder, file), os.path.join(pos_folder, file))
                    elif file[-(len(tail)+1)]=='N':
                        shutil.move(os.path.join(all_folder, file), os.path.join(neg_folder, file))
def prepare_sample_list(filelist):
    # print('i am in new')
    qc_idx = get_list_idx(filelist, 'qc')
    blk_idx_long = get_list_idx(filelist, 'blank')
    blk_idx_short = get_list_idx(filelist, 'blk')
    blk_idx = list(set(blk_idx_long+(blk_idx_short)))
    fraction_idx = []
    for i in range(len(filelist)):
        if i not in qc_idx+blk_idx:
            fraction_idx.append(i)
    # print(fraction_idx)
    # [mylist[i] for i in idx]
    # return(fraction_idx)
    # samples, qcs, blks
    return([filelist[i] for i in fraction_idx], [filelist[i] for i in qc_idx], [filelist[i] for i in blk_idx])
def get_file_list(dir, tail, with_tail = False):
    file_list = []
    for root, dirs, files in os.walk(dir):
        for file in files:
            if file.endswith(tail):
                if with_tail == True:
                    file_list.append(file)
                else:
                    file_list.append(file.split('.')[0])
    return(file_list)
def get_list_idx(lst, pattern):
    pattern = pattern.lower()
    list_lower = [x.lower() for x in lst]
    idx = []
    for i in range(0, len(list_lower)):
        if pattern in list_lower[i]:
           idx.append(i)
    return(idx)

def readin_mzml(mzml_path):
    specs = pymzml.run.Reader(mzml_path, obo_version="4.1.33",
        )

    return(specs)
def specify_column(keyword, column_name):
    score = []
    for name in column_name:
        score.append(fuzz.token_sort_ratio(keyword,name))
    return(column_name[np.argmax(score)])

def parse_file_name(filepath):
    import ntpath
    return(ntpath.basename(filepath))



def readin_peak_list(file, alighment = False):
    if alighment == True:
        df = pd.read_csv(file,
                         sep = '\t',
                         header=[4],
                         low_memory=False
                         )
    else:
        df = pd.read_csv(file,
                         sep = '\t',
                         low_memory=False
                         # header=[4]
                         )
    return(df)
def read_msp_files(msp_path, if_rt = False):

    msp_file = smart_io(msp_path, mode = 'r')
    specs = []
    print('reading in spec')
    for spec in read_one_spectrum(msp_file, include_raw=False):
        specs.append(spec)
    print('done reading in spec')
    msp = pd.DataFrame.from_dict(specs)
    # return(msp)
    # precursor_col = specify_column('')
    print('checking spectrum')
    msp_with_spectrum = msp[msp['spectrum'].map(lambda d: len(d)) > 0]
    # print('converting spectrum')
    # peaks = []
    # for index, row in tqdm(msp_with_spectrum.iterrows(), total = len(msp_with_spectrum)):
    #     peaks.append(convert_nist_to_string(row['spectrum']))
    # msp_with_spectrum['msms']=peaks
    # msp_with_spectrum = msp_with_spectrum[msp_with_spectrum['peaks'].map(lambda d: len(d)) > 0]
    # charge = []
    # for index, row in msp_with_spectrum.iterrows():
    #     if row[adduct_col][-1]=='+':
    #         charge.append(1)
    #     elif row[adduct_col][-1]=='-':
    #         charge.append(-1)
    #     else:
    #         charge.append('np.NAN')
    # msp_with_spectrum.insert(msp_with_spectrum.columns.get_loc(adduct_col)+1, 'charge', charge)
    # msp_return = msp_with_spectrum[['Name','Formula', 'Precursor_type','PrecursorMZ','InChIKey','Ion_mode', 'Instrument_type','Ionization', 'Collision_energy', 'msms']]
    return (msp_with_spectrum)


def check_missing_files(file_name_list, tail, source_dir):
    duplicate_files = []
    missing_files = []
    good_files = []
    for mix in (file_name_list):
        found_file = find_files(source_dir, f"{str(mix)+tail}")
        if len(found_file)==0:
            missing_files.append(mix)
        elif len(found_file)==1:
            good_files.append(mix)
        elif len(found_file)>1:
            duplicate_files.append(mix)
    return(good_files, duplicate_files, missing_files)
def read_in_alphapept(path_to_features, peak_purity = 0.001, ifclean = False):
    import toolsets.spectra_operations as so
    features = pd.read_csv(path_to_features)
    for index, row in features.iterrows():
        try:
            features.at[index, "peaks"]=parse_feature_peaks(row['peaks'])
        except ValueError:
            print("the scannumber is ", row['scan_number'])
            print("the problematic feature is ", path_to_features.split("/")[-1])
            features.at[index, "peaks"]=np.NAN
            # break
    features.drop_duplicates(subset = ['scan_number'], keep = "first", inplace = True, ignore_index = True)
    # the above line has been changed, it is actually the right one now!!!!
    features = num_search(features, 'ms1_precursor_intensity', 0, ">", inclusion = False)
    features = num_search(features, 'ms1_intensity_ratio', peak_purity, ">", inclusion = False)
    vc = features['charge'].value_counts().rename_axis('charge').reset_index(name='counts')
    vc.sort_values(by=['counts'], inplace=True, ascending=False)
    if vc.iloc[0]['charge']>0:
        features = string_search(features, 'charge', 1)
    elif vc.iloc[0]['charge']<0:
        features = string_search(features, 'charge', -1)
    if ifclean == True:
        peaks_cleaned = []
        for index, row in features.iterrows():
            peaks_cleaned.append(so.clean_spectrum(row['peaks'], max_mz = row['precursor_mz'],
            tolerance = 0.02, ifppm = False, noise_level = 0.00))
        features['peaks_cleaned']=peaks_cleaned
    # features_plus = string_search(features, 'charge', 1)
    # features_minus = string_search(features, 'charge', -1)
    # features['peaks'] = features.apply(parse_feature_peaks, axis = 1)
    return features

def parse_feature_peaks(peaks):
    peaks_split = peaks.replace('[','').replace(']','').split('\n')
    mass = []
    intensity = []
    for peak in peaks_split:
        peak = peak.strip()
        mass_temp = float(peak.split(" ")[0])
        intensity_temp = float(peak.split(" ")[-1])
        mass.append(mass_temp)
        intensity.append(intensity_temp)

    return(pack_spectra(mass, intensity))





from tqdm import tqdm
from rdkit import Chem
# from rdkit.Chem.rdMolDescriptors import CalcMolFormula
# from rdkit.Chem.Descriptors import ExactMolWt
from toolsets.spectra_operations import num_peaks
from toolsets.spectra_operations import spectral_entropy, truncate_spectrum
def export_library_msp(data,output_location, typeofmsms='msms_denoised', ifcollision_energy = False):
    entry = ''
    for index, row in tqdm(data.iterrows(), total = len(data)):
        if row[typeofmsms]==row[typeofmsms]:
            mol = Chem.MolFromSmiles(row['reference_smiles'])
            entry = entry + 'Name: ' + row['reference_name'] + '\n'
            entry = entry + 'InChIKey: ' + str(Chem.MolToInchiKey(mol)) + '\n'
            entry = entry + 'SMILES: ' + str(row['reference_smiles']) + '\n'
            entry = entry + 'RETENTIONTIME: ' + str(np.round(row['rt_apex'],2)) + '\n'
            entry = entry + 'reference_rt: ' + str(row['reference_rt']) + '\n'
            entry = entry +'Spectrum_type: '+'MS2'+ '\n'
            entry = entry + 'PrecursorMZ: ' + str(row['precursor_mz']) + '\n'
            # entry = entry + 'InChIKey: ' + str(row['reference_inchikey']) + '\n'

            entry = entry + 'Formula: ' + CalcMolFormula(mol) + '\n'

            entry = entry + 'ExactMass: ' + str(ExactMolWt(mol)) + '\n'
            entry = entry + 'Precursor_type: ' + row['reference_adduct'] + '\n'
            if ifcollision_energy:
                entry = entry + 'Collision_enerty: ' + str(row['Collision_energy']) + '\n'
            # entry = entry + 'RETENTIONTIME: ' + str(row['retention_time_wa']) + '\n'
            if row['reference_adduct'][-1]=='+':
                # charge = '1+'
                ionmode = 'P'
            else:
                ionmode = 'N'
            entry = entry+'Ion_mode: '+ionmode+ '\n'
            entry = (entry + 'Comment: ' +str(row['comment'])+'_ms1intensity+'+str(np.round(row['ms1_intensity'],1))
                    +'_ei_'+str(np.round(row['ei'],2)) + '\n')
            try:
                entry = entry + 'Spectrum_entropy: ' +str(spectral_entropy(truncate_msms(row[typeofmsms], row['precursor_mz']-1.6))) + '\n'
            except:
                entry = entry + 'Spectrum_entropy: ' +str(-1) + '\n'
            # entry = entry + 'Normalized_entropy: ' + str((row['normalized_entropy'])) + '\n'
            entry = entry + 'Num peaks: ' + str(num_peaks(row[typeofmsms])) + '\n'
            entry = entry + row[typeofmsms]
            # entry = entry +str(row['count'])
            entry = entry + '\n'
            entry = entry + '\n'

    #open text file
    text_file = open(output_location, "w",encoding='utf-8')
     
    #write string to file
    text_file.write(entry)
     
    #close file
    text_file.close()

def export_ms_sirius(row, output):
    # if row['Adduct'][-1]=='+':
    #         charge = '1+'
    #     else:
    #         charge = '1-'
    mass_1, intensity_1 = so.break_spectra(row['ms1'])
    pep_mass =row['Precursor m/z']
    entry = ''
    entry = entry + '>compound '+str(row['reference_file'])+'_'+str(row['PeakID'])+'\n'
    entry = entry + '>parentmass '+str((pep_mass))+'\n'
    entry = entry + '>ionization '+str((row['Adduct']))+'\n'
    entry = entry +'\n'
    entry = entry + '>collision 35' + '\n'
    entry = entry+(row['msms'])+'\n'
    entry = entry +'\n'
    entry = entry +'\n'
    entry = entry + '>ms1peaks' + '\n'
    entry = entry+(row['ms1'])+'\n'
    entry = entry +'\n'
    text_file = open(output+'.ms', "w",encoding='utf-8')
    text_file.write(entry)
    text_file.close()


def export_mgf_shortlist(row, output_dir, mode = '+'):

    entry = ''
    if mode=='+':
        charge = '1+'
    else:
        charge = '1-'
    pep_mass =row['pmz']
    # mass_1, intensity_1 = so.break_spectra(row['ms1'])
    output = os.path.join(output_dir, str(row['rank'])+'_'+str(row['pmz'])+'_'+
                          str(np.round(row['rt'],2))+'.mgf')
    # ms1
    entry = entry + 'BEGIN IONS'+'\n'
    entry = entry + 'PEPMASS='+str(pep_mass)+'\n'
    entry = entry + 'MSLEVEL=1'+'\n'
    entry = entry+'CHARGE=' + charge +'\n'
    entry = entry+(row['ms1'])+'\n'
    entry = entry + 'END IONS'+'\n'

    entry = entry +'\n'
    entry = entry + 'BEGIN IONS'+'\n'
    entry = entry + 'PEPMASS='+str(pep_mass)+'\n'
    entry = entry + 'MSLEVEL=2'+'\n'
    entry = entry+'CHARGE=' + charge +'\n'
    entry = entry+(row['msms'])+'\n'
    entry = entry + 'END IONS'+'\n'


    text_file = open(output, "w",encoding='utf-8')
    text_file.write(entry)
    text_file.close()
def export_mgf_sirius(row, output_dir):

    entry = ''
    if row['Adduct'][-1]=='+':
        charge = '1+'
    else:
        charge = '1-'
    pep_mass =row['Precursor m/z']
    # mass_1, intensity_1 = so.break_spectra(row['ms1'])
    output = os.path.join(output_dir, row['reference_file']+'+'+str(row['PeakID'])+'.mgf')
    # ms1
    entry = entry + 'BEGIN IONS'+'\n'
    entry = entry + 'PEPMASS='+str(pep_mass)+'\n'
    entry = entry + 'MSLEVEL=1'+'\n'
    entry = entry+'CHARGE=' + charge +'\n'
    entry = entry+(row['ms1'])+'\n'
    entry = entry + 'END IONS'+'\n'

    entry = entry +'\n'
    entry = entry + 'BEGIN IONS'+'\n'
    entry = entry + 'PEPMASS='+str(pep_mass)+'\n'
    entry = entry + 'MSLEVEL=2'+'\n'
    entry = entry+'CHARGE=' + charge +'\n'
    entry = entry+(row['msms'])+'\n'
    entry = entry + 'END IONS'+'\n'


    text_file = open(output, "w",encoding='utf-8')
    text_file.write(entry)
    text_file.close()
        # break

def export_mat(data,output_location, typeofms1='ms1',typeofmsms = 'msms', ifcollision_energy = True):
    entry = ''
    for index, row in data.iterrows():
        entry = entry + 'NAME: ' + row['NAME'] + '\n'
        entry = entry + 'PrecursorMZ: ' + str(row['PRECURSORMZ']) + '\n'
        entry = entry + 'PRECURSORTYPE: ' + row['Adduct'] + '\n'
        entry = entry + 'INSTRUMENTTYPE: ' +'\n' #empty line
        entry = entry + 'INSTRUMENT: ' +'\n'
        entry = entry + 'Authors: ' +'Arpana, Parker and Fanzhou'+'\n'
        entry = entry + 'License: '+'\n'
        entry = entry + 'FORMULA: ' + str(row['Formula']) + '\n'
        entry = entry + 'ONTOLOGY: ' +'\n'
        entry = entry + 'SMILES: ' +'\n'
        entry = entry + 'INCHIKEY: ' + row['InChIKey'] + '\n'
        entry = entry + 'INCHI: ' +'\n'
        entry = entry+'IONMODE: '+row['Ion_mode']+ '\n'
        if ifcollision_energy:
            entry = entry + 'Collision_enerty: ' + str(row['Collision_energy']) +'eV'+ '\n'
        entry = entry+'SPECTRUMTYPE: Centroid and composite'+ '\n'
        entry = entry + 'METABOLITENAME: ' + '\n'
        entry = entry + 'SCANNUMBER: Alignment ID '+str(row['Alignment_ID']) + '\n'
        entry = entry + 'RETENTIONTIME: ' + str(row['RETENTIONTIME']) + '\n'
        entry = entry + 'RETENTIONINDEX: N/A' +'\n'
        entry = entry + 'CCS: ' +'\n'
        entry = entry + 'INTENSITY: ' +str(row['intensity'])+'\n'
        entry = entry + '#Specific field for labeled experiment' +'\n'
        entry = entry + 'IsMarked: False' +'\n'
        entry = entry + 'Comment: ' + str(row['Comment']) + '\n'
        entry = entry + 'MSTYPE: MS1'+ '\n'
        entry = entry + 'Num Peaks: ' + str(num_peaks(row[typeofms1])) + '\n'
        mass1, intensity1= so.break_spectra(row[typeofms1])
        for i in range(0, len(mass1)):
            entry = entry + str(mass1[i]) +'\t'+ str(intensity1[i]) + '\t'+'"'+str(mass1[i])+'"'+'\n'
        entry = entry + 'MSTYPE: MS2'+ '\n'
        entry = entry + 'Num Peaks: ' + str(num_peaks(row[typeofmsms])) + '\n'
        mass2, intensity2= so.break_spectra(row[typeofmsms])
        for i in range(0, len(mass2)):
            entry = entry + str(mass2[i]) +'\t'+ str(intensity2[i]) + '\t'+'"'+str(mass2[i])+'"'+'\n'
        entry = entry + '\n'
        #
        #
        # entry = entry +'Spectrum_type: '+row['Spectrum_type']+ '\n'
        #
        # entry = entry + 'InChIKey: ' + row['InChIKey'] + '\n'
        #
        # entry = entry + 'ExactMass: ' + str(row['ExactMass']) + '\n'
        #
        #
        #
        # entry = entry + 'Num peaks: ' + str(num_peaks(row[typeofmsms])) + '\n'
        # entry = entry + row[typeofmsms]
        # # entry = entry +str(row['count'])

        # entry = entry + '\n'

    #open text file
    text_file = open(output_location, "w",encoding='utf-8')

    #write string to file
    text_file.write(entry)

    #close file
    text_file.close()

# def find_files(base, pattern):
#  this is legacy method
#     # print('i am in new')
#     score = []
#     for filename in os.listdir(base):
#         score.append(fuzz.token_sort_ratio(filename,pattern))
#     return(os.listdir(base)[np.argmax(score)])
#     '''Return list of files matching pattern in base folder.'''
#     # return [filename for filename in os.listdir(base) if re.search(pattern, filename, re.IGNORECASE)]
