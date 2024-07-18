
import os
import sys

import spectral_entropy as se

import numpy as np

import toolsets.spectra_operations as so

sys.path.append('yuanyue_code')



import pandas as pd

def _parse_information(line: bytes, spec: dict) -> int:
    """
    Parse the line in .msp file, update information in the spec
    :param line: The input line.
    :param spec: The output information collection.
    :return: The entry of the added information.
    """
    line = line.strip()
    if line:
        line = line.strip()
        lines = line.split("=")
        if len(lines) > 2 and lines[0].strip()!='TITLE':
            item = lines[0].strip()
            cont = " ".join(lines[1:])
            spec[item]=cont.strip()
            return 1
        elif lines[0].strip()=='TITLE':
            item = 'Name'
            cont = "=".join(lines[1:])
            subconts = cont.split(" ")
            spec[item] = subconts[0]
            item = 'Notes'
            spec[item] = " ".join(subconts[1:])
            return 1
        elif lines[0].strip()=='RTINSECONDS':
            item = 'RETENTIONTIME'
            cont = float(lines[1].strip())/60
            spec[item]=cont
            return 1
        elif lines[0].strip()=='PEPMASS':
            item = 'PrecursorMZ'
            spec[item]=lines[1].split(" ")[0]
            return 0
        # elif lines[0].strip()=='CHARGE':
        #     item, cont = lines
        #     item = item.strip()
        #     spec[item]=cont.strip()
        #     return 0
        else:
            item, cont = lines
            item = item.strip()
            spec[item]=cont.strip()
            return 1



def _parse_spectrum(line: str, spec: list) -> int:
    """
    Add peak data to spec
    :param line: The raw line information from .msp file
    :param spec: The spectrum will be added.
    :return: 0: success. 1: no information in this line.
    """
    line = line.strip()
    lines = line.split()
    if line.startswith("CHARGE="):
        return 1
    elif len(lines) >= 2 and lines[0]!='END':
        mz, intensity = lines[0], lines[1]
        spec.append([float(mz), float(intensity)])
        return 1
    else:
        return 0

def read_one_spectrum_mgf(fi)->dict:
    spectrum_info = {
        'spectrum': []
    }
    is_adding_information = 1
    for line in fi:
        # if line.strip() == 'END IONS':
        #     is_adding_information = 0

        if line.strip() == 'BEGIN IONS':
            is_adding_information = 1
        else:
            if is_adding_information>0:
                is_adding_information=_parse_information(line, spectrum_info)
                # is_adding_information+=1
            else:
                spec = spectrum_info['spectrum']
                r = _parse_spectrum(line, spec)
                if(r ==0):
                    yield spectrum_info
                    is_adding_information = 1
                    spectrum_info = {
                        'spectrum': []
                    }
def _export_library(data_dup,output_location, typeofmsms='msms'):
    entry = ''
    for index, row in data_dup.iterrows():
        entry = entry + 'Name: ' + row['Name'] + '\n'
        entry = entry + 'Precursor_type: ' + row['Precursor_type'] + '\n'
        entry = entry + 'Spectrum_type: ' + row['Spectrum_type'] + '\n'
        entry = entry + 'PRECURSORMZ: ' + str(row['PrecursorMZ']) + '\n'
        entry = entry + 'RETENTIONTIME: ' + str(round(row['RETENTIONTIME'],4)) + '\n'
        entry = entry + 'Comments: ' + str(row['Notes']) + '\n'
        entry = entry + 'Num peaks: ' + str(row['Num_peaks']) + '\n'
        entry = entry + (row[typeofmsms])
        # entry = entry +str(row['count'])
        entry = entry + '\n'
        entry = entry + '\n'

    #open text file
    text_file = open(output_location, "w",encoding='utf-8')

    #write string to file
    text_file.write(entry)

    #close file
    text_file.close()




def convert_mgf_to_msp(f, export_filename):
    Name = []
    Notes = []
    Precursor_type = []
    Spectrum_type = []
    PrecursorMZ = []
    Num_peaks = []
    RETENTIONTIME = []
    msms = []
    for x in read_one_spectrum_mgf(f):
        Name.append(x['Name'])
        Notes.append(x['Notes'])
        Precursor_type.append('Unknown')
        Spectrum_type.append('MS2')
        PrecursorMZ.append(x['PrecursorMZ'])
        Num_peaks.append(len(x['spectrum']))
        RETENTIONTIME.append(x['RETENTIONTIME'])
        msms.append(so.convert_nist_to_string(x['spectrum']))
    output = pd.DataFrame({'Name':Name,
                           'Notes': Notes,
                           'RETENTIONTIME':RETENTIONTIME,
                           'Precursor_type':Precursor_type,
                           'Spectrum_type':Spectrum_type,
                           'PrecursorMZ':PrecursorMZ,
                           'Num_peaks':Num_peaks,
                           'msms':msms,
                           })
    _export_library(output, export_filename, 'msms')


