import pandas as pd
import subprocess
import os
import toolsets.spectra_operations as so
import numpy as np
import random


def parallel_args(df_location, outputfoler, length, libraryname, tolerance1, tolerance2, typeofmsms):
    df_locations = df_location * length
    outputfolers = outputfoler * length
    indecies = range(0, length)
    librarynames = libraryname * length
    tolerance1s = tolerance1 * length
    tolerance2s = tolerance2 * length
    typeofmsmss = typeofmsms * length
    return (zip(df_locations, outputfolers, indecies, librarynames, tolerance1s, tolerance2s, typeofmsmss))


def parallel_args_opt(df_location, outputfoler, length, libraryname, tolerance1, tolerance2, typeofmsms, chose=10):
    df_locations = df_location * length
    outputfolers = outputfoler * length
    indecies = random.sample(range(0, length), chose)
    librarynames = libraryname * length
    tolerance1s = tolerance1 * length
    tolerance2s = tolerance2 * length
    typeofmsmss = typeofmsms * length

    return (zip(df_locations, outputfolers, indecies, librarynames, tolerance1, tolerance2, typeofmsmss), indecies)

def tree_denoising_instance(instance, output_folder, libraryname, tolerance1='40ppm', tolerance2='20ppm',
                   typeofmsms='msms'):
    mass, intensity = so.break_spectra(instance[typeofmsms])

    raw_msms = pd.DataFrame({'mass': mass, 'intensity': intensity})
    raw_msms.to_csv("data/temp_data/msms/%s%s%s.txt" % (str(instance.row_num), instance.key, libraryname),
                    header=None, sep=" ", index=False)

    if (instance['PRECURSORMZ']) <= 200:
        tolerance = '0.02Da'
    elif (instance['PRECURSORMZ']) >= 600:
        tolerance = tolerance2
    else:
        #         tolerance = '20
        tolerance = tolerance1
    if(len(mass)<60):
        peak_num_threshold = 60
    else:
        peak_num_threshold = len(mass)
    try:
        p = subprocess.run(["sirius", "-o", 'sirius_workspace/%s/%s%s%s%s' % (
            output_folder, str(instance.row_num), instance.key, libraryname, tolerance),
                            "-f", instance.Formula,
                            "-z", str(instance.PRECURSORMZ),
                            "--adduct", instance.Adduct, "-2",
                            "data/temp_data/msms/%s%s%s.txt" % (str(instance.row_num), instance.key, libraryname),
                            "config", "--MS2MassDeviation.allowedMassDeviation", tolerance,
                            "--MS2MassDeviation.standardMassDeviation", tolerance,
                            "--NoiseThresholdSettings.intensityThreshold", "0.0005",
                            "--NoiseThresholdSettings.maximalNumberOfPeaks", str(peak_num_threshold),
                            "formula"],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                           # timeout = 600
                           )
    except subprocess.TimeoutExpired:
        print("timeout!!")
def tree_denoising(df_location, output_folder, index, libraryname, tolerance1='40ppm', tolerance2='20ppm',
                   typeofmsms='msms'):
    try:
        df = pd.read_csv(df_location)
    except:
        print("No valid")
        exit()
    instance = df.iloc[index]
    mass, intensity = so.break_spectra(instance[typeofmsms])

    raw_msms = pd.DataFrame({'mass': mass, 'intensity': intensity})
    raw_msms.to_csv("temp_data/msms/%s%s%s.txt" % (str(instance.row_num), instance.key, libraryname),
                    header=None, sep=" ", index=False)

    if (instance['PRECURSORMZ']) <= 200:
        tolerance = '0.02Da'
    elif (instance['PRECURSORMZ']) >= 600:
        tolerance = tolerance2
    else:
        #         tolerance = '20
        tolerance = tolerance1
    if(len(mass)<60):
        peak_num_threshold = 60
    else:
        peak_num_threshold = len(mass)
    try:
        p = subprocess.run(["sirius", "-o", 'sirius_workspace/%s/%s%s%s%s' % (
            output_folder, str(instance.row_num), instance.key, libraryname, tolerance),
                            "-f", instance.Formula,
                            "-z", str(instance.PRECURSORMZ),
                            "--adduct", instance.Adduct, "-2",
                            "data/temp_data/msms/%s%s%s.txt" % (str(instance.row_num), instance.key, libraryname),
                            "config", "--MS2MassDeviation.allowedMassDeviation", tolerance,
                            "--MS2MassDeviation.standardMassDeviation", tolerance,
                            "--NoiseThresholdSettings.intensityThreshold", "0.0005",
                            "--NoiseThresholdSettings.maximalNumberOfPeaks", str(peak_num_threshold),
                            "formula"],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                           # timeout = 600
                           )
    except subprocess.TimeoutExpired:
        print("timeout!!")

def tree_denoising_cluster(df_location, output_folder, index, libraryname, tolerance1='40ppm', tolerance2='20ppm',
                   typeofmsms='msms'):
    try:
        df = pd.read_csv(df_location)
    except:
        print("No valid")
        exit()
    instance = df.iloc[index]
    mass, intensity = so.break_spectra(instance[typeofmsms])

    raw_msms = pd.DataFrame({'mass': mass, 'intensity': intensity})
    raw_msms.to_csv("temp_data/msms/%s%s%s.txt" % (str(instance.row_num), instance.key, libraryname),
                    header=None, sep=" ", index=False)

    if (instance['PRECURSORMZ']) <= 200:
        tolerance = '0.02Da'
    elif (instance['PRECURSORMZ']) >= 600:
        tolerance = tolerance2
    else:
        #         tolerance = '20
        tolerance = tolerance1
    if(len(mass)<60):
        peak_num_threshold = 60
    else:
        peak_num_threshold = len(mass)
    try:
        p = subprocess.run(["/share/fiehnlab/users/fzkong/fragtree_calculation/sirius/bin/sirius", "-o", 
            'sirius_workspace/%s/%s%s%s%s' % (
            output_folder, str(instance.row_num), instance.key, libraryname, tolerance),
                            "-f", instance.Formula,
                            "-z", str(instance.PRECURSORMZ),
                            "--adduct", instance.Adduct, "-2",
                            "temp_data/msms/%s%s%s.txt" % (str(instance.row_num), instance.key, libraryname),
                            "config", "--MS2MassDeviation.allowedMassDeviation", tolerance,
                            "--MS2MassDeviation.standardMassDeviation", tolerance,
                            "--NoiseThresholdSettings.intensityThreshold", "0.0005",
                            "--NoiseThresholdSettings.maximalNumberOfPeaks", str(peak_num_threshold),
                            "formula"],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                           # timeout = 600
                           )
    except subprocess.TimeoutExpired:
        print("timeout!!")

def readin_denoised_files(instance, output_folder, libraryname, tolerance1='40ppm', tolerance2='20ppm'):
    if (instance['PRECURSORMZ']) <= 200:
        tolerance = '0.02Da'
    elif (instance['PRECURSORMZ']) >= 600:
        tolerance = tolerance2
    else:
        #         tolerance = '20
        tolerance = tolerance1

    if os.path.isdir('sirius_workspace/%s/%s%s%s%s/0_unknown_/spectra/' % (
    output_folder, str(instance.row_num), instance.key, libraryname, tolerance)):
        for file in os.listdir('sirius_workspace/%s/%s%s%s%s/0_unknown_/spectra/' % (
        output_folder, str(instance.row_num), instance.key, libraryname, tolerance)):
            if file.endswith('.tsv'):
                data_denoised = pd.read_csv("sirius_workspace/%s/%s%s%s%s/0_unknown_/spectra/%s" % (
                output_folder, str(instance.row_num), instance.key, libraryname, tolerance, file),
                                            sep="\t")
                sta = pd.read_csv("sirius_workspace/%s/%s%s%s%s/formula_identifications.tsv" % (
                output_folder, str(instance.row_num), instance.key, libraryname, tolerance),
                                  sep="\t")
                break
        mz = data_denoised['mz'].apply(str).tolist()
        mz_exact = data_denoised['exactmass'].apply(str).tolist()
        intensity = data_denoised['intensity'].apply(str).tolist()
        #         os.rename(src, dst)
        return (so.pack_spectra(mz, intensity), so.pack_spectra(mz_exact, intensity), float(sta['explainedIntensity']),
                float(sta['TreeScore']))
    #         return(mz, intensity)
    else:
        return (np.NaN, np.NaN, 0, 0)


def readin_denoised_msms(instance, output_folder, libraryname, tolerance1='40ppm', tolerance2='20ppm'):
    if (instance['PRECURSORMZ']) <= 200:
        tolerance = '0.02Da'
    elif (instance['PRECURSORMZ']) >= 600:
        tolerance = tolerance2
    else:
        #         tolerance = '20
        tolerance = tolerance1
    if os.path.isdir('sirius_workspace/%s/%s%s%s%s/0_unknown_/spectra/' % (
    output_folder, str(instance.row_num), instance.key, libraryname, tolerance)):
        for file in os.listdir('sirius_workspace/%s/%s%s%s%s/0_unknown_/spectra/' % (
        output_folder, str(instance.row_num), instance.key, libraryname, tolerance)):
            if file.endswith('.tsv'):
                data_denoised = pd.read_csv("sirius_workspace/%s/%s%s%s%s/0_unknown_/spectra/%s" % (
                output_folder, str(instance.row_num), instance.key, libraryname, tolerance, file),
                                            sep="\t")
                break
        mz = data_denoised['mz'].apply(str).tolist()
        intensity = data_denoised['intensity'].apply(str).tolist()
        #         os.rename(src, dst)
        return (so.pack_spectra(mz, intensity))
    #         return(mz, intensity)
    else:
        return (np.NaN)


def readin_denoised_something(something, instance, output_folder, libraryname, tolerance1='40ppm', tolerance2='20ppm'):
    if (instance['PRECURSORMZ']) <= 200:
        tolerance = '0.02Da'
    elif (instance['PRECURSORMZ']) >= 600:
        tolerance = tolerance2
    else:
        #         tolerance = '20
        tolerance = tolerance1
    if os.path.isdir('sirius_workspace/%s/%s%s%s%s/0_unknown_/spectra/' % (
    output_folder, str(instance.row_num), instance.key, libraryname, tolerance)):
        for file in os.listdir('sirius_workspace/%s/%s%s%s%s/0_unknown_/spectra/' % (
        output_folder, str(instance.row_num), instance.key, libraryname, tolerance)):
            if file.endswith('.tsv'):
                sta = pd.read_csv("sirius_workspace/%s/%s%s%s%s/formula_identifications.tsv" % (
                output_folder, str(instance.row_num), instance.key, libraryname, tolerance),
                                  sep="\t")
                break
        return (float(sta[something]))
    else:
        return (np.NaN)

# def readin_something_pfp(instance,something, libraryname, tolerance='40ppm'):

#     if (instance['PRECURSORMZ'])<=200:
#         tolerance = '0.02Da'
#     elif (instance['PRECURSORMZ'])>=600:
#         tolerance = '20ppm'
#     else:
# #         tolerance = '20
#         tolerance = tolerance


#     if os.path.isdir('sirius_workspace/%s%s%s%s/0_unknown_/spectra/' % (str(instance.row_num),instance.key,libraryname,tolerance)):
#         for file in os.listdir('sirius_workspace/%s%s%s%s/0_unknown_/spectra/'%(str(instance.row_num),instance.key,libraryname, tolerance)):
#             if file.endswith('.tsv'):
#                 # data_denoised = pd.read_csv("sirius_workspace/%s%s%s%s/0_unknown_/spectra/%s" % (str(instance.row_num),instance.key,libraryname,tolerance, file), 
#                              # sep = "\t")
#                 sta = pd.read_csv("sirius_workspace/%s%s%s%s/formula_identifications.tsv" % (str(instance.row_num),instance.key,libraryname,tolerance), 
#                              sep = "\t")
#                 break

#         return(float(sta[something]))
# #         return(mz, intensity)
#     else:
#         return(np.NaN)
# def readin_something_hilic(instance,something, libraryname, tolerance='40ppm', typeofmsms = 'msms'):

#     if (instance['PRECURSORMZ'])<=200:
#         tolerance = '0.02Da'
#     else:
# #         tolerance = '20
#         tolerance = tolerance


#     if os.path.isdir('sirius_workspace/hilic_final/%s%s%s%s/0_unknown_/spectra/' % (str(instance.row_num),instance.key,libraryname,tolerance)):
#         for file in os.listdir('sirius_workspace/hilic_final/%s%s%s%s/0_unknown_/spectra/'%(str(instance.row_num),instance.key,libraryname, tolerance)):
#             if file.endswith('.tsv'):
#                 # data_denoised = pd.read_csv("sirius_workspace/%s%s%s%s/0_unknown_/spectra/%s" % (str(instance.row_num),instance.key,libraryname,tolerance, file), 
#                              # sep = "\t")
#                 sta = pd.read_csv("sirius_workspace/hilic_final/%s%s%s%s/formula_identifications.tsv" % (str(instance.row_num),instance.key,libraryname,tolerance), 
#                              sep = "\t")
#                 break

#         return(float(sta[something]))
# #         return(mz, intensity)
#     else:
#         return(np.NaN)
print("I am tree denoising!!!!, and i have pfp both updated!!!")
