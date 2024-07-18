#!/usr/bin/env python
# coding: utf-8

# In[6]:


import pandas as pd
import numpy as np

from matplotlib import rcParams

import toolsets.spectra_operations as so

import matplotlib.pyplot as plt

import seaborn as sns

import plotly.express as px
from toolsets.search import quick_search_values

# In[7]:



# In[169]:


# In[3]:

import ast
def head_to_tail_plot(msms1, msms2,mz_start = None, mz_end = None,pmz=None, pmz2= None,ms2_error = 0.02,
                      color1 = None, color2 = None,lower=None, upper=None, identity = False, normalize = True,
                      savepath = None, show= True, publication = False,fontsize = 12):

    if msms1 is float or msms2 is float:
        # return(np.NAN)
        return(0)
    if isinstance(msms1, str):
        msms1 = ast.literal_eval(msms1)
    if isinstance(msms2, str):
        msms2 = ast.literal_eval(msms2)
    msms1 = so.sort_spectrum(msms1)
    msms2 = so.sort_spectrum(msms2)
    if pmz is not None:
        if pmz2 is None:
            pmz2 = pmz
    print('entropy similarity is', so.entropy_identity(msms1, msms2, pmz, ms2_error = ms2_error))
    if pmz is not None and pmz2 is not None:
        msms1 = so.truncate_spectrum(msms1, pmz-1.6)
        msms2= so.truncate_spectrum(msms2, pmz2-1.6)
    mass1, intensity1 = so.break_spectra(msms1)
    intensity_nor1 = [x/np.max(intensity1)*100 for x in intensity1]

    mass2, intensity2 = so.break_spectra(msms2)
    intensity_nor2 = [x/np.max(intensity2)*100 for x in intensity2]

    # return(msms1, msms2)
    intensity_nor2=[-x for x in intensity_nor2]
        # msms1 = so.cut_msms(msms1, mz_lower = mz_start, mz_upper = mz_end)
        # msms2 = so.cut_msms(msms2, mz_lower = mz_start, mz_upper = mz_end)
    # return(msms1, msms2)
    if publication == True:
        wid = 3
        hi = 2.5
    else:
        wid = 8
        hi = 6
    fig = plt.figure(figsize = (wid, hi))#43
    plt.subplots_adjust()
    ax = fig.add_subplot()
    for i in range(len(mass1)):
        if color1 == None:
            plt.vlines(x = mass1[i], ymin = 0, ymax = intensity_nor1[i],color = 'blue')
        elif color1 != None:
            plt.vlines(x = mass1[i], ymin = 0, ymax = intensity_nor1[i],color = color1)
    if pmz != None:
        plt.vlines(x = pmz, ymin = 0, ymax = 100,color = 'grey', linestyle='dashed')
    for i in range(len(mass2)):
        if color2 ==None:
            plt.vlines(x = mass2[i], ymin = 0, ymax = intensity_nor2[i],color = 'r')
        elif color2 != None:
            plt.vlines(x = mass2[i], ymin = 0, ymax = intensity_nor2[i],color = color2)
    if pmz2 != None:
        plt.vlines(x = pmz2, ymin = -100, ymax = 0,color = 'grey', linestyle='dashed')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_xlabel(r"$m/z$")
    ax.set_ylabel(r"$Intensity\,[\%]$")
    plt.xticks(rotation='vertical')
    if(mz_start is not None and mz_end is not None):
        # print('ttt')
        ax.set_xlim(mz_start, mz_end)
    # else:
    #     ax.set_xlim(-100, 100)
    # labels = [-100, -50, 0, 50, 100]
    # ax.set_yticklabels(labels)
    # ax.set_ylim(msms2['inverse_normalized_intensity'].min(), msms1['normalized_intensity'].max())
    ax.set_ylim(-100, +100)

    plt.axhline(y=0, color='black', linestyle='-')
    start, end = ax.get_ylim()
    plt.tight_layout()
    ax.set_facecolor("none")
    # ax.yaxis.set_ticks(np.arange(start, end + 1, 10))
    ax.grid(False)

    # return(labels)

    plt.grid(True, axis="y", color='black', linestyle=':', linewidth=0.1)
    # ax.grid(False)
    # ax.set_xlim(50, upper)
    plt.tight_layout()
    if savepath != None:
        plt.savefig(savepath, dpi = 300,facecolor = 'white', edgecolor = 'none')
    if show==True:
        return(plt)
    else:
        return()


# In[10]:

from rdkit import Chem
import os
# from mimas.external.features_by_alphapept.load_mzml_data import load_mzml_data
from tqdm import tqdm
from toolsets.spectra_operations import break_spectra, pack_spectra
import matplotlib.pyplot as plt






# def faster_EIC(mix,mzml_dir, precursormz,rt_start = -1, rt_end = -1, 
#     adjusted_height = -1, vlines_location = [], ifppm = True):
#     from toolsets.parallel_functions import get_ms1_intensity_ppm, get_ms1_intensity_mda
#     filename = find_files(mzml_dir, mix+'.mzml')
#     mzml_file_path = os.path.join(mzml_dir, filename[0])
#     test_mzml = load_mzml_data(mzml_file_path)
#     # EIC_intensity = []
#     # ms1_ppm=10
#     ms1_idx = range(0, len(test_mzml['rt_list_ms1']))
#     with Pool() as pool:
#         # EIC_intensity = pool.starmap(get_ms1_intensity, zip(ms1_idx, repeat(test_mzml), repeat(precursormz)))
#         if ifppm==True:
#             EIC_intensity = pool.starmap(get_ms1_intensity_ppm, zip(ms1_idx, repeat(test_mzml), repeat(precursormz)))
#         else:
#             EIC_intensity = pool.starmap(get_ms1_intensity_mda, zip(ms1_idx, repeat(test_mzml), repeat(precursormz)))
#     # for i in tqdm(range(0, len(test_mzml['rt_list_ms1']))):
#     #     # ms1s.append(_extract_MS1(test_mzml, i))
#     #     chopped_ms1_temp = chop_msms(_extract_MS1(test_mzml, i),precursormz-0.01, precursormz+0.01)
#     #     # chopped_ms1_temp = chop_msms(_extract_MS1(test_mzml, i),precursormz*((1 - ms1_ppm * 1e-6)), precursormz*((1 + ms1_ppm * 1e-6)))
        
#     #     mass_chopped, intensity_chopped = break_spectra(chopped_ms1_temp)
#     #     EIC_intensity.append(np.sum(intensity_chopped) )
#     fig, ax = plt.subplots(
#     figsize = (9, 5)
#                       )
#     ax= sns.lineplot(x = test_mzml['rt_list_ms1'], y = EIC_intensity)
#     if rt_start != -1 and rt_end != -1:
#         ax.set_xlim(rt_start, rt_end)
#         ax.set_ylim(0, np.max(EIC_intensity)+100)
#     if adjusted_height!= -1:
#         ax.set_ylim(0, adjusted_height)
#     if len(vlines_location)>0:
#         for position in vlines_location:
#             plt.axvline(x = position, color = 'red')
#     ax.grid(False)
#     return(test_mzml['rt_list_ms1'], EIC_intensity)




# def fast_EIC(mix,mzml_dir, precursormz,rt_start = -1, rt_end = -1, vlines_location = []):
#     filename = find_files(mzml_dir, mix+'.mzml')
#     mzml_file_path = os.path.join(mzml_dir, filename[0])
#     test_mzml = load_mzml_data(mzml_file_path)
#     EIC_intensity = []
#     ms1_ppm=10
    
#     for i in tqdm(range(0, len(test_mzml['rt_list_ms1']))):
#         # ms1s.append(_extract_MS1(test_mzml, i))
#         chopped_ms1_temp = chop_msms(_extract_MS1(test_mzml, i),precursormz-0.01, precursormz+0.01)
#         # chopped_ms1_temp = chop_msms(_extract_MS1(test_mzml, i),precursormz*((1 - ms1_ppm * 1e-6)), precursormz*((1 + ms1_ppm * 1e-6)))
        
#         mass_chopped, intensity_chopped = break_spectra(chopped_ms1_temp)
#         EIC_intensity.append(np.sum(intensity_chopped) )
#     fig, ax = plt.subplots(
#     figsize = (9, 5)
#                       )
#     ax= sns.lineplot(x = test_mzml['rt_list_ms1'], y = EIC_intensity)
#     if rt_start != -1 and rt_end != -1:
#         ax.set_xlim(rt_start, rt_end)
#         ax.set_ylim(0, np.max(EIC_intensity)+100)
#     if len(vlines_location)>0:
#         for position in vlines_location:
#             plt.axvline(x = position, color = 'red')
#     ax.grid(False)
#     return(test_mzml['rt_list_ms1'], EIC_intensity)

def chop_msms(msms, lowest_allowed, highest_allowed):
    # import bisect
    import toolsets.spectra_operations as so
    # msms = so.sort_spectrum(msms)
    mass, intensity = so.break_spectra(msms)
    mass_sorted, intensity_sorted = zip(*sorted(zip(mass, intensity)))
    # upper_allowed = bisect.bisect_right(mass, highest_allowed)
    # lower_allowed = bisect.bisect_left(mass, lowest_allowed)
    # pep_index = mass.index(parention)
    index_start = np.searchsorted(mass_sorted, lowest_allowed,side = 'left')
    index_end = np.searchsorted(mass_sorted, highest_allowed,side = 'right')
    mass_choppped = mass_sorted[index_start:index_end]
    intensity_chopped = intensity_sorted[index_start:index_end]
    if len(mass_choppped)>0 and len(intensity_chopped)>0:
        return(so.pack_spectra(mass_choppped, intensity_chopped))
    else:
        return(so.pack_spectra([0],[0]))
def _extract_MS1(mzml, scan_number):
    from toolsets.spectra_operations import pack_spectra
    ms1_1 = pack_spectra(mzml['mass_list_ms1'][mzml['indices_ms1'][scan_number]:mzml['indices_ms1'][scan_number+1]],
                     mzml['int_list_ms1'][mzml['indices_ms1'][scan_number]:mzml['indices_ms1'][scan_number+1]])
    return (ms1_1)
def ms2_plot(msms_1, pmz = None, lower=None, upper=None, savepath = None, color = 'blue'):
    if type(msms_1) is list:
        msms_1 = so.convert_nist_to_string(msms_1)
    if pmz is not None:
        msms_1 = so.truncate_spectrum(msms_1, pmz-1.6)
    mass1, intensity1 = so.break_spectra(msms_1)
    mass1 = [float(x) for x in mass1]
    intensity1 = [float(x) for x in intensity1]

    if lower is not None:
        idx_left = np.searchsorted(mass1, lower, side= 'left')
    else:
        idx_left = 0
    if upper is not None:
        idx_right = np.searchsorted(mass1, upper, side = 'right')
    else:
        idx_right = len(mass1)
    mass1 = mass1[idx_left:idx_right]
    intensity1 = intensity1[idx_left:idx_right]
    normalized_intensity = [x/np.max(intensity1)*100 for x in intensity1]


    fig = plt.figure(figsize = (4, 3))
    plt.subplots_adjust()
    ax = fig.add_subplot()
    for i in range(len(mass1)):
        plt.vlines(x = mass1[i], ymin = 0, ymax = normalized_intensity[i],color = color, linewidth=2)
    if pmz != None:
        plt.vlines(x = pmz, ymin = 0, ymax = 100,color = 'grey', linestyle='dashed')
    # plt.legend()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_xlabel(r"$m/z$", fontsize = 12)
    ax.set_ylabel(r"$Intensity\,[\%]$", fontsize = 12)
    plt.xticks(rotation='vertical')
    start, end = ax.get_xlim()
    # start, end = ax.get_xlim(), 
    if(lower!=None and upper!= None):
        ax.set_xlim(lower, upper)
    ax.set_ylim(0, 100)
    plt.axhline(y=0, color='black', linestyle='-')
    start, end = ax.get_ylim()
    # ax.yaxis.set_ticks(np.arange(start, end + 1, 10))
    plt.grid(True, axis="y", color='black', linestyle=':', linewidth=0.1)
    ax.grid(False)
    ax.set_facecolor("white")
    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_color('black')
    ax.spines['right'].set_color('black')
    ax.spines['left'].set_color('black')
    # ax.set(xticklabels=[], yticklabels = [])
    fig.tight_layout()
    # fig.set(xlabel = None)
    if savepath != None:
        fig.tight_layout()
        plt.savefig(savepath, dpi = 300,facecolor = 'white', edgecolor = 'white')

    return(plt)
def ms2_overlay(msms_1=None,msms_2=None,msms_3 = None, pmz = None, savepath = None):
    #
    # if pmz is not None:
    #     msms_1 = so.truncate_spectrum(msms_1, pmz-1.6)




    fig = plt.figure(figsize = (4, 3))
    plt.subplots_adjust()
    ax = fig.add_subplot()
    if msms_1 is not None:
        mass1, intensity1 = so.break_spectra(msms_1)
        intensity1 = [x/np.max(intensity1)*100 for x in intensity1]
        for i in range(len(mass1)):

            plt.vlines(x = mass1[i], ymin = 0, ymax = intensity1[i],color = 'orange', linewidth=2)

    if msms_2 is not None:
        mass2, intensity2 = so.break_spectra(msms_2)
        intensity2 = [x/np.max(intensity2)*100 for x in intensity2]
        for i in range(len(mass2)):
            plt.vlines(x = mass2[i], ymin = 0, ymax = intensity2[i],color = 'red', linewidth=2)
    if msms_3 is not None:
        mass3, intensity3 = so.break_spectra(msms_3)
        intensity3 = [x/np.max(intensity3)*100 for x in intensity3]
        for i in range(len(mass3)):
            plt.vlines(x = mass3[i], ymin = 0, ymax = intensity3[i],color = 'blue', linewidth=2)

    if pmz != None:
        plt.vlines(x = pmz, ymin = 0, ymax = 100,color = 'grey', linestyle='dashed')
    # plt.legend()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_xlabel(r"$m/z$", fontsize = 12)
    ax.set_ylabel(r"$Intensity\,[\%]$", fontsize = 12)
    plt.xticks(rotation='vertical')
    start, end = ax.get_xlim()
    # start, end = ax.get_xlim(),
    ax.set_ylim(0, 100)
    plt.axhline(y=0, color='black', linestyle='-')
    start, end = ax.get_ylim()
    # ax.yaxis.set_ticks(np.arange(start, end + 1, 10))
    plt.grid(True, axis="y", color='black', linestyle=':', linewidth=0.1)
    ax.grid(False)
    ax.set_facecolor("white")
    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_color('black')
    ax.spines['right'].set_color('black')
    ax.spines['left'].set_color('black')
    # ax.set(xticklabels=[], yticklabels = [])
    fig.tight_layout()
    # fig.set(xlabel = None)
    if savepath != None:
        fig.tight_layout()
        plt.savefig(savepath, dpi = 300,facecolor = 'white', edgecolor = 'white')

    return(plt)
def ms_figure_simple(msms1, savepath = None, color = 'blue'):
    import toolsets.spectra_operations as so
    # msms1 = (data_temp.iloc[0]['peaks'])
    mass, intensity =so.break_spectra(msms1)
    d = {'m/z':mass, 'intensity':intensity}
    msms1 = pd.DataFrame(d)
    max_val = np.max(msms1['intensity'])
    msms1["normalized_intensity"] = msms1['intensity'] / max_val * 100.0  # normalize intensity to percent
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize = (7, 5),facecolor="white")
    # fig.figure()
    plt.subplots_adjust()
    ax = fig.subplots()
    for i in range(len(msms1['m/z'])):
        plt.vlines(x = msms1["m/z"][i], ymin = 0, ymax = msms1["normalized_intensity"][i],color = color, linewidth = 3)
    # ax.set_xlim(150, 210)
    ax.set_ylim(0, 100)
    ax.set(xticklabels=[])
    ax.set(yticklabels=[])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['bottom'].set(color = 'black', linewidth = 2)
    ax.spines['left'].set_visible(False)
    fig.tight_layout()
    ax.grid(False)
    if savepath == None:
        return(plt)
    else:
        plt.savefig(savepath, dpi = 300,transparent=True, edge_color = 'none')



def ms2_clean_noise(msms_1, msms_2, pmz1 = None, lower=None, upper=None, savepath = None, hline= None):
    mass1, intensity1 = so.break_spectra(msms_1)
    mass2, intensity2 = so.break_spectra(msms_2)
    mass1 = [float(x) for x in mass1]
    intensity1 = [float(x) for x in intensity1]
    mass2 = [float(x) for x in mass2]
    intensity2 = [float(x) for x in intensity2]
    d = {'m/z':mass1, 'intensity':intensity1}
    msms1 = pd.DataFrame(d)
    d = {'m/z':mass2, 'intensity':intensity2}
    msms2 = pd.DataFrame(d)
    max_val = np.max(intensity1+intensity2)
    msms1["normalized_intensity"] = msms1['intensity'] / max_val * 100.0  # normalize intensity to percent
    msms2["normalized_intensity"] = msms2['intensity'] / max_val * 100.0  # normalize intensity to percent
    fig = plt.figure(figsize = (4, 3))
    plt.subplots_adjust()
    ax = fig.add_subplot()
    for i in range(len(msms1['m/z'])):
        plt.vlines(x = msms1["m/z"][i], ymin = 0, ymax = msms1["normalized_intensity"][i],color = 'red', linewidth=3)
    for i in range(len(msms2['m/z'])):
        plt.vlines(x = msms2["m/z"][i], ymin = 0, ymax = msms2["normalized_intensity"][i],color = 'blue', linewidth=3)
    if pmz1 != None:
        plt.vlines(x = pmz1, ymin = 0, ymax = 100,color = 'grey', linestyle='dashed')
    # pltalegend()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if hline is not None:
        x_min, x_max = ax.get_xlim()
        # x_min = np.min(mass1+mass2)
        # x_max = pmz1
        plt.hlines(xmin = x_min, xmax = pmz1, y = hline,color = 'red', linewidth=1.5, linestyles='dashed')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_xlabel(r"$m/z$", fontsize = 12)
    ax.set_ylabel(r"$Intensity\,[\%]$", fontsize = 12)
    plt.xticks(rotation='vertical')
    start, end = ax.get_xlim()
    # start, end = ax.get_xlim(),
    if(lower!=None and upper!= None):
        ax.set_xlim(lower, upper)
    ax.set_ylim(0, 100)
    plt.axhline(y=0, color='black', linestyle='-')
    start, end = ax.get_ylim()
    # ax.yaxis.set_ticks(np.arange(start, end + 1, 10))
    plt.grid(True, axis="y", color='black', linestyle=':', linewidth=0.1)
    ax.grid(False)
    ax.set_facecolor("white")
    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_color('black')
    ax.spines['right'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.set(xticklabels=[], yticklabels = [])
    fig.tight_layout()
    # fig.set(xlabel = None)
    if savepath != None:
        plt.savefig(savepath, dpi = 300,facecolor = 'white', edgecolor = 'white')

    return(plt)
# In[17]:


