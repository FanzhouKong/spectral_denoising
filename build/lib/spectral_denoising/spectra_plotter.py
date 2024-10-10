#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np
from rdkit import Chem
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib import rcParams
from . import spectral_operations as so
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px


import ast
def head_to_tail_plot(msms1, msms2,pmz=None,mz_start = None, mz_end = None, pmz2= None,ms2_error = 0.02,title = None,
                      color1 = None, color2 = None,
                      savepath = None, show= True, publication = False,fontsize = 12):
    """
    Plots a head-to-tail comparison of two MS/MS spectra.

    Parameters:
        msms1 (np.array): First mass spectrum data in 2D np.array format. e,g. np.array([[mz1, intensity1], [mz2, intensity2], ...]).
        msms2 (np.array): Second mass spectrum data. Same as msms1.
        pmz (float or str, optional): Precursor m/z value for the first spectrum. Default is None. If given, precursors will be removed from both spectra and precursor will be shown as a grey
        dashed line in the plot.
        mz_start (float, optional): Start of the m/z range for plotting. Zoom in function. Default is None.
        mz_end (float, optional): End of the m/z range for plotting. Zoom in function. Default is None.
        pmz2 (float or str, optional): Precursor m/z value for the second spectrum. Default is None. Just in case pmz1 and pmz2 are different.
        ms2_error (float, optional): Error tolerance for m/z values. Default is 0.02.
        color1 (str, optional): Color for the first spectrum's peaks. Default is None.
        color2 (str, optional): Color for the second spectrum's peaks. Default is None.
        
        savepath (str, optional): Path to save the plot image. Default is None.
        show (bool, optional): If True, displays the plot. Default is True. Turn it off if you want to save the plot without displaying it.
        publication (bool, optional): If True, formats the plot for publication (size 3*2.5 inch for single column figure). Default is False.
        fontsize (int, optional): Font size for plot labels. Default is 12.
    Returns:
        matplotlib.pyplot or None: The plot object if show is True, otherwise None.
    """
                      
    
    if isinstance(pmz, str):
        pmz = float(pmz)
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
    print('entropy similarity is', so.entropy_similairty(msms1, msms2, pmz, ms2_error = ms2_error))
    if pmz is not None and pmz2 is not None:
        msms1 = so.truncate_spectrum(msms1, pmz-1.6)
        msms2= so.truncate_spectrum(msms2, pmz2-1.6)
    mass1, intensity1 = msms1.T[0], msms1.T[1]
    intensity_nor1 = [x/np.max(intensity1)*100 for x in intensity1]

    mass2, intensity2 = msms2.T[0], msms2.T[1]
    intensity_nor2 = [x/np.max(intensity2)*100 for x in intensity2]
    intensity_nor2=[-x for x in intensity_nor2]
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
        ax.set_xlim(mz_start, mz_end)

    ax.set_ylim(-100, +100)

    plt.axhline(y=0, color='black', linestyle='-')
    start, end = ax.get_ylim()
    plt.tight_layout()
    ax.set_facecolor("none")
    ax.grid(False)
    plt.grid(True, axis="y", color='black', linestyle=':', linewidth=0.1)
    if title != None:
        plt.title(title)
    plt.tight_layout()
    if savepath != None:
        plt.savefig(savepath, dpi = 300,facecolor = 'white', edgecolor = 'none')
    if show==True:
        return(plt)
    else:
        return()








def ms2_plot(msms_1, pmz = None, lower=None, upper=None, savepath = None, color = 'blue'):
    
    """
    Plots a single MS/MS spectrum.
    
    Parameters:
        msms_1 (numpy.ndarray): MS/MS (or MS1) spectrum in 2D np.array format. e,g. np.array([[mz1, intensity1], [mz2, intensity2], ...]).
        pmz (float, optional): Precursor m/z value. If provided, precursor will be removed from the spectrum. Default is None.
        lower (float, optional): Lower bound for m/z values to be plotted. Default is None.
        upper (float, optional): Upper bound for m/z values to be plotted. Default is None.
        savepath (str, optional): Path to save the plot image. If None, the plot will not be saved. Default is None.
        color (str, optional): Color of the spectrum lines. Default is 'blue'.
    Returns:
        matplotlib.pyplot: The plot object.
    """

    if pmz is not None:
        msms_1 = so.truncate_spectrum(msms_1, pmz-1.6)
    mass1, intensity1 = msms_1.T[0], msms_1.T[1]

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
# def ms2_overlay(msms_1=None,msms_2=None,msms_3 = None, pmz = None, savepath = None):
#     fig = plt.figure(figsize = (8, 6))
#     plt.subplots_adjust()
#     ax = fig.add_subplot()
#     if msms_1 is not None:
#         mass1, intensity1 = msms_1.T[0], msms_1.T[1]
#         intensity1 = [x/np.max(intensity1)*100 for x in intensity1]
#         for i in range(len(mass1)):

#             plt.vlines(x = mass1[i], ymin = 0, ymax = intensity1[i],color = 'orange', linewidth=2)

#     if msms_2 is not None:
#         mass2, intensity2 = msms_2.T[0], msms_2.T[1]
#         intensity2 = [x*100 for x in intensity2]
#         for i in range(len(mass2)):
#             plt.vlines(x = mass2[i], ymin = 0, ymax = intensity2[i],color = 'red', linewidth=2)
#     if msms_3 is not None:
#         mass3, intensity3 = msms_3.T[0], msms_3.T[1]
#         intensity3 = [x/np.max(intensity3)*100 for x in intensity3]
#         for i in range(len(mass3)):
#             plt.vlines(x = mass3[i], ymin = 0, ymax = intensity3[i],color = 'blue', linewidth=2)

#     if pmz != None:
#         plt.vlines(x = pmz, ymin = 0, ymax = 100,color = 'grey', linestyle='dashed')
#     # plt.legend()
#     ax.spines['top'].set_visible(False)
#     ax.spines['right'].set_visible(False)

#     ax.get_xaxis().tick_bottom()
#     ax.get_yaxis().tick_left()
#     ax.set_xlabel(r"$m/z$", fontsize = 12)
#     ax.set_ylabel(r"$Intensity\,[\%]$", fontsize = 12)
#     plt.xticks(rotation='vertical')
#     start, end = ax.get_xlim()
#     # start, end = ax.get_xlim(),
#     ax.set_ylim(0, 100)
#     plt.axhline(y=0, color='black', linestyle='-')
#     start, end = ax.get_ylim()
#     # ax.yaxis.set_ticks(np.arange(start, end + 1, 10))
#     plt.grid(True, axis="y", color='black', linestyle=':', linewidth=0.1)
#     ax.grid(False)
#     ax.set_facecolor("white")
#     ax.spines['bottom'].set_color('black')
#     ax.spines['top'].set_color('black')
#     ax.spines['right'].set_color('black')
#     ax.spines['left'].set_color('black')
#     # ax.set(xticklabels=[], yticklabels = [])
#     fig.tight_layout()
#     # fig.set(xlabel = None)
#     if savepath != None:
#         fig.tight_layout()
#         plt.savefig(savepath, dpi = 300,facecolor = 'white', edgecolor = 'white')

#     return(plt)




# def ms2_clean_noise(msms_1, msms_2, pmz1 = None, lower=None, upper=None, savepath = None, hline= None):
#     mass1, intensity1 = msms_1.T[0], msms_1.T[1]
#     mass2, intensity2 = msms2.T[0], msms2.T[1]
#     mass1 = [float(x) for x in mass1]
#     intensity1 = [float(x) for x in intensity1]
#     mass2 = [float(x) for x in mass2]
#     intensity2 = [float(x) for x in intensity2]
#     d = {'m/z':mass1, 'intensity':intensity1}
#     msms1 = pd.DataFrame(d)
#     d = {'m/z':mass2, 'intensity':intensity2}
#     msms2 = pd.DataFrame(d)
#     max_val = np.max(intensity1+intensity2)
#     msms1["normalized_intensity"] = msms1['intensity'] / max_val * 100.0  # normalize intensity to percent
#     msms2["normalized_intensity"] = msms2['intensity'] / max_val * 100.0  # normalize intensity to percent
#     fig = plt.figure(figsize = (4, 3))
#     plt.subplots_adjust()
#     ax = fig.add_subplot()
#     for i in range(len(msms1['m/z'])):
#         plt.vlines(x = msms1["m/z"][i], ymin = 0, ymax = msms1["normalized_intensity"][i],color = 'red', linewidth=3)
#     for i in range(len(msms2['m/z'])):
#         plt.vlines(x = msms2["m/z"][i], ymin = 0, ymax = msms2["normalized_intensity"][i],color = 'blue', linewidth=3)
#     if pmz1 != None:
#         plt.vlines(x = pmz1, ymin = 0, ymax = 100,color = 'grey', linestyle='dashed')
#     # pltalegend()
#     ax.spines['top'].set_visible(False)
#     ax.spines['right'].set_visible(False)
#     if hline is not None:
#         x_min, x_max = ax.get_xlim()
#         # x_min = np.min(mass1+mass2)
#         # x_max = pmz1
#         plt.hlines(xmin = x_min, xmax = pmz1, y = hline,color = 'red', linewidth=1.5, linestyles='dashed')
#     ax.get_xaxis().tick_bottom()
#     ax.get_yaxis().tick_left()
#     ax.set_xlabel(r"$m/z$", fontsize = 12)
#     ax.set_ylabel(r"$Intensity\,[\%]$", fontsize = 12)
#     plt.xticks(rotation='vertical')
#     start, end = ax.get_xlim()
#     # start, end = ax.get_xlim(),
#     if(lower!=None and upper!= None):
#         ax.set_xlim(lower, upper)
#     ax.set_ylim(0, 100)
#     plt.axhline(y=0, color='black', linestyle='-')
#     start, end = ax.get_ylim()
#     # ax.yaxis.set_ticks(np.arange(start, end + 1, 10))
#     plt.grid(True, axis="y", color='black', linestyle=':', linewidth=0.1)
#     ax.grid(False)
#     ax.set_facecolor("white")
#     ax.spines['bottom'].set_color('black')
#     ax.spines['top'].set_color('black')
#     ax.spines['right'].set_color('black')
#     ax.spines['left'].set_color('black')
#     ax.set(xticklabels=[], yticklabels = [])
#     fig.tight_layout()
#     # fig.set(xlabel = None)
#     if savepath != None:
#         plt.savefig(savepath, dpi = 300,facecolor = 'white', edgecolor = 'white')

#     return(plt)
# In[17]:


