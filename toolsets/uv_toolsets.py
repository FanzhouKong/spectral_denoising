import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import toolsets.spectra_operations as so
def read_in_uv_spec(uvpath):
    uv_df = pd.read_csv(uvpath, header=6)
    for index, row in uv_df.iterrows():
        if row['Intensity']<0:
            uv_df.loc[index, 'Intensity']=0
    uv_df['Intensity_normalized']=uv_df['Intensity']/uv_df['Intensity'].sum()
    uv_df['Intensity_standarized']=uv_df['Intensity']/uv_df['Intensity'].max()
    spec = (so.pack_spectra(uv_df['Wavelength'].tolist(), uv_df['Intensity_normalized'].tolist()))
    return (spec)
def read_in_uv(uvpath):
    uv_df = pd.read_csv(uvpath, header=6)
    for index, row in uv_df.iterrows():
        if row['Intensity']<0:
            uv_df.loc[index, 'Intensity']=0
    uv_df['Intensity_normalized']=uv_df['Intensity']/uv_df['Intensity'].sum()
    uv_df['Intensity_standarized']=uv_df['Intensity']/uv_df['Intensity'].max()
    return (uv_df)
def uv_plot_raw(uv1, intensity_col ='Intensity_normalized',save_path = None):
    fig = plt.figure(figsize = (10, 8))#43
    plt.subplots_adjust()
    ax = fig.add_subplot()
    sns.lineplot(x = uv1['Wavelength'], y = uv1[intensity_col], color='blue')
    if save_path is not None:
        plt.savefig(save_path)
def uv_stack_raw(uv1, uv2, intensity_col ='Intensity_normalized',save_path = None ):
    fig = plt.figure(figsize = (10, 8))#43
    plt.subplots_adjust()
    ax = fig.add_subplot()
    sns.lineplot(x = uv1['Wavelength'], y = uv1[intensity_col], color='blue')
    sns.lineplot(x = uv2['Wavelength'], y = uv2[intensity_col], color='red')
    if save_path is not None:
        plt.savefig(save_path)
def uv_plot(uv1, save_path = None):

    fig = plt.figure(figsize = (10, 8))#43
    plt.subplots_adjust()
    ax = fig.add_subplot()
    wl1, int1 = so.break_spectra(uv1)
    # wl2, int2 = so.break_spectra(uv2)
    sns.lineplot(x = wl1, y = int1, color='blue')
    # sns.lineplot(x = wl2, y = int2, color='red')
    ax.grid(False)
    if save_path is not None:
        plt.savefig(save_path)
def uv_stack(uv1, uv2,save_path = None):
    fig = plt.figure(figsize = (10, 8))#43
    plt.subplots_adjust()
    ax = fig.add_subplot()
    wl1, int1 = so.break_spectra(uv1)
    wl2, int2 = so.break_spectra(uv2)
    sns.lineplot(x = wl1, y = int1, color='blue')
    sns.lineplot(x = wl2, y = int2, color='red')
    ax.grid(False)
    ax.set_facecolor('none')
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path)
def uv_search(spec, uv_lib):
    score = []

    for index, row in uv_lib.iterrows():
        score.append(uv_score(spec ,row['uv_spec'] ))
    uv_lib_result = uv_lib.copy()
    uv_lib_result['score']=score
    uv_lib_result.sort_values(by = 'score', ascending=True, inplace=True)
    return(uv_lib_result)
def uv_score(uv1, uv2):
    # uv1_n = so.normalize_spectrum(uv1)
    # uv2_n = so.normalize_spectrum(uv2)
    uv1_n = so.standardize_spectra(uv1)
    uv2_n = so.standardize_spectra(uv2)
    wl1, int1 = so.break_spectra(uv1_n)
    wl2, int2 = so.break_spectra(uv2_n)
    diff = [abs(x-y) for x,y in zip(int1, int2)]
    return (np.sum(diff))