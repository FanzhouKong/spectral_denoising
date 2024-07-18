import numpy as np
from scipy.stats import ttest_ind
from toolsets.file_io import get_list_idx
import matplotlib.pyplot as plt
import seaborn as sns
import os
def biodata_prep(bio_data):
    pos_idx = get_list_idx(bio_data['mix'], 'pos')
    neg_idx = get_list_idx(bio_data['mix'], 'neg')
    pos_ctr = []
    neg_ctr = []
    for idx in pos_idx:
        pos_ctr.extend(bio_data.iloc[idx][2:4].values)
    for idx in neg_idx:
        neg_ctr.extend(bio_data.iloc[idx][2:4].values)
    ctr_idx = pos_idx+neg_idx
    bio_data_fractions = bio_data.drop(index=ctr_idx, axis = 0)
    bio_data_fractions.reset_index(inplace=True, drop=True)
    pos_p = []
    neg_p = []
    label =[]
    average_adjusted = []
    for index, row in bio_data_fractions.iterrows():
        pos_p_value = ttest_ind(pos_ctr, row[2:4].tolist(), alternative='less').pvalue
        pos_p.append(pos_p_value)
        neg_p_value = ttest_ind(neg_ctr, row[2:4].tolist(), alternative='greater').pvalue
        neg_p.append(neg_p_value)
        if neg_p_value <= 0.05:
            label.append(True)
        else:
            label.append(False)
        if row['Average'] >= np.median(neg_ctr):
            average_adjusted.append(0)
        else:
            average_adjusted.append(bio_data_fractions['Average'].max()-row['Average'])
    bio_data_fractions['average_adjusted']=average_adjusted
    bio_data_fractions['pos_p_value']=pos_p
    bio_data_fractions['neg_p_value']=neg_p
    bio_data_fractions['peak_label']=label
    return(bio_data_fractions)
def find_edge(lst, peak_idx):
    diff = np.diff(lst)
    for i in range(peak_idx, len(diff)):
        if diff[i] >=0:
            break
    right_edge= i
    for i in range(1, peak_idx):
        if diff[peak_idx-i] <=0:
            break
    left_edge= peak_idx-i+1
    return(left_edge, right_edge)

from scipy.stats import pearsonr
from tqdm import tqdm
def rank_features(alignment, bio_data_processed, mode = 'pos'):
    if mode == 'pos':
        mix_name_pos = [x+'_P' for x in bio_data_processed['mix']]
    elif mode =='neg':
        mix_name_pos = [x+'_N' for x in bio_data_processed['mix']]
    else:
        print('please provide right mode')
        exit()
    alignment_fractions = alignment[mix_name_pos]
    correlation = []
    p_values = []
    max_values = []
    for index, row in tqdm(alignment_fractions.iterrows(), total = len(alignment_fractions)):
        temp = (pearsonr(row.values, bio_data_processed['average_adjusted']))
        max_values.append(row.values.max())
        if temp[0]!=temp[0]: # this is nan
            correlation.append(0)
            p_values.append(0)
        else:
            correlation.append(temp[0])
            p_values.append(temp[1])
    correlation_result_temp = alignment.copy()
    correlation_result_temp.insert(2, 'correlation', correlation)
    correlation_result_temp.insert(3, 'max_values', max_values)
    correlation_result_temp.insert(4, 'p_values', p_values)
    correlation_result_temp = correlation_result_temp[correlation_result_temp['max_values']!=0]
    intensity_log10 = np.log10(correlation_result_temp['max_values'])
    intensity_log10 = [x/intensity_log10.max() for x in intensity_log10]
    correlation_result_temp['log_intensity_normalized']=intensity_log10
    correlation_result_temp.insert(0, 'score', correlation_result_temp['log_intensity_normalized']*correlation_result_temp['correlation'])
    # correlation_result_temp['score']=correlation_result_temp['log_intensity_normalized']*correlation_result_temp['correlation']
    correlation_result_temp.drop(['log_intensity_normalized'], inplace = True, axis = 1)
    correlation_result_temp.sort_values(by = 'score', ascending=False, inplace=True)
    correlation_result_temp.reset_index(inplace=True, drop=True)
    correlation_result_temp.insert(0, 'rank', np.arange(len(correlation_result_temp)))
    return(correlation_result_temp)

def bioacitivity_figs_batch(correlation_result_temp, bio_data_processed, top_n,save_dir, mode = 'pos'):
    for i in range(0, top_n):
        bioactivity_figs(correlation_result_temp, bio_data_processed, rank = i, save_dir = save_dir, show=False, mode = mode)
    # pass
def bioactivity_figs(correlation_result_temp, bio_data_processed, rank,save_dir = None, show = True, mode = 'pos'):
    idx = rank
    if mode == 'pos':
        mix_name_pos = [x+'_P' for x in bio_data_processed['mix']]
    elif mode =='neg':
        mix_name_pos = [x+'_N' for x in bio_data_processed['mix']]
    else:
        print('please provide the right mode name')
        exit()

    instance = correlation_result_temp.iloc[idx]
    fig, axs = plt.subplots(figsize = (10,6))
    sns.lineplot(x =  bio_data_processed['mix'], y = [x/bio_data_processed['average_adjusted'].max() for x in bio_data_processed['average_adjusted']], color = 'blue', label = 'bioactivity')
    sns.lineplot(x =  bio_data_processed['mix'], y = [x/correlation_result_temp.iloc[idx][mix_name_pos].max() for x in correlation_result_temp.iloc[idx][mix_name_pos]], color = 'orange', marker = 'D', label = 'feature_intensity')
    plt.legend(loc = 'upper right')
    axs.grid(False)
    # plt.legend()
    # plt.vlines(x = 7, ymin = 0, ymax = 1, colors='green')
    # plt.vlines(x = 10, ymin = 0, ymax = 1,colors='green')
    plt.xticks(rotation = 90)
    plt.tight_layout()
    axs.set_facecolor("none")
    if save_dir is not None:
        plt.savefig(os.path.join(save_dir, str(instance['rank'])+'_'+str(instance['pmz'])+'_'+str(np.round(instance['rt'], 3))+'.png'), dpi = 300,facecolor = 'none', edgecolor = 'none')
    if show != True:
        plt.close()
