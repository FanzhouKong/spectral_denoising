import pandas as pd
from toolsets.search import quick_search_values, quick_search_sorted, string_search
import numpy as np
from collections import Counter
import toolsets.helpers as helpers
def alignment(file_list, pl_dir):
    data_conc_all = pd.DataFrame()
    for f in file_list:
        df_temp = pd.read_csv(os.path.join(pl_dir, f+'.csv'))
        data_conc_all = pd.concat([data_conc_all, df_temp], ignore_index=True)
    data_conc_all.sort_values(by = 'precursor_mz', ascending=True, inplace=True)

    data_conc_all_working = data_conc_all.copy()
    data_conc_dereplicate = pd.DataFrame()
    iso_state_con = []
    is_dut_con = []
    is_halo_con = []
    occurance = []
    ms1_intensity = {}
    for ff in file_list:
        ms1_intensity[ff]=[]
    pmz_con = []
    rt_con = []
    msms_con = []
    ms1_con = []
    # msms_wa = []
    msms_reference = []
    while len(data_conc_all_working)>0:
        seed_feature = data_conc_all_working.iloc[np.argmax(data_conc_all_working['ms1_intensity'])]
        pmz_adj = quick_search_sorted(data_conc_all_working, 'precursor_mz', seed_feature['precursor_mz']-0.005,seed_feature['precursor_mz']+0.005 )
        pmz_rt_adj =quick_search_values(pmz_adj, 'rt_apex', seed_feature['rt_apex']-1.5/60,seed_feature['rt_apex']+1.5/60, )
        pmz_con.append(np.median(pmz_rt_adj['precursor_mz']))
        rt_average = np.median(pmz_rt_adj['rt_apex'])
        rt_con.append(rt_average)
        pmz_rt_adj['rt_offset']=abs(pmz_rt_adj['rt_apex']-rt_average)
        pmz_rt_adj.sort_values(by = 'rt_offset', ascending=True, inplace = True)
        pmz_rt_adj.drop_duplicates(subset = 'base_name', keep = 'first', inplace = True)
        occurance.append(len(pmz_rt_adj))
        # int_purity = pmz_rt_adj['ms1_intensity']*pmz_rt_adj['peak_purity']
        if pmz_rt_adj['msms'].isna().all() == True:
            msms_con.append(np.NAN)
            ms1_con.append(pmz_rt_adj.iloc[np.argmax(pmz_rt_adj['ms1_intensity'])]['ms1'])
            msms_reference.append(np.NAN)
        else:
            rep_msms_idx = np.argmax(pmz_rt_adj['ms1_intensity'])

            # agreement_counts_all = []
            # for i in range(0, len(pmz_rt_adj)):
            #     base_msms = pmz_rt_adj.iloc[i]['msms']
            #     agreements_counts = []
            #     for j in range(0, len(pmz_rt_adj)):
            #         try:
            #             agreements_counts.append(so.entropy_identity(base_msms, pmz_rt_adj.iloc[j]['msms'], pmz = seed_feature['precursor_mz']))
            #         except:
            #             agreements_counts.append(0)
            #     agreement_counts_all.append(sum(e >0.7 for e in agreements_counts))
            # if np.max(agreement_counts_all)>1:
            #     rep_msms_idx = np.argmax(agreement_counts_all)
            #
            # else:
            #     rep_msms_idx = np.argmax(pmz_rt_adj['ms1_intensity'])
            msms_reference.append(pmz_rt_adj.iloc[rep_msms_idx]['base_name'])
            ms1_con.append(pmz_rt_adj.iloc[rep_msms_idx]['ms1'])
            msms_con.append(pmz_rt_adj.iloc[rep_msms_idx]['msms'])
        # conc_row =np.zeros(len(file_list_conc))
        idx_temp = 0
        for ff in file_list:
            current_file_row = string_search(pmz_rt_adj, 'base_name', ff+'.mzML')
            if len(current_file_row)>0:
                ms1_intensity[ff].append(current_file_row.iloc[0]['ms1_intensity'])
            else:
                ms1_intensity[ff].append(0)
            # idx_temp = idx_temp+1
        # ms1_intensity.append(conc_row)
        # data_conc_dereplicate= pd.concat([data_conc_dereplicate, pd.DataFrame([pmz_rt_adj.iloc[np.argmax(int_purity)]])], ignore_index=True)
        # data_conc_dereplicate= pd.concat([data_conc_dereplicate, pd.DataFrame([pmz_rt_adj.iloc[rep_msms_idx]])], ignore_index=True)
        # msms_wa.append(so.weighted_average_spectra(pmz_rt_adj, 'msms','ms1_intensity'))
        # occurance.append(len(pmz_rt_adj['mix'].unique()))
        iso_state_con.append(Counter(pmz_rt_adj['iso_state']).most_common()[0][0])
        is_dut_con.append(Counter(pmz_rt_adj['if_dut']).most_common()[0][0])
        is_halo_con.append(Counter(pmz_rt_adj['if_halo']).most_common()[0][0])
        data_conc_all_working.drop(pmz_rt_adj.index.tolist(), inplace=True)
    # pbar.update(1)
    # break
    alignment_result = pd.DataFrame(zip(pmz_con, rt_con, occurance, iso_state_con, is_dut_con, is_halo_con, msms_con, ms1_con, msms_reference), columns=['precursor_mz', 'rt_apex', 'occurance', 'iso_state', 'is_dut','is_halo', 'msms', 'ms1', 'msms_reference_file'])
    for k in ms1_intensity.keys():
        alignment_result[k]=ms1_intensity[k]
    ms1_intensity = []
    for index, row in alignment_result.iterrows():
        ms1_intensity.append(np.median(row[file_list]))
    alignment_result.insert(2, 'ms1_intensity', ms1_intensity)
    return(alignment_result)
def find_istd_info(istd_info, feature_table):
    istd_info_working = istd_info.copy()
    name_col = helpers.specify_column('common name', istd_info.columns)
    pmz_col = helpers.specify_column('mz', istd_info.columns)
    rt_col=helpers.specify_column('rt', istd_info.columns)
    # print(rt_col)
    istd_info_working.sort_values(by = pmz_col, inplace = True, ascending=True)
    features_istd= pd.DataFrame()
    names = []
    model_mass = []
    while len(istd_info_working)>0:
        seed = istd_info_working.iloc[np.argmin(istd_info_working[rt_col])]
        #     seed = istd_info_working.loc[19]
        current_istds = quick_search_sorted(istd_info_working, pmz_col, seed[pmz_col]-0.005, seed[pmz_col]+0.005)
        current_istds = quick_search_values(current_istds, rt_col, seed[rt_col]-15/60, seed[rt_col]+15/60)
        if len(current_istds)==1:
            # print('tttt')
            _istd = quick_search_values(feature_table, 'precursor_mz', seed[pmz_col]-0.005, seed[pmz_col]+0.005)
            _istd = quick_search_values(_istd, 'rt_apex', seed[rt_col]-5/60, seed[rt_col]+5/60)
            if len(_istd)>0:
                istd_idx = np.argmax(_istd['ms1_intensity'])
                features_istd = pd.concat([features_istd, pd.DataFrame([_istd.iloc[istd_idx]])])
                names.append(seed[name_col])
                model_mass.append(seed[pmz_col])
        else:
            _istd = quick_search_values(feature_table, 'precursor_mz', seed[pmz_col]-0.005, seed[pmz_col]+0.005)
            _istd = quick_search_values(_istd, 'rt_apex', seed[rt_col]-5/60, seed[rt_col]+5/60)
            _istd = _istd[_istd['ms1_intensity']>_istd['ms1_intensity'].max()*0.5]
            _istd.sort_values(by = 'rt_apex', ascending=True, inplace=True)
            current_istds.sort_values(by = rt_col, ascending=True, inplace=True)
            if len(_istd)>= len(current_istds):# sort by order
                names.extend(current_istds[name_col])
                model_mass.extend(current_istds[pmz_col])
                features_istd = pd.concat([features_istd, _istd.iloc[0:len(current_istds)]])
            else:
                for i, j in _istd.iterrows():# try best guess
                    rt_offset = abs(j['rt_apex']-current_istds[rt_col])
                    best_guess_idx = np.argmin(rt_offset)
                    names.append(current_istds.iloc[best_guess_idx][name_col])
                    model_mass.append(current_istds.iloc[best_guess_idx][pmz_col])
                    features_istd = pd.concat([features_istd, pd.DataFrame([j])])
                # end_index = len(_istd)
                # names.extend(current_istds[name_col][0:end_index])
                # features_istd = pd.concat([features_istd, _istd.iloc[0:end_index]])
        istd_info_working.drop(current_istds.index, inplace=True)
    # return(features_istd, names)
    features_istd.insert(0,'Common name', names)
    features_istd.insert(1, 'model mass', model_mass)
    feature_return = feature_table.copy()
    feature_return.insert(2, 'compound_type', ['compound']*len(feature_return))
    feature_return.insert(3, 'annotation', ['unknown']*len(feature_return))
    feature_return.insert(4, 'model_mass',[np.NAN]*len(feature_return))

    for index, row in features_istd.iterrows():
        feature_return.loc[index, 'annotation']=row['Common name']
        feature_return.loc[index,'compound_type']='istd'
        feature_return.loc[index,'model_mass']=row['model mass']

    # features_istd.reset_index(inplace=True, drop=True)
    return (feature_return)
def find_istd_info_msdial(istd_info, feature_table):
    istd_info_working = istd_info.copy()
    name_col = helpers.specify_column('precursor_mzname', istd_info.columns)
    pmz_col = helpers.specify_column('precursor_mz', istd_info.columns)
    rt_col=helpers.specify_column('rt_apex', istd_info.columns)
    istd_info_working.sort_values(by = pmz_col, inplace = True, ascending=True)
    features_istd= pd.DataFrame()
    names = []
    height = []
    msdial_pmz_col = 'Precursor m/z'
    msdial_rt_col = 'RT (min)'
    msdial_height_col = 'Height'
    # for index, row in feature_table.iterrows():
    #     height.append(np.median(row[feature_table.columns[32:-2]]))
    # feature_table['Height']=height
    while len(istd_info_working)>0:
        seed = istd_info_working.iloc[np.argmin(istd_info_working[rt_col])]
        #     seed = istd_info_working.loc[19]
        current_istds = quick_search_sorted(istd_info_working, pmz_col, seed[pmz_col]-0.005, seed[pmz_col]+0.005)
        current_istds = quick_search_values(current_istds, rt_col, seed[rt_col]-15/60, seed[rt_col]+15/60)
        if len(current_istds)==1:
            # print('tttt')
            _istd = quick_search_values(feature_table, msdial_pmz_col, seed[pmz_col]-0.005, seed[pmz_col]+0.005)
            _istd = quick_search_values(_istd, msdial_rt_col, seed[rt_col]-15/60, seed[rt_col]+15/60)
            if len(_istd)>0:
                istd_idx = np.argmax(_istd[msdial_height_col])
                features_istd = pd.concat([features_istd, pd.DataFrame([_istd.iloc[istd_idx]])])
                names.append(seed[name_col])
        else:
            _istd = quick_search_values(feature_table, msdial_pmz_col, seed[pmz_col]-0.005, seed[pmz_col]+0.005)
            _istd = quick_search_values(_istd, msdial_rt_col, seed[rt_col]-15/60, seed[rt_col]+15/60)
            _istd = _istd[_istd[msdial_height_col]>_istd[msdial_height_col].max()*0.5]
            _istd.sort_values(by = msdial_rt_col, ascending=True, inplace=True)
            current_istds.sort_values(by = rt_col, ascending=True, inplace=True)
            if len(_istd)>= len(current_istds):
                names.extend(current_istds[name_col])
                features_istd = pd.concat([features_istd, _istd.iloc[0:len(current_istds)]])
            else:
                end_index = len(_istd)
                names.extend(current_istds[name_col][0:end_index])
                features_istd = pd.concat([features_istd, _istd.iloc[0:end_index]])
        istd_info_working.drop(current_istds.index, inplace=True)
    # return(features_istd, names)
    features_istd.insert(0,'Common name', names)
    feature_return = feature_table.copy()
    feature_return.insert(2, 'compound_type', ['compound']*len(feature_return))
    feature_return.insert(3, 'annotation', ['unknown']*len(feature_return))

    for index, row in features_istd.iterrows():
        feature_return.loc[index, 'annotation']=row['Common name']
        feature_return.loc[index,'compound_type']='istd'
    # features_istd.reset_index(inplace=True, drop=True)
    return (feature_return)
from toolsets.constants import iso_steps
iso_step_C = iso_steps['C']
iso_step_halo = iso_steps['Cl']
iso_step_H = iso_steps['H']
from scipy.stats import pearsonr
def _determine_iso_state(row, mass_sorted, intensity_sorted, index_sorted, iso_range = 2, iso_corrs_threshold = 0.75,if_debug = False):

    ion_trace_temp = flash_eic(row['precursor_mz'], mass_sorted, intensity_sorted, index_sorted)
    index_start = row['ms1_scan_range'][0]
    index_end = row['ms1_scan_range'][2]+1
    index_apex = row['ms1_scan_range'][1]
    halo_check = _check_halogen(row['precursor_mz'], mass_sorted, intensity_sorted, index_sorted, index_start, index_apex, index_end, iso_corrs_threshold = iso_corrs_threshold)
    if_dut = False
    if_halo = False
    iso_state = 0
    if halo_check[0]==True:
        if_halo = halo_check[0]
        iso_state = halo_check[1]
        mono_isotopic_mass = row['precursor_mz']-iso_state/2*iso_step_halo
        ion_trace_mono = flash_eic(mono_isotopic_mass, mass_sorted, intensity_sorted, index_sorted)
        ion_trace_dut = flash_eic(mono_isotopic_mass-iso_step_H, mass_sorted, intensity_sorted, index_sorted)
        corr, _ = pearsonr(ion_trace_temp[index_start:index_end], ion_trace_dut[index_start:index_end])
        if corr == corr and corr >iso_corrs_threshold and ion_trace_dut[index_apex]<ion_trace_mono[index_apex]:
            if_dut = True
        return(iso_state, if_dut, if_halo)

    corrs = np.zeros(2*iso_range+1)
    iso_intensity = np.zeros(2*iso_range+1)
    for step in range(-iso_range,iso_range+1):
        ion_trace_iso = flash_eic(row['precursor_mz']+(step*iso_step_C), mass_sorted, intensity_sorted, index_sorted)
        # print(row['precursor_mz']+(step*iso_step))
        lst_idx = step+iso_range
        corr, _ = pearsonr(ion_trace_temp[index_start:index_end], ion_trace_iso[index_start:index_end])
        iso_intensity[lst_idx]=(ion_trace_iso[row['ms1_scan_range'][1]])
        if corr != corr:
            corrs[lst_idx]=(0)
        else:
            corrs[lst_idx]=(corr)
    if np.argmax(iso_intensity) == 0:
        if_dut = False
    elif corrs[np.argmax(iso_intensity)-1]>iso_corrs_threshold:
        if_dut= True
    else:
        if_dut= False
    # return(iso_range, corrs, iso_intensity)
    iso_state_temp = iso_range-np.argmax(iso_intensity)
    # print(np.argmax(iso_intensity))
    # print(len(iso_intensity))
    if np.argmax(iso_intensity)+1 == len(iso_intensity) and np.max(corrs)>iso_corrs_threshold:
        iso_state = iso_state_temp
    elif np.argmax(iso_intensity)+1 == len(iso_intensity):
        iso_state = np.NAN

    elif corrs[np.argmax(iso_intensity)+1]>iso_corrs_threshold:
        iso_state = iso_state_temp
        # print(corrs[np.argmax(iso_intensity)-1],iso_corrs_threshold)



    if if_debug == False:
        return(iso_state, if_dut, if_halo)
    else:
        return(corrs, iso_intensity)
def _check_halogen(pmz, mass_sorted, intensity_sorted, index_sorted, index_start, index_apex, index_end, iso_corrs_threshold = 0.7):
    iso_step = iso_step_halo
    is_halo = False
    iso_state = 0
    ion_trace_temp = flash_eic(pmz, mass_sorted, intensity_sorted, index_sorted)
    ion_trace_C_left = flash_eic(pmz-iso_step_C, mass_sorted, intensity_sorted, index_sorted)
    ion_trace_C_right = flash_eic(pmz+iso_step_C, mass_sorted, intensity_sorted, index_sorted)
    ion_trace_halo_left = flash_eic(pmz-iso_step, mass_sorted, intensity_sorted, index_sorted)
    ion_trace_halo_right = flash_eic(pmz+iso_step, mass_sorted, intensity_sorted, index_sorted)

    corr_left, _ = pearsonr(ion_trace_temp[index_start:index_end], ion_trace_halo_left[index_start:index_end])
    corr_right, _ = pearsonr(ion_trace_temp[index_start:index_end], ion_trace_halo_right[index_start:index_end])

    if corr_left == corr_left and  corr_left>iso_corrs_threshold and ion_trace_C_left[index_apex] <ion_trace_temp[index_apex] and ion_trace_C_left[index_apex] <ion_trace_halo_left[index_apex]:
        is_halo = True
        search_left = True
    elif corr_right == corr_right and corr_right>iso_corrs_threshold and ion_trace_C_right[index_apex] <ion_trace_temp[index_apex] and ion_trace_C_right[index_apex] <ion_trace_halo_right[index_apex]:
        is_halo = True
        search_left = False
    else:
        return(is_halo, iso_state)
    if search_left == False:
        return(is_halo, 0)
    else:
        round = 1
        while(search_left==True):
            pmz = pmz-iso_step
            ion_trace_temp = flash_eic(pmz, mass_sorted, intensity_sorted, index_sorted)
            ion_trace_C_left = flash_eic(pmz-iso_step_C, mass_sorted, intensity_sorted, index_sorted)
            ion_trace_halo_left = flash_eic(pmz-iso_step, mass_sorted, intensity_sorted, index_sorted)
            corr_left, _ = pearsonr(ion_trace_temp[index_start:index_end], ion_trace_halo_left[index_start:index_end])
            if corr_left == corr_left and corr_left>iso_corrs_threshold and ion_trace_C_left[index_apex] <ion_trace_temp[index_apex] and ion_trace_C_left[index_apex] <ion_trace_halo_left[index_apex]:
                round = round+1
            else:
                search_left = False
        return(is_halo, round*2)