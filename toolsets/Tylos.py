from itertools import chain
import toolsets.raw_data_scaffold as rds
import toolsets.T_rex as trx
import pandas as pd
import numpy as np
import toolsets.spectra_operations as so
from toolsets.search import string_search, quick_search_sorted, quick_search_values
from toolsets.helpers import  find_adducts
from tqdm import tqdm
import toolsets.denoising_related_functions as drf
def libgen_function(std_list, mzml_dir, if_QE = True):
    if if_QE ==True:
        mass_error  = 0.002
    else:
        mass_error = 0.005
    adducts = find_adducts(std_list.columns)
    from itertools import  chain
    print('extracting features, and mapping by mzrt')
    matched_all = pd.DataFrame()
    for mix in tqdm(std_list['mix'].unique()):
        std_list_mix = string_search(std_list, 'mix', mix)
        masses = std_list_mix[adducts].values
        masses = np.array(list(chain.from_iterable(masses)))
        keep = masses>50
        masses = masses[keep]
        feature_targeted = find_feature_targeted(masses=masses, mix=mix, mzml_dir=mzml_dir,mass_error = 0.005, if_QE= if_QE)

        matched = feature_matching(feature_targeted, std_list_mix, adducts, mass_error = mass_error)
        matched_all = pd.concat([matched_all, matched], ignore_index=True)
    msms_denoised = []
    eis = []
    msms_raw = []
    print('denoising....')
    matched_all_predenoising = matched_all.copy()
    for index, row in tqdm(matched_all.iterrows(), total = len(matched_all)):

        entropy_temp = so.entropy_identity(row['msms'], row['msms_mf'], pmz = row['precursor_mz'])
        if entropy_temp<0.99:
            msms_d1 = drf.spectral_denoising(row['msms'], row['reference_smiles'], row['reference_adduct'], max_allowed_deviation=0.02)
            ei1 = drf.get_ei(msms_d1, row['msms'], row['precursor_mz'])
            msms_d2 = drf.spectral_denoising(row['msms_mf'], row['reference_smiles'], row['reference_adduct'])
            ei2 = drf.get_ei(msms_d2, row['msms_mf'], row['precursor_mz'])
            if ei1>ei2:
                msms_raw.append(row['msms'])
                msms_denoised.append(msms_d1)
                eis.append(ei1)
            else:
                msms_raw.append(row['msms_mf'])
                msms_denoised.append(msms_d2)
                eis.append(ei2)
        else:
            msms_d1 = drf.spectral_denoising(row['msms'], row['reference_smiles'], row['reference_adduct'])
            ei1 = drf.get_ei(msms_d1, row['msms'], row['precursor_mz'])
            msms_denoised.append(msms_d1)
            msms_raw.append(row['msms'])
            eis.append(ei1)
    matched_all.drop(columns=['msms', 'msms_pmz',
                              'msms_pmz_intensity', 'msms_idx', 'rt_offset', 'msms_mf', 'msms_mf_pmz',
                              'msms_mf_pmz_intensity', 'rt_offset_mf'], inplace=True)
    matched_all['msms_raw']=msms_raw
    matched_all['msms_denoised']=msms_denoised
    matched_all['eis']=eis
    return matched_all,matched_all_predenoising







def find_feature_targeted(masses, mix, mzml_dir, mass_error = 0.005, if_QE = True):

    mass_tolerance = mass_error
    if if_QE ==True:
        n = 2
        intensity_threshold = 30000
    else:
        n = 1
        intensity_threshold = 1000
    ms1, ms2 = rds.read_mzml(mix, mzml_dir)
    mass_sorted, intensity_sorted, index_sorted, rt_list = trx.build_index(ms1)
    pmz = []
    apex_intensity = []
    rt_apex = []
    rt_start = []
    rt_end = []
    peak_range = []
    n_scans = []
    reci_snrs = []

    for mass in (masses):
        ion_trace = trx.flash_eic(mass, mass_sorted, intensity_sorted, index_sorted)

        idx_left, idx_right = mass_sorted.searchsorted([mass-mass_tolerance, mass+mass_tolerance])
        peaks_all, raw_apex_idx_all,reci_snrs_all = trx.detect_all_peaks(ion_trace, n_neighbor=n, intensity_threshold=intensity_threshold)

        if len(peaks_all)>0:
            for p, r, a in zip(peaks_all, reci_snrs_all, raw_apex_idx_all):
                pmz_statistics = trx.guess_pmz(mass, mass_sorted,
                                               intensity_sorted, index_sorted, idx_left, idx_right, int(a), mass_error= mass_error)
                if pmz_statistics[0]==pmz_statistics[0]:
                    pmz.append(pmz_statistics[0])
                    apex_intensity.append(pmz_statistics[1])
                    rt_apex.append(trx.gaussian_estimator(tuple([int(a)-1,int(a), int(a)+1]),rt_list, ion_trace))
                    rt_start.append(rt_list[p[0]])
                    rt_end.append(rt_list[p[2]])
                    peak_range.append(p)
                    n_scans.append(p[2]-p[0])
                    reci_snrs.append(r)
    feature_targeted = pd.DataFrame(zip(pmz,
                                        rt_apex, rt_start, rt_end,
                                        apex_intensity,
                                        n_scans,peak_range,reci_snrs),
                                    columns=['precursor_mz',
                                             'rt_apex', 'rt_start', 'rt_end',
                                             'ms1_intensity',
                                             'n_scnas', 'ms1_scan_range', 'reci_snr'])
    feature_targeted_ms2 = trx.map_ms2(feature_targeted, ms2)
    return(feature_targeted_ms2)
def feature_matching(feature_targeted, std_list_mix, adducts, mass_error = 0.001,unique_identifier = None, return_raw = False):
    feature_targeted.sort_values(by = 'precursor_mz', inplace=True, ascending=True)
    mix_matched = pd.DataFrame()
    for index, row in std_list_mix.iterrows():
        compound_matched = pd.DataFrame()

        for a in adducts:
            adduct_matched = quick_search_sorted(feature_targeted, 'precursor_mz', row[a]-mass_error, row[a]+mass_error)
            if len(adduct_matched)>0:
                # adduct_matched.insert(0, 'reference_mz', row[a])
                adduct_matched.insert(1, 'reference_name', row['name'])
                adduct_matched.insert(1, 'reference_mz', row[a])
                if unique_identifier is not None:
                    adduct_matched.insert(1,unique_identifier, row[unique_identifier])
                adduct_matched.insert(2, 'reference_adduct', a)

                # adduct_matched.insert(3, 'reference_rt', row['reference_rt'])
                adduct_matched.insert(4, 'reference_smiles', row['smiles'])
                adduct_matched.insert(5, 'reference_formula', row['formula'])
                adduct_matched.insert(6, 'reference_mix', row['mix'])
                adduct_matched.insert(7, 'reference_rt', row['reference_rt'])
                compound_matched  = pd.concat([compound_matched, adduct_matched], ignore_index=True)
        if len(compound_matched)>0:
            if return_raw == False:
                compound_matched = dereplicate(compound_matched)
            mix_matched = pd.concat([mix_matched, compound_matched],ignore_index=True)
            # compound_matched = pd.concat([compound_matched, adduct_matched], ignore_index=True)




    return(mix_matched)
def dereplicate(compound_matched):
    if len(compound_matched)==0:
        return(pd.DataFrame())
    df_return = pd.DataFrame()
    guessed_rt = compound_matched.iloc[np.argmax(compound_matched['ms1_intensity'])]['rt_apex']
    if quick_search_values(compound_matched, 'rt_apex', guessed_rt-10/60, guessed_rt+10/60)['msms'].isna().all():
        return(pd.DataFrame())# must be at least 1 msms around targeted rt
    compound_matched.dropna(subset = 'msms', inplace = True)
    if len(compound_matched)==0:
        return(compound_matched)
    comment = []
    for ma in compound_matched['reference_adduct'].unique():
        # comment = []
        current_adduct = string_search(compound_matched, 'reference_adduct', ma)
        rt_matched = quick_search_values(current_adduct, 'rt_apex', guessed_rt-5/60, guessed_rt+5/60)
        if len(rt_matched)>0:
            current_adduct_left_over = current_adduct.drop(rt_matched.index.tolist())
            major = rt_matched[rt_matched['ms1_intensity']==rt_matched['ms1_intensity'].max()]
            minor = rt_matched[rt_matched['ms1_intensity']<rt_matched['ms1_intensity'].max()]
            if len(major)>1:
                major.sort_values(by = 'rt_offset', ascending=True, inplace = True)
                major = major[0:1]
            df_return = pd.concat([df_return, major], ignore_index=True)
            comment.append('Major')
            df_return = pd.concat([df_return, minor], ignore_index=True)
            comment.extend(['Minor']*len(minor))
            seed_msms = major.iloc[0]['msms']
            for i, j in current_adduct_left_over.iterrows():
                entropy_temp =so.entropy_identity(seed_msms, j['msms'], pmz = major.iloc[0]['precursor_mz'])
                if entropy_temp>0.75:
                    df_return= pd.concat([df_return, pd.DataFrame([j])], ignore_index=True)
                    comment.append('isomer')
    if len(df_return)==0:
        return (df_return)
    df_return.insert(4, 'comment', comment)
    return df_return