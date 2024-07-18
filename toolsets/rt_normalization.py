import numpy as np
import pandas as pd
from toolsets.search import quick_search_values
from tqdm import tqdm
from toolsets.file_io import readin_peak_list, get_file_list
from sklearn.ensemble import RandomForestRegressor
import os
def rt_normalization(peak_list_dir, rt_seed, istd_info, normalized_peak_list_dir, result_dir):
    all_lst = get_file_list(peak_list_dir, '.txt')
    rt_report = istd_info.copy()
    rt_report = rt_report[['compound_name', 'Precursor m/z']]
    intensity_report = istd_info.copy()
    intensity_report= intensity_report[['compound_name', 'Precursor m/z']]
    if rt_seed.endswith('.txt')==False:
        rt_seed = rt_seed+'.txt'
    seed_file =readin_peak_list(os.path.join(peak_list_dir, rt_seed), msial=False)
    # seed_file =readin_peak_list(os.path.join(peak_list_dir, rt_seed+'.txt'), msial=True)
    seed_rt, seed_intensity = find_istd(istd_info, seed_file)
    istd_found = pd.DataFrame(zip(seed_rt, istd_info['Precursor m/z']), columns=['RT_suggested','Precursor m/z'])
    # return(istd_found)
    for file in tqdm(all_lst, total = len(all_lst)):
        current_file = os.path.join(peak_list_dir, file+'.txt')
        test_file = readin_peak_list(current_file, msial=False)
        if file != rt_seed:
            raw_rt, raw_intensity = find_istd(istd_found, test_file)
            rt_report[file]=raw_rt
            intensity_report[file]=raw_intensity
            if pd.Series(raw_rt).count() >= len(istd_info)*0.8:
                seed_rt_temp = []
                raw_rt_temp = []
                for i in range(0, len(raw_rt)):
                    if raw_rt[i] == raw_rt[i]:
                        seed_rt_temp.append(seed_rt[i])
                        raw_rt_temp.append(raw_rt[i])
                rt_offset =  [seed - raw for (seed, raw) in zip(seed_rt_temp, raw_rt_temp)]
                model_temp = RandomForestRegressor()
                model_temp.fit(np.array(raw_rt_temp).reshape(-1,1), rt_offset)
                offset_pred = model_temp.predict(np.array(test_file['RT (min)']).reshape(-1,1))
                normalized_rt = [off+raw for (off, raw) in zip(offset_pred, test_file['RT (min)'])]
                test_file.insert(test_file.columns.get_loc('RT (min)'), 'RT_adjusted', normalized_rt)
            else:
                test_file.insert(test_file.columns.get_loc('RT (min)'), 'RT_adjusted', test_file['RT (min)'])
        test_file.to_csv(os.path.join(normalized_peak_list_dir, file+'.csv'), index = False)
        # return(rt_report, intensity_report)
    cv= []
    for index, row in intensity_report.iterrows():
        cv.append(row[2:].std()/row[2:].mean()*100)
    intensity_report.insert(2, 'cv', cv)
    rt_report.to_csv(os.path.join(result_dir, 'rt_report.csv'), index= False)
    intensity_report.to_csv(os.path.join(result_dir, 'intensity_resport.csv'), index= False)
    # return(rt_report, intensity_report)



def find_seeded(seeded, test_file):
    test_rt = []
    test_intensity = []
    # print('i am in new')
    for index, row in seeded.iterrows():
        istd_temp = find_feature(test_file, mz = row['pmz_protonated'], rt = row['RT'], mz_column= 'Precursor m/z', rt_column='RT_adjusted')
        if len(istd_temp)>0:
            istd_temp.sort_values(by = 'Height', ascending=False, inplace=True)
            test_rt.append(istd_temp.iloc[0]['RT_adjusted'])
            test_intensity.append(istd_temp.iloc[0]['Height'])
    return(test_rt, test_intensity)
    # pass



def find_istd(istd_info, sample_file, rt_column = 'RT (min)', return_adjusted_rt = False, echo= False):
    seed_rt = []
    seed_intensity = []
    # print('i am in new')
    for index, row in istd_info.iterrows():
        istd_temp = find_feature(sample_file, mz = row['Precursor m/z'], rt = row['RT_suggested'], mz_column= 'Precursor m/z', rt_column=rt_column)
        if len(istd_temp)>0:
            istd_temp.sort_values(by = 'Height', ascending=False, inplace=True)
            if return_adjusted_rt == False:
                seed_rt.append(istd_temp.iloc[0]['RT (min)'])
            else:
                seed_rt.append(istd_temp.iloc[0]['RT_adjusted'])
            seed_intensity.append(istd_temp.iloc[0]['Height'])
        else:
            # istd_temp = find_feature(sample_file, mz = row['Precursor m/z 2'], rt = row['RT_suggested'], mz_column= 'Precursor m/z', rt_column=rt_column)
            # if istd_temp>0:
            #     istd_temp.sort_values(by = 'Height', ascending=False, inplace=True)

            seed_rt.append(np.NAN)
            seed_intensity.append(np.NAN)
    if echo==True:
        print('the number of istd found is: ', str(pd.Series(seed_rt).count()))
        print('the percent of istd found is: ', str(np.round(pd.Series(seed_rt).count()/len(istd_info)*100,2)))
    return(seed_rt, seed_intensity)
def find_feature(feature_table, mz, rt, mz_column = 'mz', rt_column = 'rt_apex', intensity_column =None, rt_toerlance = 10):
    # print('i am nin new')
    feature_mz_search = quick_search_values(feature_table, mz_column, mz-0.005, mz+0.005, ifsorted=False)
    # return(feature_mz_search)
    feature_mzrt_search = quick_search_values(feature_mz_search, rt_column, rt-10/60, rt+10/60, ifsorted=False)
    if intensity_column != None:
        feature_mzrt_search.sort_values(by = intensity_column, inplace=True, ascending=False)
    # print(feature_mzrt_search)
    return (feature_mzrt_search)