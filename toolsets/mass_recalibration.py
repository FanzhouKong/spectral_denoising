#!/usr/bin/env python
# coding: utf-8

# In[10]:


from sklearn.linear_model import LinearRegression
import numpy as np
import pandas as pd
from tqdm import tqdm
from toolsets.spectra_operations import break_spectra, pack_spectra
from toolsets.search import string_search, num_search
def data_recalibrated(data, mix_col = 'reference_mix', peaks_col = 'peaks'):
    peaks_recalibrated = []
    data.sort_values(by = 'reference_mix', inplace = True)
    for mix in tqdm(data[mix_col].unique()):

        data_temp = string_search(data, mix_col, mix)
        peaks_recalibrated.extend(mix_recalibrate(data_temp, peaks_col = peaks_col))
    peaks_recalibrated_col = peaks_col+'_recalibrated'
    data[peaks_recalibrated_col]=peaks_recalibrated
    return(data)
def mix_recalibrate(mix, peaks_col = 'peaks', x_cols = ['precursor_mz','peak_apex_intensity', 'rt']):
    x = mix[x_cols]
    scaled_x = scalar(x)
    error = mix['reference_precursor_mz']-mix['precursor_mz']
    # model = LinearRegression()
            # model = linear_model(scalar(x), error)
    model = rf_model(scaled_x, error)
    peaks_recalibrated = []
    for index, row in mix.iterrows():
        mass_temp, intensity_temp = break_spectra(row[peaks_col]) 
        mass_recalibrated = []
        x_temp_dict = {'precursor_mz': mass_temp, 'peak_apex_intensity': [scaled_x['peak_apex_intensity'][index]]*len(mass_temp), 'retention_time': [scaled_x['rt'][index]]*len(mass_temp)}
        x_temp = pd.DataFrame.from_dict(x_temp_dict)
        error_pred = make_prediction(x_temp, model)
        mass_recalibrated = [ mass_temp[j] + error_pred[j] for j in range (len (mass_temp))]  
        peaks_recalibrated.append(pack_spectra(mass_recalibrated, intensity_temp))
    return(peaks_recalibrated)
        # print( x['ms1_precursor_intensity'][index])
        # for i in range(len(mass_temp)):
        #     x_temp = [[mass_temp[i], scaled_x['ms1_precursor_intensity'][index], scaled_x['retention_time'][index]]]
        #     mass_recalibrated.append(model.predict([[x['ms1_precursor_intensity'][index], mass_temp[i]]])[0]) 
        # peaks_recalibrated.append(pack_spectra(mass_recalibrated, intensity_temp))
    # mix['peaks_recalibrated']=peaks_recalibrated
    # return(peaks_recalibrated)
def find_diff_error(data, mix_col = 'reference_mix', cols = ['precursor_mz','peak_apex_intensity', 'rt']):
    diff = []
    # print('i am using rf method in find diff error function')
    # data_raw = data.copy()
    for mix in tqdm(data[mix_col].unique()):
        try:
            data_temp = string_search(data, mix_col, mix)
            # return(data_temp)
            x = data_temp[cols]
            # return()
            error = data_temp['reference_precursor_mz']-data_temp['precursor_mz']
            # model = linear_model(scalar(x), error)
            model = rf_model(scalar(x), error)
            error_pred = make_prediction(scalar(x), model)
            y_pred = data_temp['precursor_mz']+error_pred
            diff.extend(y_pred-data_temp['reference_precursor_mz'])
            # return(diff)

        except:
            x = data_temp[cols]
            y = data_temp['reference_precursor_mz']
            model = linear_model(x, y)
            return(x,y)
    return(diff)

def find_diff(data, mix_col = 'reference_mix', scalar = True,cols = ['ms1_precursor_intensity', 'precursor_mz']):
    diff = []
    print('i am using rf method!')
    # data_raw = data.copy()
    for mix in tqdm(data[mix_col].unique()):
        try:
            data_temp = string_search(data, mix_col, mix)
            # return(data_temp)
            x = data_temp[cols]
            # return()
            y = data_temp['reference_precursor_mz']
            if scalar ==True:
                # model = linear_model(scalar(x), y)
                model = rf_model(scalar(x),y)
                y_pred = make_prediction(scalar(x), model)
                diff.extend(y_pred-y)
            else:
                model = linear_model(x, y)
                y_pred = make_prediction(x, model)
                diff.extend(y_pred-y)
            # return(diff)

        except:
            x = data_temp[cols]
            y = data_temp['reference_precursor_mz']
            model = linear_model(x, y)
            return(x,y)
    return(diff)
def scalar(x):
    x_scaled = x.copy()
    for i in range(len(x_scaled.columns)):
        if x_scaled.columns[i]!='precursor_mz':
            x_scaled[x_scaled.columns[i]]=x_scaled[x_scaled.columns[i]]/x_scaled[x_scaled.columns[i]].max()
    return(x_scaled)
def linear_model(x, y):
    # x_scaled = scalar(x)
    model = LinearRegression()
    model.fit(x, y)
    return(model)
    # y_pred = []
from sklearn.ensemble import RandomForestRegressor
def rf_model(x,y):
    model =  RandomForestRegressor(n_estimators = 100, random_state = 0)
    model.fit(x,y)
    return(model)
  
def make_prediction(x, model): 
    x.reset_index(inplace = True, drop = True)
    y_pred = []
    for i in range(len(x)):
        x_temp = []
        for j in range(len(x.columns)):
            x_temp.append(x[x.columns[j]][i])
        # return(x_temp)

        y_pred.append(model.predict([x_temp])[0])
    return(y_pred)


# def data_recalibrate(data, save_diff = False):
#     # automatically recalibrate the whole dataset, based on all mixes
#     data_recalibrated = pd.DataFrame()
#     diff_raw = []
#     diff_recalibrated = []
#     # comments = []
#     # parent_ion = []
#     for mix in tqdm(data['reference_mix'].unique()):
#         # poly = PolynomialFeatures(degree = 3)
#         data_temp = string_search(data, "reference_mix", mix)
#         # data_temp = data.loc[data['reference_Mix label']==i]
#         x_temp = data_temp['precursor_mz'].values
#         y_temp =  data_temp['reference_precursor_mz'].values
#         diff_raw.extend((x_temp-y_temp).tolist())
#         model_temp = fit_model(x_temp, y_temp)
#         y_pred = model_temp.predict(x_temp.reshape(-1,1))
#         if len(x_temp)<=1:
            
#             data_temp['calibration']='not recalibrated_'
#         else:
#             # y_pred = model_temp.predict(poly.fit_transform(x_temp.reshape(-1,1)))
#             data_temp['calibration']='recalibrated_'
#         diff_recalibrated.extend((y_pred-y_temp).tolist())
#         # parent_ion.extend(y_pred.tolist())
#         msms_recalibrated = []
#         for n in range(len(data_temp)):

#             msms_recalibrated.append(mass_recalibrate(model_temp, data_temp.iloc[n]['peaks'] ))
#         data_temp['peaks_recalibrated']=msms_recalibrated
#         data_recalibrated = pd.concat([data_recalibrated, data_temp], ignore_index = True, axis = 0)
#     if save_diff == True:
#         data_recalibrated['diff_raw']=diff_raw
#         data_recalibrated['diff_recalibrated']=diff_recalibrated
#     # data_recalibrated['parent_ion']=parent_ion
#     # data_recalibrated['comments']=comments
#     data_recalibrated.reset_index(inplace=True, drop=True)
#     return(data_recalibrated)




# def mass_predicting(mass_measured, model):
#     # predicting mass values using constructed linear model
#     mass_cali = model.predict(np.array(mass_measured).reshape(-1,1))

    
#     mass_cali = [round(num, 6) for num in mass_cali]
# #     mass_cali = [str(integer) for integer in mass_cali]
#     return(mass_cali)
    

# from sklearn.preprocessing import PolynomialFeatures
# from sklearn.linear_model import LinearRegression
# def fit_model(x, y):
#     # fit the measured (averagemz) to actual value (precursormz) using a linear model
#     if(len(x)>1):
#         pass
#     else:
#         # if there is only 1 compound in a mix, do not change it but just return the original value
#         x = np.array([1,2,3])
#         y = np.array([1,2,3])

#         # return(model)
#         # print("i am in else")
#     # poly = PolynomialFeatures(degree = 3)
#     # X_poly = poly.fit_transform(x.reshape(-1,1))
#     # poly.fit(X_poly, y)

#     # model = LinearRegression()
#     # model.fit(X_poly, y)
#     model = LinearRegression()
#     model.fit(x.reshape(-1,1),y)
#     return (model)
    
# def mass_recalibrate(model, msms):
#     # recalibrate a single mix
#     mass, intensity = so.break_spectra(msms)
#     mass_recali = mass_predicting(mass, model)
#     msms_recali = so.pack_spectra(mass_recali, intensity)
#     return(msms_recali)
        





