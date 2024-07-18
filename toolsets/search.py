import pandas as pd
import numpy as np
import numexpr
import toolsets.helpers as helpers
def search_feature(feature, pmz=None, rt = None, mass_error =0.005, rt_error = 1.5, pmz_sorted = False):
    pmz_col = helpers.specify_column('precursor mz', feature.columns)
    rt_col =helpers.specify_column('rt apex', feature.columns)
    if pmz is not None and rt is None:
        if pmz_sorted==False:
            pmz_match = quick_search_values(feature, pmz_col, pmz-mass_error, pmz+mass_error)
        else:
            pmz_match = quick_search_sorted(feature, pmz_col, pmz-mass_error, pmz+mass_error)
        return pmz_match
    elif pmz is None and rt is not None:
        rt_match = quick_search_values(feature, rt_col, rt-rt_error/60, rt+rt_error/60)
        return rt_match
    elif pmz is not None and rt is not None:
        if pmz_sorted==False:
            pmz_match = quick_search_values(feature, pmz_col, pmz-mass_error, pmz+mass_error)
        else:
            pmz_match = quick_search_sorted(feature, pmz_col, pmz-mass_error, pmz+mass_error)
        rt_match = quick_search_values(pmz_match, rt_col, rt-rt_error/60, rt+rt_error/60)
        return rt_match
    else:
        print('both pmz and rt is None')
        return()
def search_feature_msdial(feature, pmz=None, rt = None, mass_error =0.005, rt_error = 1.5, pmz_sorted = False, if_alignment = True):
    if if_alignment:
        pmz_col = helpers.specify_column('Average Mz', feature.columns)
        # rt_col =helpers.specify_column('Average Rt(min)', feature.columns)
        rt_col =helpers.specify_column('rt_corrected', feature.columns)
    else:
        pmz_col = helpers.specify_column('Precursor m/z', feature.columns)
        rt_col =helpers.specify_column('RT (min)', feature.columns)
    # print(pmz_col, rt_col)
    if pmz is not None and rt is None:
        if pmz_sorted==False:
            pmz_match = quick_search_values(feature, pmz_col, pmz-mass_error, pmz+mass_error)
        else:
            pmz_match = quick_search_sorted(feature, pmz_col, pmz-mass_error, pmz+mass_error)
        return pmz_match
    elif pmz is None and rt is not None:
        rt_match = quick_search_values(feature, rt_col, rt-1.5/60, rt+1.5/60)
        return rt_match
    elif pmz is not None and rt is not None:
        if pmz_sorted==False:
            pmz_match = quick_search_values(feature, pmz_col, pmz-0.005, pmz+0.005)
        else:
            pmz_match = quick_search_sorted(feature, pmz_col, pmz-mass_error, pmz+mass_error)
        rt_match = quick_search_values(pmz_match, rt_col, rt-rt_error/60, rt+rt_error/60)
        return rt_match
    else:
        print('both pmz and rt is None')

        return()
def string_search(data, column_name,item, reset_index = True,reverse = False):
    if reverse == False:
        _data= data[data[column_name].to_numpy() == item]
        # return data[data[column_name].to_numpy() == item]
    else:
        _data= data[data[column_name].to_numpy() != item]
    if reset_index == True:
        _data.reset_index(inplace= True, drop = True)
    return(_data)
        # return data[data[column_name].to_numpy() != item]
def quick_search_sorted(data_raw, column_name,value_start, value_end):
    # data.sort_values(by=column_name, inplace = True)
    search_array=data_raw[column_name].to_numpy(dtype="float")
    # data = data_raw.copy()

    # index_start, index_end = search_array.searchsorted([value_start, value_end])
    index_start = np.searchsorted(search_array, value_start,side = 'left')
    index_end = np.searchsorted(search_array, value_end,side = 'right')
    return(data_raw.iloc[index_start:index_end])
def quick_search_values(data_raw, column_name,value_start, value_end):

    data_sorted = data_raw.sort_values(by=column_name)
    data_return = quick_search_sorted(data_sorted, column_name, value_start, value_end)
    # index_start = np.searchsorted(data[column_name], value_start,side = 'left')
    # index_end = np.searchsorted(data[column_name], value_end,side = 'right')
    return(data_return)

def num_search(data, column_name,number, direction, step = None,inclusion = False):
    x = data[column_name].values
    if direction == ">":
        if inclusion == False:
            return(data[numexpr.evaluate('(x > number)')])
        else:
            return(data[numexpr.evaluate('(x >= number)')])
    elif direction == '<':
        if inclusion == False:
            return(data[numexpr.evaluate('(x < number)')])
        else:
            return(data[numexpr.evaluate('(x <= number)')])

    elif direction == '==':
        return(data[numexpr.evaluate('(x == number)')])
    elif direction =='between' and step != None:
        if inclusion == False:
            temp = data[numexpr.evaluate('(x > number-step)')]
            x = temp[column_name].values
            return (temp[numexpr.evaluate('(x < number+step)')])
        else:
            temp = data[numexpr.evaluate('(x >= number-step)')]
            x = temp[column_name].values
            return (temp[numexpr.evaluate('(x <= number+step)')])
    else:
        print('the wrong method is passed')



