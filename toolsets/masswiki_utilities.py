import requests
import pandas as pd
from tqdm import tqdm
import toolsets.spectra_operations as so
def get_spectra_id(master_id_name):
    url = 'http://masswiki.us-west-2.elasticbeanstalk.com/get/masswiki_data'
    data = {"wiki_id": master_id_name}
    response = requests.post(url, json=data)
    all_spectra_id=[]
    for c in response.json()['data']['spectra_data']:
        all_spectra_id.append(c['wiki_id'])
    return all_spectra_id
def get_all_spectra(spectra_id_list):
    spectra_all = pd.DataFrame()
    for i in tqdm(spectra_id_list):
        spectra_all = pd.concat([spectra_all, get_one_spectra(i)], ignore_index=True)
    return spectra_all
def get_one_spectra( spectra_id, debug = False):
    url = 'http://masswiki.us-west-2.elasticbeanstalk.com/get/masswiki_data'
    data = {"wiki_id": spectra_id}
    response = requests.post(url, json=data)
    if debug == True:
        return response
    # return(response)
    data_dict = response.json()['data']
    data_dict['peaks']=so.convert_nist_to_string(data_dict['peaks'])
    dct = {k:[v] for k,v in data_dict.items()}
    dct = pd.DataFrame.from_dict(dct)
    if 'user_annotation' in response.json()['analysis'].keys():
        anno_dict = pd.DataFrame.from_dict(response.json()['analysis']['user_annotation'])
        # anno_dict = anno_dict[['method', 'rt', 'precurosr_mz', 'peaks', 'smiles', 'adduct', 'comments']]
        anno_dict.columns='matched_'+ anno_dict.columns
        dct = pd.concat([dct] * len(anno_dict), ignore_index=True)
        dct=pd.concat([dct,anno_dict], axis=1)
    #
    return dct