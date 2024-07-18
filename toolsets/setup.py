import subprocess
import os

import shutil
def get_dirs(workplace, mode = None, if_print = True):
    first_dirs = ['alignment_result', 'bioactivity', 'mzml', 'normalized_peak_list',
                  'peak_list', 'results', 'sirius_files', 'uv_data', 'nmr_data', 'uv_lib', 'msms_lib', 'nmr_lib']
    return_dirs = []
    if mode is not None:
        for item in first_dirs:
            if item not in ['bioactivity', 'uv_data', 'nmr_data','uv_lib', 'msms_lib', 'nmr_lib']:

                if os.path.exists(os.path.join(workplace, item, mode))== True:
                    # print('i am in true if')
                    return_dirs.append(os.path.join(workplace, item, mode))
            else:
                if os.path.exists(os.path.join(workplace, item))== True:
                    return_dirs.append(os.path.join(workplace, item))
    if if_print == True:
        for i in range(len(return_dirs)):
            print(f'at position {str(i)}: '.format(),return_dirs[i])

    return return_dirs
def set_workspace(workplace, mode = None):
    if os.path.exists(workplace)==False:
        os.makedirs(workplace)
    else:
        print("the workspace exists!")
    # print('tt')
    first_dirs = ['alignment_result', 'bioactivity', 'mzml', 'normalized_peak_list',
                  'peak_list', 'results', 'sirius_files', 'uv_data', 'nmr_data', 'uv_lib', 'msms_lib', 'nmr_lib']
    second_dir = ['pos', 'neg']
    dir_to_return = []
    for item in first_dirs:
        if item == 'msms_lib':
            try:
                os.makedirs(os.path.join(workplace, item))
            except:
                pass

            try:
                os.makedirs(os.path.join(workplace, item, 'mzml'))
            except:
                pass
            for item_sec in second_dir:
                try:
                    os.makedirs(os.path.join(workplace, item, 'mzml', item_sec))
                except:
                    pass
                try:
                    os.makedirs(os.path.join(workplace, item, 'features', item_sec))
                except:
                    pass



        if item not in ['bioactivity', 'uv_data', 'nmr_data','uv_lib', 'msms_lib', 'nmr_lib']:
            if os.path.exists(os.path.join(workplace, item))== False:
                os.makedirs(os.path.join(workplace, item))
                for item_sec in second_dir:
                    os.makedirs(os.path.join(workplace, item, item_sec))
                    # if item_sec==mode:
                        # print('i am in if')
                        # dir_to_return.append(os.path.join(workplace, item, item_sec))
        else:
            if os.path.exists(os.path.join(workplace, item))== False:
                os.makedirs(os.path.join(workplace, item))
                # dir_to_return.append(os.path.join(workplace, item))
    print('set up complete')


def split_pos_neg(all_folder):
    pos_folder = os.path.join(all_folder, 'pos')
    neg_folder = os.path.join(all_folder, 'neg')
    for folder in [pos_folder, neg_folder]:
        if os.path.exists(folder)==False:
            os.makedirs(folder)

    # file_lists = []
    for root, dirs, files in os.walk(all_folder):
        for file in files:
            # print(file)

            if file.endswith('.mzML'):
                if len(file.split('.'))==2:
                    base_name = file.split('.')[0]
                    # print(base_name)
                    if base_name[-1]=='P':
                        shutil.move(os.path.join(all_folder, file), os.path.join(pos_folder, file))
                    elif base_name[-1]=='N':
                        shutil.move(os.path.join(all_folder, file), os.path.join(neg_folder, file))