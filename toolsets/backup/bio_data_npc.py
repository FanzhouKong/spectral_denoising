def check_order(alignment, mix_name_pos):
    if len(alignment.columns)!= len(mix_name_pos):
        return(False)
    for i in range(0, len(alignment.columns)):
        if alignment.columns[i]!=mix_name_pos[i]:
            return(False)
    return(True)