import pandas as pd
from molmass import Formula
from functools import reduce
import os
import numpy as np
from fuzzywuzzy import fuzz

import re

def specify_column(df, keyword):
    pre_matched = find_matching_columns(df, [keyword])
    if len(pre_matched)==0:
        return([])
    else:
        score = []
        for name in pre_matched:
            score.append(fuzz.token_sort_ratio(keyword,name))
        return(pre_matched[np.argmax(score)])

def find_adducts(texts):
    pattern = r'\[M([+-][A-Za-z0-9]+)*\][+-]?'
    returning = []
    for text in texts:
        matches = re.findall(pattern, text)
        if matches:
            # Reconstruct full adduct patterns from captured groups
            full_matches = [text[s.start():s.end()] for s in re.finditer(pattern, text)]
            returning.append(full_matches[0])

    return returning



def find_matching_columns(df, keywords):
    matched_columns = []
    for column in df.columns:
        # Normalize column name to lowercase
        normalized_column = column.lower()
        # Check each keyword
        for keyword in keywords:
            if keyword in normalized_column:
                matched_columns.append(column)
                break  # Break to avoid adding the same column multiple times if multiple keywords match
    return matched_columns






