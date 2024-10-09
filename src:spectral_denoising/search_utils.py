import pandas as pd
import numpy as np
import numexpr

def string_search(data, column_name,item, reset_index = True,reverse = False):
    def string_search(data, column_name, item, reset_index=True, reverse=False):
        """
        Searches for rows in a DataFrame where the specified column matches (or does not match) a given item.

        Parameters:
            data (pd.DataFrame): The DataFrame to search within.
            column_name (str): The name of the column to search.
            item (str): The item to search for in the specified column.
            reset_index (bool, optional): Whether to reset the index of the resulting DataFrame. Defaults to True.
            reverse (bool, optional): If True, returns rows where the column does not match the item. Defaults to False.
        Returns:
            pd.DataFrame: A DataFrame containing the rows that match (or do not match) the search criteria.
        """

    if reverse == False:
        _data= data[data[column_name].to_numpy() == item]
    else:
        _data= data[data[column_name].to_numpy() != item]
    if reset_index == True:
        _data.reset_index(inplace= True, drop = True)
    return(_data)
        # return data[data[column_name].to_numpy() != item]
def quick_search_sorted(data_raw, column_name,value_start, value_end):
    """
    Perform a quick search on a sorted column of a DataFrame to find rows within a specified range.

    Parameters:
        data_raw (pd.DataFrame): The input DataFrame containing the data to search.
        column_name (str): The name of the column to search within.
        value_start (float): The starting value of the range.
        value_end (float): The ending value of the range.
    Returns:
        pd.DataFrame: A DataFrame containing the rows where the values in the specified column fall within the given range.
    """

    search_array=data_raw[column_name].to_numpy(dtype="float")
    index_start = np.searchsorted(search_array, value_start,side = 'left')
    index_end = np.searchsorted(search_array, value_end,side = 'right')
    return(data_raw.iloc[index_start:index_end])
def quick_search_values(data_raw, column_name,value_start, value_end):
    """
    Perform a quick search on a DataFrame to find rows where the values in a specified column fall within a given range. Basically sorting the data first
    followed by quick_search_sorted.

    Args:
        data_raw (pd.DataFrame): The raw DataFrame to search.
        column_name (str): The name of the column to search within.
        value_start (numeric): The starting value of the range.
        value_end (numeric): The ending value of the range.
    Returns:
        pd.DataFrame: A DataFrame containing rows where the values in the specified column are within the range [value_start, value_end].
    """
    data_sorted = data_raw.sort_values(by=column_name)
    data_return = quick_search_sorted(data_sorted, column_name, value_start, value_end)
    # index_start = np.searchsorted(data[column_name], value_start,side = 'left')
    # index_end = np.searchsorted(data[column_name], value_end,side = 'right')
    return(data_return)

def num_search(data, column_name,number, direction, step = None,inclusion = False):
    """
    Perform a numerical search on a specified column of a DataFrame based on given criteria.

    Parameters:
        data (pd.DataFrame): The DataFrame to search within.
        column_name (str): The name of the column to perform the search on.
        number (float or int): The reference number for the search condition.
        direction (str): The direction of the comparison. Can be one of the following: '>', '<', '==', 'between'.
        step (float or int, optional): The step value for the 'between' direction. Default is None.
        inclusion (bool, optional): Whether to include the boundary values in the comparison. Default is False.
    Returns:
        pd.DataFrame: A DataFrame containing rows that match the search criteria.
    Raises:
        ValueError: If an invalid direction is provided.
    Examples:
        >>> num_search(df, 'age', 30, '>')
        >>> num_search(df, 'age', 30, 'between', step=5, inclusion=True)
    """

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



