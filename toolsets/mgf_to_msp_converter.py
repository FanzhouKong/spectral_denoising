import os
import sys


import numpy as np

import toolsets.spectra_operations as so

import toolsets.mgf_tools as mgfc

import pandas as pd

def main():
    if len(sys.argv) <= 1:
        print("""
Usage: python mgf_containing_directory msp_output_directory.
        """)
        return
    folderpath = sys.argv[1] #this is where your mgf file locates

    for file in os.listdir(folderpath):
        if file.endswith(".mgf"):
            inputfile = (os.path.join(folderpath, file))
            outputfile =(os.path.join(sys.argv[2], file[:-3]+"msp"))
            mgfc.convert_mgf_to_msp(open(inputfile), outputfile)

if __name__ == '__main__':
    main()










