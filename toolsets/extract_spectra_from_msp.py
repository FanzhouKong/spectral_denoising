#!/usr/bin/python3
import sys
import msp_file
import filename

"""
This script will extract spectra in a specific file which contains the text indicated.
"""


def main():
    if len(sys.argv) <= 1:
        print("""
Usage: python extract_spectra_from_msp.py filename_input filename_output ion_number.
        """)
        return
    filename_in = sys.argv[1]
    filename_out = sys.argv[2]
    ion_number = int(sys.argv[3])

    fi = filename.smart_io(filename_in, "rt")
    fo = filename.smart_io(filename_out, "wt")

    for spec in msp_file.read_one_spectrum(fi, include_raw=1):
        if len(spec['spectrum']) >= ion_number:
            fo.writelines(spec["raw"])


if __name__ == '__main__':
    main()
