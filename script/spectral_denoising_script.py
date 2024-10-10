import numpy as np
import argparse
import spectral_denoising as sd

def main():
    parser = argparse.ArgumentParser()

    # Add arguments
    parser.add_argument('-i', '--input_path', type=str, required=True, help="The input msp file path")
    parser.add_argument('-o', '--output_path', type=str, help="The location of the output msp file")
    parser.add_argument('-p', '--peaks_col', type=str, help="the column name of spectra column")
    parser.add_argument('-s', '--smiles_col', type=str, help="the column name of smiles column")
    parser.add_argument('-a', '--adduct_col', type=str, help="the column name of adduct column")
    parser.add_argument('-m', '--mass_tolerance',type = float, help="the mass tolerance for spectral denoising")
    # Parse the arguments
    args = parser.parse_args()
    peaks_col = args.peaks_col if args.peaks_col else 'peaks'
    smiles_col = args.smiles_col if args.smiles_col else 'smiles'
    adduct_col = args.adduct_col if args.adduct_col else 'adduct'
    output_path = args.output_path if args.output_path else 'spectral_denoising_result.msp'
    mass_tolerance = args.mass_tolerance if args.mass_tolerance else 0.005
    query_data = sd.read_msp(args.input_path)
    query_spectra,query_smiles,query_adducts = query_data[peaks_col],query_data[smiles_col],query_data[adduct_col]

    denoised_spectra = sd.spectral_denoising_batch(query_spectra,query_smiles,query_adducts, mass_tolerance=mass_tolerance)

    query_data['peaks_denoised'] = denoised_spectra
    sd.write_to_msp(query_data, output_path, msms_col='peaks_denoised', normalize=True)

if __name__ == "__main__":
    main()