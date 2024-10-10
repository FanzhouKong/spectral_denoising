    
import argparse
import spectral_denoising as sd

def main():
    print(sd.__spec__)
    parser = argparse.ArgumentParser()

    # Add arguments
    parser.add_argument('-i', '--input_path', type=str, required=True, help="The input msp file path")
    parser.add_argument('-r', '--reference_lib_path', type=str, required=True, help="The reference msp library file path")
    parser.add_argument('-o', '--output_directory', type=str, help="The location of the output directory that msp files goes to")
    parser.add_argument('-p', '--peaks_col', type=str, help="the column name of spectra column in reference file/ query file")
    parser.add_argument('-s', '--smiles_col', type=str, help="the column name of smiles column in reference file/ query file")
    parser.add_argument('-a', '--adduct_col', type=str, help="the column name of adduct column in reference file/ query file")
    parser.add_argument('-z', '--pmz_col',type = str, help="the column name of precursor_mz column in reference file/ query file")
    parser.add_argument('-m', '--mass_tolerance',type = float, help="the mass tolerance for spectral denoising")
    parser.add_argument('-d', '--identity_search_mass_tolerance',type = float, help="the mass tolerance for identity_search")
    # Parse the arguments

    args = parser.parse_args()
    peaks_col = args.peaks_col if args.peaks_col else 'peaks'
    smiles_col = args.smiles_col if args.smiles_col else 'smiles'
    adduct_col = args.adduct_col if args.adduct_col else 'adduct'
    pmz_col = args.pmz_col if args.pmz_col else 'precursor_mz'
    identity_search_mass_tolerance = args.identity_search_mass_tolerance if args.identity_search_mass_tolerance else 0.005
    mass_tolerance = args.mass_tolerance if args.mass_tolerance else 0.005
    query_spectra= sd.read_msp(args.input_path)
    reference_library =sd.read_msp(args.reference_lib_path)
    results = sd.denoising_search_batch(query_spectra[peaks_col], query_spectra[pmz_col],reference_library,
                                        identitiy_search_mass_error=identity_search_mass_tolerance,
                                        mass_tolernace = mass_tolerance,
                                        pmz_col=pmz_col,
                                        smiles_col=smiles_col,
                                        adduct_col=adduct_col,
                                        msms_col=peaks_col
                                        )
    output_dir = args.output_directory if args.output_directory else 'denoise_search_results'
    sd.export_denoising_searches(results, output_dir, top_n=1)
    
if __name__ == "__main__":
    main()
