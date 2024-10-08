====================================
Spectral denoising: useful functions
====================================

Here we provide a list of useful functions in the project for read and write data in differnt formats. We also provide functions for manipulating, visualizing, and comparing spectra.

Read data
---------

``file_io.read_msp`` will read in all standard MSP files into a pandas.Dataframe, with basic spectra cleaning performed (sorted + zero intensity ions removed).

.. autofunction:: spectral_denoising.file_io.read_msp
    :noindex:

``file_io.read_df`` will read in a pandas.Dataframe from a csv file, and convert any stringed spectra column (m/z and intensity are separated by '\t' and new fragment are separated by '\n') to numpy array.

.. autofunction:: spectral_denoising.file_io.read_df
    :noindex:

Write data
----------

``file_io.write_to_msp`` will write a pandas.Dataframe to a MSP file with output location specified. The target spectra column to be exported to MSP file should also be specified since each dataframe file could contain multiple versions of MS/MS spectra.

.. autofunction:: spectral_denoising.file_io.write_to_msp
    :noindex:

``file_io.save_df`` will write a pandas.Dataframe to a csv file with output location specified. All spectra column will be automatically converted to string format for saving.

.. autofunction:: spectral_denoising.file_io.save_df
    :noindex:

Manipulating spectra
----------------------

``spectral_operations.break_spectrum`` will break a given np ndarray spectra into 2 np ndarrays of mass and intensities.

.. autofunction:: spectral_denoising.spectral_operations.break_spectrum
    :noindex:

``spectral_operations.pack_spectrum`` do the reverse operation of ``break_spectrum``. It will pack 2 np ndarrays of mass and intensities into a single np ndarray with shape [n, 2].

.. autofunction:: spectral_denoising.spectral_operations.pack_spectrum
    :noindex:

``spectral_operations.normalize_spectrum`` normalized the query spectrum so that the sum of intensities is 1.

.. autofunction:: spectral_denoising.spectral_operations.normalize_spectrum
    :noindex:

``spectral_operations.standardize_spectrum`` standardize the query spectrum so that the base peak has intensity of 1.

.. autofunction:: spectral_denoising.spectral_operations.standardize_spectrum
    :noindex:

``spectral_operations.sort_spectrum`` sort the query spectrum by mass.

.. autofunction:: spectral_denoising.spectral_operations.sort_spectrum
    :noindex:

``spectral_operations.remove_precursor`` will remove the precursor ion region from the query spectrum. If no pmz is provided, the function will use the max m/z in the query spectrum as the precursor ion.

.. autofunction:: spectral_denoising.spectral_operations.remove_precursor
    :noindex:

``spectral_operations.remove_zero_ions`` will remove ions with zero intensity from the query spectrum.

.. autofunction:: spectral_denoising.spectral_operations.remove_zero_ions
    :noindex:

``spectral_operations.sanitize_spectrum`` is just a wrapper function for ``spectral_operations.remove_zero_ions`` and ``spectral_operations.sort_spectrum``.

.. autofunction:: spectral_denoising.spectral_operations.sanitize_spectrum
    :noindex:

``spectral_operations.truncate_spectrum`` will truncate the query spectrum up to the max_mz provided. This is the function used in ``remove_precursor`` function.

.. autofunction:: spectral_denoising.spectral_operations.truncate_spectrum
    :noindex:

``spectral_operations.slice_spectrum`` will slice the query spectrum at the break_mz.

.. autofunction:: spectral_denoising.spectral_operations.slice_spectrum
    :noindex:

.. autofunction:: spectral_denoising.spectral_operations.pack_spectrum
    :noindex:

``spectral_operations.add_spectra`` will naively add 2 spectra together.

Formating spectra
------------------

``spectral_operations.arr_to_str`` will format the np ndarray query spectrum into a string format with m/z and intensity separated by '\t' and new fragment separated by '\n'.
Reverse of ``str_to_arr``.

.. autofunction:: spectral_denoising.spectral_operations.arr_to_str
    :noindex:

``spectral_operations.str_to_arr`` will format the stringed query spectrum with m/z and intensity separated by '\t' and new fragment separated by '\n' into a np ndarray format. 
Reverse of ``arr_to_str``.

.. autofunction:: spectral_denoising.spectral_operations.str_to_arr
    :noindex:

Comparing spectra and entropy related functions
------------------------------------------------

The package provides entropy-based similarity measures to compare 2 spectra. While spectral entropy can also be used to measure spectral information content and spectral quality.

``spectral_operations.entropy_similairty`` measures the entropy similarity between 2 spectra. If pmz is provided, the precursor region will be removed before calculating the similarity.

.. autofunction:: spectral_denoising.spectral_operations.entropy_similairty
    :noindex:

``spectral_operations.spectral_entropy`` measures the entropy of a given spectrum.

.. autofunction:: spectral_denoising.spectral_operations.spectral_entropy
    :noindex:

``spectral_operations.normalized_entropy`` will calculate the normalized entropy of the query spectrum. It is a coarse measure of the spectral quality.

.. autofunction:: spectral_denoising.spectral_operations.normalized_entropy
    :noindex:

Visualizing spectra
---------------------

``spectra_plotter.head_to_tail_plot`` will plot the head-to-tail plot of the query spectra. If pmz is given, the precursor region will be removed before plotting and precursor m/z will be marked as grey dashed line.
If savepath is provided, the plot will be saved to the location.

.. autofunction:: spectral_denoising.spectra_plotter.head_to_tail_plot
    :noindex:

``spectra_plotter.ms2_plot`` will plot the query spectrum. It takes similary parameters as ``head_to_tail_plot``.

.. autofunction:: spectral_denoising.spectra_plotter.ms2_plot
    :noindex:


Generating noise
-----------------

``noise.generate_noise`` will generate synthetic noise with m/z from 50 to pmz. The m/z follows random distribution and the intensity follows Possion distribution with lamda provided.
The parameter n is the number of noise ions to be generated.

.. autofunction:: spectral_denoising.noise.generate_noise
    :noindex:

``noise.generate_chemical_noise`` will generate synthetic chemical noise with m/z from 50 to pmz. The m/z are randomly sampled from the formula_db provided. 
The intensity is determined in the similar manner as ``generate_noise``.

.. autofunction:: spectral_denoising.noise.generate_chemical_noise
    :noindex:

``noise.add_noise`` will add noise to the query spectrum. The noise is generated by ``generate_noise`` or ``generate_chemical_noise``.

.. autofunction:: spectral_denoising.noise.add_noise
    :noindex:



