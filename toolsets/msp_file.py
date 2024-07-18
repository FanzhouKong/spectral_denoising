import logging
import numpy as np


def _parse_information(line: bytes, spec: dict) -> str:
    """
    Parse the line in .msp file, update information in the spec
    :param line: The input line.
    :param spec: The output information collection.
    :return: The entry of the added information.
    """
    line = line.strip()
    if line:
        lines = line.split(":")
        if len(lines) > 2:
            item = lines[0]
            cont = ":".join(lines[1:])
        else:
            item, cont = lines

        item = item.strip()
        spec[item] = cont.strip()
        return item
    else:
        return ""


def _parse_spectrum(line: str, spec: list) -> int:
    """
    Add peak data to spec
    :param line: The raw line information from .msp file
    :param spec: The spectrum will be added.
    :return: 0: success. 1: no information in this line.
    """
    line = line.strip()
    lines = line.split()
    if len(lines) >= 2:
        mz, intensity = lines[0], lines[1]
        spec.append([float(mz), float(intensity)])
        return 0
    else:
        return 1


def read(stream_input) -> dict:
    """
    Read information from .msp file.
    :param stream_input: a stream for input.
    :return: a dict contains a list with key 'spectra'.
        The list contains multiple dict, one dict represent a single spectrum's informaiton.
    """

    exp = {
        'spectra': []
    }

    for spectrum_info in read_one_spectrum(stream_input):
        exp['spectra'].append(spectrum_info)
    return exp


def read_one_spectrum(fi, include_raw=0) -> dict:
    """
    Read one spectrum from .msp file.
    :param stream_input: a stream for input.
    :param include_raw: whether output raw spectrum or not.
    :return: a dict contains one spectrum information.
    """

    spectrum_info = {
        'spectrum': []
    }
    is_adding_information = 1
    raw = []
    for line in fi:
        if not isinstance(line, str):
            line = line.decode()

        if include_raw:
            raw.append(line)

        if is_adding_information:
            item = _parse_information(line, spectrum_info)
            if item.startswith("Num") and item.lower() == "num peaks":
                spectrum_info[item] = int(spectrum_info[item])
                is_adding_information = 0
                peak_num = spectrum_info[item]
        else:
            spec = spectrum_info['spectrum']
            _parse_spectrum(line, spec)
            if len(spec) == peak_num:
                if include_raw:
                    raw.append("\n")
                    spectrum_info["raw"] = "".join(raw)
                    raw = []

                yield spectrum_info

                # Preparing adding the next one.
                is_adding_information = 1
                spectrum_info = {
                    'spectrum': []
                }


def write(experiment: dict) -> bool:
    pass
