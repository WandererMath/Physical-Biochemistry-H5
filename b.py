from pyteomics import mass
from pyteomics import mgf # function that will help read in the .msg file
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import requests

mgf_file = "exp.mgf"
mgf_pre="pre.mgf"



spectra_list_pre = []


with mgf.read(mgf_pre) as spectra:
    for spectrum in spectra:
        spectra_list_pre.append(spectrum)
mgf_pre="pre.mgf"



spectra_list_pre = []


with mgf.read(mgf_pre) as spectra:
    for spectrum in spectra:
        #breakpoint()
        spectra_list_pre.append(spectrum)

def compute_fragment_ions(peptide, charge=1):
    """Compute theoretical b-ion and y-ion masses for a peptide."""
    b_ions = []
    y_ions = []

    for i in range(1, len(peptide)):
        b_ions.append(mass.calculate_mass(sequence=peptide[:i], ion_type='b', charge=charge))
        y_ions.append(mass.calculate_mass(sequence=peptide[i:], ion_type='y', charge=charge))

    return b_ions, y_ions



def find_closest_ion(mz_value, b_ions, y_ions, tolerance=0.5):
    """Find the closest b-ion or y-ion to the given m/z value."""
    all_ions = {'b': b_ions, 'y': y_ions}
    closest_ion = None
    min_diff = float("inf")

    for ion_type, ion_list in all_ions.items():
        for idx, ion_mz in enumerate(ion_list):
            diff = abs(ion_mz - mz_value)
            if diff < min_diff and diff <= tolerance:
                min_diff = diff
                closest_ion = f"{ion_type}{idx+1}"

    return closest_ion
for i, spectrum in enumerate(spectra_list_pre):
    seq=spectrum['params']['title'][:-3]
    b_ions, y_ions = compute_fragment_ions(seq)
    experimental_intensity = spectrum['intensity array']

    # Find the most intense peak
    max_intensity_index = np.argmax(experimental_intensity)
    experimental_mz = spectrum['m/z array']
    most_intense_mz = experimental_mz[max_intensity_index]
    most_intense_ion = find_closest_ion(most_intense_mz, b_ions, y_ions)

    print(f"Spectrum {i+1} : Most intense peak at {most_intense_mz:.4f} m/z "
          f" → Identified as {most_intense_ion}")