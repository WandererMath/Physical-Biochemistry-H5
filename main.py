from pyteomics import mgf # function that will help read in the .msg file
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import requests


mgf_file = "exp.mgf"
mgf_pre="pre.mgf"


# read the .mgf file with Pyteomics
spectra_list = []
spectra_list_pre = []
with mgf.read(mgf_file) as spectra:
    for spectrum in spectra:
        spectra_list.append(spectrum)

with mgf.read(mgf_pre) as spectra:
    for spectrum in spectra:
        spectra_list_pre.append(spectrum)



k=0
# get experimental mz and intensity data
for i, _ in enumerate(spectra_list):
    for j, _ in enumerate(spectra_list_pre):
        k+=1
        #breakpoint()
        experimental_mz = spectra_list[i]['m/z array']
        experimental_intensity = spectra_list[i]['intensity array']
        experimental_intensiry_normalized = experimental_intensity / np.max(experimental_intensity) # normalize intensities to max intensity = 1

        # get predicted mz and intensity data
        predicted_mz = spectra_list_pre[j]['m/z array']
        predicted_intensity = spectra_list_pre[j]['intensity array'] # no normalization needed because the predicted spectrum is already normalized

        fig, ax = plt.subplots(dpi = 150)

        # plot experimental data
        ax.vlines(experimental_mz, 0, experimental_intensiry_normalized, color='cornflowerblue', label='experimental')

        # plot predicted data
        ax.vlines(predicted_mz, 0, -1*predicted_intensity, color='indianred', label='predicted') # multiply intensities by -1 so they plot below

        # add line at y = 0
        ax.plot([0, 1500], [0, 0], color='k')

        # set axes limits
        ax.set_xlim(100, 1300) # change axis as needed
        ax.set_ylim(-1.05, 1.05)

        # label axes
        ax.set_xlabel('m/z')
        ax.set_ylabel('Normalized intensity')

        # show legend
        ax.legend(loc = 'lower left')

        # add title
        ax.set_title(f'EIGINAYGHR +2 Butterfly Plot experiement {i+1} vs predicted {j+1}')

        plt.tight_layout()

        plt.savefig(f"{k}.pdf")
        plt.clf()




