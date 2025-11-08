from constants import *
import numpy as np
from MECA0029_Group16_1 import nat_freqs
import pandas as pd
from globalMassStiffMatrices import create_globalMass_and_globalStiffness
from scipy.sparse.linalg import eigsh


M_global, K_global = create_globalMass_and_globalStiffness()

eigvals, eigvecs = eigsh(K_global, k=desiredFreqNb,
                         M=M_global, sigma=0.0, which='LM')

nat_freqs = np.sqrt(eigvals)/(2*np.pi)

"""Sort natural frequencies and mode shapes"""
idx = np.argsort(nat_freqs)
nat_freqs = nat_freqs[idx]

freq_index = np.arange(1, len(nat_freqs) + 1)

data = {
    'Frequency Index': freq_index,
    'Natural Frequency (Hz)': nat_freqs
}
df = pd.DataFrame(data)

filename = f"convergencePt1/freq{elemPerBar}EPB.csv"
df.to_csv(filename, index=False)
