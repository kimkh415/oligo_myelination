import scrublet as scr
import scipy.io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import numpy as np

import sys

args=sys.argv

matrix=args[1]
savename=args[2]  # prefix

print("Load Matrix!")
counts_matrix = scipy.io.mmread(matrix).T.tocsc()
print(counts_matrix.shape)

print("Make Scrublet Object")
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.05)
print("Run Scrublet")
doub_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2)

scrub.plot_histogram()
plt.savefig(savename+"_hist.png")

print("Save!")
np.savetxt(savename+"_score.txt",doub_scores)
np.savetxt(savename+"_pred.txt",predicted_doublets)

