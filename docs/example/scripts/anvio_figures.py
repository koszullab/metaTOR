#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import metator.io as mio
import metator.figures as mtf
import pandas as pd
import sys
from os.path import join 

outdir = sys.argv[1]
contigs_data_file = sys.argv[2]
checkm_dir = sys.argv[3]
threshold = 500000

contigs_data = pd.read_csv(contigs_data_file, sep="\t")

checkm_file = join(checkm_dir, "checkM_results_complete_2.txt")
checkm_taxonomy_file = join(checkm_dir, "1.txt")

checkm_summary = mio.read_results_checkm(checkm_file, checkm_taxonomy_file)

mtf.plot_figures(outdir, contigs_data, checkm_summary, threshold)
