import sys
import os
import pandas as pd

filenames =sys.argv
os.system("cd ../")
os.system("mkdir results")
counts = pd.read_csv(filenames[1], sep='\t', index_col=0)["NumReads"]

for f in filenames[2:]:
    counts = pd.concat([counts, pd.read_csv(f, sep='\t', index_col=0)["NumReads"]], axis=1)
counts = round(counts).astype(int)
counts.columns = list(map(lambda f: f.split("/")[1].rstrip("_quant"), filenames[1:]))
counts.to_csv("results/counts.tsv", sep='\t')
