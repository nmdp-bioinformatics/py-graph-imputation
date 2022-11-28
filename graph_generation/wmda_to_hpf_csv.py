import os
import pandas as pd

wmda_orig = pd.read_csv("data/wmda/freqs.txt", sep=";", names=["hap", "freq"])

# in case we want to change to read a gz version
# wmda_orig = pd.read_csv('data/wmda/freqs.txt.gz', sep=';', compression='gzip', names = ['hap', 'freq'])

project_dir = "../../"

wmda_orig["pop"] = "CAU"
wmda = wmda_orig[["hap", "pop", "freq"]]

os.makedirs("output", exist_ok=True)
pop_ratio_dir = project_dir + "imputation/graph_generation/output/pop_ratio.txt"
pop_ratio_file = open(pop_ratio_dir, "w")
pop_ratio_file.write("{},{},{}\n".format("CAU", 20 / (1e-5), 1))
wmda.to_csv("output/hpf.csv", index=False)
