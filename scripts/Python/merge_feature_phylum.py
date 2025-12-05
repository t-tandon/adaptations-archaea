import pandas as pd

# load gtdb metadata
meta = pd.read_csv("archaea_representatives.tsv", sep="\t")

#clean accession
meta["clean_acc"] = meta["accession"].str.replace(r"^(GB_|RS_)", "", regex=True)
meta["clean_acc"] = meta["clean_acc"].str.replace(r"\.\d+$", "", regex=True)

# Keep only relevant 
meta = meta[["clean_acc", "phylum", "class"]].rename(columns={"clean_acc": "accession"})

#load feature df
df = pd.read_csv("../feature_set.tsv", sep = "\t")

#merge
df = df.merge(meta, on = "accession", how = "left")

#save df
df.to_csv("feature_set_phylum_class.tsv", sep = "\t", index=False)
