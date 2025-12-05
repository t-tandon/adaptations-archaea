import pandas as pd

#load table
df = pd.read_csv("GTDB-Archaea-domain-GTDB-rep-metadata.tsv.complete.c90", delimiter="\t")

#filter by completeness and contamination 
filtered = df[
    (df["checkm_completeness"] > 90) &
    (df["checkm_contamination"] < 5) ] 

#choose one representative genome per genus
representatives = (filtered.sort_values("checkm_completeness", ascending=False).drop_duplicates(subset="genus", keep="first")) 

#save to file
representatives.to_csv("archaea_representatives.tsv", sep="\t", index=False)

#extract accessions
accessions = representatives["accession"].str.replace("RS_", "").str.replace("GB_", "")
accessions.to_csv("accessions.txt", index=False, header=False)
