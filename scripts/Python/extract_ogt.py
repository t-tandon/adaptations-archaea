import os
import pandas as pd
from glob import glob

#path to file
input_dir = "genomespot" 
output_file = "genomespot_combined.tsv"

files = glob(os.path.join(input_dir, "*.genomespot.predictions.tsv"))

records = []
for file in files:
    accession = os.path.basename(file).split(".")[0]  # strips file name, gives only accession
    df = pd.read_csv(file, sep="\t")
    # convert the 'value' column into a dictionary keyed by 'target'
    data = df.set_index("target")["value"].to_dict()
    data["accession"] = accession
    records.append(data)

# combine into a single df
combined_df = pd.DataFrame(records)

# reorder columns (accession first)
combined_df = combined_df.set_index('accession').reset_index()

combined_df.head()
