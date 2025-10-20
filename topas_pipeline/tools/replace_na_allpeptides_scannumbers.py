import sys
import re
import argparse

import numpy as np
import pandas as pd



for folder in meta_input_df['mq_txt_folder']:
    matches = re.findall(r'Batch([A-Za-z]*\d+)', folder)
            

def main(argv):

    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        "--config",
        dest="config",
        required=True,
        help="Absolute path to configuration file.",
    )



    allpep = pd.read_csv(folder + '/allPeptides.txt', sep='\t')
    if allpep['Min scan number'].isna().any() or allpep['Max scan number'].isna().any():

        msms = pd.read_csv(folder + '/msms.txt', sep='\t')
        msms_max_scans = msms.groupby('Raw file')['Precursor full scan number'].max()

        allpep.loc[allpep['Max scan number'].isna(), 'Max scan number'] = allpep.loc[allpep['Max scan number'].isna(), 'Raw file'].map(msms_max_scans)
        allpep.loc[allpep['Min scan number'].isna(), 'Min scan number'] = 1
        allpep.to_csv(folder + '/allPeptides.txt', sep='\t', index=False)



"""
python3 -m topas_pipeline.tools.replace_na_allpeptides_scannumbers -c config_patients.json
"""
if __name__ == "__main__":
    main(sys.argv[1:])
