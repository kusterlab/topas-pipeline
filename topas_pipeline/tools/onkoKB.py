import requests
import pandas as pd
from pathlib import Path
import re
import logging 
logger = logging.getLogger(__name__)
logging.basicConfig(filename='onkoKB.log', level=logging.INFO)

"""
idea: to make a pre-computed csv file where the rows are patients/samples and columns are protein names
      the values are string annotation for oncogenicity of SNVs

input: a csv file with the rows as samples/patients and columns as protein names
        the values should be as cnv:***_snv:***p.E23G*_fusion:******

USAGES: python scriptname.py input.csv
"""

AUTHENTICATION_KEY="" # this token should be taken from onkoKB API
# the key for the authentication can be retrieved after registering in the OncoKB portal 

def get_data_from_the_ONKOKB_api(gene_name,alteration,AUTHENTICATION_KEY=AUTHENTICATION_KEY,data_type='variantSummary'):
    """
    this function retrives the oncoKB annoation from oncoKB API as string 
    :gene_name: is the symbol gene name
    :alteration is the SNV mutation like R989H
    :data_type is the field you like to extract
    """
    alteration = str(alteration).replace('p.','')
    url = f'https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange?hugoSymbol={gene_name}&alteration={alteration}'
    headers = {'accept': 'application/json', "Authorization": AUTHENTICATION_KEY}
    final = requests.get(url, headers=headers).json()[data_type]
    logger.info(f'{final} for ## {alteration} ## on {gene_name}')
    return final



def load_genomics_table(genomics_path="") -> pd.DataFrame:
    genomics_df = pd.read_csv(Path(genomics_path))
    genomics_df = genomics_df.loc[:, ~genomics_df.columns.str.contains("^Unnamed")]
    return genomics_df


def get_all_snv_ONKOKB_per_snvname(X,protein_name,pattern_alteration = r'p.[A-Z][0-9]+[A-Z]'):
    """
    This function is the wrapper around the function for extracting the SNVs 
    :X: should be like cnv:CNN_snv:C_T_exonic_20_p.P266S_fusion:n.d	 if it is multi they shold be seprated by ;
    """
    try:
        
        all = X.split(';')
        all_alteration = [x.split('snv:')[-1].split('fusion:')[0] for x in all]

        finallist = []
        all_matches = [re.finditer(pattern_alteration,x)  for x in all_alteration]
        for match in all_matches:
            for x in match:
                finallist.append(x.group()) 
        finallist = list(set(finallist))

        if len(finallist) > 0:
            return (';').join([get_data_from_the_ONKOKB_api(protein_name,x) for x in finallist ])
        else:
            return ''
    except Exception as err:
        #print(f'something wrong{err}')
        return ''


if __name__ == "__main__":
    import sys
    genomocics_df = load_genomics_table(genomics_path=sys.argv[1])
    df = genomocics_df.copy()
    df.set_index('Sample name',inplace=True)
    for j in range(len(df.columns)):
        gene_name = df.columns[j]
        for i in range(len(df)):
            df.iloc[i,j] = get_all_snv_ONKOKB_per_snvname(df.iloc[i,j],gene_name)
    df.to_csv("onkokb2portal.csv")
