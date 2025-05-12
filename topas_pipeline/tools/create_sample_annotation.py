"""

Description:
Script to create input sample annotation file for pipeline from metadata file. 
The location of MQ searches are given by Cohort column info and is either hard coded to be Sarcoma (for patients) or Workflow_Test (for models) based on Batch name.
Reference channels are added for each batch as TMT channel 10 and 11.

Usage:
    
    python tools/create_sample_annotation.py Metadata_file.xlsx

"""

import sys
import re
import datetime
import pandas as pd



def main(argv):
    metadata_file = argv[0]
    metadata = pd.read_excel(metadata_file, engine='openpyxl')
    # metadata = metadata[metadata['Paper_extv2'] == 'yes']
    metadata = metadata[metadata['QC'] != 'exclude']
    # metadata = metadata[['Sample name', 'Program', 'code_oncotree', 'tissue_topology', 'Batch_No', 'TMT_Channel', 'QC']]
    metadata = metadata.rename(columns={'Batch_No': 'Batch Name', 'TMT_Channel': 'TMT Channel', 'code_oncotree': 'Entity', 'tissue_topology': 'Histologic subtype'})

    new_ref_rows = []

    # Iterate through each unique Batch Name
    for batch in metadata['Batch Name'].unique():
        for channel in [10, 11]:
            # Create a new row as a dictionary
            new_row = {
                'Batch Name': batch,
                'Sample name': f"ref_channel_{channel}_batch{batch}",
                'TMT Channel': channel,
                'QC': 'passed',
                'is_reference': True
            }
            # Append the new row to the list
            new_ref_rows.append(new_row)

    # Convert the list of new rows into a DataFrame
    new_rows_df = pd.DataFrame(new_ref_rows)

    # Concatenate the new rows DataFrame with the original DataFrame
    metadata = pd.concat([metadata, new_rows_df], ignore_index=True)

    metadata['Cohort'] = 'Sarcoma'

    def update_cohort(row):
        batch_name = str(row['Batch Name'])  # Convert the value to string
        # Check if the batch name contains both letters and numbers
        if re.search(r'[A-Za-z]', batch_name) and re.search(r'\d', batch_name):
            return "Workflow_Test"
        return row['Cohort']

    metadata['Cohort'] = metadata.apply(update_cohort, axis=1)
    today = datetime.datetime.today().date()
    metadata.to_csv(f'sample_annotation_{today}.csv', index=False)


if __name__ == '__main__':
    main(sys.argv[1:])
