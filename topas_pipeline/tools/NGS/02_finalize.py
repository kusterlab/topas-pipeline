import pandas as pd
from pathlib import Path



def _aggregator(x):
    try: 
        return ';'.join(sorted(set(x)))
    except:
        return None


def main(
    path_to_folder =  '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/CJ/MASTER_genomics/2025.07.29/MASTER_genome',
    final_file = 'master_patients2portal.csv',
    mapping_file = 'MASTER_mapping_table_AS.xlsx'):
    """
    Process and aggregate genomic data (CNV, SNV, fusion) for patient-protein combinations and 
    generate a pivoted NGS annotation matrix for downstream portal use.

    This function:
    - Reads input CSV files for SNVs, CNVs, and fusions
    - Constructs per-patient, per-protein annotations for each modality
    - Aggregates annotations across modalities into a single 'NGS' column
    - Pivots the data into a wide-format table (patients x proteins)
    - Maps sample names using an external Excel mapping file
    - Saves the final pivoted data to a CSV file

    Parameters
    ----------
    path_to_folder : str, optional
        Path to the folder containing input and output files. Defaults to a pre-defined project path.

    final_file : str, optional
        File name for the final pivoted output CSV. Defaults to 'master_patients2portal.csv'.

    mapping_file : str, optional
        Excel file containing the mapping between internal sample IDs and final patient labels.
        Must contain columns 'NGS_mappingID' and 'Sample name'.

    Returns
    -------
    None
        Writes the processed pivot table to `final_file` in the `path_to_folder`.

    Input Files Required (inside `path_to_folder`)
    ---------------------------------------------
    - 'fusion_df.csv' : Fusion events per patient and gene
    - 'cnv_df.csv' : Copy number variants per patient and gene
    - 'snv_df.csv' : Single nucleotide variants per patient and gene
    - mapping_file : Excel sheet with patient ID mapping

    Output Files (written into `path_to_folder`)
    --------------------------------------------
    - 'fusion_df_processed.csv' : Processed fusion annotations
    - 'cnv_df_processed.csv' : Processed CNV annotations
    - 'snv_df_processed.csv' : Processed SNV annotations
    - `final_file` : Final pivoted annotation matrix for upload

    Notes
    -----
    - Assumes existence of a `_aggregator` function to deduplicate and join multiple annotations per group.
    - CNV and SNV annotations are normalized with flags like 'LOH', 'functional', 'clinical'.
    - Fusion confidence types are annotated and grouped similarly.

    Example
    -------
    >>> main(
    ...     path_to_folder='/path/to/data',
    ...     final_file='upload_ready.csv',
    ...     mapping_file='sample_mapping.xlsx'
    ... )

    """

    fusion_df = pd.read_csv(Path(path_to_folder) / Path('fusion_df.csv'))
    cnv_df = pd.read_csv(Path(path_to_folder) / Path('cnv_df.csv'))
    snv_df = pd.read_csv(Path(path_to_folder) / Path('snv_df.csv'))



    print('going with fusions')
    fusion_df['fustion_annotation'] = fusion_df['Type'].astype(str) + '(' + fusion_df['Confidence'].astype(str) + ')'
    fusion_df['fustion_annotation'] = fusion_df.groupby(['patient', 'protein_name'])['fustion_annotation'].transform(_aggregator)
    fusion_df = fusion_df[['fustion_annotation','patient','protein_name']].drop_duplicates(subset=['patient','protein_name'])
    fusion_df.to_csv(Path(path_to_folder) / Path('fusion_df_processed.csv'),index=None)

    print('going with CNVs')
    cnv_df['LOH'][cnv_df.LOH] = 'LOH'
    cnv_df['LOH'][cnv_df.LOH != 'LOH'] = ''
    cnv_df['cnv_annotation'] = cnv_df['Type'].astype(str) + '(' + cnv_df['LOH'].astype(str) + ')'
    deduped_df = cnv_df[['patient', 'protein_name', 'cnv_annotation']].drop_duplicates()
    cnv_df_agg = (
        deduped_df
        .groupby(['patient', 'protein_name'], as_index=False)
        .agg({'cnv_annotation': lambda x: _aggregator(x)})
    )
    cnv_df_agg = cnv_df_agg[['patient', 'protein_name', 'cnv_annotation']].drop_duplicates()
    cnv_df_agg.to_csv(Path(path_to_folder) / Path('cnv_df_processed.csv'),index=None)

    print('processing SNVs')
    snv_df['IsFun'][snv_df.IsFun] = 'functional'
    snv_df['IsFun'][snv_df.IsFun != 'functional'] = ''
    snv_df['IsFunOnTarget'][snv_df.IsFunOnTarget] = 'Targetfunctional'
    snv_df['IsFunOnTarget'][snv_df.IsFunOnTarget != 'Targetfunctional'] = ''
    snv_df['PassedClinicalWorkflow'][snv_df.PassedClinicalWorkflow] = 'clinical'
    snv_df['PassedClinicalWorkflow'][snv_df.PassedClinicalWorkflow != 'clinical'] = ''
    snv_df['snv_annotation'] = snv_df['CanonicalTranscript'].astype(str) + '(' + snv_df['IsFun'].astype(str) + '_' + snv_df['IsFunOnTarget'].astype(str)  + '_' + snv_df['PassedClinicalWorkflow'].astype(str)  + ')'
    snv_df['snv_annotation'] = snv_df.groupby(['patient', 'protein_name'])['snv_annotation'].transform(_aggregator)
    snv_df = snv_df[['snv_annotation','patient','protein_name']].drop_duplicates(subset=['patient','protein_name'])
    snv_df.to_csv(Path(path_to_folder) / Path('snv_df_processed.csv'),index=None)

    snv_df = pd.read_csv(Path(path_to_folder) / Path('snv_df_processed.csv'))
    cnv_df = pd.read_csv(Path(path_to_folder) / Path('cnv_df_processed.csv'))
    fusion_df = pd.read_csv(Path(path_to_folder) / Path('fusion_df_processed.csv'))

    cnv_snv = cnv_df.merge(snv_df,on=['patient','protein_name'],how='outer')
    final_df = cnv_snv.merge(fusion_df,on=['patient','protein_name'],how='outer')
    final_df['cnv_annotation'] = final_df['cnv_annotation'].fillna('n.d') 
    final_df['cnv_annotation'] = 'cnv:' + final_df['cnv_annotation']

    final_df['snv_annotation'] = final_df['snv_annotation'].fillna('n.d') 
    final_df['snv_annotation'] = 'snv:' + final_df['snv_annotation']

    final_df['fustion_annotation'] = final_df['fustion_annotation'].fillna('n.d') 
    final_df['fustion_annotation'] = 'fusion:' + final_df['fustion_annotation']
    final_df['NGS'] = final_df['cnv_annotation'] + '_' + final_df['snv_annotation'] + '_' + final_df['fustion_annotation']
    df_pivot = final_df[['protein_name','patient','NGS']].pivot(index='patient', columns='protein_name', values='NGS')

    mapping_df = pd.read_excel(Path(path_to_folder) / Path(mapping_file))
    mapping_df = mapping_df.dropna(subset='NGS_mappingID')
    mapping_dic = mapping_df[['Sample name','NGS_mappingID']].set_index('NGS_mappingID').to_dict()['Sample name']
    df_pivot = df_pivot[df_pivot.index.isin(mapping_df.NGS_mappingID.unique().tolist())]
    df_pivot['Sample name'] = df_pivot.index.map(mapping_dic)
    df_pivot = df_pivot.set_index('Sample name')
    df_pivot.to_csv(Path(path_to_folder) / Path(final_file))




if __name__ == '__main__':

    main(
    path_to_folder =  '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/CJ/MASTER_genomics/2025.07.29/MASTER_genome',
    final_file = 'master_patients2portal.csv',
    mapping_file = 'MASTER_mapping_table_AS.xlsx'
    )

    main(
    path_to_folder =  '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/CJ/MASTER_genomics/2025.07.29_CHDM',
    final_file = 'chdm_patients2portal.csv',
    mapping_file = 'Chordoma_mapping_table_AS.xlsx')
