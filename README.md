# TOPAS-pipeline

Automated (phospho)proteomics processing pipeline for large patient cohorts providing cohort as well as patient-specific insights. 

The results of the pipeline can be explored on the web-based TOPAS portal: https://github.com/kusterlab/TOPAS-portal.git
A public instance of the portal can be found here: https://topas-portal.kusterlab.org/

## Supported inputs

- MaxQuant (TMT, LFQ)
- SIMSI-Transfer (TMT)

## Runtime overview

Runtimes exclude processing time of SIMSI-Transfer.

| Dataset | #channels | #samples | #cores | runtime (h) | max memory (GB) |
| --- | --- | --- | --- | --- | --- |
| CPTAC UCEC | 170 | 153 | 8 | 1.2 | 30 |
| CPTAC BRCA | 170 | 153 | 8 | 1.4 | 30 |
| CPTAC LUAD | 250 | 225 | 8 | 2 | 40 |
| MTB cohort | 2068 | 1284 | 8 | 13 | 150 |

## Configuration

The pipeline needs configurations input from a JSON file format. 
Examples can be found in `config.json` (full version) and `config_minimal.json` (only required configs).
Both relative and absolute paths are allowed.

Input parameters:

| Parameter | Required | Description | Example | Default |
| --- | --- | --- | --- | --- |
| results_folder | **yes** | Path to the folder where results will be written. | `"results/example_run"` | |
| sample_annotation | **yes** | Path to the sample annotation file (CSV). | `"example/annotation.csv"` | |
| metadata_annotation | **yes** | Path to the metadata annotation file (Excel). | `"example/METADATA_UCEC.xlsx"` | |
| raw_file_folders | **yes** | List of raw file folders for proteomics and phosphoproteomics data. | `["example/raw_fp", "example/raw_pp"]` | |
| data_types | | List of data types to process: "fp" for proteome and "pp" for phosphoproteome. | `["fp", "pp"]` | `["fp", "pp"]` |
| slack_webhook_url | | URL for the Slack webhook. | `""` | `""` |
| **simsi** |  |  |  |  |
| run_simsi | | Boolean indicating whether to run SIMSI analysis. | `true` | `true` |
| simsi_folder | **yes** | Path to the folder for writing SIMSI-Transfer results. | `"results/SIMSI"` | N/A |
| tmt_ms_level | | MS level for TMT quantification. | `"ms2"` | `"ms2"` |
| stringencies | | Stringency value for MaRaCluster. | `10` | `10` |
| tmt_requantify | | Boolean indicating whether to requantify TMT data. | `false` | `false` |
| maximum_pep | | Maximum posterior error probability in percent for peptide ID propagation. | `1` | `1` |
| num_threads | | Number of threads to use for SIMSI-Transfer. | `8` | `8` |
| **preprocessing** |  |  |  |  |
| raw_data_location | **yes** | Path to the folder containing MaxQuant search result folders. | `"example/CPTAC_searches"` | N/A |
| fasta_file | **yes** | Path to the FASTA file for protein sequences. | `"example/uniprot_human.fasta"` | N/A |
| picked_fdr | | False discovery rate threshold for protein groups. | `0.01` | `0.01` |
| fdr_num_threads | | Number of threads to use in MaxLFQ computation. | `8` | `8` |
| imputation | | Perform data imputation within batch on phosphoproteome level. | `true` | `true` |
| debug | | Run in debug mode. | `false` | `false` |
| run_lfq | | Input is from LFQ experiments. | `false` | `false` |
| normalize_to_reference | | Normalize channel intensities to the reference channel. | `false` | `false` |
| **clinic_proc** |  |  |  |  |
| pspFastaFile | **yes** | Path to the PSP FASTA file. | `"PSP_annotations/Phosphosite_seq.fasta"` | N/A |
| pspKinaseSubstrateFile | **yes** | Path to the PSP kinase-substrate dataset. | `"PSP_annotations/Kinase_Substrate_Dataset"` | N/A |
| pspAnnotationFile | **yes** | Path to the PSP phosphorylation site dataset. | `"PSP_annotations/Phosphorylation_site_dataset"` | N/A |
| pspRegulatoryFile | **yes** | Path to the PSP regulatory sites file. | `"PSP_annotations/Regulatory_sites"` | N/A |
| prot_baskets | **yes** | Path to the annotation file for TOPAS scores and proteins of interest. | `"TOPASscores_POI_AS_250307.xlsx"` | N/A |
| extra_kinase_annot | | Path to the annotation file with custom kinase-substrate relations. | `""` | `""` |
| **report** |  |  |  |  |
| samples_for_report | | Which samples to include in the report. | `"all"` | `"all"` |
| **portal** |  |  |  |  |
| update | | Automatically update the TOPAS portal once the run has finished | `false` | `false` |
| cohort | | Specifies the cohort name that should be updated in the TOPAS portal. | `""` | `""` |
| url | | URL of the TOPAS portal. | `""` | `""` |
| config | | Configuration file for the TOPAS portal. | `""` | `""` |


### Install webhook for slack (optional)

If you want the pipeline to post update messages (finished runs, error messages) to 
your slack channel, follow these steps:

1. Create a new slack app here: https://api.slack.com/apps?new_app=1, use the `from scratch` option.
2. Select your slack workspace and pick an appropriate name for the app, e.g. `topas-pipeline`.
3. Navigate to `Incoming webhooks` in the left menu.
4. Set `Activate Incoming webhooks` to `On` if this was not already the case.
5. Click on `Add New Webhook to Workspace` at the bottom of the page.
6. Select the channel you want to post messages in.
7. Copy the generated `Webhook URL` to your config file as the `slack_webhook_url` property.

Source: https://api.slack.com/messaging/webhooks


## Running the pipeline

### With Docker (recommended)

Requirements:

- git
- docker
- make


1. Clone this repository
    ```
    git clone https://github.com/kusterlab/TOPAS-pipeline.git
    ```
2. Build the docker image
    ```
    make build
    ```
3. Create a config file named `config_patients.json` in the repository with your configurations (see section `Configuration`) and run:
    ```
    make docker_all 
    ```
   You can use a custom config file (works only with relative paths) and adjust the memory and cores (default: 300GB, 8 cores):
    ```
    CONFIG_FILE=./path/to/config.json MEMORY_LIMIT=300gb CPU_LIMIT=16 make docker_all
    ```

### With conda and poetry

Requirements:

- git
- Python ">=3.9, <=3.11"
- poetry
- make
- conda

1. Create environment and install required packages from poetry.lock file:
    ```
    conda create --name topas-pipeline python=3.9.12
    conda activate topas-pipeline
    ```
2. Clone this repository
    ```
    git clone https://github.com/kusterlab/TOPAS-pipeline.git
    ```
3. Install dependencies and start a poetry shell
    ```
    poetry install
    poetry shell
    ```
4. Adjust the file paths in `config_minimal.json` and run:
    ```
    make all
    ```


Note that it is also possible to run individual pipeline modules, e.g.:
```
# run simsi
python -m topas_pipeline.simsi -c config.json

# run whole pipeline following simsi
python -m topas_pipeline.main -c config.json 

# run clinical annotation
python -m topas_pipeline.clinical_annotation -c config.json
```


## Example

An example of the project folder setup and configuration file can be found in the `/example` folder.
Check the ReadMe in the `/example` folder for details.

## Integration tests

If problems arise with running the pipeline, check if the integration tests pass.
```
pytest ./tests/integration_tests/test_simsi.py
pytest ./tests/integration_tests/test_picked_group.py
pytest ./tests/integration_tests/test_clinical_tools.py
```


## Pipeline result files

The pipeline creates a folder with multiple output files:

| Output file | Description | Used on portal |
| --- | --- | ---  |
| config.json | Copy of the configuration file in JSON format used for this pipeline run as described above. |  |
| sample_annot.tsv | Saved copy of current version of sample annotation/metadata given as input |   |
| sample_annot_filtered.tsv | Subset of sample annotation/metadata after filtering out QC failed samples |   |
| meta_input_file_{data_type}.tsv | Location per batch of search folder input, raw files and TMT correction factor file |   |
| {data_type}_qc_numbers.csv | Per sample count of peptides, median intensities and summed intensities |   |
| {data_type}_qc_batch_wise.csv | Per batch median and summed intensities |   |
| {data_type}_in_batch_correction_factors.csv | Per sample correction factors for in-batch median centering |   |
| {data_type}_ms1_correction_factors.csv | Per batch correction factors for MS1 median centering |   |
| evidence.txt | Combined output from SIMSI-Transfer (same format as MQ evidence file) |   |
| pickedGeneGroups.txt | Output from Picked Protein Group FDR (using gene-level) |   |
| pickedGeneGroups_with_quant.txt | Groups from Picked Protein Group FDR (using gene-level) containing quant using MaxLFQ algorithm  |   |
| preprocessed_{data_type}.csv | Data matrix with patients as columns and normalized abundances of proteins or phosphopeptides as rows |  |
| preprocessed_{data_type}_with_ref.csv | Same as preprocessed_{data_type}.csv but also containing the QC channels |   |
| annot_{data_type}.csv | Same as preprocessed_{data_type}.csv but with clinical annotations |   |
| annot_{data_type}_with_ref.csv | Same as preprocessed_{data_type}_with_ref.csv but with clinical annotations |   |
| {data_type}_measures_rank.tsv | Data matrix with patients as columns and in-cohort rank of proteins or phosphopeptides as rows |   |
| {data_type}_measures_fc.tsv | Same as {data_type}_measures_rank.tsv but with fold changes  |   |
| {data_type}_measures_z.tsv | Same as {data_type}_measures_rank.tsv but with z-scores |   |
| {data_type}_measures_p.tsv | Same as {data_type}_measures_rank.tsv but with p-values derived from the z-scores |   |
| basket_scores_4th_gen_zscored.tsv | Data matrix with patients as columns and TOPAS scores as rows |   |
| subbasket_scores_{topas_rtk}.tsv | Data matrix with patients as columns and TOPAS subscores as rows for each TOPAS RTK |   |
| kinase_results/kinase_scores.tsv | Data matrix with patients as columns and TOPAS substrate phosphorylation scores as rows |   |
| protein_results/protein_scores.tsv | Data matrix with patients as columns and TOPAS protein phosphorylation scores as rows |   |
| Reports/{patient_id}.xlsx | Patient-specific reports |   |
| Pipeline_log.txt | Log messages printed by the pipeline |   |
|  |  |   |
