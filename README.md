# TOPAS-pipeline

Automated (phospho)proteomics processing pipeline for large patient cohorts providing cohort as well as patient-specific insights. 

The results of the pipeline can be explored on the web-based TOPAS portal: https://github.com/kusterlab/TOPAS-portal.git
A public instance of the portal can be found here: https://topas-portal.kusterlab.org/

## Supported inputs

- MaxQuant (TMT, LFQ)
- SIMSI-Transfer (TMT)

For a pipeline run on ~1,000 samples runtime is ... and requires xxx gb memory

For a pipeline run on xxx samples runtime is ... and requires xxx gb memory

## Configuration

The pipeline needs configurations input from a JSON file format. 
Examples can be found in `config.json` (full version) and `config_minimal.json` (only required configs).
Both relative and absolute paths are allowed.

Input parameters:

| Parameter | Description | Required |
| --- | --- | ---  |
| results_folder |  | yes |
| sample_annotation |  | yes |
| metadata_annotation |  | yes |
| simsi_folder |  | yes |
| tmt_ms_level |  |  |
| stringencies |  |  |
| tmt_requantify |  |  |
| maximum_pep |  |  |
| num_threads |  |  |
| raw_data_location |  | yes |
| picked_fdr |  |  |
| fdr_num_threads |  |  |
| imputation |  |  |
| run_simsi |  |  |
| debug |  |  |
| run_lfq |  |  |
| normalize_to_reference |  |  |
| fasta_file |  |  |
| pspFastaFile |  |  |
| pspKinaseSubstrateFile |  |  |
| pspAnnotationFile |  |  |
| pspRegulatoryFile |  |  |
| prot_baskets |  |  |
| samples_for_report |  |  |
| data_types |  |  |
| update |  |  |
| cohort |  |  |
| url |  |  |
| config |  |  |
| slack_webhook_url |  |  |
|  |  |   |

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
3. Adjust the file paths in `config_minimal.json` and run:
    ```
    make docker_all 
    ```
   You can use a custom config file and adjust the memory and cores (default: 300GB, 8 cores):
    ```
    CONFIG_FILE=/path/to/config.json MEMORY_LIMIT=300gb CPU_LIMIT=16 make docker_all
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
python -m topas_pipeline.clinical_process -c config.json
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
| config.json |  |  |
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
| preprocessed_{data_type}.csv | i want to delete |   |
| preprocessed_{data_type}_with_ref.csv | i want to make it standard and remove _with_ref |   |
| annot_{data_type}.csv |  |   |
| annot_{data_type}_with_ref.csv | same story as above plus remove intensities? |   |
| topas_annot_dict_{data_type}.json  |  |   |
| poi_annot_dict.json |  |   |
| {data_type}_measures_rank.tsv |  |   |
| {data_type}_measures_fc.tsv |  |   |
| {data_type}_measures_z.tsv |  |   |
| {data_type}_measures_p.tsv |  |   |
| subbasket_scores_{subbasket}.tsv |  |   |
| basket_scores_4th_gen_zscored.tsv |  |   |
| Pipeline_log.txt |  |   |
|  |  |   |
