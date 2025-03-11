# TOPAS-pipeline

Automated (phospho)-proteomics processing pipeline for large patient cohorts providing cohort as well as patient-specific insights. 

After analysis can be explored on the web-based portal: https://github.com/kusterlab/TOPAS-portal.git

A public instance of the portal can be found here: https://topas-portal.kusterlab.org/master_mtb_portal/ 

Pipeline can run on both TMT and label-free data. 


## Installation and running


### Requirements

The TOPAS pipeline requires Python ">=3.9, <=3.11". 
OS: Ubuntu >=18

The pipeline needs configurations input from a JSON file format. Examples can be found in `config.json` (full version) and `config_minimal_example.json` (only required configs)

Minimum cores and memory?

For a pipeline run on ~1,000 samples runtime is ... and requires xxx gb memory

For a pipeline run on xxx samples runtime is ... and requires xxx gb memory


### Running pipeline with Docker (recommended)

TOPAS-pipeline can be installed from this repository:

```
git clone https://github.com/kusterlab/TOPAS-pipeline.git
pip install .
```

For running with default configurations:
adjust the file paths in `config_minimal_example.json` and run:

```
make docker_all 
```

If you want to run it not using the local config file or if you want to set a different memory limit (default: 100gb), this can be done as following:

#### #TODO --> what default values do we want for memory and cpu?
```
CONFIG_FILE=/path/to/config.json MEMORY_LIMIT=300gb CPU_LIMIT=16 make all
```

### How to run pipeline with conda and poetry

Create environment and install required packages from poetry.lock file:

```
conda create --name topas-pipeline python=3.9.12
conda activate topas-pipeline
poetry install
```

To run the pipeline, adjust the file paths in `config_minimal_example.json` and run:
```
make all
```

Alternatively, to run pipeline module by module:
```
# run simsi only
python -m topas_pipeline.simsi -c config.json

# run whole pipeline following simsi
python -m topas_pipeline.main -c config.json 

# run only clinical annotation
python -m topas_pipeline.clinical_process -c config.json
```


### Setup example and Integration tests

An example of the used project folder setup is presented in the `/example` folder. A `test_config.json` is found which is used for integration tests found in `/tests/integration_tests`. After cloning and if problems arise, it is recommended to first test pipeline integration by running the tests.
```
pytest ./tests/integration_tests/test_simsi.py
pytest ./tests/integration_tests/test_picked_group.py
pytest ./tests/integration_tests/test_clinical_tools.py

```

#### # TODO: --> try to add PSP files and see if it would work


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



## Configurations explained

Input parameters for running the pipeline can be adjusted in config.json file.

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



## Output files explained


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
