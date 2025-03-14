# Example

This example applies the pipeline to the first two batches of the CPTAC UCEC cohort.


## 1. Download input files

To run this example, some files need to be downloaded:

1. Unzip the fasta file in `uniprot_proteome.zip` inside the current directory.
2. Download the raw files for Batch 1 and 2 of the UCEC cohort from the CPTAC data portal and put them in `raw_fp/<batch_name>` and `raw_pp/<batch_name>` respectively, e.g. `raw_fp/01CPTAC_UCEC_Proteome_PNNL_20170922_raw` should contain all raw files for the proteome of Batch 1.
    - Proteome: https://proteomic.datacommons.cancer.gov/pdc/study/PDC000125
    - Phosphoproteome: https://proteomic.datacommons.cancer.gov/pdc/study/PDC000126
3. Download the MaxQuant search results for Batch 1 and 2 from our PRIDE repository (https://www.ebi.ac.uk/pride/archive/projects/PXD061316) and put them in `CPTAC_searches/UCEC/<batch_name>`, e.g. `CPTAC_searches/UCEC/Batch01_FP_CPTAC_UCEC`. Remove the prefixes from the MaxQuant output files, e.g. `UCEC_Batch1_Proteome_evidence.txt` => `evidence.txt`.
4. Download the `Kinase_Substrate_Dataset`, `Phosphorylation_site_dataset`, `Phosphosite_seq.fasta` and `Regulatory_sites` files from PhosphositePlus: https://www.phosphosite.org/staticDownloads (registration needed).

After this, your directory structure should look like this:

```
example/
|-- CPTAC_searches
|   `-- UCEC
|       |-- Batch01_FP_CPTAC_UCEC
|       |   |-- combined
|       |   |   `-- txt
|       |   |       |-- allPeptides.txt
|       |   |       |-- evidence.txt
|       |   |       |-- msms.txt
|       |   |       |-- msmsScans.txt
|       |   |       |-- summary.txt
|       |-- Batch01_PP_CPTAC_UCEC
|       |   |-- combined
|       |   |   `-- txt
|       |   |       |-- allPeptides.txt
|       |   |       |-- evidence.txt
|       |   |       |-- msms.txt
|       |   |       |-- msmsScans.txt
|       |   |       |-- summary.txt
|       |-- Batch02_FP_CPTAC_UCEC
|       |   |-- combined
|       |   |   `-- txt
|       |   |       |-- allPeptides.txt
|       |   |       |-- evidence.txt
|       |   |       |-- msms.txt
|       |   |       |-- msmsScans.txt
|       |   |       |-- summary.txt
|       `-- Batch02_PP_CPTAC_UCEC
|           |-- combined
|           |   `-- txt
|           |       |-- allPeptides.txt
|           |       |-- evidence.txt
|           |       |-- msms.txt
|           |       |-- msmsScans.txt
|           |       |-- summary.txt
|-- METADATA_UCEC.xlsx
|-- PSP_annotations
|   |-- Kinase_Substrate_Dataset
|   |-- Phosphorylation_site_dataset
|   |-- Phosphosite_seq.fasta
|   `-- Regulatory_sites
|-- ReadMe.md
|-- TOPASscores_POI_AS_250307.xlsx
|-- cptac_annotation_subset.csv
|-- raw_fp
|   |-- 01CPTAC_UCEC_Proteome_PNNL_20170922_raw
|   |   |-- 01CPTAC_UCEC_W_PNNL_20170922_B1S1_f01.raw
|   |   |-- 01CPTAC_UCEC_W_PNNL_20170922_B1S1_f02.raw
|   |   |-- 01CPTAC_UCEC_W_PNNL_20170922_B1S1_f03.raw
|   |   |-- 01CPTAC_UCEC_W_PNNL_20170922_B1S1_f04.raw
|   |   |-- 01CPTAC_UCEC_W_PNNL_20170922_B1S1_f05.raw
|   |   |-- 01CPTAC_UCEC_W_PNNL_20170922_B1S1_f06.raw
|   |   |-- 01CPTAC_UCEC_W_PNNL_20170922_B1S1_f07.raw
|   |   |-- 01CPTAC_UCEC_W_PNNL_20170922_B1S1_f08.raw
|   |   |-- 01CPTAC_UCEC_W_PNNL_20170922_B1S1_f09.raw
|   |   |-- 01CPTAC_UCEC_W_PNNL_20170922_B1S1_f10.raw
|   |   |-- 01CPTAC_UCEC_W_PNNL_20170922_B1S1_f11.raw
|   |   |-- 01CPTAC_UCEC_W_PNNL_20170922_B1S1_f12.raw
|   |   |-- 01CPTAC_UCEC_W_PNNL_20170922_B1S1_f13.raw
|   |   |-- 01CPTAC_UCEC_W_PNNL_20170922_B1S1_f14.raw
|   |   |-- 01CPTAC_UCEC_W_PNNL_20170922_B1S1_f15.raw
|   |   |-- 01CPTAC_UCEC_W_PNNL_20170922_B1S1_f16.raw
|   |   |-- 01CPTAC_UCEC_W_PNNL_20170922_B1S1_f17.raw
|   |   |-- 01CPTAC_UCEC_W_PNNL_20170922_B1S1_f18.raw
|   |   |-- 01CPTAC_UCEC_W_PNNL_20170922_B1S1_f19.raw
|   |   |-- 01CPTAC_UCEC_W_PNNL_20170922_B1S1_f20.raw
|   |   |-- 01CPTAC_UCEC_W_PNNL_20170922_B1S1_f21.raw
|   |   |-- 01CPTAC_UCEC_W_PNNL_20170922_B1S1_f22.raw
|   |   |-- 01CPTAC_UCEC_W_PNNL_20170922_B1S1_f23.raw
|   |   |-- 01CPTAC_UCEC_W_PNNL_20170922_B1S1_f24.raw
|   `-- 02CPTAC_UCEC_Proteome_PNNL_20170922_raw
|       |-- 02CPTAC_UCEC_W_PNNL_20170922_B1S2_f01.raw
|       |-- 02CPTAC_UCEC_W_PNNL_20170922_B1S2_f02.raw
|       |-- 02CPTAC_UCEC_W_PNNL_20170922_B1S2_f03.raw
|       |-- 02CPTAC_UCEC_W_PNNL_20170922_B1S2_f04.raw
|       |-- 02CPTAC_UCEC_W_PNNL_20170922_B1S2_f05.raw
|       |-- 02CPTAC_UCEC_W_PNNL_20170922_B1S2_f06.raw
|       |-- 02CPTAC_UCEC_W_PNNL_20170922_B1S2_f07.raw
|       |-- 02CPTAC_UCEC_W_PNNL_20170922_B1S2_f08.raw
|       |-- 02CPTAC_UCEC_W_PNNL_20170922_B1S2_f09.raw
|       |-- 02CPTAC_UCEC_W_PNNL_20170922_B1S2_f10.raw
|       |-- 02CPTAC_UCEC_W_PNNL_20170922_B1S2_f11.raw
|       |-- 02CPTAC_UCEC_W_PNNL_20170922_B1S2_f12.raw
|       |-- 02CPTAC_UCEC_W_PNNL_20170922_B1S2_f13.raw
|       |-- 02CPTAC_UCEC_W_PNNL_20170922_B1S2_f14.raw
|       |-- 02CPTAC_UCEC_W_PNNL_20170922_B1S2_f15.raw
|       |-- 02CPTAC_UCEC_W_PNNL_20170922_B1S2_f16.raw
|       |-- 02CPTAC_UCEC_W_PNNL_20170922_B1S2_f17.raw
|       |-- 02CPTAC_UCEC_W_PNNL_20170922_B1S2_f18.raw
|       |-- 02CPTAC_UCEC_W_PNNL_20170922_B1S2_f19.raw
|       |-- 02CPTAC_UCEC_W_PNNL_20170922_B1S2_f20.raw
|       |-- 02CPTAC_UCEC_W_PNNL_20170922_B1S2_f21.raw
|       |-- 02CPTAC_UCEC_W_PNNL_20170922_B1S2_f22.raw
|       |-- 02CPTAC_UCEC_W_PNNL_20170922_B1S2_f23.raw
|       `-- 02CPTAC_UCEC_W_PNNL_20170922_B1S2_f24.raw
|-- raw_pp
|   |-- 01CPTAC_UCEC_Phosphoproteome_PNNL_20170922_raw
|   |   |-- 01CPTAC_UCEC_P_PNNL_20170922_B1S1_f01.raw
|   |   |-- 01CPTAC_UCEC_P_PNNL_20170922_B1S1_f02.raw
|   |   |-- 01CPTAC_UCEC_P_PNNL_20170922_B1S1_f03.raw
|   |   |-- 01CPTAC_UCEC_P_PNNL_20170922_B1S1_f04.raw
|   |   |-- 01CPTAC_UCEC_P_PNNL_20170922_B1S1_f05.raw
|   |   |-- 01CPTAC_UCEC_P_PNNL_20170922_B1S1_f06.raw
|   |   |-- 01CPTAC_UCEC_P_PNNL_20170922_B1S1_f07.raw
|   |   |-- 01CPTAC_UCEC_P_PNNL_20170922_B1S1_f08.raw
|   |   |-- 01CPTAC_UCEC_P_PNNL_20170922_B1S1_f09.raw
|   |   |-- 01CPTAC_UCEC_P_PNNL_20170922_B1S1_f10.raw
|   |   |-- 01CPTAC_UCEC_P_PNNL_20170922_B1S1_f11.raw
|   |   `-- 01CPTAC_UCEC_P_PNNL_20170922_B1S1_f12.raw
|   `-- 02CPTAC_UCEC_Phosphoproteome_PNNL_20170922_raw
|       |-- 02CPTAC_UCEC_P_PNNL_20170922_B1S2_f01.raw
|       |-- 02CPTAC_UCEC_P_PNNL_20170922_B1S2_f02.raw
|       |-- 02CPTAC_UCEC_P_PNNL_20170922_B1S2_f03.raw
|       |-- 02CPTAC_UCEC_P_PNNL_20170922_B1S2_f04.raw
|       |-- 02CPTAC_UCEC_P_PNNL_20170922_B1S2_f05.raw
|       |-- 02CPTAC_UCEC_P_PNNL_20170922_B1S2_f06.raw
|       |-- 02CPTAC_UCEC_P_PNNL_20170922_B1S2_f07.raw
|       |-- 02CPTAC_UCEC_P_PNNL_20170922_B1S2_f08.raw
|       |-- 02CPTAC_UCEC_P_PNNL_20170922_B1S2_f09.raw
|       |-- 02CPTAC_UCEC_P_PNNL_20170922_B1S2_f10.raw
|       |-- 02CPTAC_UCEC_P_PNNL_20170922_B1S2_f11.raw
|       `-- 02CPTAC_UCEC_P_PNNL_20170922_B1S2_f12.raw
`-- uniprot_proteome_up000005640_03112020_cdkn2a_isoforms.fasta
```

## 2. Run the pipeline

Make sure the pipeline is installed as outlined in the main ReadMe.

### Docker

Go to the main directory of this repository and run:
```
CONFIG_FILE=config.json make docker_all
```

### Local environment

Go to the main directory of this repository and run:
```
CONFIG_FILE=config.json make all
```