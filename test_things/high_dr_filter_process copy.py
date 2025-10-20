import json
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from topas_pipeline import preprocess
from topas_pipeline import config

configs = config.load(("./config_patients.json"))

# run preprocessing module

# okay now next thing is printing shape instead and letting it run through and get results.. hopefully it then works



# we have to do it on a branch to make it work.... we can already check if it did not work for fp after grouping


# if it fails again then make it import already existing data to save time...

preprocess.preprocess_raw(
    results_folder=configs.results_folder,
    sample_annotation_file=configs.sample_annotation,
    metadata_annotation=configs.metadata_annotation,
    run_simsi=configs.simsi.run_simsi,
    simsi_folder=configs.simsi.simsi_folder,
    preprocessing_config=configs.preprocessing,
    data_types=configs.data_types,
)
