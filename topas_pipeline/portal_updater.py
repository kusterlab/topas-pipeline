import json
import urllib.request

from . import config
"""
For Automatic updating of the portal with the TOPAS-pipeline generated results 

USAGE:
python portal_updater.py -c configs.json
"""


def main(results_folder: str, sample_annotation: str, metadata_annotation: str, config_portal: config.Portal):
    if not config_portal.update:
        return

    new_portal_paths = {
        "sample_annotation_path": sample_annotation,
        "patient_annotation_path": metadata_annotation,
        "report_directory": results_folder,
    }
    # updating the portal config file
    for path_name, path in new_portal_paths.items():
        config_updater(config_portal.config, config_portal.cohort, path_name, path)
    endpoint_runner(config_portal.url, config_portal.cohort)


def config_reader(config_path: str):
    with open(config_path, "r") as f:
        config = json.load(f)
    return config


def config_writer(config, config_path: str):
    jsonFile = open(config_path, "w+")
    jsonFile.write(json.dumps(config, indent=4))
    jsonFile.close()


def config_updater(config_path: str, cohort: str, key: str, value: str):
    config = config_reader(config_path)
    config[key][cohort] = value
    # Save our changes to JSON file
    config_writer(config, config_path)


def endpoint_runner(portal_url, cohort, endpoint: str = "reload"):
    adress = f"{portal_url}/{endpoint}/{cohort}"
    with urllib.request.urlopen(adress) as url:
        data = url.read()
        # json.loads(data)


if __name__ == "__main__":
    import sys
    import json
    import argparse
    from . import config

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--config", required=True, help="Absolute path to configuration file."
    )

    args = parser.parse_args(sys.argv[1:])
    configs = config.load(args.config)
    main(configs)
