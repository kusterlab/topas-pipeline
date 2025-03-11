import json


def load(config_file):
    with open(config_file, "r") as inp:
        configs = json.load(inp)

    # TODO: replace this with a default configuration file
    if "data_types" not in configs:
        configs["data_types"] = ["fp", "pp"]

    if "run_lfq" not in configs["preprocessing"]:
        configs["preprocessing"]["run_lfq"] = False
    
    if "normalize_to_reference" not in configs["preprocessing"]:
        configs["preprocessing"]["normalize_to_reference"] = False

    if "num_threads" not in configs["simsi"]:
        configs["simsi"]["num_threads"] = 4

    if "portal" not in configs:
        configs["portal"] = {"update": 0}

    if "slack_webhook_url" not in configs:
        configs["slack_webhook_url"] = ""

    return configs
