import json

from pathlib import Path


def load(config_file):
    with open(Path(__file__).parent.parent / "config_defaults.json", "r") as inp:
        config_defaults = json.load(inp)

    with open(config_file, "r") as inp:
        configs = json.load(inp)
    
    configs = merge_dicts(config_defaults, configs)

    return configs


def merge_dicts(defaults, custom):
    for key, value in custom.items():
        if isinstance(value, dict) and key in defaults and isinstance(defaults[key], dict):
            defaults[key] = merge_dicts(defaults[key], value)  # Recursive merge
        else:
            defaults[key] = value  # Overwrite default value
    return defaults