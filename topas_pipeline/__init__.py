import sys
import logging.handlers
import time
from pathlib import Path

# get version number
__version__ = "0.0.0"
try:
    from importlib.metadata import version, PackageNotFoundError

    try:
        __version__ = version(__name__)
    except PackageNotFoundError:
        pass
except ImportError:
    from pkg_resources import get_distribution, DistributionNotFound

    try:
        __version__ = get_distribution(__name__).version
    except DistributionNotFound:
        pass

__git_commit_hash__ = "unknown-commit"
git_commit_hash_file = Path(__file__).resolve().parent.parent / "hash.file"
if git_commit_hash_file.exists():
    with open(git_commit_hash_file, "r") as file:
        __git_commit_hash__ = (
            file.readline().strip()
        )  # Reads the first line and removes trailing whitespace

__copyright__ = """Copyright (c) 2021-2024 Cecilia Jensen, Firas Hamood, Amirhossein Sakhteman & Matthew The. All rights reserved.
Written by:
- Cecilia Jensen (cecilia.jensen@tum.de)
- Firas Hamood (firas.hamood@tum.de)
- Amirhossein Sakhteman (amirhossein.sakhteman@tum.de)
- Matthew The (matthew.the@tum.de)
at the Chair of Proteomics and Bioanalytics at the Technical University of Munich."""

CONSOLE_LOG_LEVEL = logging.INFO
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
if len(logger.handlers) == 0:
    # formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(name)s::%(funcName)s %(message)s")
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    # formatter.converter = time.gmtime

    # add console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(CONSOLE_LOG_LEVEL)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # add error handler
    error_handler = logging.StreamHandler()
    error_handler.setLevel(logging.ERROR)
    error_handler.setFormatter(formatter)
    logger.addHandler(error_handler)
else:
    logger.info("Logger already initizalized. Resuming normal operation.")
