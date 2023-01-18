#!/usr/bin/env python3
"""
consolidate_md.py
====================================
Consolidates Zarr metadata
"""

import logging
import fire
from zarr import consolidate_metadata
from pathlib import Path


def main(file_in: str) -> None:
    """Function to consolidate the metadata of a Zarr file

    Args:
        file_in (str): Path to Zarr file
    """
    log = logging.getLogger()
    log.setLevel(logging.INFO)
    stem = Path(file_in).stem
    logging.info(stem)
    consolidated = consolidate_metadata(file_in)
    logging.info(consolidated.info)


if __name__ == "__main__":
    fire.Fire(main)
