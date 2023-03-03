#!/usr/bin/env python3
"""
process_xenium.py
====================================
Processes Xenium output using spatialdata_io package
"""

import fire
from process_h5ad import h5ad_to_zarr
import spatialdata_io as sio
from spatialdata_io._constants._constants import XeniumKeys
from glob import glob
from pathlib import Path


def xenium_to_zarr(
    path: str,
    stem: str,
    save_h5ad: bool = False,
) -> str:
    """Main function to convert Xenium data in a folder to Zarr

    Args:
        path (str): Path to a SpaceRanger output directory
        dataset_id (str): Prefix of the Xenium experiment
        stem (str): Prefix for the output Zarr filename
        save_h5ad (bool, optional): If the AnnData object should also be written to an h5ad file. Defaults to False.

    Returns:
        str: Output Zarr filename
    """
    spec_path = glob(f"{p}/*_{XeniumKeys.XENIUM_SPECS}")[0]
    dataset_id = Path(spec_path).name.replace(XeniumKeys.XENIUM_SPECS, "")
    sdata = sio.xenium(path, dataset_id)
    if save_h5ad:
        sdata.table.write_h5ad(f"{stem}.h5ad")
    return h5ad_to_zarr(adata=sdata.table, stem=stem, **kwargs)


if __name__ == "__main__":
    fire.Fire(xenium_to_zarr)
