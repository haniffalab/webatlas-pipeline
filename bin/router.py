#!/usr/bin/env python3

import fire

from process_h5ad import h5ad_to_zarr
from process_molecules import tsv_to_json

def main(**args):
    file_type = args.pop('type')
    file_path = args.pop('file')
    
    if file_type == 'h5ad':
        h5ad_to_zarr(file=file_path, **args)
    if file_type == 'molecules':
        tsv_to_json(file=file_path, **args)

if __name__ == "__main__":
    fire.Fire(main)