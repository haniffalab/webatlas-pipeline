#!/usr/bin/env python3

import fire

from process_h5ad import h5ad_to_zarr
from process_molecules import tsv_to_json

def main(**args):
    file_type = args.pop('type')
    file_path = args.pop('file')
    stem = args.pop('stem')

    out_file = None
    
    if file_type == 'h5ad':
        out_file = h5ad_to_zarr(file=file_path, stem=stem, **args)
    if file_type == 'molecules':
        out_file = tsv_to_json(file=file_path, stem=stem, **args)
    
    return out_file

if __name__ == "__main__":
    fire.Fire(main)