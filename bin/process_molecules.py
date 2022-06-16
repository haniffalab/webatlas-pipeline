#!/usr/bin/env python3

import os
import csv
import json
import fire

SUFFIX = 'molecules.json'

def tsv_to_json(
    file,
    stem,
    has_header=True,
    gene_col_name='Name',
    x_col_name='x_int',
    y_col_name='y_int',
    delimiter=',',
    gene_col_idx=None,
    x_col_idx=None,
    y_col_idx=None
    ):


    with open(file) as f:
        reader = csv.reader(f, delimiter=delimiter)

        if has_header:
            header = next(reader, None)
            try:
                if gene_col_idx is None:
                    gene_col_idx = header.index(gene_col_name)
                if x_col_idx is None or y_col_idx is None:
                    x_col_idx = header.index(x_col_name)
                    y_col_idx = header.index(y_col_name)
            except ValueError as e:
                print(f"Column name(s), ({gene_col_name}, {x_col_name}, {y_col_name}) not in header")
                return e

        molecules_json = {}
        for row in reader:
            try:
                gene = row[gene_col_idx]
                molecules_json.setdefault(gene, []).append([float(row[x_col_idx]), float(row[y_col_idx])])
            except ValueError as e:
                print(e)

    json_file = f"{stem}_{SUFFIX}"

    with open(json_file, 'w') as out_file:
        json.dump(molecules_json, out_file)
    
    return json_file


if __name__ == "__main__":
    fire.Fire(tsv_to_json)
