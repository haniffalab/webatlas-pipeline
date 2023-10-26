#!/usr/bin/env python3
"""
process_molecules.py
====================================
Processes molecules files
"""

import os
import csv
import json
import fire
from constants.suffixes import MOLECULES_JSON_SUFFIX


def tsv_to_json(
    path: str,
    stem: str,
    has_header: bool = True,
    gene_col_name: str = "Name",
    x_col_name: str = "x_int",
    y_col_name: str = "y_int",
    delimiter: str = "\t",
    x_scale: float = 1.0,
    y_scale: float = 1.0,
    x_offset: float = 0.0,
    y_offset: float = 0.0,
    gene_col_idx: int = None,
    x_col_idx: int = None,
    y_col_idx: int = None,
    filter_col_name: str = None,
    filter_col_idx: int = None,
    filter_col_value: str = None,
) -> str:
    """This function loads a TSV/CSV file containing gene names, X and Y coordinates
    and writes them to a JSON file supported by Vitessce

    Args:
        path (str): Path to tsv/csv file
        stem (str): Prefix for output JSON file
        has_header (bool, optional): If input file contains a header row. Defaults to True.
        gene_col_name (str, optional): Column header name where gene names are stored. Defaults to "Name".
        x_col_name (str, optional): Column header name where `X` coordinates are stored. Defaults to "x_int".
        y_col_name (str, optional): Column header name where `Y` coordinates are stored. Defaults to "y_int".
        delimiter (str, optional): Input file delimiter. Defaults to "\t".
        x_scale (float, optional): Scale to multiply `X` coordinates by. Defaults to 1.0.
        y_scale (float, optional): Scale to multiply `Y` coordinates by. Defaults to 1.0.
        x_offset (float, optional): Offset to add to `X` coordinates. Defaults to 0.0.
        y_offset (float, optional): Offset to add to `Y` coordinates. Defaults to 0.0.
        gene_col_idx (int, optional): Column index where gene names are stored if header is not present. Defaults to None.
        x_col_idx (int, optional): Column index where `X` coordinates are stored if header is not present. Defaults to None.
        y_col_idx (int, optional): Column index where `Y` coordinates are stored if header is not present. Defaults to None.
        filter_col_name (str, optional): Column header name storing values to filter data. Defaults to None.
        filter_col_idx (int, optional): Column index storing values to filter data if header is not present. Defaults to None.
        filter_col_value (str, optional): Value expected in filter column.
            If a row has a different value it will not be written to output file. Defaults to None.

    Raises:
        SystemExit: If any column header name is not in the header row.
        e: If coordinate values cannot be parsed to float

    Returns:
        str: Output JSON filename
    """

    with open(path) as f:
        reader = csv.reader(f, delimiter=delimiter)

        if has_header:
            header = next(reader, None)
            try:
                if gene_col_idx is None:
                    gene_col_idx = header.index(gene_col_name)
                if x_col_idx is None or y_col_idx is None:
                    x_col_idx = header.index(x_col_name)
                    y_col_idx = header.index(y_col_name)
                if filter_col_name is not None:
                    filter_col_idx = header.index(filter_col_name)
            except ValueError as e:
                raise SystemExit(
                    f"Column name(s), ({gene_col_name}, {x_col_name}, {y_col_name}) not in header"
                )

        molecules_json = {}
        for row in reader:
            try:
                if filter_col_idx is not None:
                    if row[filter_col_idx] != filter_col_value:
                        continue
                gene = row[gene_col_idx]
                molecules_json.setdefault(gene, []).append(
                    [
                        (float(row[x_col_idx]) * x_scale) + x_offset,
                        (float(row[y_col_idx]) * y_scale) + y_offset,
                    ]
                )
            except ValueError as e:
                # TODO: add message
                raise e

    json_file = f"{stem}-{MOLECULES_JSON_SUFFIX}"
    json_file = (
        f"{stem}-{MOLECULES_JSON_SUFFIX}"
        if not stem.endswith("-" + os.path.splitext(MOLECULES_JSON_SUFFIX)[0])
        else f"{stem}.{os.path.splitext(MOLECULES_JSON_SUFFIX)[1]}"
    )

    with open(json_file, "w") as out_file:
        json.dump(molecules_json, out_file)

    return json_file


if __name__ == "__main__":
    fire.Fire(tsv_to_json)
