from anndata import read_zarr
import numpy as np
import json
import sys
import time
import os
import webatlas2.utils as utils

def process(project_annotations_path,
            rnaseq_expressions_path,
            zarr_dir):

    if not rnaseq_expressions_path or not zarr_dir:
        print("Please provide the output file name for RNASeq expressions json, and the directory containing the RNAseq anndata.zarr")
        sys.exit(1)

    dot_plot_entity_types = \
        utils.get_project_annotation(project_annotations_path, "dot_plot_entity_types").split(",")
    rnaseq_plot_entity_type = \
        utils.get_project_annotation(project_annotations_path, "rnaseq_plot_entity_type")
    rnaseq_plot_entity_type_obs_col = \
        utils.get_project_annotation(project_annotations_path, "{}_obs_col".format(rnaseq_plot_entity_type))
    rnaseq_plot_entities = \
        utils.get_project_annotation(project_annotations_path, "rnaseq_plot_entities")
    scrnaseq_zarr = \
        utils.get_project_annotation(project_annotations_path, "scrnaseq_zarr")
    if scrnaseq_zarr is None:
        print("No scrnaseq_zarr file found - exiting")
        sys.exit(0)
    if rnaseq_plot_entities is not None:
        rnaseq_plot_entities = rnaseq_plot_entities.split(",")
    zarr =  os.path.join(zarr_dir, scrnaseq_zarr)
    # E.g. entity_type = 'gene', feature = 'A2M', annot = 'Hypoxic'
    entity_type2feature2expressions = {}

    try:
        start = time.time()
        o = read_zarr(zarr)
        for entity_type in dot_plot_entity_types:
            print("Processing entity_type: {} ...".format(entity_type), flush=True)
            if entity_type not in entity_type2feature2expressions:
                entity_type2feature2expressions[entity_type] = {}
            feature_type = utils.get_project_annotation(project_annotations_path, entity_type)
            features = None
            for col_name in utils.feature_type_colnames_alternatives:
                if col_name in o.var:
                    feature_type_series = o.var[col_name]
                    features = list(o.var.index[feature_type_series == feature_type])
                    break
            if features is None:
                print("WARNING: Did not find entity_type in {} - skipping".format(entity_type, zarr), flush=True)
                continue
            if len(features) > 0 and rnaseq_plot_entities is not None:
                for annot in rnaseq_plot_entities:
                    print("Processing annot: {} ...".format(annot))
                    # Filter o.X (cols: features, rows: cells) by feature_type and annot
                    row_mask = o.obs[rnaseq_plot_entity_type_obs_col] == annot
                    col_mask = o.var[col_name] == feature_type
                    filtered_x = o.X[:, col_mask]
                    filtered_x = filtered_x[row_mask, :]
                    # Mean across all cells annotated with annot, per feature
                    annot_mean = np.mean(filtered_x, axis=0)
                    # min_max normalisation (default for scanpy)
                    annot_mean = (annot_mean - annot_mean.min())/(annot_mean.max() - annot_mean.min())
                    mean_expressions = [round(x, 2) for x in annot_mean.tolist()]
                    # Fraction of cells with expression  > 0
                    annot_fraction = np.sum(filtered_x > 0, axis=0) / filtered_x.shape[0]
                    annot_fractions = [round(x, 2) for x in annot_fraction.tolist()]
                    for idx, feature in enumerate(features):
                        expression = mean_expressions[idx]
                        fraction = annot_fractions[idx]
                        if feature not in entity_type2feature2expressions[entity_type]:
                            entity_type2feature2expressions[entity_type][feature] = []
                        entity_type2feature2expressions[entity_type][feature].append((expression, fraction))
        end = time.time()
        print("Duration: {}s ".format(round(end - start, 0)), flush=True)
    except Exception as e:
        print("ERROR: there was an error '{}' reading zarr {} - exiting".format(e, zarr))
        sys.exit(1)
    with open(rnaseq_expressions_path, 'w') as f:
        f.write(json.dumps(entity_type2feature2expressions))

