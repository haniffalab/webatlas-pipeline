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

    continuous_entity_types = \
        utils.get_project_annotation(project_annotations_path, "continuous_entity_types").split(",")
    rnaseq_plot_entity_type = \
        utils.get_project_annotation(project_annotations_path, "rnaseq_plot_entity_type")
    rnaseq_plot_entity_type_obs_col = \
        utils.get_project_annotation(project_annotations_path, "{}_obs_col".format(rnaseq_plot_entity_type))
    rnaseq_plot_entities = \
        utils.get_project_annotation(project_annotations_path, "rnaseq_plot_entities").split(",")
    scrnaseq_zarr = \
        utils.get_project_annotation(project_annotations_path, "scrnaseq_zarr")
    zarr =  os.path.join(zarr_dir, scrnaseq_zarr)
    # E.g. entity_type = 'gene', feature = 'A2M', annot = 'Hypoxic'
    entity_type2feature2expressions = {}

    try:
        start = time.time()
        o = read_zarr(zarr)
        for entity_type in continuous_entity_types:
            print("Processing {} ...".format(entity_type))
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
                print("WARNING: Did not find entity_type in {} - skipping".format(entity_type, zarr))
                continue
            if len(features) > 0:
                for annot in rnaseq_plot_entities:
                    print("Processing {} ...".format(annot))
                    # Filter o.X by feature_type and annot
                    row_mask = o.obs[rnaseq_plot_entity_type_obs_col] == annot
                    col_mask = o.var['feature_types'] == feature_type
                    filtered_x = o.X[:, col_mask]
                    filtered_x = filtered_x[row_mask, :]
                    mean_expressions = np.mean(filtered_x, axis=0)
                    mean_expressions = [round(x, 1) for x in mean_expressions.tolist()]
                    for idx, feature in enumerate(features):
                        expression = mean_expressions[idx]
                        if feature not in entity_type2feature2expressions[entity_type]:
                            entity_type2feature2expressions[entity_type][feature] = []
                        entity_type2feature2expressions[entity_type][feature].append(expression)
        end = time.time()
        print("Duration: {}s ".format(round(end - start, 0)))
    except Exception as e:
        print("ERROR: there was an error '{}' reading zarr {} - exiting".format(e, zarr))
        sys.exit(1)
    with open(rnaseq_expressions_path, 'w') as f:
        f.write(json.dumps(entity_type2feature2expressions))

