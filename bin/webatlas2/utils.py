import csv

# Map entity types used in WebAtlas2.0 UI to the o.var.feature_types values in -anndata.zarr files
spatialxy_colnames_alternatives = ['spatial', 'X_spatial']
feature_type_colnames_alternatives = ['feature_type', 'feature_types']

def get_project_annotation(project_annotations_path, data_key) -> str:
    value = None
    with open(project_annotations_path, 'r') as csvfile:
        csvreader = csv.reader(csvfile, delimiter='\t')
        # skip header
        next(csvreader)
        for row in csvreader:
            key = row[0]
            if key == data_key:
                value = row[1]
                break
    return value