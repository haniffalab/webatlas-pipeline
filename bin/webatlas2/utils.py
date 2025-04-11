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

def get_img_name2scaling_factors(section_annotations_path) -> dict:
    img_name2scaling_factors = {}
    with open(section_annotations_path, 'r') as csvfile:
        csvreader = csv.reader(csvfile, delimiter='\t')
        # skip header
        next(csvreader)
        for row in csvreader:
            img_name = row[0]
            # row[1] - title
            # row[2] - description
            scaling_factor_x = row[3]
            scaling_factor_y = row[4]
            img_name2scaling_factors[img_name] = \
                [float(scaling_factor_x), float(scaling_factor_y)]
    return img_name2scaling_factors

def get_image_name2entity_type2visium_intensity_cutoff(section_annotations_path) -> dict:
    image_name2entity_type2visium_intensity_cutoff = {}
    with open(section_annotations_path, 'r') as csvfile:
        csvreader = csv.reader(csvfile, delimiter='\t')
        # skip header
        next(csvreader)
        for row in csvreader:
            img_name = row[0]
            visium_intensity_cutoffs = row[5]
            for token in visium_intensity_cutoffs.split(","):
                arr = token.split(":")
                entity_type = arr[0]
                visium_intensity_cutoff = arr[1]
                if img_name not in image_name2entity_type2visium_intensity_cutoff:
                    image_name2entity_type2visium_intensity_cutoff[img_name] = {}
                image_name2entity_type2visium_intensity_cutoff[img_name][entity_type] = float(visium_intensity_cutoff)
    return image_name2entity_type2visium_intensity_cutoff