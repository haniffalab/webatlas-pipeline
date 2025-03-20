from anndata import read_zarr
import pandas as pd
import json
import sys
import time
import csv
import webatlas2.utils as utils

def remove_zeros_rows_cols(df):
    # Drop rows with all zeros
    df = df.loc[(df!=0).any(axis=1)]
    # Drop columns with all zeros
    df = df.loc[:, (df != 0).any(axis=0)]
    return df

def process(project_annotations_path,
            section_annotations_path,
            feature_coordinates_path,
            anndata_zarrs):

    if not section_annotations_path or not feature_coordinates_path or len(anndata_zarrs) < 1:
        print("Please provide the name of the file containing Visium intensity mins, " +
              "the output file name for feature coordinates json, and at least one *-anndata.zarr path")
        sys.exit(1)

    entity_type2img_name2feature2xy_coords_intensity_list = {}
    entity_type2img_name2feature2stat2intensity = {}
    image_name2entity_type2visium_intensity_cutoff = {}
    # Read in visium_intensity_cutoffs
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

    entity_type2feature2min_max_intensity = {}

    continuous_entity_types = \
        utils.get_project_annotation(project_annotations_path, "continuous_entity_types").split(",")
    for entity_type in continuous_entity_types:
        entity_type2feature2min_max_intensity[entity_type] = {}

    # Iterate over zarr files, for each applying the corresponding visium_intensity_cutoff in order to
    # obtain coordinates of features to be displayed over the thumbnails
    for zarr_dir in anndata_zarrs:
        print("Processing {}".format(zarr_dir))
        zarr_fname = zarr_dir.split("/")[-1]
        img_name = "{}.jpeg".format(zarr_fname.replace("-anndata.zarr",""))
        if img_name not in image_name2entity_type2visium_intensity_cutoff:
            print("WARNING: visium_intensity_cutoff missing for {} - skipping".format(img_name))
            continue
            # sys.exit(1)
        try:
            o = read_zarr(zarr_dir)
            for entity_type in continuous_entity_types:
                if entity_type not in image_name2entity_type2visium_intensity_cutoff[img_name]:
                    print("WARNING: visium_intensity_cutoff missing for section: {} - entity_type: {} - defaulting to 0".format(img_name, entity_type))
                    # It may be that some features are missing in certain sections - if that's the case, defaulting the cutoff to 0 is a no-op
                    image_name2entity_type2visium_intensity_cutoff[img_name][entity_type] = 0
                # Initialise data structure for entity_type - img_name
                if entity_type not in entity_type2img_name2feature2xy_coords_intensity_list:
                    entity_type2img_name2feature2xy_coords_intensity_list[entity_type] = {}
                    entity_type2img_name2feature2stat2intensity[entity_type] = {}
                if img_name not in entity_type2img_name2feature2xy_coords_intensity_list[entity_type]:
                    entity_type2img_name2feature2xy_coords_intensity_list[entity_type][img_name] = {}
                    entity_type2img_name2feature2stat2intensity[entity_type][img_name] = {}
                feature2xy_coords_intensity_list = \
                    entity_type2img_name2feature2xy_coords_intensity_list[entity_type][img_name]
                feature2stat2intensity = entity_type2img_name2feature2stat2intensity[entity_type][img_name]

                feature_type = utils.get_project_annotation(project_annotations_path, entity_type)
                features = None
                feature_type_series = None
                for col_name in utils.feature_type_colnames_alternatives:
                    if col_name in o.var:
                        feature_type_series = o.var[col_name]
                        features = list(o.var.index[feature_type_series==feature_type])
                        break
                if features is None:
                    print(
                        "ERROR: none of the feature_type col name alternatives: {} where found in o.var for {}".format(
                            ", ".join(utils.feature_type_colnames_alternatives), zarr_dir))
                    sys.exit(1)
                if len(features) > 0:
                    filtered_x = o.X.T[feature_type_series==feature_type].copy()
                    visium_intensity_cutoff = image_name2entity_type2visium_intensity_cutoff[img_name][entity_type]
                    # Initialise feature2xy_coords_intensity_list for all features
                    for feature in features:
                        if feature not in entity_type2feature2min_max_intensity[entity_type]:
                            # sys.maxsize in python3 corresponds to sys.maxint in python2
                            # see: https://docs.python.org/3/whatsnew/3.0.html#integers
                            entity_type2feature2min_max_intensity[entity_type][feature] = [sys.maxsize, 0]
                        if feature not in feature2xy_coords_intensity_list:
                            feature2xy_coords_intensity_list[feature] = []
                    # Having recorded all features now apply visium_intensity_cutoff
                    filtered_x[filtered_x < visium_intensity_cutoff] = 0
                    barcodes = None
                    barcode_colnames_alternatives = ['label_id', 'cell_id', 'spot_id']
                    for col_name in barcode_colnames_alternatives:
                        if o.obs.columns.isin([col_name]).any():
                            barcodes = list(o.obs[col_name])
                            break
                    if barcodes is None:
                        print("ERROR: none of the barcode col name alternatives: {} where found in o.obs for {}".format(
                              ", ".join(barcode_colnames_alternatives), zarr_dir))
                        sys.exit(1)
                    spatial_xy = None
                    for col_name in utils.spatialxy_colnames_alternatives:
                        if col_name in o.obsm.keys():
                            spatial_xy = list(o.obsm[col_name])
                            break
                    if spatial_xy is None:
                        print("ERROR: none of the spatial_xy col name alternatives: {} where found in o.obsm for {}".format(
                              ", ".join(utils.spatialxy_colnames_alternatives), zarr_dir))
                        sys.exit(1)
                    df = pd.DataFrame(data=filtered_x, index=features, columns=barcodes)
                    df = remove_zeros_rows_cols(df)
                    dict = df.to_dict()
                    start = time.time()
                    feature2total_intensity = {}
                    for barcode in dict:
                        idx = barcodes.index(barcode)
                        if idx % 1000 == 0:
                            end = time.time()
                            print("{} {}s so far".format(idx, round(end - start, 0)))
                        for feature in dict[barcode]:
                            intensity = dict[barcode][feature]
                            if feature not in feature2total_intensity:
                                feature2total_intensity[feature] = intensity
                            else:
                                feature2total_intensity[feature] += intensity
                            if intensity > 0:
                                current_min_intensity = entity_type2feature2min_max_intensity[entity_type][feature][0]
                                current_max_intensity = entity_type2feature2min_max_intensity[entity_type][feature][1]
                                if intensity > current_max_intensity:
                                    entity_type2feature2min_max_intensity[entity_type][feature] = [current_min_intensity, intensity]
                                current_min_intensity = entity_type2feature2min_max_intensity[entity_type][feature][0]
                                current_max_intensity = entity_type2feature2min_max_intensity[entity_type][feature][1]
                                if intensity < current_min_intensity:
                                    entity_type2feature2min_max_intensity[entity_type][feature] = [intensity, current_max_intensity]
                                xy = spatial_xy[idx]
                                x = int(xy[0].astype(object))
                                y = int(xy[1].astype(object))
                                intensity = round(intensity, 2)
                                feature2xy_coords_intensity_list[feature].append((x, y, intensity))
                                if feature not in feature2stat2intensity:
                                    feature2stat2intensity[feature] = {}
                                if 'max' not in feature2stat2intensity[feature] or intensity > feature2stat2intensity[feature]['max']:
                                    feature2stat2intensity[feature]['max'] = intensity
                    for feature in feature2total_intensity:
                        feature2stat2intensity[feature]['avg'] =  int(feature2total_intensity[feature] / len(barcodes))
        except Exception as e:
            print("ERROR: there was an error '{}' reading zarr {} - exiting".format(e, zarr_dir))
            sys.exit(1)

    for entity_type in entity_type2feature2min_max_intensity:
        if entity_type in entity_type2img_name2feature2xy_coords_intensity_list:
            for img_name in entity_type2img_name2feature2xy_coords_intensity_list[entity_type]:
                for feature in entity_type2feature2min_max_intensity[entity_type]:
                    # For a given feature, as the first element of the array of coordinates-intensities
                    # store a tuple of:
                    # 1. the minimum intensity across all sections
                    # 2. the maximum intensity across all sections
                    # 3. the maximum intensity in the section corresponding to img_name (or 1. if feature is not expressed at all in this section)
                    # 4. ditto but with the average intensity
                    # Storing 1 and 2 is so that a given feature the intensity colours shown across all thumbnails are
                    # comparable visually. Storing 3 enables thumbnails to be sorted in the UI - by the highest expression
                    # (of the selected feature) first.
                    # 3 and 4 are shown on a stacked plot when the user selects a feature in the UI
                    min_max = [int(m) for m in entity_type2feature2min_max_intensity[entity_type][feature]]
                    if feature in entity_type2img_name2feature2stat2intensity[entity_type][img_name]:
                        max_intensity_in_section = int(entity_type2img_name2feature2stat2intensity[entity_type][img_name][feature]['max'])
                        avg_intensity_in_section = int(entity_type2img_name2feature2stat2intensity[entity_type][img_name][feature]['avg'])
                    else:
                        if min_max[0] == sys.maxsize:
                            # Feature is not expressed in any section
                            min_max[0] = 0
                        max_intensity_in_section = min_max[0]
                        avg_intensity_in_section = min_max[0]
                    stats = min_max + [max_intensity_in_section, avg_intensity_in_section]

                    if feature in entity_type2img_name2feature2xy_coords_intensity_list[entity_type][img_name]:
                        if len(entity_type2img_name2feature2xy_coords_intensity_list[entity_type][img_name][feature]) == 0:
                            # if feature has no expressions above the minimum for a given section, both min and max
                            # should be minimum_intensity cutoff
                            stats[1] = stats[0]
                        entity_type2img_name2feature2xy_coords_intensity_list[entity_type][img_name][feature].insert(0, stats)

    with open(feature_coordinates_path, 'w') as f:
        f.write(json.dumps(entity_type2img_name2feature2xy_coords_intensity_list))

