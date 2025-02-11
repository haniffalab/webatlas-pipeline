import os
from anndata import read_zarr
import json
import sys
import csv
import webatlas2.utils as utils

def get_entity2children(output_path: str, project_annotations_path: str, entity_type: str):
    entity2children = {}
    hierarchy_tsv_fname = utils.get_project_annotation(project_annotations_path, entity_type)
    if hierarchy_tsv_fname is not None:
        hierarchy_tsv_fpath = os.path.join(os.path.dirname(project_annotations_path), hierarchy_tsv_fname)
        with open(hierarchy_tsv_fpath, 'r') as csvfile:
            csvreader = csv.reader(csvfile, delimiter='\t')
            # skip header
            next(csvreader)
            for row in csvreader:
                parents = []
                for entity in row:
                    if entity:
                        if entity not in entity2children:
                            entity2children[entity] = set([])
                        for parent in parents:
                            entity2children[parent].add(entity)
                        parents.append(entity)
    return entity2children

def process(output_path, project_annotations_path, zarr_anndata_paths):
    if not output_path or not project_annotations_path or not zarr_anndata_paths:
        print("Please provide the output path, "
              " the file path for project_annotations.tsv"
              " and at least one *-anndata.zarr path (input")
        sys.exit(1)

    categorical_entity_types = \
        utils.get_project_annotation(project_annotations_path, "categorical_entity_types").split(",")
    for entity_type in categorical_entity_types:
        print("About to retrieve hierarchical entity coordinates for entity_type: {} ..".format(entity_type))

        # NB. The hierarchical entity value coordinates are stored at the lowest level of the hierarchy - hence
        # for each hierarchical entity value we need to aggregate coordinates from all its children. Therefore:

        # 1. Collect in entity children all recursive children of each hierarchical entity value
        entity2children = get_entity2children(output_path, project_annotations_path, entity_type)

        # 1a. Collect all entities from entity2children into all_entities - used for QC below
        all_entities = set([])
        for key, values in entity2children.items():
            all_entities.add(key)
            for val in values:
                all_entities.add(val)

        if len(all_entities) > 0:
            # 2. Collect coordinates at the lowest level of the hierarchy
            img_name2entity2xy_coords_list = {}
            for zarr_dir in zarr_anndata_paths:
                print("Processing " + zarr_dir + "..")
                zarr_fname = zarr_dir.split("/")[-1]
                img_name = "{}.jpeg".format(zarr_fname.replace("-anndata.zarr",""))
                img_name2entity2xy_coords_list[img_name] = {}
                o = read_zarr(zarr_dir)

                spatial_xy = None
                for col_name in utils.spatialxy_colnames_alternatives:
                    if col_name in o.obsm.keys():
                        spatial_xy = list(o.obsm[col_name])
                        break
                if spatial_xy is None:
                    print("ERROR: none of the spatial_xy col name alternatives: {} where found in o.obsm for {}".format(
                        ", ".join(utils.spatialxy_colnames_alternatives), zarr_dir))
                    sys.exit(1)

                obsColName = utils.get_project_annotation(project_annotations_path, "{}_obs_col".format(entity_type))
                hierarchical_entity_values = o.obs[obsColName].values.tolist()
                for idx, entity in enumerate(hierarchical_entity_values):
                    if entity not in img_name2entity2xy_coords_list[img_name]:
                        img_name2entity2xy_coords_list[img_name][entity] = []
                    xy = spatial_xy[idx]
                    img_name2entity2xy_coords_list[img_name][entity].append((xy[0].astype(object), xy[1].astype(object)))

            # 3. Now propagate the coordinates of the children to all their parents recursively
            for img_name in img_name2entity2xy_coords_list:
                for entity in entity2children:
                    for child in entity2children[entity]:
                        if child in img_name2entity2xy_coords_list[img_name]:
                            if entity not in img_name2entity2xy_coords_list[img_name]:
                                # This could be because 1. hierarchical entity value is misspelled in the data compared
                                # to that in the hierarchy for entity_type, or
                                # 2. Since coordinates in the data are included only at the lowest level of the
                                # hierarchy, coordinates for the parents won't be there - hence we need to aggregate
                                # them here - so that users can search by hierarchical entity values at all levels of the hierarchy.
                                # Note that 1. above is reported in 3a. section below.
                                img_name2entity2xy_coords_list[img_name][entity] = []
                            img_name2entity2xy_coords_list[img_name][entity] += img_name2entity2xy_coords_list[img_name][child]

                # 3a. Report hierarchical entity values in img_name2entity2xy_coords_list[img_name] that are not found in all_entities
                # (could be because of misspelling)
                missing = []
                for entity in img_name2entity2xy_coords_list[img_name]:
                    if entity not in all_entities:
                        missing.append(str(entity))
                if missing:
                    print("WARNING: The following hierarchical entity values for {} are missing in the hierarchy for {}:".format(img_name, entity_type))
                    print("\n".join(sorted(missing)))

            # 4. Now output the coordinates into out_fpath
            output_fpath = "{}/{}{}".format(output_path, entity_type, "_coordinates.json")
            with open(output_fpath, 'w') as f:
                f.write(json.dumps(img_name2entity2xy_coords_list))

