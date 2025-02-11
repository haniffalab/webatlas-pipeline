import os
import sys
import re
from ome_zarr.io import parse_url
from ome_zarr.reader import Node, Plate, Reader, Well
import tifffile as tiff
import numpy as np
from anndata import read_zarr
import webatlas2.utils as utils

RAW_ZARR_REGEX='raw.zarr'
ANNDATA_ZARR_REGEX='anndata.zarr'

def process(output_dir, in_zarr_fname, level, project_annotations_path):

    if not output_dir or not in_zarr_fname or not level:
        print("Please provide the output dir, the zarr file name and the level of the thumbnail to be extracted")
        sys.exit(1)

    out_fname = in_zarr_fname.split("/")[-1].replace("-{}".format(RAW_ZARR_REGEX), ".jpeg")
    reader = Reader(parse_url(str("{}/0".format(in_zarr_fname))))()
    image_node = list(reader)[0]
    num_levels_available = len(image_node.data)
    # DEBUG print("Number of levels available: {}".format(num_levels_available))
    if level > num_levels_available:
        print("ERROR: The requested level: {} does not exist in the {}".format(level, in_zarr_fname))
        sys.exit(1)
    nd_arr = image_node.data[level]
    if nd_arr.dtype == ">u2":
        # Pre-emptive conversion to np.uint8 - images in Xenium raw zarr files appear to have
        # dtype=>u2 (> == big endian)
        nd_arr = nd_arr.astype(np.uint8)
    output_fpath = os.path.join(output_dir, out_fname)
    tiff.imwrite(output_fpath, nd_arr, compression='jpeg')
    # Retrieve title
    # TODO: Review if the section title extraction logic works for all projects
    m = re.search(r'^.+?\-(.*)\.jpeg$', out_fname)
    section_title = m.group(1)
    # Retrieve min Visium intensity from ANNDATA_ZARR_REGEX
    anndata_zarr_fname = in_zarr_fname.replace("-{}".format(RAW_ZARR_REGEX), "-{}".format(ANNDATA_ZARR_REGEX))
    o = read_zarr(anndata_zarr_fname)
    visium_intensity_cutoffs = ""
    continuous_entity_types = \
        utils.get_project_annotation(project_annotations_path, "continuous_entity_types").split(",")
    for entity_type in continuous_entity_types:
        feature_type = utils.get_project_annotation(project_annotations_path, entity_type)
        for col_name in utils.feature_type_colnames_alternatives:
            if col_name in o.var:
                if o.var[o.var[col_name]== feature_type].shape[0] > 0:
                    # features of feature_type are present in o.var
                    hist, bins = np.histogram(o.X.T[o.var[col_name] == feature_type], 100)
                    # DEBUG print(out_fname, hist[0:5], bins[4])
                    # The Visium intensity cutoff corresponds to >=0.05 chance of encountering it
                    visium_intensity_cutoff = bins[4]
                    prefix = ""
                    if visium_intensity_cutoffs != "":
                        prefix = ","
                    visium_intensity_cutoffs = \
                        visium_intensity_cutoffs + prefix + entity_type + ":" + str(visium_intensity_cutoff)
                break
    return (section_title, visium_intensity_cutoffs)