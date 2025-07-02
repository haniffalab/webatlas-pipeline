import json
import urllib.request, urllib.parse
import sys
import os
from pathlib import Path

raw_zarr_regex='raw.zarr'

def process(in_zarr_fname, level):

    if not in_zarr_fname or not level:
        print("Please provide name of the zarr file and the level of the thumbnail " +
              "for which the scaling factors are to be extracted")
        sys.exit(1)

    x_scaling_factor = ""
    y_scaling_factor = ""
    fpath = os.path.join(in_zarr_fname, "0", ".zattrs")
    f = Path(fpath)
    if f.is_file():
        r = urllib.request.urlopen("file:{}".format(fpath))
        data = json.load(r)
        out_fpath = in_zarr_fname.replace("-{}".format(raw_zarr_regex), ".jpeg")
        out_fname = out_fpath.split("/")[-1]
        pixel_xy_sizes = []
        for el in data['multiscales'][0]['datasets']:
            if 'coordinateTransformations' in el:
                l = el['coordinateTransformations'][0]['scale']
                pixel_xy_sizes.append((l[4], l[3]))
        if pixel_xy_sizes:
            base_level = 0
            x_scaling_factor = pixel_xy_sizes[base_level][0] / pixel_xy_sizes[level][0]
            y_scaling_factor = pixel_xy_sizes[base_level][1] / pixel_xy_sizes[level][1]
    return (str(x_scaling_factor), str(y_scaling_factor))


