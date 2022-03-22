#!/usr/bin/env python3

# from tifffile import TiffFile
# from apeer_ometiff_library.omexmlClass import OMEXML
# from numcodecs import Zlib
# import zarr
# import numpy as np

import argparse
# from pathlib import Path
import json
# import urllib

# from tile_zarr_base import tile_zarr

# DEFAULT_COMPRESSOR = Zlib(level=1)
# DEFAULT_TILE_SIZE = 512


# def should_be_pyramid(shape):
#     # Generate image pyramid if any dim >= 4096
#     return bool((np.log2(shape) >= 12).any())


# class OmeTiffReader:
#     def __init__(self, img_path):
#         self.tiff = TiffFile(str(img_path))
#         self.metadata = OMEXML(self.tiff.ome_metadata)
#         self.img_path = img_path
#         self.base_series = self.tiff.series[0]
#         self.dim_order = self.base_series.axes
#         self.shape = tuple(dim for dim in self.base_series.shape if dim > 1)
#         self.dtype = self.base_series.dtype

#     def get_raster_dimensions(self):
#         image = self.metadata.image()
#         pixels = image.Pixels
#         channel_names = [
#             str(pixels.Channel(i).Name) for i in range(pixels.channel_count)
#         ]
#         return [
#             {"field": "channel", "type": "nominal", "values": channel_names},
#             {"field": "y", "type": "quantitative", "values": None},
#             {"field": "x", "type": "quantitative", "values": None},
#         ]

#     def to_zarr(
#         self,
#         output_path,
#         tile_size,
#         is_pyramid_base=False,
#         compressor=DEFAULT_COMPRESSOR,
#     ):
#         arr_kwargs = {
#             "chunks": (1, tile_size, tile_size),
#             "compressor": compressor,
#             "shape": self.shape,
#             "dtype": self.dtype,
#         }

#         if is_pyramid_base:
#             group = zarr.open(str(output_path))
#             z = group.create("0", **arr_kwargs)
#         else:
#             z = zarr.open(str(output_path), **arr_kwargs)

#         for idx, page in enumerate(self.base_series.pages):
#             z[idx] = page.asarray()



def write_raster_json():
    image_json = {
        "name": "ISS Human Brain Sample",
        "version": "1.0.0",
        "description": "",
        "public": True,
        "datasets": [
            {
            "uid": "iss-human-brain-simple",
            "name": "iss-human-brain-simple",
            "description": "In-situ Sequencing Human Brain Sample. Raster cell segmentation [bioformats2raw v0.2.6], with additional Vitessce components",
            "files": [
                {
                "type": "raster",
                "fileType": "raster.json",
                "options": {
                    "renderLayers": [
                    "Microscopy Image",
                    "Cells Segmentations"
                    ],
                    "schemaVersion": "0.0.2",
                    "images": [
                    {
                        "name": "Cells Segmentations",
                        "url": "https://a04fcc815aa920b9c7e028bb79f7c2db29d0682c.cog.sanger.ac.uk/39bfafda600ff69122887bce04f4efb88f767caa/various_label_formats/out_opt_flow_registered_label_expanded_0.2.6",
                        "type": "zarr",
                        "metadata": {
                        "isBitmask": True,
                        "dimensions": [
                            {
                            "field": "channel",
                            "type": "nominal",
                            "values": [
                                "Cells"
                            ]
                            },
                            {
                            "field": "y",
                            "type": "quantitative",
                            "values": None
                            },
                            {
                            "field": "x",
                            "type": "quantitative",
                            "values": None
                            }
                        ],
                        "isPyramid": True,
                        "transform": {
                            "translate": {
                            "y": 0,
                            "x": 0
                            },
                            "scale": 1
                        }
                        }
                    },                
                    {
                        "name": "Microscopy Image",
                        "url": "https://a04fcc815aa920b9c7e028bb79f7c2db29d0682c.cog.sanger.ac.uk/da4b9237bacccdf19c0760cab7aec4a8359010b0",
                        "type": "zarr",
                        "metadata": {
                        "dimensions": [
                            {
                            "field": "channel",
                            "type": "nominal",
                            "values": [
                                "c01 DAPI",
                                "c01 Alexa 488",
                                "c01 Atto 425",
                                "c01 Alexa 568",
                                "c01 Alexa 647",
                                "c02 DAPI",
                                "c02 Alexa 488",
                                "c02 Atto 425",
                                "c02 Alexa 568",
                                "c02 Alexa 647",
                                "c03 DAPI",
                                "c03 Alexa 488",
                                "c03 Atto 425",
                                "c03 Alexa 568",
                                "c03 Alexa 647",
                                "c04 DAPI",
                                "c04 Alexa 488",
                                "c04 Atto 425",
                                "c04 Alexa 568",
                                "c04 Alexa 647",
                                "c05 DAPI",
                                "c05 Alexa 488",
                                "c05 Atto 425",
                                "c05 Alexa 568",
                                "c05 Alexa 647",
                                "c06 DAPI",
                                "c06 Alexa 488",
                                "c06 Atto 425",
                                "c06 Alexa 568",
                                "c06 Alexa 647",
                                "c07 DAPI",
                                "c07 Alexa 488",
                                "c07 Atto 425",
                                "c07 Alexa 568",
                                "c07 Alexa 647"
                            ]
                            },
                            {
                            "field": "y",
                            "type": "quantitative",
                            "values": None
                            },
                            {
                            "field": "x",
                            "type": "quantitative",
                            "values": None
                            }
                        ],
                        "isPyramid": True,
                        "transform": {
                            "translate": {
                            "y": 0,
                            "x": 0
                            },
                            "scale": 1
                        }
                        }
                    }
                    ]
                }
                }        
            ]
            }
        ],
        "initStrategy": "auto",
        "layout": [
            {
            "component": "spatial",
            "x": 0,
            "y": 0,
            "w": 9,
            "h": 6
            },
            {
            "component": "layerController",
            "x": 9,
            "y": 0,
            "w": 3,
            "h": 6
            }
        ]
    }
    out_file = open("config.json", "w")
    json.dump(image_json, out_file, indent=2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert OME-TIFF image to full resolution zarr"
    )
    parser.add_argument(
        "--outdir", required=True, help="Path to input OME-TIFF"
    )

    # args = parser.parse_args()

    # img_path = Path(args.input_tiff)
    # zarr_path = Path(args.output_zarr)
    # full_dest_url = urllib.parse.urljoin(args.dest_url, zarr_path.name)

    # reader = OmeTiffReader(img_path)

    # is_pyramid_base = should_be_pyramid(reader.shape)



    # write_raster_json(
    #     json_file="",
    #     url=full_dest_url,
    #     name=args.image_name,
    #     dimensions=reader.get_raster_dimensions(),
    #     is_pyramid=is_pyramid_base,
    # )
    write_raster_json()    