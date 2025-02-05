from __future__ import annotations

import logging
import os

import fire
import pyvips
import tifffile as tf

logging.getLogger().setLevel(logging.INFO)


def multitif2tif(path: str, prefix: str, outdir: str = "."):
    logging.info(f"Reading path {path}")
    channels = [
        os.path.join(path, x)
        for x in os.listdir(path)
        if x.endswith(".tiff") or x.endswith(".tif")
    ]
    logging.info(f"Found {len(channels)} channels")

    c_imgs = [pyvips.Image.new_from_array(tf.imread(x, is_ome=False)) for x in channels]

    full_img = pyvips.Image.arrayjoin(c_imgs, across=1)
    full_img = full_img.copy()
    full_img.set_type(pyvips.GValue.gint_type, "page-height", c_imgs[0].height)

    xml_channels = "".join(
        [
            f"""<Channel ID="Channel:0:{x}" SamplesPerPixel="1" Name="{os.path.basename(c)}"><LightPath/></Channel>"""
            for x, c in enumerate(channels)
        ]
    )

    full_img.set_type(
        pyvips.GValue.gstr_type,
        "image-description",
        " ".join(
            f"""<?xml version="1.0" encoding="UTF-8"?>
        <OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06"
            xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
            xsi:schemaLocation="http://www.openmicroscopy.org/Schemas/OME/2016-06 http://www.openmicroscopy.org/Schemas/OME/2016-06/ome.xsd">
            <Image ID="Image:0">
                <Pixels DimensionOrder="XYCZT"
                        ID="Pixels:0"
                        SizeC="{len(c_imgs)}"
                        SizeT="1"
                        SizeX="{c_imgs[0].width}"
                        SizeY="{c_imgs[0].height}"
                        SizeZ="1"
                        Type="uint16">
                        {xml_channels}
                </Pixels>
            </Image>
        </OME>""".split()
        ),
    )

    out_filename = prefix + os.path.basename(os.path.dirname(path)) + ".tif"
    out_path = os.path.join(outdir, out_filename)

    logging.info(f"Writing file to {out_path}")
    full_img.tiffsave(out_path, bigtiff=True)


if __name__ == "__main__":
    fire.Fire(multitif2tif)
