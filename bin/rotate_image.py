#!/usr/bin/env python3
"""
rotate_image.py
====================================
Rotate tif images
"""

from __future__ import annotations
import fire
import logging
import pyvips
from typing import T
import xml.etree.ElementTree as ET

logging.getLogger().setLevel(logging.INFO)


def rotate_image(
    input_path: str, output_path: str, degrees: T.Literal[90, 180, 270]
) -> None:
    if degrees not in [90, 180, 270]:
        raise ValueError("degrees must be 90, 180, or 270")

    image = pyvips.Image.tiffload(input_path, n=-1)
    # No sequential access as rotating would cause issues

    logging.info(
        f"""Input image:
        Width: {image.width}
        Height: {image.height}
        Pages: {image.get_n_pages()}
        Page height: {image.get_page_height()}
        """
    )

    pages = image.pagesplit()
    for i in range(len(pages)):
        if degrees == 90:
            pages[i] = pages[i].rot270()
        elif degrees == 180:
            pages[i] = pages[i].rot180()
        elif degrees == 270:
            pages[i] = pages[i].rot90()

    first_page, *rest = pages
    rotated_image = first_page.pagejoin(rest)

    logging.info(
        f"""Rotated image:
        Width: {rotated_image.width}
        Height: {rotated_image.height}
        Pages: {rotated_image.get_n_pages()}
        Page height: {rotated_image.get_page_height()}
        """
    )

    try:
        image_description = image.get("image-description")
        if image_description:
            NS = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}
            ET.register_namespace("", NS["ome"])
            root = ET.fromstring(image_description)
            pixels = root.find("./*/ome:Pixels", NS)
            if pixels is not None:
                pixels.set("SizeX", str(rotated_image.width))
                pixels.set("SizeY", str(rotated_image.height))
            rotated_image.set_type(
                pyvips.GValue.gstr_type,
                "image-description",
                ET.tostring(root, encoding="unicode", xml_declaration=True),
            )
    except:
        pass

    rotated_image.tiffsave(output_path, bigtiff=True)


if __name__ == "__main__":
    fire.Fire(rotate_image)
