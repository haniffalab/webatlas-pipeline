#!/usr/bin/env python3

import json
import fire
from xml.etree import ElementTree as ET


def get_image_basic_metadata(xml_path):
    NS = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}

    ome_metadata = ET.parse(xml_path)
    dimOrder = ome_metadata.find("./*/ome:Pixels", NS).attrib["DimensionOrder"]
    X = ome_metadata.find("./*/ome:Pixels", NS).attrib["SizeX"]
    Y = ome_metadata.find("./*/ome:Pixels", NS).attrib["SizeY"]
    Z = ome_metadata.find("./*/ome:Pixels", NS).attrib["SizeZ"]
    C = ome_metadata.find("./*/ome:Pixels", NS).attrib["SizeC"]
    T = ome_metadata.find("./*/ome:Pixels", NS).attrib["SizeT"]
    channels = [
        channel.attrib["Name"]
        for channel in ome_metadata.findall("./**/ome:Channel", NS)
        if "Name" in channel.attrib
    ]
    return dimOrder, channels, X, Y, Z, C, T


def main(xml_path):
    dimOrder, channel_names, X, Y, Z, C, T = get_image_basic_metadata(xml_path)
    md = {
        "dimOrder": dimOrder,
        "channel_names": channel_names,
        "X": X,
        "Y": Y,
        "Z": Z,
        "C": C,
        "T": T,
    }
    return json.dumps(md)


if __name__ == "__main__":
    fire.Fire(main)
