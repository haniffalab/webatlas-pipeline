#!/usr/bin/env python3

import json
import fire
from xml.etree import ElementTree as ET


def get_metadata(xml_path: str) -> str:
    """Function that parses an OME XML file
    and dumps basic metadata as a JSON formatted str

    Args:
        xml_path (str): Path to OME XML file

    Returns:
        str: JSON formatted metadata
    """
    NS = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}

    ome_metadata = ET.parse(xml_path)
    dimOrder = ome_metadata.find("./*/ome:Pixels", NS).attrib["DimensionOrder"]
    X = ome_metadata.find("./*/ome:Pixels", NS).attrib["SizeX"]
    Y = ome_metadata.find("./*/ome:Pixels", NS).attrib["SizeY"]
    Z = ome_metadata.find("./*/ome:Pixels", NS).attrib["SizeZ"]
    C = ome_metadata.find("./*/ome:Pixels", NS).attrib["SizeC"]
    T = ome_metadata.find("./*/ome:Pixels", NS).attrib["SizeT"]

    channel_names = [
        channel.attrib["Name"]
        for channel in ome_metadata.findall("./**/ome:Channel", NS)
        if "Name" in channel.attrib
    ]

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
    fire.Fire(get_metadata)
