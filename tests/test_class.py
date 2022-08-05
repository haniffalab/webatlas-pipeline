import pytest

import os
import csv
import json
from xml.etree.ElementTree import ElementTree

import xmlschema
import zarr
import numpy as np
import anndata as ad

from bin.process_h5ad import h5ad_to_zarr
from bin.process_molecules import tsv_to_json
from bin.consolidate_md import main as consolidate_md
from bin.router import main as router
from bin.ome_zarr_metadata import main as ome_zarr_metadata
from bin.generate_label import main as generate_label
# from bin.build_config import write_json

class TestClass:

    @pytest.fixture(scope='class')
    def anndata_h5ad_file(self, tmp_path_factory):
        adata = ad.AnnData(np.array([[100.0]*3]*3), dtype=float)
        adata.uns['spatial'] = {'anndata': {'scalefactors': {'spot_diameter_fullres': 10}}}
        adata.obsm['spatial'] = np.array([[20.0,20.0],[40.0,40.0],[60.0,60.0]])
        fn = tmp_path_factory.mktemp("data") / "anndata.h5ad"
        adata.write_h5ad(fn)
        return fn
    
    @pytest.fixture(scope='class')
    def molecules_tsv_file(self, tmp_path_factory):
        fn = tmp_path_factory.mktemp("data") / "molecules.tsv"
        with open(fn, 'w') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter='\t')
            csvwriter.writerow(["Name", "x_int", "y_int"])
            csvwriter.writerow(["background", "0", "0"])
        return fn
    
    @pytest.fixture(scope='class')
    def zarr_file(self, tmp_path_factory):
        fn = tmp_path_factory.mktemp("data") / "dummy.zarr"
        z = zarr.open(
            fn, mode='w', shape=(100, 100),
            chunks=(10, 10), dtype='i4')
        return fn
    
    @pytest.fixture(scope='class')
    def ome_xml_file(self, tmp_path_factory):
        xs = xmlschema.XMLSchema('https://www.openmicroscopy.org/Schemas/OME/2016-06/ome.xsd')
        xml_data = {'@xmlns': 'http://www.openmicroscopy.org/Schemas/OME/2016-06',
            '@xmlns:xsi': 'http://www.w3.org/2001/XMLSchema-instance',
            '@xsi:schemaLocation': 'http://www.openmicroscopy.org/Schemas/OME/2016-06 http://www.openmicroscopy.org/Schemas/OME/2016-06/ome.xsd',
            'Image': [{'@ID': 'Image:0',
                'Pixels': {
                    '@DimensionOrder': 'XYZCT', '@ID': 'Pixels:0', '@Type': 'uint32',
                    '@SizeC': 1, '@SizeT': 1, '@SizeX': 100, '@SizeY': 100, '@SizeZ': 1,
                    'Channel': [{'@ID': 'Channel:0:0', '@Name': 'Channel_0'}], 'TiffData': []}}]}
        xml = xmlschema.from_json(json.dumps(xml_data), schema=xs, preserve_root=True,
            namespaces={"":"http://www.openmicroscopy.org/Schemas/OME/2016-06"}, path="OME")
        xs = xmlschema.XMLSchema("https://www.openmicroscopy.org/Schemas/OME/2016-06/ome.xsd")
        fn = tmp_path_factory.mktemp("data") / "METADATA.ome.xml"
        ElementTree(xml).write(fn)
        return fn
    
    def test_h5ad_to_zarr(self, monkeypatch, anndata_h5ad_file):
        monkeypatch.chdir(os.path.dirname(anndata_h5ad_file))
        stem = "test"
        out_file = h5ad_to_zarr(anndata_h5ad_file, stem)
        assert os.path.exists(out_file)
        assert out_file == stem + "_anndata.zarr"
    
    def test_tsv_to_json(self, monkeypatch, molecules_tsv_file):
        monkeypatch.chdir(os.path.dirname(molecules_tsv_file))
        stem = "test"
        out_file = tsv_to_json(molecules_tsv_file, stem)
        assert os.path.exists(out_file)
        assert out_file == stem + "_molecules.json"
    
    def test_zarr(self, zarr_file):
        md = consolidate_md(zarr_file)
        assert os.path.exists(os.path.join(zarr_file, '.zmetadata'))
    
    def test_route_h5ad(self, monkeypatch, anndata_h5ad_file):
        monkeypatch.chdir(os.path.dirname(anndata_h5ad_file))
        stem = "test"
        out_file = router("h5ad", anndata_h5ad_file, stem, {})
        assert os.path.exists(out_file)
        assert out_file == stem + "_anndata.zarr"
    
    def test_route_molecules(self, monkeypatch, molecules_tsv_file):
        monkeypatch.chdir(os.path.dirname(molecules_tsv_file))
        stem = "test"
        out_file = router("molecules", molecules_tsv_file, stem, {})
        assert os.path.exists(out_file)
        assert out_file == stem + "_molecules.json"

    def test_ome_metadata(self, monkeypatch, ome_xml_file):
        monkeypatch.chdir(os.path.dirname(ome_xml_file))
        stem = "test"
        out_json = json.loads(ome_zarr_metadata(ome_xml_file))
        md = {
            "dimOrder": "XYZCT", "channel_names": ["Channel_0"],
            "X": "100", "Y": "100", "Z": "1", "C": "1", "T": "1"
        }
        assert out_json == md
    
    def test_generate_label(self, monkeypatch, anndata_h5ad_file):
        monkeypatch.chdir(os.path.dirname(anndata_h5ad_file))
        stem = "test"
        ome_md = {"X": 100, "Y": 100}
        generate_label(stem, ome_md, anndata_h5ad_file)
        assert os.path.exists(stem + ".tif")
    
    # def test_build_config(self, request, tmp_path_factory):
    #     input_dir = os.path.dirname(request.fspath) # to get files in tests dir
    #     fn = tmp_path_factory.mktemp("data")
    #     for config_type in ["minimal", "simple", "advanced", "custom"]:
    #         input_file = os.path.join(input_dir, "input", config_type + "_config.json")
    #         with open(input_file) as jsonfile:
    #             input_json = json.load(jsonfile)
    #         input_json["outdir"] = fn
    #         expected_file = os.path.join(input_dir, "expected_output", config_type + "_test_config.json")
    #         with open(expected_file) as jsonfile:
    #             expected_json = json.load(jsonfile)
    #         write_json(**input_json)
    #         output_file = os.path.join(fn, config_type + "_test_config.json")
    #         with open(output_file) as jsonfile:
    #             output_json = json.load(jsonfile)
    #         assert output_json == expected_json
