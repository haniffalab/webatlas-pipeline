import os

from bin.build_config_multimodal import write_json as write_json_multimodal


class TestClass:
    def test_build_config_multimodal(self):
        input = {
            "integrated": True,
            "project": "test_project",
            "extended_features": "celltype",
            "url": "http://localhost/",
            "config_filename_suffix": "config.json",
            "datasets": {
                "iss_dataset": {
                    "file_paths": ["test_project-iss_dataset-anndata.zarr"],
                    "images": {
                        "raw": [
                            {
                                "path": "/path/to/raw-image.zarr",
                                "md": {
                                    "dimOrder": "XYZT",
                                    "channel_names": ["Channel_1"],
                                    "X": 10,
                                    "Y": 10,
                                    "Z": 1,
                                    "C": 1,
                                    "T": 0,
                                },
                            }
                        ],
                        "label": [
                            {
                                "path": "/path/to/label-image.zarr",
                                "md": {
                                    "dimOrder": "XYZT",
                                    "channel_names": ["Channel_1"],
                                    "X": 10,
                                    "Y": 10,
                                    "Z": 1,
                                    "C": 1,
                                    "T": 0,
                                },
                            }
                        ],
                    },
                    "options": {
                        "matrix": "X",
                        "factors": ["obs/sample"],
                        "mappings": {"obsm/X_umap": [0, 1]},
                        "sets": ["obs/cluster"],
                        "spatial": {"xy": "obsm/spatial"},
                    },
                    "obs_type": "cell",
                    "is_spatial": True,
                },
                "scrnaseq_dataset": {
                    "file_paths": ["test_project-scrnaseq_dataset-anndata.zarr"],
                    "options": {
                        "matrix": "X",
                        "factors": ["obs/sample"],
                        "mappings": {"obsm/X_umap": [0, 1]},
                        "spatial": {"xy": "obsm/spatial"},
                    },
                    "obs_type": "cell",
                    "is_spatial": False,
                },
            },
        }

        write_json_multimodal(**input)

        assert os.path.exists(
            f"{input['project']}-multimodal-{input['config_filename_suffix']}"
        )
