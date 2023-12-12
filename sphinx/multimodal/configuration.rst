.. _multimodal_configuration:

########################
Multimodal configuration
########################

After running the main conversion pipeline you can populate the required YAML parameters file to run the multimodal integration pipeline.

.. _multimodal_parameters_file:

***************
Parameters file
***************

The parameters file looks like this:

.. code-block:: yaml

    outdir: "/path/to/output/"

    url: http://localhost:3000/
    project: my_project
    title: "My Project"

    data:
      -
        dataset: scrnaseq
        obs_type: cell
        anndata: /path/to/main/output/scrnaseq-anndata.zarr
        offset: 0
        is_spatial: false
        vitessce_options:
          spatial:
            xy: obsm/spatial
          mappings:
            obsm/X_umap: [0,1]
          matrix: X
      -
        dataset: visium
        obs_type: spot
        anndata: /path/to/main/output/visium-anndata.zarr
        offset: 1000000
        is_spatial: true
        raw_image: /path/to/main/output/visium-raw.zarr
        label_image: /path/to/main/output/visium-label.zarr
        vitessce_options:
          spatial:
            xy: obsm/spatial
          matrix: X

In contrast to the main conversion pipeline's parameters file, this file includes a single `project` to which multiple `datasets` belong.

Each ``dataset`` block defines the name of the dataset and paths to the converted data and image files (if any).

Each ``dataset`` also requires a set of ``vitessce_options`` that specify the location of certain data (spatial coordinates, embeddings, expression matrix, etc.) within the AnnData object that is processed/generated.
This follows the same structure as in the :ref:`main pipeline's vitessce_options <vitessce_options>`.

Additionally, each ``dataset`` requires:

* ``obs_type``, a string indicating the type of observation of the dataset. For example, "cell" or "spot".
* ``offset``, an integer offset to add to the dataset's ID's so they don't clash with the other datasets.
* ``is_spatial``, a boolean indicating whether the dataset contains spatial information and has associated image files (raw and/or label images)

Given that raw images are only read but not modified the pipeline does not generate new output from them.
In order for the output directory (defined by ``outdir``) to contain all necessary files that need to be served for the web application to consume,
by default, the pipeline copies the raw images to the output directory (unless a file with the same name already exists in the output directory).
This process can take a long time depending on the size of the image.
You may want to manually copy or move the image or serve it from its own directory separate from the rest of the output.
The default copying can be disabled by setting ``copy_raw: false`` as a project-wide parameter (at the same level as ``outdir``, ``project``, etc).
For example,

.. code-block:: yaml

    outdir: "/path/to/output/"
    url: http://localhost:3000/
    project: my_project
    title: "My Project"
    copy_raw: false


With additional features
========================

Using the above example parameters file to run the multimodal integration pipeline will run the reindexing and intersection steps.
To perform the concatenation of additional features (like celltypes) to visualise them as continuous values, some extra parameters need to be added.

As a project-wide parameter (at the same level as ``outdir``, ``project``, etc.):

* ``extend_feature_name``, the name of the additional feature. For example, "celltype"

And at a ``dataset`` level:

* ``extend_feature``, the location of the additional feature information.
  This can be either the path to a *cell2location* output file, or the location within the AnnData object where the feature is stored as a categorical within ``obs``.
  For example, ``/path/to/c2l.h5ad`` containing predicted continuous values, or ``obs/celltype`` containing categoricals.

The full parameters file will then look like this

.. code-block:: yaml

    outdir: "/path/to/output/"

    url: http://localhost:3000/
    project: my_project
    title: "My Project"

    extend_feature_name: celltype

    data:
      -
        dataset: scrnaseq
        obs_type: cell
        anndata: /path/to/main/output/scrnaseq-anndata.zarr
        extend_feature: obs/celltype
        offset: 0
        is_spatial: false
        vitessce_options:
          spatial:
            xy: obsm/spatial
          mappings:
            obsm/X_umap: [0,1]
          matrix: X
      -
        dataset: visium
        obs_type: spot
        anndata: /path/to/main/output/visium-anndata.zarr
        extend_feature: /path/to/c2l.h5ad
        offset: 1000000
        is_spatial: true
        raw_image: /path/to/main/output/visium-raw.zarr
        label_image: /path/to/main/output/visium-label.zarr
        vitessce_options:
          spatial:
            xy: obsm/spatial
          matrix: X

With these parameters the multimodal integration pipeline will concatenate the expression matrix with the additional feature values so both can be queried and visualised across datasets within the same portal.

In the case of providing a *cell2location* output file, you can further configure ``extend_feature`` with arguments for how the file should be processed.
Instead of only setting the path to the file you would need to define ``extend_feature`` as a map containing ``path`` and optional ``args``.

.. code-block:: yaml

    extend_feature_name: celltype
    data:
      -
        dataset: visium
        obs_type: spot
        anndata: /path/to/main/output/visium-anndata.zarr
        extend_feature: 
          path: /path/to/c2l.h5ad
          args:
            sample: ["library_id", "sample_1"] # tuple containing the obs column name and value to filter the object. By default the object is not filtered.
            q: "q05_cell_abundance_w_sf" # matrix in obsm to use. Defaults to "q05_cell_abundance_w_sf".
            sort_index: "index_column" # column in the AnnData object that contains an index that matches the index in cell2location.
            sort: True # can be set to False to skip ordering the cell2location matrix but observations might not match in order between files. Defaults to True.

For example, ``sample`` can be used when a *cell2location* output file contains predictions for multiple samples.
Setting ``sample`` to filter the output file enables the pipeline to obtain the appropriate prediction matrix for the data being processed,
without having to split the *cell2location* output file for each sample. Otherwise, if a file with multiple sample prediction is input
it will not match the number of observations of the AnnData object and the process will throw an error.

``q`` can be set to use a different prediction matrix from the *cell2location* output file.
It defaults to ``"q05_cell_abundance_w_sf"``

``sort`` and ``sort_index`` can be used to define how a *cell2location* output file matches the AnnData object.
By default the pipeline will try to ensure the order of observations between the prediction matrix and AnnData object match
so values are correctly concatenated.
The pipeline will attempt to order the prediction matrix given the index of the AnnData object 
(or the original index if the main pipeline re-indexed it).
However you can override the observations column of the AnnData object that contains the index that the prediction matrix should match using ``sort_index``.
``sort`` can be set to ``False`` to disable any re-ordering. If disabled, the prediction matrix would be concatenated as-is into the AnnData object
without checking if observations' IDs match.
