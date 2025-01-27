# Map entity types used in WebAtlas2.0 UI to the o.var.feature_types values in -anndata.zarr files
# TODO: Check that feature type is stored in the same way in o.var.feature_type for all projects (for: GBM, ASD)
entity_type2_feature_type = {
    "gene" : "gene_expression",
    "cell_type_lvl7" : "cell_type_lvl7",
    "cell_type_lvl3" : "cell_type_lvl3",
    "annotation" : "annotations"
}
spatialxy_colnames_alternatives = ['spatial', 'X_spatial']
feature_type_colnames_alternatives = ['feature_type', 'feature_types']