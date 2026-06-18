workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // Preserve the requested final output filename across all steps.
    | map { id, state ->
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
    // -- step 1: log_normalize (normalize_total + log1p) --
    | log_normalize.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "output": "workflow_output"
      ],
      args: ["output_layer": "log_normalized"],
      toState: ["input": "output"]
    )
    // -- step 2: highly_variable_genes --
    | highly_variable_genes.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "output": "workflow_output"
      ],
      args: ["input_layer": "log_normalized"],
      toState: ["input": "output"]
    )
    // -- step 3: pca --
    | pca.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "output": "workflow_output"
      ],
      args: [
        "layer": "log_normalized",
        "obsm_output": "X_pca"
      ],
      toState: ["input": "output"]
    )
    // -- step 4: neighbors_leiden_umap (neighbors + leiden + umap) --
    | neighbors_leiden_umap.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "output": "workflow_output"
      ],
      args: [
        "obsm_input": "X_pca",
        "uns_neighbors": "neighbors",
        "obsp_neighbor_distances": "distances",
        "obsp_neighbor_connectivities": "connectivities",
        "obs_cluster": "leiden",
        "leiden_resolution": [1.0],
        "obsm_umap": "X_umap"
      ],
      toState: ["input": "output"]
    )
    // -- step 5: spatial_autocorr (Moran's I on the spatial graph) --
    | spatial_autocorr.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "output": "workflow_output"
      ],
      args: [
        "obsp_connectivities": "spatial_connectivities",
        "mode": "moran"
      ],
      toState: ["input": "output"]
    )
    | setState(["output": "input"])

  emit:
  output_ch
}
