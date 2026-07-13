workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // Preserve the requested final output filename across the step.
    | map { id, state ->
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
    // -- GPU variant (rapids-singlecell) --
    | neighbors_leiden_umap_gpu.run(
      runIf: { id, state -> state.device_type == "gpu" },
      fromState: [
        "input": "input",
        "obsm_input": "obsm_input",
        "modality": "modality",
        "output": "workflow_output",
        "uns_neighbors": "uns_neighbors",
        "obsp_neighbor_distances": "obsp_neighbor_distances",
        "obsp_neighbor_connectivities": "obsp_neighbor_connectivities",
        "obs_cluster": "obs_cluster",
        "leiden_resolution": "leiden_resolution",
        "obsm_umap": "obsm_umap"
      ],
      toState: ["input": "output"]
    )
    // -- CPU variant (scanpy) --
    | neighbors_leiden_umap_cpu.run(
      runIf: { id, state -> state.device_type == "cpu" },
      fromState: [
        "input": "input",
        "obsm_input": "obsm_input",
        "modality": "modality",
        "output": "workflow_output",
        "uns_neighbors": "uns_neighbors",
        "obsp_neighbor_distances": "obsp_neighbor_distances",
        "obsp_neighbor_connectivities": "obsp_neighbor_connectivities",
        "obs_cluster": "obs_cluster",
        "leiden_resolution": "leiden_resolution",
        "obsm_umap": "obsm_umap"
      ],
      toState: ["input": "output"]
    )
    | setState(["output": "input"])

  emit:
  output_ch
}
