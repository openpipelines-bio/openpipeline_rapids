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
    | bbknn_gpu.run(
      runIf: { id, state -> state.device_type == "gpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "obsm_input": "obsm_input",
        "batch_key": "batch_key",
        "output": "workflow_output",
        "output_compression": "output_compression",
        "uns_output": "uns_output",
        "obsp_distances": "obsp_distances",
        "obsp_connectivities": "obsp_connectivities",
        "overwrite": "overwrite",
        "neighbors_within_batch": "neighbors_within_batch",
        "n_pcs": "n_pcs",
        "metric": "metric",
        "algorithm": "algorithm",
        "trim": "trim",
        "random_state": "random_state"
      ],
      toState: ["input": "output"]
    )
    // -- CPU variant (openpipeline) --
    | bbknn_cpu.run(
      runIf: { id, state -> state.device_type == "cpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "obsm_input": "obsm_input",
        "obs_batch": "batch_key",
        "output": "workflow_output",
        "output_compression": "output_compression",
        "uns_output": "uns_output",
        "obsp_distances": "obsp_distances",
        "obsp_connectivities": "obsp_connectivities",
        "n_neighbors_within_batch": "neighbors_within_batch",
        "n_pcs": "n_pcs",
        "n_trim": "trim"
      ],
      toState: ["input": "output"]
    )
    | setState(["output": "input"])

  emit:
  output_ch
}
