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
    | neighbors_gpu.run(
      runIf: { id, state -> state.device_type == "gpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "obsm_input": "obsm_input",
        "output": "workflow_output",
        "output_compression": "output_compression",
        "uns_output": "uns_output",
        "obsp_distances": "obsp_distances",
        "obsp_connectivities": "obsp_connectivities",
        "num_neighbors": "num_neighbors",
        "metric": "metric",
        "n_pcs": "n_pcs",
        "algorithm": "algorithm",
        "method": "method",
        "random_state": "random_state"
      ],
      toState: ["input": "output"]
    )
    // -- CPU variant (openpipeline) --
    | neighbors_cpu.run(
      runIf: { id, state -> state.device_type == "cpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "obsm_input": "obsm_input",
        "output": "workflow_output",
        "output_compression": "output_compression",
        "uns_output": "uns_output",
        "obsp_distances": "obsp_distances",
        "obsp_connectivities": "obsp_connectivities",
        "num_neighbors": "num_neighbors",
        "metric": "metric",
        "seed": "random_state"
      ],
      toState: ["input": "output"]
    )
    | setState(["output": "input"])

  emit:
  output_ch
}
