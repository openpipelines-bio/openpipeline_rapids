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
    | leiden_gpu.run(
      runIf: { id, state -> state.device_type == "gpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "obsp_connectivities": "obsp_connectivities",
        "output": "workflow_output",
        "obsm_name": "obsm_name",
        "output_compression": "output_compression",
        "resolution": "resolution",
        "n_iterations": "n_iterations",
        "random_state": "random_state",
        "theta": "theta",
        "use_weights": "use_weights"
      ],
      toState: ["input": "output"]
    )
    // -- CPU variant (openpipeline) --
    | leiden_cpu.run(
      runIf: { id, state -> state.device_type == "cpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "obsp_connectivities": "obsp_connectivities",
        "output": "workflow_output",
        "obsm_name": "obsm_name",
        "output_compression": "output_compression",
        "resolution": "resolution",
        "n_iterations": "n_iterations",
        "seed": "random_state",
        "flavor": "flavor"
      ],
      toState: ["input": "output"]
    )
    | setState(["output": "input"])

  emit:
  output_ch
}
