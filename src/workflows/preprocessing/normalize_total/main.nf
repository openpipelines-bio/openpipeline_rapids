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
    | normalize_total_gpu.run(
      runIf: { id, state -> state.device == "gpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "input_layer": "input_layer",
        "output": "workflow_output",
        "output_layer": "output_layer",
        "output_compression": "output_compression",
        "target_sum": "target_sum",
        "exclude_highly_expressed": "exclude_highly_expressed",
        "max_fraction": "max_fraction"
      ],
      toState: ["input": "output"]
    )
    // -- CPU variant (openpipeline) --
    | normalize_total_cpu.run(
      runIf: { id, state -> state.device == "cpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "input_layer": "input_layer",
        "output": "workflow_output",
        "output_layer": "output_layer",
        "output_compression": "output_compression",
        "target_sum": "target_sum",
        "exclude_highly_expressed": "exclude_highly_expressed"
      ],
      toState: ["input": "output"]
    )
    | setState(["output": "input"])

  emit:
  output_ch
}
