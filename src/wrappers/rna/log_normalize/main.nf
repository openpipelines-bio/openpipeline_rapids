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
    | log_normalize_gpu.run(
      runIf: { id, state -> state.device_type == "gpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "layer": "layer",
        "target_sum": "target_sum",
        "output_layer": "output_layer",
        "output": "workflow_output"
      ],
      toState: ["input": "output"]
    )
    // -- CPU variant (scanpy) --
    | log_normalize_cpu.run(
      runIf: { id, state -> state.device_type == "cpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "layer": "layer",
        "target_sum": "target_sum",
        "output_layer": "output_layer",
        "output": "workflow_output"
      ],
      toState: ["input": "output"]
    )
    | setState(["output": "input"])

  emit:
  output_ch
}
