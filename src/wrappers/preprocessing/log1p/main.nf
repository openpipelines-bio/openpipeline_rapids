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
    | log1p_gpu.run(
      runIf: { id, state -> state.device_type == "gpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "input_layer": "input_layer",
        "input_obsm": "input_obsm",
        "output": "workflow_output",
        "output_layer": "output_layer",
        "output_obsm": "output_obsm",
        "output_compression": "output_compression",
        "base": "base"
      ],
      toState: ["input": "output"]
    )
    // -- CPU variant (openpipeline) --
    | log1p_cpu.run(
      runIf: { id, state -> state.device_type == "cpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "input_layer": "input_layer",
        "output": "workflow_output",
        "output_layer": "output_layer",
        "output_compression": "output_compression",
        "base": "base"
      ],
      toState: ["input": "output"]
    )
    | setState(["output": "input"])

  emit:
  output_ch
}
