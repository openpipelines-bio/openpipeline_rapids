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
    | scale_gpu.run(
      runIf: { id, state -> state.device_type == "gpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "input_layer": "input_layer",
        "output": "workflow_output",
        "output_layer": "output_layer",
        "output_compression": "output_compression",
        "zero_center": "zero_center",
        "max_value": "max_value"
      ],
      toState: ["input": "output"]
    )
    // -- CPU variant (openpipeline) --
    // The CPU scale exposes --zero_center as a boolean_false: its PRESENCE
    // disables centering. The workflow arg --zero_center is a boolean where
    // true means "do center", so only set the CPU flag when it is false.
    | scale_cpu.run(
      runIf: { id, state -> state.device_type == "cpu" },
      fromState: { id, state ->
        def m = [
          "input": state.input,
          "modality": state.modality,
          "input_layer": state.input_layer,
          "output": state.workflow_output,
          "output_layer": state.output_layer,
          "output_compression": state.output_compression,
          "max_value": state.max_value
        ]
        if (state.zero_center == false) { m.zero_center = true }
        m
      },
      toState: ["input": "output"]
    )
    | setState(["output": "input"])

  emit:
  output_ch
}
