workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // -- GPU variant (rapids-singlecell) --
    | regress_out_gpu.run(
      runIf: { id, state -> state.device_type == "gpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "input_layer": "input_layer",
        "obs_keys": "obs_keys",
        "output_layer": "output_layer",
        "output_compression": "output_compression",
        "batchsize": "batchsize"
      ],
      toState: ["output": "output"]
    )
    // -- CPU variant (openpipeline) --
    | regress_out_cpu.run(
      runIf: { id, state -> state.device_type == "cpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "input_layer": "input_layer",
        "obs_keys": "obs_keys",
        "output_layer": "output_layer",
        "output_compression": "output_compression"
      ],
      toState: ["output": "output"]
    )
    | setState(["output"])

  emit:
  output_ch
}
