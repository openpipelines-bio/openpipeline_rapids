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
    | pca_gpu.run(
      runIf: { id, state -> state.device == "gpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "layer": "layer",
        "var_input": "var_input",
        "num_components": "num_components",
        "chunked": "chunked",
        "chunk_size": "chunk_size",
        "random_state": "random_state",
        "output": "workflow_output",
        "output_compression": "output_compression",
        "obsm_output": "obsm_output",
        "varm_output": "varm_output",
        "uns_output": "uns_output",
        "overwrite": "overwrite"
      ],
      toState: ["input": "output"]
    )
    // -- CPU variant (openpipeline) --
    | pca_cpu.run(
      runIf: { id, state -> state.device == "cpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "layer": "layer",
        "var_input": "var_input",
        "num_components": "num_components",
        "chunked": "chunked",
        "chunk_size": "chunk_size",
        "seed": "random_state",
        "output": "workflow_output",
        "output_compression": "output_compression",
        "obsm_output": "obsm_output",
        "varm_output": "varm_output",
        "uns_output": "uns_output",
        "overwrite": "overwrite"
      ],
      toState: ["input": "output"]
    )
    | setState(["output": "input"])

  emit:
  output_ch
}
