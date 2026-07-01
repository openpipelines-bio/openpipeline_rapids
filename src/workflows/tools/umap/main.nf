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
    | umap_gpu.run(
      runIf: { id, state -> state.device == "gpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "uns_neighbors": "uns_neighbors",
        "output": "workflow_output",
        "obsm_output": "obsm_output",
        "output_compression": "output_compression",
        "min_dist": "min_dist",
        "spread": "spread",
        "num_components": "num_components",
        "max_iter": "max_iter",
        "alpha": "alpha",
        "negative_sample_rate": "negative_sample_rate",
        "init_pos": "init_pos",
        "random_state": "random_state"
      ],
      toState: ["input": "output"]
    )
    // -- CPU variant (openpipeline) --
    | umap_cpu.run(
      runIf: { id, state -> state.device == "cpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "uns_neighbors": "uns_neighbors",
        "output": "workflow_output",
        "obsm_output": "obsm_output",
        "output_compression": "output_compression",
        "min_dist": "min_dist",
        "spread": "spread",
        "num_components": "num_components",
        "max_iter": "max_iter",
        "alpha": "alpha",
        "negative_sample_rate": "negative_sample_rate",
        "init_pos": "init_pos",
        "gamma": "gamma"
      ],
      toState: ["input": "output"]
    )
    | setState(["output": "input"])

  emit:
  output_ch
}
