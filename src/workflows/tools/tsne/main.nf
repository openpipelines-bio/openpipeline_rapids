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
    | tsne_gpu.run(
      runIf: { id, state -> state.device == "gpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "obsm_input": "obsm_input",
        "n_pcs": "n_pcs",
        "output": "workflow_output",
        "obsm_output": "obsm_output",
        "output_compression": "output_compression",
        "overwrite": "overwrite",
        "perplexity": "perplexity",
        "early_exaggeration": "early_exaggeration",
        "learning_rate": "learning_rate",
        "metric": "metric",
        "method": "method"
      ],
      toState: ["input": "output"]
    )
    // -- CPU variant (openpipeline) --
    | tsne_cpu.run(
      runIf: { id, state -> state.device == "cpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "use_rep": "obsm_input",
        "n_pcs": "n_pcs",
        "output": "workflow_output",
        "obsm_output": "obsm_output",
        "output_compression": "output_compression",
        "perplexity": "perplexity",
        "early_exaggeration": "early_exaggeration",
        "learning_rate": "learning_rate",
        "metric": "metric",
        "min_dist": "min_dist",
        "random_state": "random_state"
      ],
      toState: ["input": "output"]
    )
    | setState(["output": "input"])

  emit:
  output_ch
}
