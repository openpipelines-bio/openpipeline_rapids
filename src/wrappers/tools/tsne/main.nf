workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // -- GPU variant (rapids-singlecell) --
    | tsne_gpu.run(
      runIf: { id, state -> state.device_type == "gpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "obsm_input": "obsm_input",
        "n_pcs": "n_pcs",
        "obsm_output": "obsm_output",
        "output_compression": "output_compression",
        "overwrite": "overwrite",
        "perplexity": "perplexity",
        "early_exaggeration": "early_exaggeration",
        "learning_rate": "learning_rate",
        "metric": "metric",
        "method": "method"
      ],
      toState: ["output": "output"]
    )
    // -- CPU variant (openpipeline) --
    | tsne_cpu.run(
      runIf: { id, state -> state.device_type == "cpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "use_rep": "obsm_input",
        "n_pcs": "n_pcs",
        "obsm_output": "obsm_output",
        "output_compression": "output_compression",
        "perplexity": "perplexity",
        "early_exaggeration": "early_exaggeration",
        "learning_rate": "learning_rate",
        "metric": "metric",
        "min_dist": "min_dist",
        "random_state": "random_state"
      ],
      toState: ["output": "output"]
    )
    | setState(["output"])

  emit:
  output_ch
}
