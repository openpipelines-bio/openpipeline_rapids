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
    | harmony_integrate_gpu.run(
      runIf: { id, state -> state.device == "gpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "obsm_input": "obsm_input",
        "obs_covariates": "obs_covariates",
        "output": "workflow_output",
        "obsm_output": "obsm_output",
        "output_compression": "output_compression",
        "overwrite": "overwrite",
        "theta": "theta",
        "flavor": "flavor",
        "n_clusters": "n_clusters",
        "max_iter_harmony": "max_iter_harmony",
        "random_state": "random_state"
      ],
      toState: ["input": "output"]
    )
    // -- CPU variant (openpipeline) --
    | harmony_integrate_cpu.run(
      runIf: { id, state -> state.device == "cpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "obsm_input": "obsm_input",
        "obs_covariates": "obs_covariates",
        "output": "workflow_output",
        "obsm_output": "obsm_output",
        "output_compression": "output_compression",
        "theta": "theta"
      ],
      toState: ["input": "output"]
    )
    | setState(["output": "input"])

  emit:
  output_ch
}
