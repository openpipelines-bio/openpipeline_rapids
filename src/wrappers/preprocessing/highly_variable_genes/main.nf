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
    | highly_variable_genes_gpu.run(
      runIf: { id, state -> state.device_type == "gpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "input_layer": "input_layer",
        "output": "workflow_output",
        "output_compression": "output_compression",
        "var_name_filter": "var_name_filter",
        "varm_name": "varm_name",
        "flavor": "flavor",
        "n_top_features": "n_top_features",
        "min_mean": "min_mean",
        "max_mean": "max_mean",
        "min_disp": "min_disp",
        "max_disp": "max_disp",
        "span": "span",
        "n_bins": "n_bins",
        "theta": "theta",
        "clip": "clip",
        "chunksize": "chunksize",
        "n_samples": "n_samples",
        "obs_batch_key": "obs_batch_key",
        "check_values": "check_values"
      ],
      toState: ["input": "output"]
    )
    // -- CPU variant (openpipeline) --
    | highly_variable_genes_cpu.run(
      runIf: { id, state -> state.device_type == "cpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "layer": "input_layer",
        "output": "workflow_output",
        "output_compression": "output_compression",
        "var_name_filter": "var_name_filter",
        "varm_name": "varm_name",
        "flavor": "flavor",
        "n_top_features": "n_top_features",
        "min_mean": "min_mean",
        "max_mean": "max_mean",
        "min_disp": "min_disp",
        "max_disp": "max_disp",
        "span": "span",
        "n_bins": "n_bins",
        "obs_batch_key": "obs_batch_key",
        "var_input": "var_input",
        "features_to_exclude": "features_to_exclude"
      ],
      toState: ["input": "output"]
    )
    | setState(["output": "input"])

  emit:
  output_ch
}
