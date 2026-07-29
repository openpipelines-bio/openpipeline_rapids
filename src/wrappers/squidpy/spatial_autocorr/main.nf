workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // -- GPU variant (rapids-singlecell) --
    | spatial_autocorr_gpu.run(
      runIf: { id, state -> state.device_type == "gpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "layer": "layer",
        "obsp_connectivities": "obsp_connectivities",
        "input_genes": "input_genes",
        "output_compression": "output_compression",
        "mode": "mode",
        "n_perms": "n_perms",
        "transformation": "transformation",
        "two_tailed": "two_tailed",
        "corr_method": "corr_method",
        "use_sparse": "use_sparse",
        "output": "output"
      ],
      toState: ["output": "output"]
    )
    // -- CPU variant (squidpy, openpipeline_spatial) --
    | spatial_autocorr_cpu.run(
      runIf: { id, state -> state.device_type == "cpu" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "layer": "layer",
        "obsp_neighborhood_graph": "obsp_connectivities",
        "input_genes": "input_genes",
        "mode": "mode",
        "n_perms": "n_perms",
        "use_all_genes": "use_all_genes",
        "attr": "attr",
        "use_raw": "use_raw",
        "output": "output"
      ],
      toState: ["output": "output"]
    )
    | setState(["output"])

  emit:
  output_ch
}
