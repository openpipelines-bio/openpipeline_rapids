workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // Normalize raw counts straight into the final output layer, leaving the
    // source matrix (.X or --layer) untouched.
    | normalize_total.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "input_layer": "layer",
        "target_sum": "target_sum",
        "output_layer": "output_layer"
      ],
      toState: [
        "input": "output"
      ]
    )
    // Log-transform the normalized counts in place in the output layer.
    | log1p.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "input_layer": "output_layer",
        "output_compression": "output_compression"
      ],
      toState: [
        "output": "output"
      ]
    )
    | setState(["output"])

  emit:
  output_ch
}
