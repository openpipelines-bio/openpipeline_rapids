nextflow.enable.dsl=2

params.rootDir = params.rootDir ?: projectDir + "/../../../.."

include { neighbors_leiden_umap } from params.rootDir + "/target/nextflow/workflows/multiomics/neighbors_leiden_umap/main.nf"

params.resources_test = params.resources_test ?: params.rootDir + "/resources_test/pbmc_1k_protein_v3"

workflow test_wf {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "simple_execution_test",
        input: resources_test.resolve("pbmc_1k_protein_v3_mms.h5mu"),
        obsm_input: "X_pca",
        uns_neighbors: "neighbors",
        obsp_neighbor_distances: "distances",
        obsp_neighbor_connectivities: "connectivities",
        obs_cluster: "leiden",
        leiden_resolution: [1.0],
        obsm_umap: "X_umap"
      ]
    ])
    | map { state -> [state.id, state] }
    | neighbors_leiden_umap
    | view { output ->
      assert output.size() == 2 : "Outputs should contain two elements; [id, state]"

      // check id
      def id = output[0]
      assert id == "simple_execution_test"

      // check output
      def state = output[1]
      assert state instanceof Map : "State should be a map. Found: ${state}"
      assert state.containsKey("output") : "Output should contain key 'output'."
      assert state.output.isFile() : "'output' should be a file."
      assert state.output.toString().endsWith(".h5mu") : "Output file should end with '.h5mu'. Found: ${state.output}"

      "Output: $output"
    }
    | toSortedList({ a, b -> a[0] <=> b[0] })
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain 1 event"
      assert output_list.collect{ it[0] } == ["simple_execution_test"]
    }
}
