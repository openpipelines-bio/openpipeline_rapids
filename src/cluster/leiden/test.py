import sys
import pytest
import mudata as mu

## VIASH START
meta = {
    "name": "leiden",
    "resources_dir": "resources_test/",
    "config": "src/cluster/leiden/config.vsh.yaml",
    "executable": "target/docker/cluster/leiden/leiden",
}
## VIASH END

input = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"


def test_run(run_component, random_h5mu_path):
    output = random_h5mu_path()
    cmd_pars = [
        "--input",
        input,
        "--output",
        output,
        "--output_compression",
        "gzip",
        "--resolution",
        "1.0",
    ]
    run_component(cmd_pars)

    assert output.is_file(), "No output was created."

    mu_output = mu.read_h5mu(output)
    assert "rna" in mu_output.mod, 'Output should contain data.mod["rna"].'

    rna_out = mu_output.mod["rna"]
    assert "leiden" in rna_out.obsm, "Output should contain .obsm['leiden']."
    leiden = rna_out.obsm["leiden"]
    assert "1.0" in leiden.columns, (
        "Default run should contain a column for resolution 1.0."
    )
    assert leiden.shape[0] == rna_out.n_obs, "Output should have one row per cell."


def test_resolution(run_component, random_h5mu_path):
    """Multiple resolutions should produce a column per resolution and
    different resolutions should yield different cluster assignments."""
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            input,
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--resolution",
            "1.0;0.25",
        ]
    )

    rna_out = mu.read_h5mu(output).mod["rna"]
    leiden = rna_out.obsm["leiden"]
    assert "1.0" in leiden.columns, "Output should contain resolution 1.0."
    assert "0.25" in leiden.columns, "Output should contain resolution 0.25."
    assert not leiden["1.0"].equals(leiden["0.25"]), (
        "Different resolutions should produce different cluster assignments."
    )


def test_obsm_name(run_component, random_h5mu_path):
    """A custom --obsm_name should be used as the key in .obsm."""
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            input,
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--resolution",
            "1.0",
            "--obsm_name",
            "my_leiden",
        ]
    )

    rna_out = mu.read_h5mu(output).mod["rna"]
    assert "my_leiden" in rna_out.obsm, "Output should contain .obsm['my_leiden']."
    assert "leiden" not in rna_out.obsm, (
        "Default key 'leiden' should not be present when --obsm_name is set."
    )


def test_obsp_connectivities(run_component, random_h5mu_path):
    """A custom --obsp_connectivities should be picked up as the adjacency."""
    mu_orig = mu.read_h5mu(input)
    mu_orig.mod["rna"].obsp["custom_connectivities"] = (
        mu_orig.mod["rna"].obsp["connectivities"].copy()
    )
    del mu_orig.mod["rna"].obsp["connectivities"]
    input_custom = random_h5mu_path()
    mu_orig.write(input_custom)

    output = random_h5mu_path()
    run_component(
        [
            "--input",
            str(input_custom),
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--resolution",
            "1.0",
            "--obsp_connectivities",
            "custom_connectivities",
        ]
    )

    rna_out = mu.read_h5mu(output).mod["rna"]
    assert "leiden" in rna_out.obsm, (
        "Output should contain .obsm['leiden'] when using a custom adjacency."
    )


def test_seed(run_component, random_h5mu_path):
    """Same seed should produce identical cluster assignments across runs."""
    output_a = random_h5mu_path()
    output_b = random_h5mu_path()

    for output in (output_a, output_b):
        run_component(
            [
                "--input",
                input,
                "--output",
                output,
                "--output_compression",
                "gzip",
                "--resolution",
                "1.0",
                "--seed",
                "42",
            ]
        )

    leiden_a = mu.read_h5mu(output_a).mod["rna"].obsm["leiden"]
    leiden_b = mu.read_h5mu(output_b).mod["rna"].obsm["leiden"]
    assert leiden_a.equals(leiden_b), (
        "Identical --seed should produce identical cluster assignments."
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
