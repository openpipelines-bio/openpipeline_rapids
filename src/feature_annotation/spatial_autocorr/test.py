import sys
import pytest
import mudata as mu
import pandas as pd
import scanpy as sc

## VIASH START
meta = {
    "name": "spatial_autocorr",
    "resources_dir": "resources_test/",
    "config": "src/feature_annotation/spatial_autocorr/config.vsh.yaml",
    "executable": "target/docker/feature_annotation/spatial_autocorr/spatial_autocorr",
}
## VIASH END

input = f"{meta['resources_dir']}/xenium/xenium_tiny.qc.neighbors.h5mu"


@pytest.fixture
def filtered_input(random_h5mu_path):
    """rsc.gr.spatial_autocorr's default float32 reduction yields nan/inf
    on constant or low-variance genes, which the xenium test panel has
    several of. Drop genes detected in fewer than 3 cells before passing
    to the component — matches typical pre-spatial-autocorr filtering."""
    mu_in = mu.read_h5mu(input)
    rna = mu_in.mod["rna"]
    sc.pp.filter_genes(rna, min_cells=3)
    path = random_h5mu_path()
    mu_in.write(path)
    return path


def test_run(run_component, random_h5mu_path, filtered_input):
    output = random_h5mu_path()
    cmd_pars = [
        "--input",
        filtered_input,
        "--output",
        output,
        "--output_compression",
        "gzip",
        "--mode",
        "moran",
        "--n_perms",
        "10",
    ]
    run_component(cmd_pars)

    assert output.is_file(), "No output was created."

    mu_output = mu.read_h5mu(output)
    assert "rna" in mu_output.mod, 'Output should contain data.mod["rna"].'
    rna_out = mu_output.mod["rna"]

    assert "moranI" in rna_out.uns, "moranI key missing from .uns."
    df = rna_out.uns["moranI"]
    assert isinstance(df, pd.DataFrame), "moranI should be a DataFrame."
    assert not df.empty, "moranI DataFrame is empty."
    assert "I" in df.columns, "Expected statistic column 'I' in moranI."
    # Moran's I lies in roughly [-1, 1].
    assert df["I"].max() <= 1.0, "Moran's I should not exceed 1."
    assert df["I"].min() >= -1.0, "Moran's I should not fall below -1."


def test_geary(run_component, random_h5mu_path, filtered_input):
    """Geary's C mode should write a `gearyC` DataFrame to .uns."""
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            filtered_input,
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--mode",
            "geary",
            "--n_perms",
            "10",
        ]
    )

    rna_out = mu.read_h5mu(output).mod["rna"]
    assert "gearyC" in rna_out.uns, "gearyC key missing from .uns."
    df = rna_out.uns["gearyC"]
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert "C" in df.columns, "Expected statistic column 'C' in gearyC."


def test_genes_subset(run_component, random_h5mu_path, filtered_input):
    """Passing --genes should restrict the analysis to that subset."""
    mu_input = mu.read_h5mu(filtered_input)
    genes = list(mu_input.mod["rna"].var_names[:5])

    output = random_h5mu_path()
    run_component(
        [
            "--input",
            filtered_input,
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--mode",
            "moran",
            "--n_perms",
            "10",
            "--genes",
            ",".join(genes),
        ]
    )

    df = mu.read_h5mu(output).mod["rna"].uns["moranI"]
    assert len(df) == len(genes), (
        f"Expected results for {len(genes)} genes, got {len(df)}."
    )
    for g in genes:
        assert g in df.index, f"Gene {g} missing from output index."


def test_obsp_connectivities(run_component, random_h5mu_path, filtered_input):
    """A custom --obsp_connectivities should be picked up as the adjacency."""
    mu_orig = mu.read_h5mu(filtered_input)
    mu_orig.mod["rna"].obsp["custom_connectivities"] = (
        mu_orig.mod["rna"].obsp["spatial_connectivities"].copy()
    )
    del mu_orig.mod["rna"].obsp["spatial_connectivities"]
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
            "--mode",
            "moran",
            "--n_perms",
            "10",
            "--obsp_connectivities",
            "custom_connectivities",
        ]
    )

    rna_out = mu.read_h5mu(output).mod["rna"]
    assert "moranI" in rna_out.uns, (
        "Output should contain .uns['moranI'] when using a custom adjacency."
    )


def test_corr_method_disabled(run_component, random_h5mu_path, filtered_input):
    """Passing an empty --corr_method should leave p-values uncorrected.

    The exact column naming for corrected p-values depends on the
    rapids-singlecell version, so we only assert that the run succeeds and
    that no `*_fdr_bh` column is present in the output.
    """
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            filtered_input,
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--mode",
            "moran",
            "--n_perms",
            "10",
            "--corr_method",
            "",
        ]
    )

    df = mu.read_h5mu(output).mod["rna"].uns["moranI"]
    assert not any("fdr_bh" in c for c in df.columns), (
        "Disabling --corr_method should not produce FDR-corrected columns."
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
