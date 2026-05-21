import sys
import pytest
import mudata as mu
import numpy as np
import scanpy as sc
from openpipeline_testutils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    "name": "find_neighbors",
    "resources_dir": "resources_test/",
    "config": "src/neighbors/find_neighbors/config.vsh.yaml",
    "executable": "target/docker/neighbors/find_neighbors/find_neighbors",
}
## VIASH END

input = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"


def _prepare_input_with_pca(random_h5mu_path):
    """Read the raw input, compute a PCA on log-normalized counts so that
    .obsm['X_pca'] exists, and write it back to a fresh h5mu path."""
    mdata = mu.read_h5mu(input)
    adata = mdata.mod["rna"]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, n_comps=20)
    path = random_h5mu_path()
    mdata.write(path)
    return path


def test_run(run_component, random_h5mu_path):
    input_with_pca = _prepare_input_with_pca(random_h5mu_path)
    output = random_h5mu_path()
    cmd_pars = [
        "--input",
        str(input_with_pca),
        "--output",
        output,
        "--output_compression",
        "gzip",
    ]
    run_component(cmd_pars)

    assert output.is_file(), "No output was created."

    mu_input = mu.read_h5mu(input_with_pca)
    mu_output = mu.read_h5mu(output)

    assert "rna" in mu_output.mod, 'Output should contain data.mod["rna"].'

    rna_in = mu_input.mod["rna"]
    rna_out = mu_output.mod["rna"]

    assert rna_in.shape == rna_out.shape, "Should have same shape as before"

    # Default output keys
    assert "neighbors" in rna_out.uns, "Output should have .uns['neighbors']"
    assert "distances" in rna_out.obsp, "Output should have .obsp['distances']"
    assert "connectivities" in rna_out.obsp, (
        "Output should have .obsp['connectivities']"
    )

    # Input should not have those neighbor outputs
    assert "neighbors" not in rna_in.uns
    assert "distances" not in rna_in.obsp
    assert "connectivities" not in rna_in.obsp

    # Each row of the KNN distance matrix should have num_neighbors-1 non-zero
    # entries by default (15 neighbors, self distance is zero).
    distances = rna_out.obsp["distances"]
    nnz_per_row = (distances != 0).sum(axis=1)
    assert int(np.max(nnz_per_row)) <= 15

    # Copy over the new slots so the rest of the object can be compared
    rna_in.uns["neighbors"] = rna_out.uns["neighbors"]
    rna_in.obsp["distances"] = rna_out.obsp["distances"]
    rna_in.obsp["connectivities"] = rna_out.obsp["connectivities"]
    assert_annotation_objects_equal(mu_input, mu_output)


def test_uns_output(run_component, random_h5mu_path):
    """When --uns_output is set to a custom key, neighbor metadata is stored
    under that key and obsp keys are namespaced with that prefix."""
    input_with_pca = _prepare_input_with_pca(random_h5mu_path)
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            str(input_with_pca),
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--uns_output",
            "foo_neigh",
        ]
    )

    rna_out = mu.read_h5mu(output).mod["rna"]

    assert "foo_neigh" in rna_out.uns, "Output should have .uns['foo_neigh']"
    assert "foo_neigh_distances" in rna_out.obsp, (
        "Output should have .obsp['foo_neigh_distances']"
    )
    assert "foo_neigh_connectivities" in rna_out.obsp, (
        "Output should have .obsp['foo_neigh_connectivities']"
    )
    assert "neighbors" not in rna_out.uns
    assert "distances" not in rna_out.obsp
    assert "connectivities" not in rna_out.obsp


def test_num_neighbors(run_component, random_h5mu_path):
    """A smaller num_neighbors should produce a sparser KNN graph."""
    input_with_pca = _prepare_input_with_pca(random_h5mu_path)
    output_default = random_h5mu_path()
    output_small = random_h5mu_path()

    run_component(
        [
            "--input",
            str(input_with_pca),
            "--output",
            output_default,
            "--output_compression",
            "gzip",
        ]
    )
    run_component(
        [
            "--input",
            str(input_with_pca),
            "--output",
            output_small,
            "--output_compression",
            "gzip",
            "--num_neighbors",
            "5",
        ]
    )

    rna_default = mu.read_h5mu(output_default).mod["rna"]
    rna_small = mu.read_h5mu(output_small).mod["rna"]

    nnz_default = rna_default.obsp["distances"].nnz
    nnz_small = rna_small.obsp["distances"].nnz
    assert nnz_small < nnz_default, (
        "Smaller num_neighbors should produce a sparser distance matrix"
    )
    assert rna_small.uns["neighbors"]["params"]["n_neighbors"] == 5


def test_obsm_input(run_component, random_h5mu_path):
    """Passing a non-default obsm slot should be used as the neighbor input."""
    mdata = mu.read_h5mu(input)
    adata = mdata.mod["rna"]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, n_comps=20)
    # Move the PCA into a custom obsm key
    adata.obsm["X_custom"] = adata.obsm.pop("X_pca")
    input_with_custom = random_h5mu_path()
    mdata.write(input_with_custom)

    output = random_h5mu_path()
    run_component(
        [
            "--input",
            str(input_with_custom),
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--obsm_input",
            "X_custom",
        ]
    )

    rna_out = mu.read_h5mu(output).mod["rna"]
    assert rna_out.uns["neighbors"]["params"]["use_rep"] == "X_custom"
    assert "distances" in rna_out.obsp
    assert "connectivities" in rna_out.obsp


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
