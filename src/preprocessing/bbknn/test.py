import sys
import pytest
import subprocess
import mudata as mu
import pandas as pd

## VIASH START
meta = {
    "name": "bbknn",
    "resources_dir": "resources_test/",
    "config": "src/preprocessing/bbknn/config.vsh.yaml",
    "executable": "target/docker/preprocessing/bbknn/bbknn",
}
## VIASH END

input = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"


def _prepare_input_with_batch(random_h5mu_path):
    """Read the input (which already has .obsm['X_pca']) and add a categorical
    .obs['batch'] column alternating between 'a' and 'b', then write it back to
    a fresh h5mu path. The stale neighbor-graph slots carried over from the
    upstream pipeline are dropped so bbknn can write its default output slots
    without --overwrite."""
    mdata = mu.read_h5mu(input)
    adata = mdata.mod["rna"]
    labels = ["a", "b"] * (adata.n_obs // 2 + 1)
    adata.obs["batch"] = pd.Categorical(labels[: adata.n_obs])
    adata.uns.pop("neighbors", None)
    adata.obsp.pop("distances", None)
    adata.obsp.pop("connectivities", None)
    path = random_h5mu_path()
    mdata.write(path)
    return path


def test_run(run_component, random_h5mu_path):
    input_with_batch = _prepare_input_with_batch(random_h5mu_path)
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            str(input_with_batch),
            "--output",
            output,
            "--batch_key",
            "batch",
            "--output_compression",
            "gzip",
        ]
    )

    assert output.is_file(), "No output was created."

    mu_output = mu.read_h5mu(output)

    assert "rna" in mu_output.mod, 'Output should contain data.mod["rna"].'

    rna_out = mu_output.mod["rna"]

    # Default output keys
    assert "neighbors" in rna_out.uns, "Output should have .uns['neighbors']"
    assert "distances" in rna_out.obsp, "Output should have .obsp['distances']"
    assert "connectivities" in rna_out.obsp, (
        "Output should have .obsp['connectivities']"
    )

    distances = rna_out.obsp["distances"]
    connectivities = rna_out.obsp["connectivities"]
    n_obs = rna_out.n_obs
    assert distances.shape == (n_obs, n_obs)
    assert connectivities.shape == (n_obs, n_obs)


def test_custom_output_slots(run_component, random_h5mu_path):
    """Custom --uns_output / --obsp_distances / --obsp_connectivities slots are
    stored under exactly the requested keys and the default keys are absent."""
    input_with_batch = _prepare_input_with_batch(random_h5mu_path)
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            str(input_with_batch),
            "--output",
            output,
            "--batch_key",
            "batch",
            "--output_compression",
            "gzip",
            "--uns_output",
            "foo_neigh",
            "--obsp_distances",
            "foo_dist",
            "--obsp_connectivities",
            "foo_conn",
        ]
    )

    rna_out = mu.read_h5mu(output).mod["rna"]

    assert "foo_neigh" in rna_out.uns, "Output should have .uns['foo_neigh']"
    assert "foo_dist" in rna_out.obsp, "Output should have .obsp['foo_dist']"
    assert "foo_conn" in rna_out.obsp, "Output should have .obsp['foo_conn']"
    # The metadata should point at the custom obsp slots
    assert rna_out.uns["foo_neigh"]["distances_key"] == "foo_dist"
    assert rna_out.uns["foo_neigh"]["connectivities_key"] == "foo_conn"
    assert "neighbors" not in rna_out.uns
    assert "distances" not in rna_out.obsp
    assert "connectivities" not in rna_out.obsp


def test_neighbors_within_batch(run_component, random_h5mu_path):
    """A smaller neighbors_within_batch should produce a sparser KNN graph."""
    input_with_batch = _prepare_input_with_batch(random_h5mu_path)
    output_default = random_h5mu_path()
    output_small = random_h5mu_path()

    run_component(
        [
            "--input",
            str(input_with_batch),
            "--output",
            output_default,
            "--batch_key",
            "batch",
            "--output_compression",
            "gzip",
        ]
    )
    run_component(
        [
            "--input",
            str(input_with_batch),
            "--output",
            output_small,
            "--batch_key",
            "batch",
            "--output_compression",
            "gzip",
            "--neighbors_within_batch",
            "2",
        ]
    )

    rna_default = mu.read_h5mu(output_default).mod["rna"]
    rna_small = mu.read_h5mu(output_small).mod["rna"]

    nnz_default = rna_default.obsp["distances"].nnz
    nnz_small = rna_small.obsp["distances"].nnz
    assert nnz_small < nnz_default, (
        "Smaller neighbors_within_batch should produce a sparser distance matrix"
    )


def test_raise_if_obsm_input_missing(run_component, random_h5mu_path):
    """A missing --obsm_input slot should raise an error."""
    input_with_batch = _prepare_input_with_batch(random_h5mu_path)
    output = random_h5mu_path()
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                str(input_with_batch),
                "--output",
                output,
                "--batch_key",
                "batch",
                "--obsm_input",
                "X_nonexistent",
            ]
        )
    assert "obsm slot 'X_nonexistent' not found" in err.value.stdout.decode("utf-8")


def test_raise_if_batch_key_missing(run_component, random_h5mu_path):
    """A missing --batch_key column should raise an error."""
    input_with_batch = _prepare_input_with_batch(random_h5mu_path)
    output = random_h5mu_path()
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                str(input_with_batch),
                "--output",
                output,
                "--batch_key",
                "nonexistent_batch",
            ]
        )
    assert "batch_key 'nonexistent_batch' not found" in err.value.stdout.decode("utf-8")


def test_overwrite_existing_slot(run_component, random_h5mu_path):
    """Writing into existing .uns/.obsp slots should fail without --overwrite
    and succeed with it."""
    input_with_batch = _prepare_input_with_batch(random_h5mu_path)
    first = random_h5mu_path()
    run_component(
        [
            "--input",
            str(input_with_batch),
            "--output",
            first,
            "--batch_key",
            "batch",
            "--output_compression",
            "gzip",
        ]
    )

    second = random_h5mu_path()
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                str(first),
                "--output",
                second,
                "--batch_key",
                "batch",
                "--output_compression",
                "gzip",
            ]
        )
    assert "but field already exists" in err.value.stdout.decode("utf-8")

    run_component(
        [
            "--input",
            str(first),
            "--output",
            second,
            "--batch_key",
            "batch",
            "--output_compression",
            "gzip",
            "--overwrite",
            "true",
        ]
    )
    assert "neighbors" in mu.read_h5mu(second).mod["rna"].uns


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
