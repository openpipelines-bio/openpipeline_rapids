import sys
import subprocess
import re
import pytest
import mudata as mu
import numpy as np
import scanpy as sc
from openpipeline_testutils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    "name": "pca",
    "resources_dir": "resources_test/",
    "config": "src/preprocessing/pca/config.vsh.yaml",
    "executable": "target/docker/preprocessing/pca/pca",
}
## VIASH END

input = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"


@pytest.fixture
def clean_input(random_h5mu_path):
    """Prepare the PBMC test data for rapids-singlecell PCA.

    The "_mms" input is the PBMC dataset after the standard openpipeline RNA
    single-sample -> multi-sample -> dimensionality-reduction pipeline, so it
    is normalized/log-transformed and already carries PCA slots. Here we filter
    out zero-expression genes (rsc PCA refuses them) and drop the existing
    default PCA output slots so PCA can write into them without --overwrite."""
    mu_in = mu.read_h5mu(input)
    rna = mu_in.mod["rna"]
    sc.pp.filter_genes(rna, min_counts=1)
    rna.obsm.pop("X_pca", None)
    rna.varm.pop("pca_loadings", None)
    rna.uns.pop("pca_variance", None)
    path = random_h5mu_path()
    mu_in.write(path)
    return path


def test_run(run_component, random_h5mu_path, clean_input):
    output = random_h5mu_path()
    cmd_pars = [
        "--input",
        clean_input,
        "--output",
        output,
        "--output_compression",
        "gzip",
        "--num_components",
        "26",
    ]
    run_component(cmd_pars)

    assert output.is_file(), "No output was created."

    mu_input = mu.read_h5mu(clean_input)
    mu_output = mu.read_h5mu(output)

    assert "rna" in mu_output.mod, 'Output should contain data.mod["rna"].'

    rna_in = mu_input.mod["rna"]
    rna_out = mu_output.mod["rna"]

    assert rna_in.shape == rna_out.shape, "Should have same shape as before"
    assert rna_out.obsm["X_pca"].shape == (rna_in.n_obs, 26), (
        "obsm['X_pca'] should hold the PCA embedding with the requested "
        "number of components."
    )
    assert rna_out.varm["pca_loadings"].shape == (rna_in.n_vars, 26), (
        "varm['pca_loadings'] should hold the PC loadings."
    )
    assert "pca_variance" in rna_out.uns
    assert "variance" in rna_out.uns["pca_variance"]
    assert "variance_ratio" in rna_out.uns["pca_variance"]
    assert not np.array_equal(
        rna_out.uns["pca_variance"]["variance"],
        rna_out.uns["pca_variance"]["variance_ratio"],
    )

    # Copy over the PCA outputs so that the rest of the objects can be compared
    mu_input["rna"].obsm["X_pca"] = mu_output["rna"].obsm["X_pca"]
    mu_input["rna"].varm["pca_loadings"] = mu_output["rna"].varm["pca_loadings"]
    mu_input["rna"].uns["pca_variance"] = mu_output["rna"].uns["pca_variance"]
    assert_annotation_objects_equal(mu_input, mu_output)


def test_custom_output_slots(run_component, random_h5mu_path, clean_input):
    """Custom obsm/varm/uns output keys should be honored."""
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            clean_input,
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--num_components",
            "26",
            "--obsm_output",
            "X_foo",
            "--varm_output",
            "foo_loadings",
            "--uns_output",
            "foo_variance",
        ]
    )

    rna_out = mu.read_h5mu(output).mod["rna"]
    assert rna_out.obsm["X_foo"].shape == (rna_out.n_obs, 26)
    assert rna_out.varm["foo_loadings"].shape == (rna_out.n_vars, 26)
    assert "foo_variance" in rna_out.uns
    assert "X_pca" not in rna_out.obsm
    assert "pca_loadings" not in rna_out.varm
    assert "pca_variance" not in rna_out.uns


def test_input_layer(run_component, random_h5mu_path, clean_input):
    """PCA should run on the requested layer when --layer is provided."""
    mu_orig = mu.read_h5mu(clean_input)
    mu_orig.mod["rna"].layers["counts"] = mu_orig.mod["rna"].X.copy()
    input_with_layer = random_h5mu_path()
    mu_orig.write(input_with_layer)

    output = random_h5mu_path()
    run_component(
        [
            "--input",
            str(input_with_layer),
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--num_components",
            "26",
            "--layer",
            "counts",
        ]
    )

    rna_out = mu.read_h5mu(output).mod["rna"]
    assert rna_out.obsm["X_pca"].shape == (rna_out.n_obs, 26)
    assert rna_out.varm["pca_loadings"].shape == (rna_out.n_vars, 26)


def test_var_input(run_component, random_h5mu_path, clean_input):
    """When --var_input is provided, only genes flagged True should
    contribute to the loadings (other rows must be zero)."""
    mu_orig = mu.read_h5mu(clean_input)
    n_vars = mu_orig.mod["rna"].n_vars
    mask = np.zeros(n_vars, dtype=bool)
    mask[: n_vars // 2] = True
    mu_orig.mod["rna"].var["filter_with_hvg"] = mask
    input_with_mask = random_h5mu_path()
    mu_orig.write(input_with_mask)

    output = random_h5mu_path()
    run_component(
        [
            "--input",
            str(input_with_mask),
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--num_components",
            "26",
            "--var_input",
            "filter_with_hvg",
        ]
    )

    rna_out = mu.read_h5mu(output).mod["rna"]
    loadings = rna_out.varm["pca_loadings"]
    assert loadings.shape == (n_vars, 26)
    assert np.all(loadings[~mask] == 0), "Loadings for masked-out genes should be zero."
    assert not np.all(loadings[mask] == 0), (
        "Loadings for selected genes should not all be zero."
    )


def test_var_input_missing_column_raises(run_component, random_h5mu_path):
    output = random_h5mu_path()
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input,
                "--output",
                output,
                "--output_compression",
                "gzip",
                "--num_components",
                "26",
                "--var_input",
                "does_not_exist",
            ]
        )
    assert re.search(
        r"Requested to use \.var column does_not_exist as a selection of "
        r"genes to run the PCA on, but the column is not available",
        err.value.stdout.decode("utf-8"),
    )


def test_missing_layer_raises(run_component, random_h5mu_path):
    output = random_h5mu_path()
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input,
                "--output",
                output,
                "--output_compression",
                "gzip",
                "--num_components",
                "26",
                "--layer",
                "does_not_exist",
            ]
        )
    assert re.search(
        r"does_not_exist was not found in modality rna",
        err.value.stdout.decode("utf-8"),
    )


def test_chunked_requires_chunk_size(run_component, random_h5mu_path):
    output = random_h5mu_path()
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input,
                "--output",
                output,
                "--output_compression",
                "gzip",
                "--num_components",
                "26",
                "--chunked",
            ]
        )
    assert re.search(
        r"Requested to perform an incremental PCA \('chunked'\), but the "
        r"chunk size is not set",
        err.value.stdout.decode("utf-8"),
    )


def test_chunk_size_smaller_than_num_components_raises(run_component, random_h5mu_path):
    output = random_h5mu_path()
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input,
                "--output",
                output,
                "--output_compression",
                "gzip",
                "--num_components",
                "26",
                "--chunked",
                "--chunk_size",
                "25",
            ]
        )
    assert re.search(
        r"The requested chunk size \(25\) must not be smaller than the "
        r"number of components \(26\)",
        err.value.stdout.decode("utf-8"),
    )


def test_overwrite_existing_slot(run_component, random_h5mu_path, clean_input):
    """Re-running PCA into an existing slot should fail without --overwrite
    and succeed with it."""
    first = random_h5mu_path()
    run_component(
        [
            "--input",
            clean_input,
            "--output",
            first,
            "--output_compression",
            "gzip",
            "--num_components",
            "26",
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
                "--output_compression",
                "gzip",
                "--num_components",
                "26",
            ]
        )
    assert re.search(
        r"but field already exists",
        err.value.stdout.decode("utf-8"),
    )

    # With --overwrite it should succeed
    run_component(
        [
            "--input",
            str(first),
            "--output",
            second,
            "--output_compression",
            "gzip",
            "--num_components",
            "26",
            "--overwrite",
        ]
    )
    rna_out = mu.read_h5mu(second).mod["rna"]
    assert rna_out.obsm["X_pca"].shape == (rna_out.n_obs, 26)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
