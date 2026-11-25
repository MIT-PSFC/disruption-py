The [`XarrayDataConnection`][disruption_py.inout.xr.XarrayDataConnection] class reads shot data from [Xarray](https://xarray.dev/) stores, such as Zarr or NetCDF files, either local or hosted on S3.
It implements the same [`DataConnection`][disruption_py.inout.base.DataConnection] interface as [`MDSConnection`][disruption_py.inout.mds.MDSConnection], so physics methods use it through `params.get_data(...)` without knowing which backend is in use.

## How it is selected { .doc .doc-heading }

[`get_process_connection`][disruption_py.workflow.get_process_connection] picks the backend from the tokamak's configuration:
an `inout.mds` section selects MDSplus, otherwise an `inout.xarray` section selects Xarray.
Currently, MAST is the only device configured for Xarray:

```toml
# disruption_py/machine/mast/config.toml
[mast.inout.xarray]
endpoint_url = "https://s3.echo.stfc.ac.uk"
file_ext = "zarr"
file_path = "mast/level2/shots"
```

Each shot is opened as an Xarray `DataTree` from `{endpoint_url}/{file_path}/{shot_id}.{file_ext}`.
Files with the `zarr` extension are read with the Zarr engine, and any other extension with the NetCDF4 engine.

To read from a local copy of the data instead, override these values in your [user configuration file](../../installation/configuration.md), for example:

```toml
# ~/.config/disruption-py/user.toml
[mast.inout.xarray]
endpoint_url = ""
file_path = "/path/to/local/shots"
```

## Data paths { .doc .doc-heading }

Variables are addressed by their path within the `DataTree`, for example `"summary/ip"`.
The optional `group` argument is prepended to paths that do not already contain a `/`, so `get_data("ip", group="summary")` is equivalent to `get_data("summary/ip")`.
Pass `return_xarray=True` to get an `xr.DataArray`, with its coordinates, instead of a NumPy array.

## Xarray Connection Reference { .doc .doc-heading }

::: disruption_py.inout.xr
    handler: python
    options:
        heading_level: 3
        show_root_heading: false
        show_root_toc_entry: false
        filters: ["!^_[^_]"]
        members_order: "source"
