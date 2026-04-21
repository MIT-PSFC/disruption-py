DisruptionPy defines abstract base classes for data connections, allowing physics methods to work with any backend (MDSplus, Xarray/Zarr, etc.) through a unified interface.

## Connection Architecture { .doc .doc-heading }

[`DataConnection`][disruption_py.inout.base.DataConnection] is the per-shot data access interface. Each instance is bound to a single shot and provides methods for retrieving data and dimensions. Concrete implementations include [`MDSConnection`][disruption_py.inout.mds.MDSConnection] for MDSplus and `XarrayDataConnection` for Xarray/Zarr.

[`ProcessConnection`][disruption_py.inout.base.ProcessConnection] is the per-process factory that creates `DataConnection` instances for each shot. This separation ensures multiprocessing safety -- each process holds its own `ProcessConnection`, which produces independent `DataConnection` objects per shot.

## Base Class Reference { .doc .doc-heading }

::: disruption_py.inout.base
    handler: python
    options:
        heading_level: 3
        show_root_heading: false
        show_root_toc_entry: false
        filters: ["!^_[^_]"]
        members_order: "source"
