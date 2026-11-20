DisruptionPy uses logbook sql databases for convenience when retrieving data from MDSplus. Users may also use DisruptionPy to directly retrieve data from the logbook database's `disruption_warning` tables.

## The `disruption_warning` table { .doc .doc-heading }
The `disruption_warning` sql tables for C-Mod and DIII-D contain important disruption parameters for a large number of shots.
Each dataset is described on its device page:

- [C-Mod `disruption_warning` dataset](../../devices/cmod/disruption_warning.md)
- [DIII-D `disruption_warning` dataset](../../devices/d3d/disruption_warning.md)

## Retrieving data from the SQL database { .doc .doc-heading }
Here is an example retrieving data from `disruption_warning` or `disruptions` table
```python
--8<--
examples/sql.py
--8<--
```

## Database Class Reference { .doc .doc-heading }

::: disruption_py.inout.sql
    handler: python
	options:
	  heading_level: 3
	  show_root_heading: false
	  show_root_toc_entry: false
	  filters: ["!^_[^_]"]
        
