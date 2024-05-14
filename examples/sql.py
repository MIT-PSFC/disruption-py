#!/usr/bin/env python3

"""
execute a few meaningful queries to test DB connection.
"""

import os
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.mappings.tokamak_helpers import get_tokamak_from_environment, get_tokamak_database

queries = [
    "select count(distinct shot) from disruption_warning",
    "select count(distinct shot) from disruption_warning"
    + " where shot not in (select shot from disruptions)",
    "select count(distinct shot) from disruption_warning"
    + " where shot in (select shot from disruptions)",
    "select count(distinct shot) from disruptions",
]
tokamak = get_tokamak_from_environment()
db = get_tokamak_database(tokamak)

if tokamak is Tokamak.D3D:
    vals = [13245, 8055, 5190, 24219]
elif tokamak is Tokamak.CMOD:
    vals = [10435, 6640, 3795, 13785]
else:
    raise ValueError(f"Tokamak {tokamak} not supported for this example")

print(f"Initialized DB: {db.user}@{db.host}/{db.db_name}")

while queries:

    query = queries.pop(0)
    print(">", query.strip(" "))

    out = db.query(query)
    print("=", out.shape)

    print(out.iloc[0] if out.shape[0] == 1 else out)
    if vals:
        assert out.iloc[0, 0] == vals.pop(0)

    if queries:
        print()
        continue

    if not __debug__ or "PYTEST_CURRENT_TEST" in os.environ:
        break

    try:
        query = input("\n> ")
        if query:
            queries += [query]
    except (EOFError, KeyboardInterrupt):
        print()
        break
