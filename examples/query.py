#!/usr/bin/env python3

"""
execute a few meaningful queries to test DB connection.
"""

import os
from disruption_py.databases import CModDatabase, D3DDatabase

queries = [
    "select count(distinct shot) from disruption_warning",
    "select count(distinct shot) from disruption_warning"
    + " where shot not in (select shot from disruptions)",
    "select count(distinct shot) from disruption_warning"
    + " where shot in (select shot from disruptions)",
    "select count(distinct shot) from disruptions",
]

if os.path.exists("/fusion/projects/disruption_warning"):
    db = D3DDatabase.default()
else:
    db = CModDatabase.default()
print(f"Initialized DB: {db.user}@{db.host}/{db.db_name}")

while queries:

    query = queries.pop(0)
    print(">", query.strip(" "))

    out = db.query(query)
    print("=", out.shape)

    print(out.iloc[0] if out.shape[0] == 1 else out)

    if queries:
        print()
        continue

    try:
        query = input("\n> ")
        if query:
            queries += [query]
    except (EOFError, KeyboardInterrupt):
        print()
        break
