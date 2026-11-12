#!/usr/bin/env python3

"""
Demonstrate the implications of the efit nickname cascade (PR #559).

Run on the CMOD cluster (mfews04):

    uv run python scripts/nickname_cascade_validation.py

Greg's ask on PR #559: show what a cascade like ``[efit18, analysis, efit21]``
actually does versus the single-tree settings, and assert on the output so the
behaviour is pinned, not just eyeballed. The "primary difference between those
efit trees" is the number of timeslices, so we report the per-shot row count and
assert the equality relationships.

Shots (CMOD):
  * 1150805012 - disruptive; has efit18, analysis, and efit21
  * 1150805013 - non-disruptive; has analysis and efit21, but NO efit18

Builds on the snippet @yumouwei posted in the PR thread.
"""

import numpy as np

from disruption_py.settings import RetrievalSettings
from disruption_py.workflow import get_shots_data

DISRUPTIVE = 1150805012  # has efit18, analysis, efit21
NON_DISRUPTIVE = 1150805013  # has analysis, efit21; NO efit18
SHOTLIST = [DISRUPTIVE, NON_DISRUPTIVE]


def fetch(efit_nickname_setting):
    """Run get_shots_data for one efit setting and return the dataset."""
    settings = RetrievalSettings(
        efit_nickname_setting=efit_nickname_setting,
        run_columns=["li"],
        time_setting="efit",
    )
    return get_shots_data(
        shotlist_setting=SHOTLIST,
        retrieval_settings=settings,
    )


def rows_for(ds, shot):
    """Return the slice of the flat (idx,) dataset belonging to one shot.

    The default 'dataset' output concatenates every shot along ``idx`` with
    ``shot``/``time`` as non-dimension coords, so we mask on shot rather than
    use ``.sel(shot=...)`` (shot is not a dimension coordinate).
    """
    mask = ds["shot"].values == shot
    return ds.isel(idx=np.where(mask)[0])


def length(ds, shot):
    """Number of timeslices retrieved for a shot (0 if the shot is absent)."""
    return int((ds["shot"].values == shot).sum())


def report(label, ds):
    """Print the per-shot row count for one setting."""
    print(f"\n=== {label} ===")
    for shot in SHOTLIST:
        n = length(ds, shot)
        present = "present" if n else "ABSENT"
        print(f"  shot {shot}: {n:>5} rows ({present})")


def main():
    """Fetch each setting, report row counts, and assert the cascade behaviour."""
    # Cascade: prefer efit18, fall back to analysis, then efit21.
    cascade = fetch(["efit18", "analysis", "efit21"])
    efit18 = fetch("efit18")
    analysis = fetch("analysis")
    efit21 = fetch("efit21")

    report("cascade [efit18, analysis, efit21]", cascade)
    report("efit18 only", efit18)
    report("analysis only", analysis)
    report("efit21 only", efit21)

    print("\n=== Assertions ===")

    # Timeslice count is the primary difference between efit trees: the cascade
    # picks efit18, so it matches efit18's length, not analysis's.
    assert length(efit18, DISRUPTIVE) != length(analysis, DISRUPTIVE)
    assert length(cascade, DISRUPTIVE) == length(efit18, DISRUPTIVE)
    print(
        f"  OK  len(efit18[{DISRUPTIVE}])={length(efit18, DISRUPTIVE)}"
        f" != len(analysis[{DISRUPTIVE}])={length(analysis, DISRUPTIVE)} (trees differ in length)"
    )

    # Disruptive shot has efit18 -> cascade picks efit18, bit-identical to efit18-only.
    assert rows_for(cascade, DISRUPTIVE).equals(rows_for(efit18, DISRUPTIVE))
    print(f"  OK  cascade[{DISRUPTIVE}] == efit18[{DISRUPTIVE}] (first item wins)")

    # ...and efit18 differs from analysis for that shot (distinct equilibria).
    assert not rows_for(cascade, DISRUPTIVE).equals(rows_for(analysis, DISRUPTIVE))
    print(f"  OK  cascade[{DISRUPTIVE}] != analysis[{DISRUPTIVE}] (trees differ)")

    # Non-disruptive shot lacks efit18 entirely.
    assert length(efit18, NON_DISRUPTIVE) == 0
    print(f"  OK  efit18[{NON_DISRUPTIVE}] absent (no efit18 tree)")

    # ...so the cascade falls through to analysis: same length, same data.
    assert length(cascade, NON_DISRUPTIVE) == length(analysis, NON_DISRUPTIVE)
    assert rows_for(cascade, NON_DISRUPTIVE).equals(rows_for(analysis, NON_DISRUPTIVE))
    print(
        f"  OK  cascade[{NON_DISRUPTIVE}] == analysis[{NON_DISRUPTIVE}] (fell through)"
    )

    print("\nAll assertions passed.")


if __name__ == "__main__":
    main()
