#!/usr/bin/env python3

import argparse
from multiprocessing import Pool
import numpy as np
import pandas as pd
import sys
import warnings

from disruption_py.handlers.cmod_handler import CModHandler
from disruption_py.utils.constants import MAX_PROCESSES
import MDSplus as mds
from MDSplus.tree import TreeNode

NODENAMES_SHOT_LIST_PATH = "disruption_py/data/compare_nodenames_default.txt"


class NodeNameComparer:
    """Compare two node names in MDSPlus to determine which one you should rely
    on for fetching data and which one is an alias that may have only been created
    for newer shots.

    Terminology in this class:
    - A 'ref' node is one that contains data or an expression, and not merely a
    pointer to another node.
    - A 'new' node is an alias pointing to another node.

    E.g. creating a comparison table for `ali` as ref and `li` as new shows `ali`
    has data for all shots and `li` is an alias for only a fraction of shots because
    it was added only to new shots.
    """

    def __init__(self, shots: list[int], ref: str, new: str = None, num_processes=1):
        """
        Params:
            shots: List of shot ids. If None, then the default is 20 old shots
                and 20 new shots
            ref: Referent node name
            new: (optional) New, alias, node name
            num_processes: number of processes with which to open MDSPlus trees
        """
        if ref is None:
            raise ValueError("ref must be specified")

        if isinstance(shots, list):
            self.shots = shots
        else:
            self.shots = NodeNameComparer.get_default_shot_list()
        self.ref = ref
        self.new = new
        self.columns = ["values", "alias", "empty"]
        self.parent_name = r"\efit_aeqdsk"
        self.num_processes = min([len(self.shots), num_processes, MAX_PROCESSES])

    @staticmethod
    def get_default_shot_list(use_default_list=True) -> list[int]:
        """Return 20 new and 20 old shots ids."""
        with open(NODENAMES_SHOT_LIST_PATH, "r") as f:
            return [int(s) for s in f.readlines()]

    @staticmethod
    def is_alias(node: TreeNode) -> bool:
        """Return True if the node points to another record that has data."""
        try:
            node.getRecord().getNid()
            return True
        except (mds.mdsExceptions.TreeNODATA, AttributeError):
            pass
        return False

    @staticmethod
    def get_node(node_name: str, parent: TreeNode) -> TreeNode:
        """Return the MDSPlus node off of the parent node."""
        try:
            return parent.getNode(node_name)
        except mds.mdsExceptions.TreeNNF:
            return None

    def compare_names_one_shot(self, shot: int) -> tuple[int, int]:
        """Return the indices of the comparison table representing a relationship
        between ref and node."""
        try:
            tree = mds.Tree("analysis", shot)
        except mds.mdsExceptions.TreeFOPENR:
            return None
        parent = tree.getNode(self.parent_name)
        ref_node, new_node = None, None

        ref_node = NodeNameComparer.get_node(self.ref, parent)
        new_node = NodeNameComparer.get_node(self.new, parent)

        # ref has no data
        if ref_node is None:
            warnings.warn(f"ref node '{self.ref}' not found, shot: {shot}")
            # new has no data
            if new_node is None:
                return (2, 2)
            # new is an alias
            elif NodeNameComparer.is_alias(new_node):
                return (2, 1)
            # new has data
            else:
                return (0, 2)
        # ref is an alias
        elif NodeNameComparer.is_alias(ref_node):
            if new_node is None:
                warnings.warn(
                    f"ref '{self.ref}' is an alias, but new '{self.new}' has no "
                    + f"data, shot: {shot}"
                )
                return (2, 1)
            elif NodeNameComparer.is_alias(new_node):
                raise ValueError(
                    f"No data found: ref '{self.ref}' is an alias and new '{self.new}'"
                    + f" is an alias, shot: {shot}"
                )
            else:
                # ref is an alias for new data
                if ref_node.getRecord().getNid() == new_node.getNid():
                    return (0, 1)
                # ref is an alias, new is data, but ref does not point to new
                else:
                    raise ValueError(
                        f"No connection: ref '{self.ref}' does not point to new "
                        + f"'{self.new}', shot: {shot}"
                    )
        # ref has data
        else:
            if new_node is None:
                return (2, 0)
            elif NodeNameComparer.is_alias(new_node):
                # new is an alias for ref
                if new_node.getRecord().getNid() == ref_node.getNid():
                    return (1, 0)
                # new is an alias, but not for ref
                else:
                    raise ValueError(
                        f"No connection: new '{self.new}' does not point to ref "
                        + f"'{self.ref}', shot: {shot}"
                    )
            # new has data
            else:
                raise ValueError(
                    f"No connection: new '{self.new}' and ref '{self.ref}' both "
                    + f"have data, shot: {shot}"
                )

    def get_comparison_table(self) -> np.ndarray:
        """
        Create a Numpy array for the comparison matrix.

        The columns and rows, for ref and new respectively, are
        | values | alias | empty |
        """
        comparison_table = np.zeros(shape=(3, 3))
        with Pool(self.num_processes) as p:
            indices = p.map(self.compare_names_one_shot, self.shots)
            for idx in indices:
                if idx is None:
                    continue
                comparison_table[idx] += 1

        return comparison_table

    def get_pretty_table(self) -> pd.DataFrame:
        """Convert values in table to percents, add column and row labels."""
        comparison_table = self.get_comparison_table()
        comparison_table = np.vectorize(
            lambda item: f"{item / len(self.shots) * 100:.0f}% {item:.0f}"
        )(comparison_table)
        comparison_table = pd.DataFrame(comparison_table)
        comparison_table.columns = self.columns
        comparison_table.index = self.columns
        return comparison_table


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--shot-ids",
        nargs="*",
        action="store",
        default=None,
        help="Shot number(s) to test, txt file(s) containing the shots, or pipe "
        + "shot numbers from stdin. If not specified, use 20 new and 20 old shots.",
    )
    parser.add_argument(
        "--ref",
        type=str,
        action="store",
        required=True,
        help="Referent node name",
    )
    parser.add_argument(
        "--new",
        type=str,
        action="store",
        required=True,
        help="New (alias) node name",
    )
    parser.add_argument(
        "--num-processes",
        type=int,
        action="store",
        default=1,
        help="Number of processes in which to open MDSPlus trees (default=1)",
    )

    args = parser.parse_args()
    shot_ids = None

    # Grab from stdin
    if len(args.shot_ids) == 0:
        if not sys.stdin.isatty():
            shot_ids = [int(s.strip()) for s in sys.stdin.readlines()]
        else:
            shot_ids = None
    # Parse list of shot ids
    elif args.shot_ids[0].isdigit():
        shot_ids = [int(s) for s in args.shot_ids]
    # Parse list of txt files
    elif args.shot_ids[0].endswith(".txt"):
        shot_ids = []
        for file in args.shot_ids:
            with open(file, "r") as f:
                shot_ids += [int(s) for s in f.readlines()]
    else:
        shot_ids = None

    compare = NodeNameComparer(shot_ids, args.ref, args.new, args.num_processes)
    table = compare.get_pretty_table()
    print(f"ref -> \t {compare.parent_name}:{args.ref}")
    print(table)
    print(f"{compare.parent_name}:{args.new}")
