#!/usr/bin/env python3

import argparse
import json
from multiprocessing import Pool
import re
import sys
import warnings

import numpy as np
import pandas as pd

import MDSplus as mds
from MDSplus.tree import TreeNode
from MDSplus.compound import Signal
from MDSplus.mdsExceptions import TreeNODATA, TreeFOPENR


class NodeNameComparer:
    """Compare two node names in MDSPlus to determine which one you should rely
    on for fetching data and which one is an alias that may have only been created
    for newer shots.

    Terminology in this class:
    - A 'ref' node is the one that we want to exist for as many shots as possible
    so that we can count on using it for all shots.
    - A 'new' node is the node that may not exist for all shots. It could be that
    the new node is an alias that was created for new shots for convenience.

    E.g. creating a comparison table for `ali` as ref and `li` as new shows `ali`
    has data for all shots and `li` is an alias for only a fraction of shots because
    it was added only to new shots.
    """

    def __init__(self, shots: list[int], ref: str, new: str = None, num_processes=1, expression:str=None):
        """
        Params:
            shots: List of shot ids. If None, then the default is 20 old shots
                and 20 new shots
            ref: Referent node name
            new: (optional) New, alias, node name
            num_processes: number of processes with which to open MDSPlus trees
            expression: str, the expected expression (if applicable)
        """
        if ref is None:
            raise ValueError("ref must be specified")

        if isinstance(shots, list):
            self.shots = shots
        else:
            self.shots = [
                941130004,
                941130006,
                941130007,
                941130009,
                941130010,
                1010626011,
                1010626012,
                1010627003,
                1010627006,
                1010627007,
            ]
        self.ref = ref
        self.new = new

        self.labels = ["values", "alias", "expression", "empty"]
        self.VALUES_IDX = 0
        self.ALIAS_IDX = 1
        self.EXPRESSION_IDX = 2
        self.EMPTY_IDX = 3

        self.parent_name = r"\efit_aeqdsk"
        self.num_processes = min(len(self.shots), num_processes)

        self.expression = expression

    @staticmethod
    def is_alias(record: Signal) -> bool:
        """Return True if the node points to another record that has data."""
        try:
            record.getNid()
            return True
        except (TreeNODATA, TreeFOPENR, AttributeError):
            pass
        return False

    def is_expression(self, record: Signal) -> bool:
        """Return True if the node is an expression made of other nodes."""
        operators = ["/", "*", "+", "-"]
        node_str = str(record)
        match = re.search(r"Build_Signal\(([^,]+),", node_str)
        expression = match.group(1) if match is not None else node_str 
        for op in operators:
            if op in expression:
                if self.expression != expression:
                    warnings.warn(f"expression does not match expected '{expression}' != '{self.expression}'")
                return True
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

        ref_record, new_record = None, None

        ref_node = NodeNameComparer.get_node(f"{self.parent_name}:{self.ref}", tree)
        if ref_node is not None:
            try:
                ref_record = ref_node.getRecord()
            except (TreeNODATA, TreeFOPENR, AttributeError):
                pass

        new_node = NodeNameComparer.get_node(f"{self.parent_name}:{self.new}", tree)
        if new_node is not None:
            try:
                new_record = new_node.getRecord()
            except (TreeNODATA, TreeFOPENR, AttributeError):
                pass

        # ref has no data
        if ref_node is None:
            warnings.warn(f"ref node '{self.ref}' not found, shot: {shot}")
            # new has no data
            if new_node is None:
                return (self.EMPTY_IDX, self.EMPTY_IDX)
            # new is an alias
            elif NodeNameComparer.is_alias(new_record):
                return (self.ALIAS_IDX, self.EMPTY_IDX)
            # new is an expression
            elif self.is_expression(new_record):
                return (self.EXPRESSION_IDX, self.EMPTY_IDX)
            # new has data
            else:
                return (self.VALUES_IDX, self.EMPTY_IDX)
        # ref is an alias
        elif NodeNameComparer.is_alias(ref_record):
            if new_node is None:
                warnings.warn(
                    f"ref '{self.ref}' is an alias, but new '{self.new}' has no "
                    + f"data, shot: {shot}"
                )
                return (self.EMPTY_IDX, self.ALIAS_IDX)
            elif NodeNameComparer.is_alias(new_record):
                warnings.warn(
                    f"No data found: ref '{self.ref}' is an alias and new '{self.new}'"
                    + f" is an alias, shot: {shot}"
                )
            elif self.is_expression(new_record):
                return (self.EXPRESSION_IDX, self.ALIAS_IDX)
            else:
                # ref is an alias for new data
                if ref_record.getNid() == new_node.getNid():
                    return (self.VALUES_IDX, self.ALIAS_IDX)
                # ref is an alias, new is data, but ref does not point to new
                else:
                    warnings.warn(
                        f"No connection: ref '{self.ref}' does not point to new "
                        + f"'{self.new}', shot: {shot}"
                    )
        # ref is expression
        elif self.is_expression(ref_record):
            if new_node is None:
                return (self.EMPTY_IDX, self.EXPRESSION_IDX)
            elif NodeNameComparer.is_alias(new_record):
                # new is an alias for ref
                if new_record.getNid() == ref_node.getNid():
                    return (self.ALIAS_IDX, self.EXPRESSION_IDX)
                # new is an alias, but not for ref
                else:
                    warnings.warn(
                        f"No connection: new '{self.new}' does not point to ref "
                        + f"'{self.ref}', shot: {shot}"
                    )
            # new has data
            else:
                return (self.VALUES_IDX, self.EXPRESSION_IDX)
        # ref has data
        else:
            if new_node is None:
                return (self.EMPTY_IDX, self.VALUES_IDX)
            elif NodeNameComparer.is_alias(new_record):
                # new is an alias for ref
                if new_record.getNid() == ref_node.getNid():
                    return (self.ALIAS_IDX, self.VALUES_IDX)
                # new is an alias, but not for ref
                else:
                    warnings.warn(
                        f"No connection: new '{self.new}' does not point to ref "
                        + f"'{self.ref}', shot: {shot}"
                    )
            elif self.is_expression(new_record):
                return (self.EXPRESSION_IDX, self.VALUES_IDX)
            # new has data
            else:
                warnings.warn(
                    f"No connection: new '{self.new}' and ref '{self.ref}' both "
                    + f"have data, shot: {shot}"
                )

    def get_comparison_table(self) -> np.ndarray:
        """
        Create a Numpy array for the comparison matrix and saves a json of the shots.

        The columns and rows, for ref and new respectively, are
        | values | alias | expression | empty |
        """
        json_keys = {
            (self.VALUES_IDX, self.VALUES_IDX): "vv",
            (self.VALUES_IDX, self.ALIAS_IDX): "va",
            (self.VALUES_IDX, self.EXPRESSION_IDX): "vx",
            (self.VALUES_IDX, self.EMPTY_IDX): "ve",

            (self.ALIAS_IDX, self.VALUES_IDX): "av",
            (self.ALIAS_IDX, self.ALIAS_IDX): "aa",
            (self.ALIAS_IDX, self.EXPRESSION_IDX): "ax",
            (self.ALIAS_IDX, self.EMPTY_IDX): "ae",
            
            (self.EXPRESSION_IDX, self.VALUES_IDX): "xv",
            (self.EXPRESSION_IDX, self.ALIAS_IDX): "xa",
            (self.EXPRESSION_IDX, self.EXPRESSION_IDX): "xx",
            (self.EXPRESSION_IDX, self.EMPTY_IDX): "xe",
            
            (self.EMPTY_IDX, self.VALUES_IDX): "ev",
            (self.EMPTY_IDX, self.ALIAS_IDX): "ea",
            (self.EMPTY_IDX, self.EXPRESSION_IDX): "ex",
            (self.EMPTY_IDX, self.EMPTY_IDX): "ee",
        }
        shot_log = {}
        comparison_table = np.zeros(shape=(len(self.labels), len(self.labels)))
        with Pool(self.num_processes) as p:
            table_indices = p.map(self.compare_names_one_shot, self.shots)
            for i in range(len(table_indices)):
                idx = table_indices[i]
                if idx is None:
                    continue
                comparison_table[idx] += 1

                if json_keys[idx] not in shot_log:
                    shot_log[json_keys[idx]] = []
                shot_log[json_keys[idx]].append(self.shots[i])

        with open(f"{self.ref}-{self.new}.json", "w") as f:
            f.write(json.dumps(shot_log))
        return comparison_table

    def get_pretty_table(self) -> pd.DataFrame:
        """Convert values in table to percents, add column and row labels."""
        comparison_table = self.get_comparison_table()
        comparison_table = np.vectorize(
            lambda item: f"{item:.0f} ({item / len(self.shots) * 100:.0f}%)"
        )(comparison_table)
        comparison_table = pd.DataFrame(comparison_table)
        comparison_table.columns = self.labels
        comparison_table.index = self.labels
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

    parser.add_argument(
        "--expr",
        type=str,
        action="store",
        help="Expected expression, e.g. '\\EFIT_AEQDSK:AOUT / 100.'",
    )

    args = parser.parse_args()
    shot_ids = None

    # Grab from stdin
    if args.shot_ids is None and not sys.stdin.isatty():
        shot_ids = [int(s.strip()) for s in sys.stdin.read().split()]
    # Parse list of shot ids
    elif args.shot_ids and len(args.shot_ids) > 0 and args.shot_ids[0].isdigit():
        shot_ids = [int(s) for s in args.shot_ids]
    # Parse txt files
    elif args.shot_ids and args.shot_ids[0].endswith(".txt"):
        with open(args.shot_ids[0], "r") as f:
            shot_ids = [int(s) for s in f.readlines()]

    compare = NodeNameComparer(shot_ids, args.ref, args.new, args.num_processes, args.expr)
    table = compare.get_pretty_table()
    print("      ref:" + f"{compare.parent_name}:{args.ref}")
    print(table)
    print(f"{compare.parent_name}:{args.new}")
