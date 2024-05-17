from disruption_py.mdsplus_integration.mds_connection import MDSConnection
from disruption_py.shots.helpers.cached_method_props import CachedMethodProps
from dataclasses import dataclass, field

from typing import Callable, Dict, List, Set


class MethodOptimizer:

    @dataclass
    class MethodDependecy:
        dependencies: Set[str] = field(default_factory=set)
        dependent_on: Set[str] = field(default_factory=set)

    def __init__(
        self,
        mds_conn: MDSConnection,
        parameter_cached_method_props: List[CachedMethodProps],
        all_cached_method_props: List[CachedMethodProps],
        pre_cached_method_names: List[str],
    ):
        self._mds_conn = mds_conn
        self._cached_method_props_by_name: Dict[str, CachedMethodProps] = {
            method.name: method for method in all_cached_method_props
        }

        # compute method dependency graph
        graph: dict[str, MethodOptimizer.MethodDependecy] = {}
        visited_method_names = {
            parameter_method.name for parameter_method in parameter_cached_method_props
        }
        parameter_cached_method_props_stack = parameter_cached_method_props.copy()
        while len(parameter_cached_method_props_stack) != 0:
            cached_method_props: CachedMethodProps = (
                parameter_cached_method_props_stack.pop()
            )
            # if cached_method.name pre cached we don't want to run it or run any of its contained cached methods
            if cached_method_props.name in pre_cached_method_names:
                continue
            graph.setdefault(
                cached_method_props.name, MethodOptimizer.MethodDependecy()
            )
            # contained cached methods is a list of names of used cached methods
            for contained_cached_method_name in cached_method_props.get_param_value(
                "contained_cached_methods", []
            ):
                if (
                    contained_cached_method_name
                    not in self._cached_method_props_by_name
                ):
                    continue
                if contained_cached_method_name not in visited_method_names:
                    visited_method_names.add(contained_cached_method_name)
                    parameter_cached_method_props_stack.append(
                        self._cached_method_props_by_name[contained_cached_method_name]
                    )
                # get method dependencies for contained cached method, creating one if it does not already exist
                contained_cached_method_dependency = graph.setdefault(
                    contained_cached_method_name, MethodOptimizer.MethodDependecy()
                )
                # add in contained cached method that cached method is dependent on it
                contained_cached_method_dependency.dependencies.add(
                    cached_method_props.name
                )
                # add that cached method is dependent on contained cached method
                # (contained cached method must complete before it can run)
                graph[cached_method_props.name].dependent_on.add(
                    contained_cached_method_name
                )
        self._method_dependencies = graph

        # compute used trees by all methods
        tree_remaining_count = {}
        for cached_method_name in set(self._method_dependencies.keys()) | getattr(
            self, "_methods_in_progress", set()
        ):
            for tree_name in self._get_method_used_trees(cached_method_name):
                # count initialized to 0, as represents remaining count after chosen method completed
                cur_count = tree_remaining_count.get(tree_name, 0)
                tree_remaining_count[tree_name] = cur_count + 1
        self._tree_remaining_count = tree_remaining_count

    def next_method(self) -> CachedMethodProps:
        if len(self._method_dependencies) == 0:
            return None
        # get methods that are not dependent on another cached method
        # guranteed to never be empty unless there is a cycle
        methods_to_consider = [
            cached_method_name
            for cached_method_name in self._method_dependencies
            if len(self._method_dependencies[cached_method_name].dependent_on) == 0
        ]
        if len(methods_to_consider) == 0:
            raise Exception("Cycle detected in cached method dependencies")

        method_tree_counts = []
        open_tree_names = self._mds_conn.open_trees
        for cached_method_name in methods_to_consider:
            used_trees = self._get_method_used_trees(cached_method_name)

            # estimated remaining open tree count after this method has been run
            remaining_open_tree_count = sum(
                self._tree_remaining_count[tree_name] for tree_name in used_trees
            ) - len(used_trees)
            close_estimate = len(
                [
                    tree_name
                    for tree_name in used_trees
                    if self._tree_remaining_count[tree_name] == 1
                    and tree_name in open_tree_names
                ]
            )
            open_estiamte = len(
                [
                    tree_name
                    for tree_name in used_trees
                    if self._tree_remaining_count[tree_name] != 1
                    and tree_name not in open_tree_names
                ]
            )
            method_tree_counts.append(
                (
                    open_estiamte - close_estimate,
                    remaining_open_tree_count,
                    cached_method_name,
                )
            )
        next_method_name = min(method_tree_counts)[2]

        # update self._method_dependencies for method removal
        self._method_dependencies.pop(next_method_name)

        next_method = self._cached_method_props_by_name[next_method_name]

        return next_method

    def run_methods_sync(self, method_executor: Callable[[CachedMethodProps], None]):
        while True:
            # get next method
            next_method: CachedMethodProps = self.next_method()
            if next_method is None:
                break

            # run method
            method_executor(next_method)

            # handle method completion
            for method_to_check in self._method_dependencies:
                # remove method that has compeleted from all methods dependent_on
                self._method_dependencies[method_to_check].dependent_on.discard(
                    next_method.name
                )

            # update self._tree_remaining_count
            for tree_name in self._get_method_used_trees(next_method.name):
                self._tree_remaining_count[tree_name] -= 1

                # close trees that are no longer needed
                if self._tree_remaining_count[tree_name] <= 0:
                    self._tree_remaining_count.pop(tree_name)
                    self._mds_conn.close_tree(tree_name)

    def _get_method_used_trees(self, method_name):
        """
        Get all used trees for a given method.

        Specifies unique name as must treat nicknames as there real tree name.
        """
        used_trees = self._cached_method_props_by_name[method_name].get_param_value(
            "used_trees", []
        )
        return [self._mds_conn.tree_name(used_tree) for used_tree in used_trees]
