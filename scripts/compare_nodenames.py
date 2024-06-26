import MDSplus as mds
from pprint import pprint

shot = 1150805022

tree = mds.Tree("cmod", shot)

efit_cols = {
    "beta_n": r'\efit_aeqdsk:betan',
    "beta_p": r'\efit_aeqdsk:betap',
    "kappa": r'\efit_aeqdsk:eout',
    "li": r'\efit_aeqdsk:ali',
    "upper_gap": r'\efit_aeqdsk:otop',
    "lower_gap": r'\efit_aeqdsk:obott',
    "q0": r'\efit_aeqdsk:qqmagx',
    "qstar": r'\efit_aeqdsk:qsta',
    "q95": r'\efit_aeqdsk:qpsib',
    "v_loop_efit": r'\efit_aeqdsk:vloopt',
    "Wmhd": r'\efit_aeqdsk:wplasm',
    "ssep": r'\efit_aeqdsk:ssep',
    "n_over_ncrit": r'\efit_aeqdsk:xnnc',
    "tritop": r'\efit_aeqdsk:doutu',
    "tribot":  r'\efit_aeqdsk:doutl',
    "a_minor": r'\efit_aeqdsk:aout',
    "R0":r'\efit_aeqdsk:rmagx',
    "chisq":r'\efit_aeqdsk:tsaisq',
    "area":r'\efit_aeqdsk:areao',
    }

def get_original(node_name: str, parent=tree.getNode("mhd.analysis.efit.results.a_eqdsk")):
    """
    Get the original column name for an alias in MDSPlus.

    E.g. after 2000 people started using `li` to refer to `ali`.
    `get_original('li', parent)` --> `'ali'`
    
    Params:
        node_name: str, node name to get the original of
        parent: TreeNode, MDSPlus parent node under which is the original
            and alias. E.g. 
    Returns:
        the original name or expression if node_name exists in the parent, otherwise False.
    """
    try:
        node = parent.getNode(node_name)
        node_name_is_alias = type(node.getRecord()) != mds.compound.Signal
        if node_name_is_alias:
            return node.getRecord()
        else:
            return node_name
    except mds.mdsExceptions.TreeNNF:
        return False

all_node_names = list(efit_cols.keys()) + list(efit_cols.values())

alias_to_original = {}
for node_name in all_node_names:
    alias_to_original[node_name] = get_original(node_name)
    
pprint(alias_to_original)
