from transmission_models.models import *
from transmission_models.models.didelot_unsampled import *
from transmission_models.utils import *



def tree_to_newick(g, root=None):
    if root is None:
        roots = list(filter(lambda p: p[1] == 0, g.in_degree()))
        assert 1 == len(roots)
        root = roots[0][0]
    subgs = []
    # print("root ",root)
    for child in g[root]:
        if len(g[child]) > 0:
            subgs.append(tree_to_newick(g, root=child))
        else:
            subgs.append(child)
#     print("------",len(subgs),subgs)
    L = []

    for h in subgs:
        if isinstance(h,str):
            L.append(h)
        else:
            L.append(str(h)+":{}".format(h.t_inf-root.t_inf))
    return "("+",".join(L)+"){}:{}".format(str(root),abs(root.t_inf))
