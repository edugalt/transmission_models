from transmission_models.models import *
from transmission_models.models.didelot_unsampled import *
from transmission_models.utils import *

import scipy.special as sc
import scipy.stats as st

import numpy as np


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


def pdf_in_between(model,Dt,t):
    return st.beta(a=model.k_inf,b=model.k_inf,scale=Dt).pdf(t)
def sample_in_between(model,Dt):
    return st.beta(a=model.k_inf,b=model.k_inf,scale=Dt).rvs()
