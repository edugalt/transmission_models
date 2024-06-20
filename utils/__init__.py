# from transmission_models.models import *
from transmission_models.models.didelot_unsampled import *
# from transmission_models.utils import *

import scipy.special as sc
import scipy.stats as st

import numpy as np


def tree_to_newick(g,lengths=True, root=None):
    if root is None:
        roots = list(filter(lambda p: p[1] == 0, g.in_degree()))
        assert 1 == len(roots)
        root = roots[0][0]
    subgs = []
    # print("root ",root)
    for child in g[root]:
        if len(g[child]) > 0:
            subgs.append(tree_to_newick(g, root=child,lengths=lengths))
        else:
            subgs.append(child)
#     print("------",len(subgs),subgs)
    L = []

    for h in subgs:
        if isinstance(h,str):
            L.append(h)
        else:
            if lengths:
                L.append(str(h)+":{}".format(h.t_inf-root.t_inf))
            else:
                L.append(str(h))
    if lengths:
        return "("+",".join(L)+"){}:{}".format(str(root),abs(root.t_inf))
    else:
        # print("------","("+",".join(L)+"){}".format(str(root)))
        return "("+",".join(L)+"){}".format(str(root))


def pdf_in_between(model,Dt,t):
    return st.beta(a=model.k_inf,b=model.k_inf,scale=Dt).pdf(t)
def sample_in_between(model,Dt):
    return st.beta(a=model.k_inf,b=model.k_inf,scale=Dt).rvs()
def random_combination(iterable, r=1):
    "Random selection from itertools.combinations(iterable, r)"
    pool = tuple(iterable)
    n = len(pool)
    indices = sorted(sample(range(n), r))
    return tuple(pool[i] for i in indices)


def hierarchy_pos(G, root=None, width=1., vert_gap=0.2, vert_loc=0, xcenter=0.5):
    '''
    From Joel's answer at https://stackoverflow.com/a/29597209/2966723.
    Licensed under Creative Commons Attribution-Share Alike

    If the graph is a tree this will return the positions to plot this in a
    hierarchical layout.

    G: the graph (must be a tree)

    root: the root node of current branch
    - if the tree is directed and this is not given,
      the root will be found and used
    - if the tree is directed and this is given, then
      the positions will be just for the descendants of this node.
    - if the tree is undirected and not given,
      then a random choice will be used.

    width: horizontal space allocated for this branch - avoids overlap with other branches

    vert_gap: gap between levels of hierarchy

    vert_loc: vertical location of root

    xcenter: horizontal location of root
    '''
    if not nx.is_tree(G):
        raise TypeError('cannot use hierarchy_pos on a graph that is not a tree')

    if root is None:
        if isinstance(G, nx.DiGraph):
            root = next(iter(nx.topological_sort(G)))  # allows back compatibility with nx version 1.11
        else:
            root = random.choice(list(G.nodes))

    def _hierarchy_pos(G, root, width=1., vert_gap=0.2, vert_loc=0, xcenter=0.5, pos=None, parent=None):
        '''
        see hierarchy_pos docstring for most arguments

        pos: a dict saying where all nodes go if they have been assigned
        parent: parent of this branch. - only affects it if non-directed

        '''

        if pos is None:
            pos = {root: (xcenter, vert_loc)}
        else:
            pos[root] = (xcenter, vert_loc)
        children = list(G.neighbors(root))
        if not isinstance(G, nx.DiGraph) and parent is not None:
            children.remove(parent)
        if len(children) != 0:
            dx = width / len(children)
            nextx = xcenter - width / 2 - dx / 2
            for child in children:
                nextx += dx
                pos = _hierarchy_pos(G, child, width=dx, vert_gap=vert_gap,
                                     vert_loc=vert_loc - vert_gap, xcenter=nextx,
                                     pos=pos, parent=root)
        return pos

    return _hierarchy_pos(G, root, width, vert_gap, vert_loc, xcenter)


def hierarchy_pos_times(G, root=None, width=1., vert_gap=0.2, vert_loc=0, xcenter=0.5):
    '''
    From Joel's answer at https://stackoverflow.com/a/29597209/2966723.
    Licensed under Creative Commons Attribution-Share Alike

    If the graph is a tree this will return the positions to plot this in a
    hierarchical layout.

    G: the graph (must be a tree)

    root: the root node of current branch
    - if the tree is directed and this is not given,
      the root will be found and used
    - if the tree is directed and this is given, then
      the positions will be just for the descendants of this node.
    - if the tree is undirected and not given,
      then a random choice will be used.

    width: horizontal space allocated for this branch - avoids overlap with other branches

    vert_gap: gap between levels of hierarchy

    vert_loc: vertical location of root

    xcenter: horizontal location of root
    '''
    if not nx.is_tree(G):
        raise TypeError('cannot use hierarchy_pos on a graph that is not a tree')

    if root is None:
        if isinstance(G, nx.DiGraph):
            root = next(iter(nx.topological_sort(G)))  # allows back compatibility with nx version 1.11
        else:
            root = random.choice(list(G.nodes))

    def _hierarchy_pos(G, root, width=1., vert_gap=0.2, vert_loc=0, xcenter=0.5, pos=None, parent=None):
        '''
        see hierarchy_pos docstring for most arguments

        pos: a dict saying where all nodes go if they have been assigned
        parent: parent of this branch. - only affects it if non-directed

        '''

        if pos is None:
            pos = {root: (xcenter, vert_loc)}
        else:
            pos[root] = (xcenter, vert_loc)
        children = list(G.neighbors(root))
        if not isinstance(G, nx.DiGraph) and parent is not None:
            children.remove(parent)
        if len(children) != 0:
            dx = width / len(children)
            nextx = xcenter - width / 2 - dx / 2
            for child in children:
                nextx += dx
                pos = _hierarchy_pos(G, child, width=dx, vert_gap=vert_gap,
                                     vert_loc=(child.t_inf), xcenter=nextx,
                                     pos=pos, parent=root)
        return pos

    return _hierarchy_pos(G, root, width, vert_gap, vert_loc, xcenter)



def plot_transmision_network(T,nodes_labels=False,pos = None, highlighted_nodes = None, ax=None,to_frame=False, title=None, filename=None, show=True):
    if pos is None:
        pos = graphviz_layout(T, prog="dot")
    colors = ["red" if not h.sampled else "blue" for h in T]
    ColorLegend = {'tested': 2,'no tested': 1}

    if highlighted_nodes is not None:
        edgecolors = ["black" if h not in highlighted_nodes else "green" for h in T]
        linewidths = [0.0 if h not in highlighted_nodes else 3.0 for h in T]
    else:
        edgecolors = colors
        linewidths = 1

    if ax is None:
        ax = plt.gca()

    if title is not None:
        ax.set_title(title)

    nx.draw(T, pos,with_labels=nodes_labels,node_color=colors, edgecolors=edgecolors, linewidths=linewidths,ax=ax)
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='tested',markerfacecolor='b', markersize=15),
        Line2D([0], [0], marker='o', color='w', label='no tested',markerfacecolor='r', markersize=15),
    ]

    if ax is None:
        plt.legend(handles=legend_elements, loc='upper right')
    else:
        ax.legend(handles=legend_elements, loc='upper right')


    if filename is not None:
        plt.savefig(filename)

    if to_frame:
        plt.savefig('temp.png')
        plt.close()  # Close the plot to prevent display
        # Read the saved image and return it
        return imageio.imread('temp.png')
    if show:
        plt.show()
