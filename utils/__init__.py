# from transmission_models.models import *
from transmission_models.models.didelot_unsampled import *
# from transmission_models.utils import *
from transmission_models import host

import scipy.special as sc
import scipy.stats as st

import numpy as np
import json

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

def search_firsts_sampled_siblings(host, T):
    """
    Search the firsts sampled siblings of a host in the transmission tree.

    Parameters:
    -----------
        host: Host
            The host to search the siblings from
        T: nx.DiGraph
            The transmission tree where the hot belongs to.
    Returns:
    --------
        sampled_hosts: list of Hosts
            The list of the firsts sampled siblings of the host

    """
    sampled_hosts = []
    for h in T.successors(host):
        if h.sampled:
            sampled_hosts.append(h)
        else:
            sampled_hosts += search_firsts_sampled_siblings(h, T)

    return sampled_hosts


def search_first_sampled_parent(host, T, root):
    """
    Search the first sampled parent of a host in the transmission tree.

    Parameters:
    -----------
        host: Host
            The host to search the parent from. If the host is the root of the tree, it returns None
        T: nx.DiGraph
            The transmission tree where the hot belongs to.
        root: Host
            The root of the transmission tree
    Returns:
    --------
        parent: Host
            The first sampled parent of the host
    """

    if host == root:
        return None

    parent = next(T.predecessors(host))

    if not parent.sampled:
        return search_first_sampled_parent(parent, T, root)
    else:
        return parent

def Delta_log_gamma(Dt_ini,Dt_end,k,theta):
    """
    Compute the log likelihood of the gamma distribution for the time between two events.

    Parameters:
    -----------
        Dt_ini: float
            Initial time
        Dt_end: float
            End time
    Returns:
    --------
        float
            Difference of the log likelihood of the gamma distribution
    """
    return (k-1)*np.log(Dt_end/Dt_ini) - ((Dt_end-Dt_ini)/theta)

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
    ColorLegend = {'tested': 2,'not tested': 1}

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
        Line2D([0], [0], marker='o', color='w', label='not tested',markerfacecolor='r', markersize=15),
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


def tree_to_dict(model, h):
    # print("host",h)
    try:
        dict_tree = {"name": str(h),
                     "index": int(h),
                     "Infection time": h.t_inf,
                     "Sampled": h.sampled,
                     }
    except TypeError as e:
        print("error",h, str(h))
        raise e
    # print(h.dict_attributes,h.dict_attributes!={})
    if h.dict_attributes != {}:
        dict_tree["Attributes"] = h.dict_attributes
    if h.sampled:
        try:
            dict_tree["Sampling time"] = h.t_sample
        except AttributeError:
            print("error", str(h), int(h), h.sampled)
    if h==model.root_host:
        dict_tree["root"] = True
    else:
        dict_tree["root"] = False
    if model.out_degree(h) > 0:
        dict_tree["children"] = [tree_to_dict(model, h2) for h2 in model.T[h]]

    return dict_tree

def cast_types(value, types_map):
    """
    recurse into value and cast any np.int64 to int

    fix: TypeError: Object of type int64 is not JSON serializable

    import numpy as np
    import json
    data = [np.int64(123)]
    data = cast_types(data, [
        (np.int64, int),
        (np.float64, float),
    ])
    data_json = json.dumps(data)
    data_json == "[123]"

    https://stackoverflow.com/a/75552723/10440128
    """
    if isinstance(value, dict):
        # cast types of dict keys and values
        return {cast_types(k, types_map): cast_types(v, types_map) for k, v in value.items()}
    if isinstance(value, list):
        # cast types of list values
        return [cast_types(v, types_map) for v in value]
    for f, t in types_map:
        if isinstance(value, f):
            return t(value) # cast type of value
    return value # keep type of value

def tree_to_json(model, filename):
    dict_tree = tree_to_dict(model,model.root_host)
    dict_model = {"log_likelihood": model.log_likelihood,
                 "tree": dict_tree}

    # Genetic prior
    if model.genetic_prior is not None:
        dict_model["genetic_prior"] = model.genetic_log_prior

    dict_model["parameters"] = {"sampling_params":{},
                                "offspring_params":{},
                                "infection_params":{}
                                }

    dict_model["parameters"]["sampling_params"]["pi"] = model.pi
    dict_model["parameters"]["sampling_params"]["k_samp"] = model.k_samp
    dict_model["parameters"]["sampling_params"]["theta_samp"] = model.theta_samp

    dict_model["parameters"]["offspring_params"]["r"] = model.r  # rate of infection
    dict_model["parameters"]["offspring_params"]["p_inf"] = model.p_inf  # probability of infection
    dict_model["parameters"]["offspring_params"]["R"] = model.R  # Reproduction number

    dict_model["parameters"]["infection_params"]["k_inf"] = model.k_inf #
    dict_model["parameters"]["infection_params"]["theta_inf"] = model.theta_inf



    # Convert and write JSON object to file
    with open(filename, "w") as outfile:
        json.dump(cast_types(dict_model, [
            (np.int64, int),
            (np.float64, float),
        ]), outfile)


def get_host_from_dict(dict_tree):
    if dict_tree["Sampled"]:
        Host = host(dict_tree["name"], dict_tree["index"], t_sample=dict_tree["Sampling time"],
                         t_inf=dict_tree["Infection time"])
        if "Attributes" in dict_tree:
            for k in dict_tree["Attributes"]:
                Host.dict_attributes[k] = dict_tree["Attributes"][k]
    else:
        Host = host(dict_tree["name"], dict_tree["index"], t_inf=dict_tree["Infection time"])
    return Host


def read_tree_dict(dict_tree, h1=None, edge_list=[]):
    if h1 is None: h1 = get_host_from_dict(dict_tree)

    if "children" in dict_tree:
        for d2 in dict_tree["children"]:
            h2 = get_host_from_dict(d2)

            # print(h2,d2['Infection time'],h2.t_inf)
            edge_list.append((h1, h2))
            read_tree_dict(d2, h1=h2, edge_list=edge_list)

    return edge_list

def json_to_tree(filename,sampling_params = None,offspring_params = None,infection_params = None):
    with open(filename) as json_data:
        dict_tree = json.load(json_data)
        json_data.close()

    if sampling_params is not None:
        sampling_params = dict_tree["parameters"]["sampling_params"]
    if offspring_params is not None:
        offspring_params = dict_tree["parameters"]["offspring_params"]
    if infection_params is not None:
        infection_params = dict_tree["parameters"]["infection_params"]

    model = didelot_unsampled(sampling_params, offspring_params, infection_params)
    model.log_likelihood = dict_tree["log_likelihood"]

    return model
