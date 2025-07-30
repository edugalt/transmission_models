"""
Utilities Module.

This module contains utility functions for tree manipulation, visualization,
and data conversion in transmission network analysis.

Main Functions
--------------
tree_to_newick : Convert transmission tree to Newick format
search_firsts_sampled_siblings : Find first sampled siblings in tree
search_first_sampled_parent : Find first sampled parent in tree
plot_transmision_network : Visualize transmission network
tree_to_json : Convert tree to JSON format
json_to_tree : Convert JSON to tree format
tree_slicing_step : Tree topology manipulation functions

Visualization
-------------
hierarchy_pos : Generate hierarchical layout positions
hierarchy_pos_times : Generate time-based hierarchical layout
plot_transmision_network : Plot transmission network with various options

Data Conversion
---------------
tree_to_newick : Convert to Newick format for phylogenetic software
tree_to_json : Convert to JSON for data storage
json_to_tree : Convert from JSON back to tree structure
"""

# from transmission_models.classes import *
from transmission_models.classes.didelot_unsampled import *
# from transmission_models.utils import *
from transmission_models.classes import host

# Import topology movement functions
from .topology_movements import *

import scipy.special as sc
import scipy.stats as st

import numpy as np
import json

def tree_to_newick(g, lengths=True, root=None):
    """
    Convert a transmission tree to Newick format for phylogenetic software.

    Parameters
    ----------
    g : networkx.DiGraph
        The transmission tree as a directed graph.
    lengths : bool, optional
        Whether to include branch lengths in the Newick string. Default is True.
    root : node, optional
        The root node of the tree. If None, the root will be inferred.

    Returns
    -------
    str
        The Newick string representation of the tree.
    """
    if root is None:
        roots = list(filter(lambda p: p[1] == 0, g.in_degree()))
        assert 1 == len(roots)
        root = roots[0][0]
    subgs = []
    for child in g[root]:
        if len(g[child]) > 0:
            subgs.append(tree_to_newick(g, root=child, lengths=lengths))
        else:
            subgs.append(child)
    L = []
    for h in subgs:
        if isinstance(h, str):
            L.append(h)
        else:
            if lengths:
                L.append(str(h) + ":{}".format(h.t_inf - root.t_inf))
            else:
                L.append(str(h))
    if lengths:
        return "(" + ",".join(L) + "){}:{}".format(str(root), abs(root.t_inf))
    else:
        return "(" + ",".join(L) + "){}".format(str(root))


def pdf_in_between(model, Dt, t):
    """
    Compute the probability density function (PDF) value for the infection time between two events
    using a beta distribution parameterized by the model's infection parameters.

    Parameters
    ----------
    model : object
        The model object containing infection parameters (expects attributes k_inf).
    Dt : float
        The scale parameter (duration between events).
    t : float
        The time at which to evaluate the PDF.

    Returns
    -------
    float
        The PDF value at time t.
    """
    return st.beta(a=model.k_inf, b=model.k_inf, scale=Dt).pdf(t)

def sample_in_between(model, Dt):
    """
    Sample a random infection time between two events using a beta distribution
    parameterized by the model's infection parameters.

    Parameters
    ----------
    model : object
        The model object containing infection parameters (expects attributes k_inf).
    Dt : float
        The scale parameter (duration between events).

    Returns
    -------
    float
        A random sample from the beta distribution.
    """
    return st.beta(a=model.k_inf, b=model.k_inf, scale=Dt).rvs()

def random_combination(iterable, r=1):
    """
    Randomly select a combination of r elements from the given iterable.

    Parameters
    ----------
    iterable : iterable
        The input iterable to select elements from.
    r : int, optional
        The number of elements to select. Default is 1.

    Returns
    -------
    tuple
        A tuple containing r randomly selected elements from the iterable.

    Notes
    -----
    This function is equivalent to a random selection from itertools.combinations(iterable, r).
    """
    pool = tuple(iterable)
    n = len(pool)
    indices = sorted(sample(range(n), r))
    return tuple(pool[i] for i in indices)

def search_firsts_sampled_siblings(host, T):
    """
    Search the firsts sampled siblings of a host in the transmission tree.

    Parameters
    ----------
    host : Host
        The host to search the siblings from.
    T : nx.DiGraph
        The transmission tree where the host belongs to.

    Returns
    -------
    list
        The list of the firsts sampled siblings of the host.
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

    Parameters
    ----------
    host : Host
        The host to search the parent from. If the host is the root of the tree,
        it returns None.
    T : nx.DiGraph
        The transmission tree where the host belongs to.
    root : Host
        The root of the transmission tree.

    Returns
    -------
    Host or None
        The first sampled parent of the host, or None if host is the root.
    """

    if host == root:
        return None

    parent = next(T.predecessors(host))

    if not parent.sampled:
        return search_first_sampled_parent(parent, T, root)
    else:
        return parent

def Delta_log_gamma(Dt_ini, Dt_end, k, theta):
    """
    Compute the log likelihood of the gamma distribution for the time between two events.

    Parameters
    ----------
    Dt_ini : float
        Initial time.
    Dt_end : float
        End time.
    k : float
        Shape parameter of the gamma distribution.
    theta : float
        Scale parameter of the gamma distribution.

    Returns
    -------
    float
        Difference of the log likelihood of the gamma distribution.
    """
    return (k-1)*np.log(Dt_end/Dt_ini) - ((Dt_end-Dt_ini)/theta)

def hierarchy_pos(G, root=None, width=1., vert_gap=0.2, vert_loc=0, xcenter=0.5):
    """
    Compute hierarchical layout positions for a tree graph.

    Parameters
    ----------
    G : networkx.Graph
        The graph (must be a tree).
    root : node, optional
        The root node of the current branch. If None, the root will be found automatically.
    width : float, optional
        Horizontal space allocated for this branch. Default is 1.0.
    vert_gap : float, optional
        Gap between levels of hierarchy. Default is 0.2.
    vert_loc : float, optional
        Vertical location of root. Default is 0.
    xcenter : float, optional
        Horizontal location of root. Default is 0.5.

    Returns
    -------
    dict
        A dictionary of positions keyed by node.

    Notes
    -----
    This function is adapted from Joel's answer at https://stackoverflow.com/a/29597209/2966723.
    Licensed under Creative Commons Attribution-Share Alike.
    """
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
    """
    Compute hierarchical layout positions for a tree graph, using time as vertical position.

    Parameters
    ----------
    G : networkx.Graph
        The graph (must be a tree).
    root : node, optional
        The root node of the current branch. If None, the root will be found automatically.
    width : float, optional
        Horizontal space allocated for this branch. Default is 1.0.
    vert_gap : float, optional
        Gap between levels of hierarchy. Default is 0.2.
    vert_loc : float, optional
        Vertical location of root. Default is 0.
    xcenter : float, optional
        Horizontal location of root. Default is 0.5.

    Returns
    -------
    dict
        A dictionary of positions keyed by node.

    Notes
    -----
    This function is adapted from Joel's answer at https://stackoverflow.com/a/29597209/2966723.
    Licensed under Creative Commons Attribution-Share Alike.
    """
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
                                     vert_loc=(pos[root][1]+(child.t_inf-root.t_inf)), xcenter=nextx,
                                     pos=pos, parent=root)
        return pos

    return _hierarchy_pos(G, root, width, vert_gap, vert_loc, xcenter)



def plot_transmision_network(T, nodes_labels=False, pos=None, highlighted_nodes=None, ax=None, to_frame=False, title=None, filename=None, show=True):
    """
    Visualize a transmission network using matplotlib and networkx.

    Parameters
    ----------
    T : networkx.DiGraph
        The transmission network to plot.
    nodes_labels : bool, optional
        Whether to display node labels. Default is False.
    pos : dict, optional
        Node positions for layout. If None, uses graphviz_layout. Default is None.
    highlighted_nodes : list, optional
        List of nodes to highlight. Default is None.
    ax : matplotlib.axes.Axes, optional
        The axes to plot on. If None, uses current axes. Default is None.
    to_frame : bool, optional
        If True, saves the plot to a temporary image and returns it as an array. Default is False.
    title : str, optional
        Title for the plot. Default is None.
    filename : str, optional
        If provided, saves the plot to this file. Default is None.
    show : bool, optional
        Whether to display the plot. Default is True.

    Returns
    -------
    image : ndarray or None
        The image array if to_frame is True, otherwise None.
    """
    if pos is None:
        pos = graphviz_layout(T, prog="dot")
    colors = ["red" if not h.sampled else "blue" for h in T]
    ColorLegend = {'tested': 2, 'not tested': 1}

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

    nx.draw(T, pos, with_labels=nodes_labels, node_color=colors, edgecolors=edgecolors, linewidths=linewidths, ax=ax)
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='tested', markerfacecolor='b', markersize=15),
        Line2D([0], [0], marker='o', color='w', label='not tested', markerfacecolor='r', markersize=15),
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
    """
    Convert a host and its descendants to a nested dictionary suitable for JSON export.

    Parameters
    ----------
    model : object
        The transmission model containing the tree structure.
    h : host
        The host node to convert.

    Returns
    -------
    dict
        A nested dictionary representing the host and its descendants.
    """
    try:
        dict_tree = {"name": str(h),
                     "index": int(h),
                     "Infection time": h.t_inf,
                     "Sampled": h.sampled,
                     }
    except TypeError as e:
        print("error", h, str(h))
        raise e
    if h.dict_attributes != {}:
        dict_tree["Attributes"] = h.dict_attributes
    if h.sampled:
        try:
            dict_tree["Sampling time"] = h.t_sample
        except AttributeError:
            print("error", str(h), int(h), h.sampled)
    if h == model.root_host:
        dict_tree["root"] = True
    else:
        dict_tree["root"] = False
    if model.out_degree(h) > 0:
        dict_tree["children"] = [tree_to_dict(model, h2) for h2 in model.T[h]]

    return dict_tree

def cast_types(value, types_map):
    """
    Recursively cast types in a nested data structure.

    This function recursively traverses a nested data structure (dict, list)
    and casts any values that match the types in types_map to their target types.
    Useful for fixing JSON serialization issues with numpy types.

    Parameters
    ----------
    value : any
        The value to cast. Can be a dict, list, or any other type.
    types_map : list of tuples
        List of (from_type, to_type) tuples specifying type conversions.

    Returns
    -------
    any
        The value with types cast according to types_map.

    Examples
    --------
    >>> import numpy as np
    >>> import json
    >>> data = [np.int64(123)]
    >>> data = cast_types(data, [(np.int64, int), (np.float64, float)])
    >>> data_json = json.dumps(data)
    >>> data_json == "[123]"
    True

    Notes
    -----
    This function is useful for fixing "TypeError: Object of type int64 is not
    JSON serializable" errors when working with numpy arrays and JSON.
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
    """
    Save a transmission model and its tree to a JSON file.

    Parameters
    ----------
    model : object
        The transmission model to export.
    filename : str
        The path to the output JSON file.
    """
    dict_tree = tree_to_dict(model, model.root_host)
    dict_model = {"log_likelihood": model.log_likelihood,
                  "tree": dict_tree}

    # Genetic prior
    if model.genetic_prior is not None:
        dict_model["genetic_prior"] = model.genetic_log_prior

    dict_model["parameters"] = {"sampling_params": {},
                                 "offspring_params": {},
                                 "infection_params": {}
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
    """
    Create a host object from a dictionary representation (as used in JSON trees).

    Parameters
    ----------
    dict_tree : dict
        The dictionary representing a host (from JSON).

    Returns
    -------
    host
        The reconstructed host object.
    """
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
    """
    Recursively read a tree dictionary and extract edges as (parent, child) tuples.

    Parameters
    ----------
    dict_tree : dict
        The dictionary representing the tree (from JSON).
    h1 : host, optional
        The parent host node. If None, will be created from dict_tree.
    edge_list : list, optional
        The list to append edges to. Default is an empty list.

    Returns
    -------
    list
        A list of (parent, child) edge tuples.
    """
    if h1 is None: h1 = get_host_from_dict(dict_tree)

    if "children" in dict_tree:
        for d2 in dict_tree["children"]:
            h2 = get_host_from_dict(d2)

            # print(h2,d2['Infection time'],h2.t_inf)
            edge_list.append((h1, h2))
            read_tree_dict(d2, h1=h2, edge_list=edge_list)

    return edge_list

def json_to_tree(filename, sampling_params=None, offspring_params=None, infection_params=None):
    """
    Load a transmission model from a JSON file and reconstruct the model object.

    Parameters
    ----------
    filename : str
        Path to the JSON file.
    sampling_params : dict, optional
        Sampling parameters to override those in the file. Default is None.
    offspring_params : dict, optional
        Offspring parameters to override those in the file. Default is None.
    infection_params : dict, optional
        Infection parameters to override those in the file. Default is None.

    Returns
    -------
    didelot_unsampled
        The reconstructed transmission model.
    """
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
    
    if 'tree' in dict_tree:
        edge_list = read_tree_dict(dict_tree['tree'])
        model.T = nx.DiGraph()
        model.T.add_edges_from(edge_list)
        
        # Find root host
        roots = [h for h in model.T if model.T.in_degree(h) == 0]
        if roots:
            model.root_host = roots[0]
    
    return model


def build_infection_based_network(model, host_list):
    """
    Generate a transmission tree network given a list of sampled hosts.
    
    This function creates a transmission tree from the dataset. It uses 
    the model's sampling and infection parameters
    to construct a plausible initial transmission network.

    For each host, we get a number of infected hosts and then we toss a 
    coin to each host to see if they are connected given the infection time.

    At the end, we add a virtual root host to connect all disconnected components.
    
    Parameters
    ----------
    model : didelot_unsampled
        The transmission model with sampling and infection parameters.
    host_list : list
        List of host objects representing the sampled data.
        
    Returns
    -------
    model: The updated model with the generated transmission tree
        
    Notes
    -----
    This function implements the algorithm described in the notebook for
    generating initial transmission networks. It creates a directed graph
    representing the transmission tree and adds a virtual root host to
    connect all disconnected components.
    """
    import networkx as nx
    
    host_list = sorted(host_list, key=lambda host_obj: host_obj.t_inf)

    # Generate a transmission tree with our dataset to initialize the MCMC
    edge_list = []
    K_dict = {h: model.samp_offspring() for h in host_list}
    
    for i, h1 in enumerate(host_list[::-1]):
        P = 0
        event = False
        if i < 20:
            s = 0
        else:
            s = i
        K = model.samp_offspring()
        if K > 0:
            for h2 in host_list[s::-1]:
                if int(h1) == int(h2):
                    continue
                if K_dict[h2] == 0:
                    continue
                if h1.t_inf > h2.t_inf:
                    P_new = model.pdf_infection(h1.t_inf - h2.t_inf)
                    
                    if P < P_new:
                        P = P_new
                        h2_infec = h2
                        K_dict[h2] -= 1
                        event = True
            if event:
                edge_list.append((h2_infec, h1))
    
    # Create the transmission tree
    T = nx.DiGraph()
    T.add_edges_from(edge_list)
    
    # Find roots and create virtual root host
    roots = [h for h in T if T.in_degree(h) == 0]
    if roots:
        t_min = min(roots, key=lambda h: h.t_inf).t_inf
    else:
        # If no roots found, use minimum infection time from host_list
        t_min = min(h.t_inf for h in host_list)
    
    root_host = host("Virtual_host", -1, t_inf=t_min - 1.5 * model.Delta_crit)
    
    # Connect roots to virtual host
    for h in roots:
        T.add_edge(root_host, h)
    
    # Connect any disconnected hosts to virtual host
    for h in [h for h in host_list if h not in T]:
        T.add_edge(root_host, h)
    
    return T

