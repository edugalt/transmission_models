from random import choice,randint,random,sample,choices
from scipy.special import gamma as GAMMA
from scipy.stats import nbinom, gamma, binom, expon, poisson
import numpy as np
from transmission_models import utils
from itertools import combinations
import networkx as nx


def check_attribute_sampling(host, attribute):
    """
    Check if a host has a sampled attribute.

    Parameters
    ----------
    host : object
        The host to check.
    attribute : str
        The attribute to check for sampling.

    Returns
    -------
    bool
        True if the attribute is sampled, False otherwise.
    """
    return  host.sampled and np.isnan(attribute[int(host),int(host)])

def search_partial_sampled_siblings(host, T):
    """
    Search for partially sampled siblings of a host in the transmission tree.

    Parameters
    ----------
    host : object
        The host to search siblings for.
    T : object
        The transmission tree.

    Returns
    -------
    list
        List of partially sampled siblings.
    """

    sampled_hosts = []
    for h in T.successors(host):
        if check_attribute_sampling(host, T.distance_matrix):
            sampled_hosts.append(h)
        else:
            sampled_hosts += search_partial_sampled_siblings(host, T)

    return sampled_hosts


def search_partial_sampled_parent(host, T, root):
    """
    Search for the partially sampled parent of a host in the transmission tree.

    Parameters
    ----------
    host : object
        The host to search the parent for.
    T : object
        The transmission tree.
    root : object
        The root of the transmission tree.

    Returns
    -------
    object or None
        The partially sampled parent, or None if not found.
    """

    for h in T.predecessors(host):
        if check_attribute_sampling(host, T.distance_matrix):
            return h
        else:
            return search_partial_sampled_parent(host, T, root)
