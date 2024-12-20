from random import choice,randint,random,sample,choices
from scipy.special import gamma as GAMMA
from scipy.stats import nbinom, gamma, binom, expon, poisson
import numpy as np
import transmission_models.utils as utils
from itertools import combinations
import networkx as nx


def check_attribute_sampling(host,distance_matrix):
    """
    Check if the host has been sampled with an attribute
    Args:
        host:
        distance_matrix:

    Returns:

    """
    return  host.sampled and np.isnan(distance_matrix[int(host),int(host)])

def search_partial_sampled_siblings(selected_host, T, distance_matrix):
    """
    Search for the siblings of a host who have been sampled with an attribute given the distance matrix
    Args:
        T:
        selected_host:
        distance_matrix:

    Returns:

    """

    sampled_hosts = []
    for h in T.successors(selected_host):
        if check_attribute_sampling(selected_host,distance_matrix):
            sampled_hosts.append(h)
        else:
            sampled_hosts += search_partial_sampled_siblings(selected_host, T, distance_matrix)

    return sampled_hosts


def search_partial_sampled_parent(selected_host, T, distance_matrix):
    """
    Search for the ancestor of the selected_host who has been sampled with an attribute given the distance matrix
    Args:
        T:
        selected_host:
        distance_matrix:

    Returns:

    """

    for h in T.predecessors(selected_host):
        if check_attribute_sampling(selected_host,distance_matrix):
            return h
        else:
            return search_partial_sampled_parent(selected_host, T, distance_matrix)
