from random import random, sample

import networkx as nx
import numpy as np


def tree_slicing_to_offspring(model, selected_host=None, forced=False, verbose=False):
    """
    Slices a node reconnecting it with its grandparent. It passes from a chain model to a offspring model for the
    selected_host, its parent and grandparent.

    Parameters:
    -----------
        model: transmission_models.models.didelot_unsampled.didelotUnsampled
            model with the transmission network to apply the transformation
        selected_host: host object. Default None
            host to be sliced. If None, it is randomly selected
        forced: bool. Default False
            If True, the movement is forced because the other is not possible
        verbose: bool. Default False
            If True, it prints information about the movement

    Returns:
    --------
        T_new: nx.DiGraph
            New transmission network with the moves applied
        gg: float
            Ratio of proposals
        selected_host: host object
            Host sliced
        parent: host object
            Parent of the selected_host
        grandparent: host object
            Grandparent of the selected_host
        to_chain: bool
            If True, the proposal was to reconnect selected_host to be connected with one of its brother.
            Else, it was reconnected to its grandparent

    """
    candidates = [h for h in model.T.nodes()
                  if h != model.root_host and model.root_host not in model.T.predecessors(h)]
    # Selecting node and its relatives
    if selected_host is None:
        # If there are no candidates, we slice to sibling and take into account the new ratio of proposals
        if len(candidates) == 0:
            T_new, gg, selected_host, parent, grandparent, to_chain = tree_slicing_to_chain(model, forced=True,
                                                                                            verbose=verbose)
            # gg = 2*gg
            return T_new, gg, selected_host, parent, grandparent, to_chain

        selected_host = sample(candidates, 1)[0]

    #     if selected_host is None:
    #         parent = model.root_host
    #         selected_host = list(model.T.successors(model.root_host))[0]
    #         while selected_host == model.root_host or parent == model.root_host:
    #             selected_host = sample(list(model.T.nodes()),1)[0]
    #             parent = list(model.T.predecessors(selected_host))[0]
    #     print("after",int(selected_host))
    parent = list(model.T.predecessors(selected_host))[0]
    grandparent = list(model.T.predecessors(parent))[0]

    # candidates to slice to sibling in new network
    N_new = len([h for h in model.T if h != model.root_host and
                 model.out_degree(model.parent(h)) != 1
                 ])
    # if verbose: print("before N_new",N_new,"N_candidates_to_chain",model.N_candidates_to_chain,model.N_candidates_to_chain==N_new)
    #Checking if parent and grandparent are candidates in the new network for slicing to chain
    if model.T.out_degree(parent) == 1:
        N_new += 1
        model.N_candidates_to_chain += 1
    elif model.T.out_degree(parent) == 2:
        N_new -= 1
        model.N_candidates_to_chain -= 1
    if model.T.out_degree(grandparent) == 1:
        N_new += 1
        model.N_candidates_to_chain += 1

    # if verbose: print("N_new",N_new,"N_candidates_to_chain",model.N_candidates_to_chain,model.N_candidates_to_chain==N_new)
    # Ratio of proposal
    gg = len(candidates) / (N_new * model.T.out_degree(grandparent))

    # Creating new network
    T_new = nx.DiGraph(model.T)
    T_new.remove_edge(parent, selected_host)
    T_new.add_edge(grandparent, selected_host)

    # if verbose: print("len(T_new)-T_new.out_degree(model.root_host) -1 == 0:",len(T_new)-T_new.out_degree(model.root_host) -1 == 0)
    # if verbose: print("forced",forced , "model.N_candidates_to_chain_old==0", model.N_candidates_to_chain_old == 0)
    # if verbose: print("forced or ddd",forced or model.N_candidates_to_chain_old == 0)
    # if verbose: print("N_new==0",N_new==0)

    # Correcting ratio of proposals (gg) and number of candidates (N_new)
    # If this movement is forced because the other is not possible,
    # we need to take into account the new ratio of proposals
    if verbose:
        print(
            f"len(T_new)-T_new.out_degree(model.root_host) -1 == 0:{len(T_new) - T_new.out_degree(model.root_host) - 1 == 0}")
        print(f"forced: {forced}, model.N_candidates_to_chain_old==0: {model.N_candidates_to_chain_old == 0}")
    if len(T_new) - T_new.out_degree(model.root_host) - 1 == 0:
        gg = 2 * gg
    # If this movement is forced because the other is not possible, we need to take into account the new ratio of proposals
    if forced or model.N_candidates_to_chain_old == 0:
        gg = 0.5 * gg
    #     LL_new = model.log_likelihood_transmission(T_new)
    if verbose:
        print(f"slicing node to be parent: Selected host: {selected_host}, Parent: {parent}, Grandparent:{grandparent}")
        print(
            f"\tgg: {gg}, N_new: {N_new}, N_new2: {model.N_candidates_to_chain}, k_out_grandparent: {model.out_degree(grandparent)}, Num candidates: {len(candidates)}, Num candidates old: {model.N_candidates_to_chain_old}")
    return T_new, gg, selected_host, parent, grandparent, False


def tree_slicing_to_chain(model, selected_host=None, selected_sibling=None, forced=False, verbose=False):
    """
    Slices a node reconnecting it with one of its sibling. It passes from a offspring model to a chain model for the selected_host, its parent and the choosen sibling.

    Parameters:
    -----------
     model : transmission_models.models.didelot_unsampled.didelot_unsampled
        model with the transmission network to apply the transformation
     selected_host: host object. Default None
        host to be sliced. If None, it is randomly selected
     selected_sibling: host object. Default None
        sibling to connect the selected_host. If None, it is randomly selected
     forced: bool. Default False
        If True, the movement is forced because the other is not possible
     verbose: bool. Default False
        If True, it prints information about the movement

    Returns:
    --------
    T_new: nx.DiGraph
        New transmission network with the moves applied
    gg: float
        Ratio of proposals
    selected_host: host object
        Host sliced
    parent: host object
        Parent of the selected_host
    selected_sibling: host object
        Sibling now connected to the selected_host
    to_chain: bool
        If True, the proposal was to reconnect selected_host to be connected with one of its brother.
        Else, it was reconnected to its grandparent
    """
    candidates = [h for h in model.T.nodes()
                  if h != model.root_host and
                  (model.T.out_degree(list(model.T.predecessors(h))[0]) >= 2 or False)
                  ]

    # Selecting node and its relatives
    if selected_host is None:
        # If there are no candidates, we slice to sibling and take into account the new ratio of proposals
        if len(candidates) == 0:
            T_new, gg, selected_host, parent, grandparent, to_chain = tree_slicing_to_offspring(model, forced=True,
                                                                                                verbose=verbose)
            # gg = 2*gg
            return T_new, gg, selected_host, parent, grandparent, to_chain
        try:
            selected_host = sample(candidates, 1)[0]
        except:
            raise ValueError("Merda!!! {}".format(candidates))

    parent = list(model.T.predecessors(selected_host))[0]

    siblings = list(model.T.successors(parent))
    if selected_sibling is None:
        siblings.remove(selected_host)
        selected_sibling = sample(siblings, 1)[0]

    # Number of candidates of the new network
    N_new = len(model.T) - model.T.out_degree(model.root_host) - 1

    if selected_host in model.T.successors(model.root_host):
        N_new += 1

    # Ratio of proposal
    gg = (len(candidates) * (model.out_degree(parent) - 1)) / N_new

    # Creating new network
    T_new = nx.DiGraph(model.T)
    T_new.remove_edge(parent, selected_host)
    T_new.add_edge(selected_sibling, selected_host)

    # Correcting ratio of proposals (gg) and number of candidates (N_new)
    # If this movement is forced because the other is not possible,
    #  we need to take into account the new ratio of proposals

    # Checking new parent
    N_new2 = len(candidates)
    if T_new.out_degree(selected_sibling) == 1:
        N_new2 -= 1
        model.N_candidates_to_chain -= 1
    elif T_new.out_degree(selected_sibling) == 2:
        N_new2 += 1
        model.N_candidates_to_chain += 1
    # Checinkg old parent
    if T_new.out_degree(parent) == 0:
        N_new2 += 1
        model.N_candidates_to_chain += 1
    elif T_new.out_degree(parent) == 1:
        N_new2 -= 1
        model.N_candidates_to_chain -= 1

    candidates2 = [h for h in T_new.nodes()
                   if h != model.root_host and
                   (T_new.out_degree(list(T_new.predecessors(h))[0]) >= 2 or False)
                   ]

    if forced or len(model.T) - 1 - model.out_degree(model.root_host) == 0:
        gg = 0.5 * gg
    if len(candidates2) == 0:
        gg = 2 * gg
    if verbose:
        print(
            f"slicing node to be sibling: Selected host: {selected_host}, Parent: {parent}, Sibling:{selected_sibling}")
        print(f"\tgg: {gg}, N_new: {N_new}, k_out_parent:{model.out_degree(parent)}, Num candidates: {len(candidates)}")
    return T_new, gg, selected_host, parent, selected_sibling, True


def tree_slicing_step(model, verbose=False):
    """
    Performs a tree slicing step in the transmission tree. Can be either to parent or sibling with equal probability.

    Parameters
    ----------
    model: transmission_models.models.didelot_unsampled.didelotUnsampled
        Transmission model

    verbose: bool. Default False
        If True, it prints information about the proposal

    """

    # L_old = model.get_log_likelihood_transmission()

    if random() > 0.5:
        if verbose:
            print(f"Slicing to chain")
        T_new, gg, selected_host, h_a, h_b, to_chain = tree_slicing_to_chain(model)
    else:
        if verbose:
            print(f"Slicing to offspring")
        T_new, gg, selected_host, h_a, h_b, to_chain = tree_slicing_to_offspring(model)

    # If to_chain and a link is impossible to exists, automatically reject the proposal
    if (h_b.t_inf > selected_host.t_inf) and to_chain:
        if verbose:
            print(f"Impossible infection proposed!!: From {h_b} to {selected_host}")
        model.N_candidates_to_chain = model.N_candidates_to_chain_old
        return T_new, gg, 0, 0, selected_host, False


    if to_chain:
        Delta = model.log_likelihood_hosts_list([selected_host, h_a, h_b],
                                                T_new) - model.log_likelihood_hosts_list(
            [selected_host, h_b, h_a], model.T)
    else:
        Delta = model.log_likelihood_hosts_list([selected_host, h_a, h_b],
                                                T_new) - model.log_likelihood_hosts_list(
            [selected_host, h_b, h_a], model.T)

    # L_new = model.log_likelihood_transmission_tree(T_new)
    pp = np.exp(Delta)
    # pic_pala = L_new-L_old
    # print(f"Delta: {Delta}, A pico y pala: {L_new-L_old}, error {np.abs(Delta-pic_pala)}")

    P = gg * pp

    # Metropolis-Hastings algorithm
    if P > 1:
        accepted = True
        if verbose:
            print(f"\t-- Slicing accepted")
        model.T = T_new
        model.log_likelihood += Delta
        # L_old = L_new
        # model.get_newick()
        model.N_candidates_to_chain_old = model.N_candidates_to_chain

    else:
        if random() < P:
            accepted = True
            if verbose:
                print(f"\t-- Slicing accepted")
            model.T = T_new
            model.log_likelihood += Delta
            # L_old = L_new
            # model.get_newick()
            model.N_candidates_to_chain_old = model.N_candidates_to_chain
        else:
            accepted = False
            if verbose:
                print(f"\t-- Slicing rejected")
            model.N_candidates_to_chain = model.N_candidates_to_chain_old
    return T_new, gg, pp, P, selected_host, accepted
