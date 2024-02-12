import random

import networkx as nx
import numpy as np


def tree_slicing_to_parent(model, selected_host = None, parent = None, forced=False, verbose=False):
    candidates = [h for h in model.T.nodes()
                  if h != model.root_host and model.root_host not in model.T.predecessors(h)]
    #Selecting node and its relatives
    if selected_host is None:
        #If there are no candidates, we slice to sibling and take into account the new ratio of proposals
        if len(candidates) == 0:
            T_new,gg,selected_host,parent,grandparent = tree_slicing_to_sibling(model, forced=True, verbose=verbose)
            # gg = 2*gg
            return T_new,gg,selected_host,parent,grandparent
        try:
            selected_host = random.sample(candidates, 1)[0]
        except:
            raise ValueError("Merda!!! {}".format(candidates))
#     if selected_host is None:
#         parent = model.root_host
#         selected_host = list(model.T.successors(model.root_host))[0]
#         while selected_host == model.root_host or parent == model.root_host:
#             selected_host = sample(list(model.T.nodes()),1)[0]
#             parent = list(model.T.predecessors(selected_host))[0]
#     print("after",int(selected_host))
    if parent is None:
        parent = list(model.T.predecessors(selected_host))[0]
    grandparent = list(model.T.predecessors(parent))[0]

    #candidates to slice to sibling in new network
    N_new = len([h for h in model.T if h != model.root_host and
                 model.out_degree(model.parent(h)) != 1
              ])
    # if verbose: print("before N_new",N_new,"N_candidates_to_chain",model.N_candidates_to_chain,model.N_candidates_to_chain==N_new)
    ##Checking if parent and grandparent are candidates in the new network for slicing to chain
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
    #Ratio of proposal
    gg = len(candidates)/(N_new*model.T.out_degree(grandparent))


    #Creating new network
    T_new = nx.DiGraph(model.T)
    T_new.remove_edge(parent,selected_host)
    T_new.add_edge(grandparent,selected_host)

    # if verbose: print("len(T_new)-T_new.out_degree(model.root_host) -1 == 0:",len(T_new)-T_new.out_degree(model.root_host) -1 == 0)
    # if verbose: print("forced",forced , "model.N_candidates_to_chain_old==0", model.N_candidates_to_chain_old == 0)
    # if verbose: print("forced or ddd",forced or model.N_candidates_to_chain_old == 0)
    # if verbose: print("N_new==0",N_new==0)

    # Correcting ratio of proposals (gg) and number of candidates (N_new)
    ## If this movement is forced because the other is not possible, we need to take into account the new ratio of proposals
    if len(T_new)-T_new.out_degree(model.root_host) -1 == 0:
        gg = 2*gg
    #If this movement is forced because the other is not possible, we need to take into account the new ratio of proposals
    if forced or model.N_candidates_to_chain_old == 0:
        gg = 0.5*gg
#     LL_new = model.log_likelihood_transmission(T_new)
    if verbose:
        print(f"slicing node to be parent: Selected host: {selected_host}, Parent: {parent}, Grandparent:{grandparent}")
        print(f"\tgg: {gg}, N_new: {N_new}, N_new2: {model.N_candidates_to_chain}, k_out_grandparent: {model.out_degree(grandparent)}, Num candidates: {len(candidates)}")
    return T_new,gg,selected_host,parent,grandparent


def tree_slicing_to_sibling(model, selected_host=None, parent=None, selected_sibling=None, forced=False, verbose=False):
    """

    Parameters
    ----------
    forced
    model : transmission_models.models.didelot_unsampled.didelot_unsampled
    """
    candidates = [h for h in model.T.nodes()
                  if h != model.root_host and
                  (model.T.out_degree(list(model.T.predecessors(h))[0]) >= 2 or False)
                  ]

    # Selecting node and its relatives
    if selected_host is None:
        #If there are no candidates, we slice to sibling and take into account the new ratio of proposals
        if len(candidates) == 0:
            T_new,gg,selected_host,parent,grandparent =  tree_slicing_to_parent(model, forced=True, verbose=verbose)
            # gg = 2*gg
            return T_new,gg,selected_host,parent,grandparent
        try:
            selected_host = random.sample(candidates, 1)[0]
        except:
            raise ValueError("Merda!!! {}".format(candidates))

    if parent is None:
        parent = list(model.T.predecessors(selected_host))[0]

    siblings = list(model.T.successors(parent))
    if selected_sibling is None:
        siblings.remove(selected_host)
        selected_sibling = random.sample(siblings, 1)[0]


    # Number of candidates of the new network
    N_new = len(model.T) - model.T.out_degree(model.root_host)-1

    if selected_host in model.T.successors(model.root_host):
        N_new += 1

    # Ratio of proposal
    gg =  (len(candidates)*(model.out_degree(parent)-1)) / N_new



    #Creating new network
    T_new = nx.DiGraph(model.T)
    T_new.remove_edge(parent, selected_host)
    T_new.add_edge(selected_sibling, selected_host)

    # Correcting ratio of proposals (gg) and number of candidates (N_new)
    ## If this movement is forced because the other is not possible, we need to take into account the new ratio of proposals

    #Checking new parent
    N_new2 = len(candidates)
    if T_new.out_degree(selected_sibling) == 1:
        N_new2 -= 1
        model.N_candidates_to_chain -= 1
    elif T_new.out_degree(selected_sibling) == 2:
        N_new2 += 1
        model.N_candidates_to_chain += 1
    #Checinkg old parent
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

    if forced or len(model.T)-1-model.out_degree(model.root_host) == 0:
        gg = 0.5 * gg
    elif len(candidates2)==0:
        gg = 2 * gg
    if verbose:
        print(f"slicing node to be sibling: Selected host: {selected_host}, Parent: {parent}, Sibling:{selected_sibling} ")
        print(f"\tgg: {gg}, N_new: {N_new}, k_out_parent:{model.out_degree(parent)}, Num candidates: {len(candidates)}")
    return T_new, gg, selected_host, parent, selected_sibling

def tree_slicing_step(model):
    """
    Performs a tree slicing step in the transmission tree. Can be either to parent or sibling with equal probability.

    Parameters
    ----------
    model: transmission_models.models.didelot_unsampled.didelotUnsampled
        Transmission model

    """

    L_old = model.get_log_likelihood_transmission()

    if random.random() > 0.5:
        T_new, gg, selected_host, h_a, h_b = tree_slicing_to_sibling(model)

    else:
        T_new, gg, selected_host, h_a, h_b = tree_slicing_to_parent(model)

    L_new = model.log_likelihood_transmission_tree(T_new)
    pp = np.exp(L_new-L_old)
    P = gg*pp



    #Metropolis-Hastings algorithm
    if P>1:
        model.T = T_new
        model.log_likelihood = L_new
        L_old = L_new
        model.get_newick()
        model.N_candidates_to_chain_old = model.N_candidates_to_chain

    else:
        if random.random()<P:
            model.T = T_new
            model.log_likelihood = L_new
            L_old = L_new
            model.get_newick()
            model.N_candidates_to_chain_old = model.N_candidates_to_chain
