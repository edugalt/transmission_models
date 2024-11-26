import numpy as np
import transmission_models.utils as utils
from transmission_models.utils import hierarchy_pos,hierarchy_pos_times,plot_transmision_network,tree_to_newick,search_firsts_sampled_siblings
from transmission_models.models import didelot_unsampled as du
from transmission_models.models.topology_movements import *


class MCMC():
    def __init__(self, model, P_rewire, P_add_remove, P_t_shift):
        """
        Initializes a new instance of the MCMC class.

        Parameters:
        - model (transmission_tree): The transmission tree model to sample from.
        - P_rewire (float): The probability of rewiring a transmission tree.
        - P_add_remove (float): The probability of adding or removing an unsampled host in the transmission tree.
        - P_t_shift (float): The probability of shifting the infection time of the host in the transmission tree.
        """

        self.model = model
        self.P_rewire = P_rewire
        self.P_add_remove = P_add_remove
        self.P_t_shift = P_t_shift

    def MCMC_iteration(self,verbose=False):
        """
        Performs an MCMC iteration on the transmission tree model.

        """
        # Randomly select a move
        move = np.random.choice(["rewire", "add_remove", "time_shift"], p=[self.P_rewire, self.P_add_remove, self.P_t_shift])


        if move or len(model.T) == 2:

            ptype = "add_remove"

            T_new, gg, pp, P, added, accepted, DL = model.add_remove_step()
            pp = P / gg
            # print(len(model.T), accepted, P)

        elif rnd_type < 2 / 3:
            ## SLICES
            ptype = "slices"
            T_new, gg, pp, P, selected_host, accepted, DL = tree_slicing_step(model, verbose=0)
            # model.get_log_likelihood_transmission()


        else:
            ptype = "time_shift"
            if verbose:
                print(f"\t-- Time shift")
            # selected_host = choice(model.unsampled_hosts)

            len_new = len(model.T)
            # nn = tree_to_newick(model.T, lengths=False, root=model.root_host)
            t_inf_new, gg, pp, P, selected_host, accepted, DL = self.model.infection_time_from_infection_model_step(
                verbose=verbose)
            # L_off = model.get_log_likelihood_transmission()
            # print(gg,pp,P,selected_host)
            # model.get_newick(lengths=False)
        # print(f"{itt}\t{len(model.T)}\t{len(model.unsampled_hosts)}\t{model.log_likelihood}\t{model.genetic_log_prior}\t{np.exp(model.log_likelihood)}\t{P}\t{ptype}\t{accepted}")
