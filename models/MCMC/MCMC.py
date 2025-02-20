import numpy as np
import transmission_models.utils as utils
from transmission_models.utils import hierarchy_pos,hierarchy_pos_times,plot_transmision_network,tree_to_newick,search_firsts_sampled_siblings
from transmission_models.models import didelot_unsampled as du
from transmission_models.models.topology_movements import tree_slicing_step


class MCMC():
    def __init__(self, model, P_rewire=1/3, P_add_remove=1/3, P_t_shift=1/3, P_add=0.5, P_rewire_add=0.5,P_offspring_add=0.5,P_to_offspring=0.5):
        """
        Initializes a new instance of the MCMC class.

        Parameters:
        - model (transmission_tree): The transmission tree model to sample from.
        - P_rewire (float): The probability of rewiring a transmission tree.
        - P_add_remove (float): The probability of adding or removing an unsampled host in the transmission tree.
        - P_t_shift (float): The probability of shifting the infection time of the host in the transmission tree.
        - P_add (float): The probability of adding a new host to the transmission tree once the add/remove have been proposed.
        - P_rewire_add (float): The probability of rewiring the new unsampled host once the add have been proposed.
        - P_offspring_add (float): The probability that the new unsampled host is an offspring once the add and rewire have been proposed.
        """

        self.model = model
        self.P_rewire = P_rewire
        self.P_add_remove = P_add_remove
        self.P_t_shift = P_t_shift
        self.P_add = P_add
        self.P_rewire_add = P_rewire_add
        self.P_offspring_add = P_offspring_add
        self.P_to_offspring = P_to_offspring

    def MCMC_iteration(self,verbose=False):
        """
        Performs an MCMC iteration on the transmission tree model.

        Parameters:
        -----------
            verbose: bool
                Whether to print the progress of the MCMC iteration.

        Returns:
        - move (str): The type of move proposed.
        - gg (float): The ratio of proposals probabilities.
        - pp (float): The ratio of posteriors probabilities.
        - P (float): The Acceptance probability.
        - accepted (bool): Whether the move was accepted.
        - DL (float): The difference in log likelihood.
        - selected_host (Host): The host selected for the move.

        """

        self.model.get_N_candidates_to_chain()
        self.model.get_candidates_to_chain()
        self.model.N_candidates_to_chain = len(self.model.candidates_to_chain)
        self.model.N_candidates_to_chain_old = len(self.model.candidates_to_chain)

        # Randomly select a move
        move = np.random.choice(["rewire", "add_remove", "time_shift"], p=[self.P_rewire, self.P_add_remove, self.P_t_shift])


        if move == "add_remove":
            if verbose:
                print(f"\t-- Add or remove")

            T_new, gg, pp, P, added, accepted, DL = self.model.add_remove_step(P_add=self.P_add,P_rewiring=self.P_rewire_add,P_off=self.P_offspring_add,verbose=verbose)
            # pp = P / gg
        elif move == "rewire":
            if verbose:
                print(f"\t-- Rewiring")
            ## SLIDES
            T_new, gg, pp, P, selected_host, accepted, DL = tree_slicing_step(self.model,P_to_offspring = self.P_to_offspring , verbose=verbose)


        elif move == "time_shift":
            if verbose:
                print(f"\t-- Time shift")
            t_inf_new, gg, pp, P, selected_host, accepted, DL = self.model.infection_time_from_infection_model_step(verbose=verbose)

        return move,gg,pp,P,accepted,DL
