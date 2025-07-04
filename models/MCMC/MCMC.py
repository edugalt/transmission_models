import numpy as np
import transmission_models.utils as utils
from transmission_models.utils import hierarchy_pos,hierarchy_pos_times,plot_transmision_network,tree_to_newick,search_firsts_sampled_siblings
from transmission_models.models import didelot_unsampled as du
from transmission_models.models.topology_movements import tree_slicing_step


class MCMC():
    """
    Markov Chain Monte Carlo sampler for transmission tree inference.

    This class implements MCMC sampling algorithms for transmission network
    inference using various proposal mechanisms.

    Parameters
    ----------
    model : didelot_unsampled
        The transmission tree model to sample from.
    P_rewire : float, optional
        The probability of rewiring a transmission tree. Default is 1/3.
    P_add_remove : float, optional
        The probability of adding or removing an unsampled host in the
        transmission tree. Default is 1/3.
    P_t_shift : float, optional
        The probability of shifting the infection time of the host in the
        transmission tree. Default is 1/3.
    P_add : float, optional
        The probability of adding a new host to the transmission tree once
        the add/remove have been proposed. Default is 0.5.
    P_rewire_add : float, optional
        The probability of rewiring the new unsampled host once the add
        have been proposed. Default is 0.5.
    P_offspring_add : float, optional
        The probability that the new unsampled host is an offspring once
        the add and rewire have been proposed. Default is 0.5.
    P_to_offspring : float, optional
        The probability of moving to offspring model during rewiring.
        Default is 0.5.

    Attributes
    ----------
    model : didelot_unsampled
        The transmission model being sampled.
    P_rewire : float
        Probability of rewiring moves.
    P_add_remove : float
        Probability of add/remove moves.
    P_t_shift : float
        Probability of time shift moves.
    P_add : float
        Probability of adding vs removing hosts.
    P_rewire_add : float
        Probability of rewiring added hosts.
    P_offspring_add : float
        Probability of offspring vs chain model for added hosts.
    P_to_offspring : float
        Probability of moving to offspring model.
    """

    def __init__(self, model, P_rewire=1/3, P_add_remove=1/3, P_t_shift=1/3, P_add=0.5, P_rewire_add=0.5, P_offspring_add=0.5, P_to_offspring=0.5):
        """
        Initialize the MCMC sampler.

        Parameters
        ----------
        model : didelot_unsampled
            The transmission tree model to sample from.
        P_rewire : float, optional
            The probability of rewiring a transmission tree. Default is 1/3.
        P_add_remove : float, optional
            The probability of adding or removing an unsampled host in the
            transmission tree. Default is 1/3.
        P_t_shift : float, optional
            The probability of shifting the infection time of the host in the
            transmission tree. Default is 1/3.
        P_add : float, optional
            The probability of adding a new host to the transmission tree once
            the add/remove have been proposed. Default is 0.5.
        P_rewire_add : float, optional
            The probability of rewiring the new unsampled host once the add
            have been proposed. Default is 0.5.
        P_offspring_add : float, optional
            The probability that the new unsampled host is an offspring once
            the add and rewire have been proposed. Default is 0.5.
        P_to_offspring : float, optional
            The probability of moving to offspring model during rewiring.
            Default is 0.5.
        """

        self.model = model
        self.P_rewire = P_rewire
        self.P_add_remove = P_add_remove
        self.P_t_shift = P_t_shift
        self.P_add = P_add
        self.P_rewire_add = P_rewire_add
        self.P_offspring_add = P_offspring_add
        self.P_to_offspring = P_to_offspring

    def MCMC_iteration(self, verbose=False):
        """
        Perform an MCMC iteration on the transmission tree model.

        Parameters
        ----------
        verbose : bool, optional
            Whether to print the progress of the MCMC iteration. Default is False.

        Returns
        -------
        tuple
            A tuple containing:

            - move : str
                The type of move proposed ('rewire', 'add_remove', or 'time_shift').
            - gg : float
                The ratio of proposal probabilities.
            - pp : float
                The ratio of posterior probabilities.
            - P : float
                The acceptance probability.
            - accepted : bool
                Whether the move was accepted.
            - DL : float
                The difference in log likelihood.

        Notes
        -----
        The function operates as follows:

        1. Selects a move type at random.
        2. Performs the move and computes acceptance probability.
        3. Returns move details and acceptance status.
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
