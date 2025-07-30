import json
from itertools import combinations
from math import factorial

from transmission_models import *
from transmission_models import utils
from transmission_models.classes.host import *
# from transmission_models.utils import tree_to_newick
from transmission_models.utils.topology_movements import *

from random import random, randint,choice
from scipy.stats import nbinom, gamma, binom, expon, norm
from scipy.special import gamma as GAMMA
import networkx as nx
from networkx.exception import NetworkXError

from transmission_models.classes.genetic_prior import genetic_prior_tree
from transmission_models.classes.location_prior import same_location_prior_tree


# from ..utils import tree_to_newick


class didelot_unsampled():
    """
    Didelot unsampled transmission model.

    This class implements the Didelot et al. (2017) framework for transmission
    tree inference with unsampled hosts. It provides methods for building
    transmission networks, computing likelihoods, and performing MCMC sampling.

    The model incorporates three main components:
    1. Sampling model: Gamma distribution for sampling times
    2. Offspring model: Negative binomial distribution for offspring number
    3. Infection model: Gamma distribution for infection times

    Parameters
    ----------
    sampling_params : dict
        Parameters for the sampling model containing:
        - pi : float, sampling probability
        - k_samp : float, shape parameter for gamma distribution
        - theta_samp : float, scale parameter for gamma distribution
    offspring_params : dict
        Parameters for the offspring model containing:
        - r : float, rate of infection
        - p_inf : float, probability of infection
    infection_params : dict
        Parameters for the infection model containing:
        - k_inf : float, shape parameter for gamma distribution
        - theta_inf : float, scale parameter for gamma distribution

    Attributes
    ----------
    T : networkx.DiGraph
        The transmission tree.
    host_dict : dict
        Dictionary mapping host IDs to host objects.
    log_likelihood : float
        Current log likelihood of the model.
    genetic_prior : genetic_prior_tree, optional
        Prior for genetic data.
    same_location_prior : same_location_prior_tree, optional
        Prior for location data.

    References
    ----------
    Didelot, X., Gardy, J., & Colijn, C. (2017). Bayesian inference of
    transmission chains using timing of events, contact and genetic data.
    PLoS computational biology, 13(4), e1005496.
    """

    def __init__(self, sampling_params, offspring_params, infection_params, T=None):
        """
        Initialize the Didelot unsampled transmission model.

        Parameters
        ----------
        sampling_params : dict
            Parameters for the sampling model containing:
            - pi : float, sampling probability
            - k_samp : float, shape parameter for gamma distribution
            - theta_samp : float, scale parameter for gamma distribution
        offspring_params : dict
            Parameters for the offspring model containing:
            - r : float, rate of infection
            - p_inf : float, probability of infection
        infection_params : dict
            Parameters for the infection model containing:
            - k_inf : float, shape parameter for gamma distribution
            - theta_inf : float, scale parameter for gamma distribution
        T : networkx.DiGraph, optional
            The transmission tree. If provided, the model will be initialized
            with this tree. Default is None.

        Raises
        ------
        KeyError
            If any required parameter is missing from the input dictionaries.
        """
        # Reading parameters
        try:
            self.pi = sampling_params["pi"]
            self.k_samp = sampling_params["k_samp"]
            self.theta_samp = sampling_params["theta_samp"]
        except KeyError as e:
            raise

        try:
            self.r = offspring_params["r"]  # rate of infection
            self.p_inf = offspring_params["p_inf"]  # probability of infection
            self.R = self.r * (1 - self.p_inf) / self.p_inf  # Reproduction number
        except KeyError as e:
            raise

        try:
            self.k_inf = infection_params["k_inf"]  #
            self.theta_inf = infection_params["theta_inf"]
        except KeyError as e:
            raise


        # Distributions functions of the models
        self.dist_sampling = gamma(a=self.k_samp, scale=self.theta_samp)
        self.pdf_sampling = lambda t: self.dist_sampling.pdf(t)
        self.samp_sampling = lambda: self.dist_sampling.rvs()

        self.dist_infection = gamma(a=self.k_inf, scale=self.theta_inf)
        self.pdf_infection = lambda t: self.dist_infection.pdf(t)
        self.samp_infection = lambda: self.dist_infection.rvs()

        self.dist_offspring = nbinom(self.r, self.p_inf)
        self.pmf_offspring = lambda k: self.dist_offspring.pmf(k)
        self.samp_offspring = lambda: self.dist_offspring.rvs()

        #Distribution in between

        self.pdf_infection_in_between = lambda t,Dt: (pdf_in_between(self,Dt,t))


        self.T = T
        self.G = nx.DiGraph()
        self.host_dict = {}

        self.log_likelihood = 0
        self.likelihood = 0

        self.infection_likelihood = 0
        self.sampling_likelihood = 0
        self.offspring_likelihood = 0
        self.infection_log_likelihood = 0
        self.sampling_log_likelihood = 0
        self.offspring_log_likelihood = 0

        self.newick = None
        self.root_host = None
        self.unsampled_hosts = []
        self.N_candidates_to_chain = 0
        self.N_candidates_to_chain_old = 0
        self.candidates_to_chain = []
        self.roots_subtrees = []


        self.genetic_log_prior = 0
        self.genetic_prior = None

        self.same_location_log_prior = 0
        self.same_location_prior = None


        self.Delta_crit = 4/((1-self.pi)*self.pmf_offspring(1)*((self.theta_inf)**(-self.k_inf))/(GAMMA(self.k_inf)))**(1/(self.k_inf-1))

        # def generate_networks(self):

    @property
    def T(self):
        return self._T
    
    @T.setter
    def T(self, value):
        self._T = value
        if self.T is None:
            return
        return value
    
    def set_T(self, T):
        if T is None:
            raise ValueError("T is None!!")
        self.T = T       
        # Find root host
        roots = [h for h in self.T if self.T.in_degree(h) == 0]
        if roots:
            self.root_host = roots[0]
        # Get unsampled hosts
        self.unsampled_hosts = self.get_unsampled_hosts()

        # Get candidates to chain
        self.candidates_to_chain = self.get_candidates_to_chain()

    def samp_t_inf_between(self, h1, h2):
        """
        Sample a time of infection between two hosts.

        Uses a rejection sampling method to sample the time of infection of
        the infected host using the chain model from Didelot et al. 2017.

        Parameters
        ----------
        h1 : host
            Infector host.
        h2 : host
            Infected host.

        Returns
        -------
        float
            Time of infection of the host infected by h1 and the infector of h2.

        Notes
        -----
        This method implements the rejection sampling algorithm described in
        Didelot et al. (2017) for sampling infection times in transmission chains.
        """
        # Dt = abs(h1 - h2)

        Dt = h2.t_inf-h1.t_inf
        return sample_in_between(self,Dt)

    def add_root(self, t_sampl, id="0", genetic_data=[], t_inf=0, t_sample=None):
        """
        Add the root host to the transmission tree.

        Parameters
        ----------
        t_sampl : float
            Sampling time of the root host.
        id : str, optional
            Identifier for the root host. Default is "0".
        genetic_data : list, optional
            Genetic data for the root host. Default is empty list.
        t_inf : float, optional
            Infection time of the root host. Default is 0.
        t_sample : float, optional
            Sampling time of the root host. Default is None.

        Returns
        -------
        host
            The root host object.
        """
        self.root_host = host(id, 0, genetic_data, t_inf, t_sampl)
        return self.root_host

    def successors(self, host):
        """
        Get the successors (children) of a given host in the transmission tree.

        Parameters
        ----------
        host : host
            The host node whose successors are to be returned.

        Returns
        -------
        iterator
            An iterator over the successors of the host.
        """
        return self.T.successors(host)

    def parent(self, host):
        """
        Get the parent (infector) of a given host in the transmission tree.

        Parameters
        ----------
        host : host
            The host node whose parent is to be returned.

        Returns
        -------
        host
            The parent host object.
        """
        return list(self.T.predecessors(host))[0]

    def out_degree(self, host):
        """
        Get the out-degree (number of children) of a host in the transmission tree.

        Parameters
        ----------
        host : host
            The host node whose out-degree is to be returned.

        Returns
        -------
        int
            The out-degree of the host.
        """
        return self.T.out_degree(host)

    def choose_successors(self, host, k=1):
        """
        Choose k unique successors of a given host.

        Parameters
        ----------
        host : host
            Host whose successors will be chosen.
        k : int, optional
            Number of successors to choose. Default is 1.

        Returns
        -------
        list
            List of k randomly chosen successors of the host.
        """
        return sample(list(self.successors(host)), k)

    def compute_Delta_loc_prior(self, T_new):
        """
        Compute the change in the location prior log-likelihood for a new tree.

        Parameters
        ----------
        T_new : networkx.DiGraph
            The new transmission tree.

        Returns
        -------
        tuple
            (Delta log prior, new log prior, old log prior, old correction log-likelihood)
        """
        LP_sloc_top_old = self.same_location_prior.correction_LL
        LP_sloc_old = self.same_location_prior.log_prior
        LP_sloc_new = self.same_location_prior.log_prior_T(T_new)
        DL_prior_same_location = LP_sloc_new - LP_sloc_old

        return DL_prior_same_location, LP_sloc_new, LP_sloc_old, LP_sloc_top_old

    def get_candidates_to_chain(self):
        """
        Get the list of candidate hosts for chain moves in the transmission tree.

        Returns
        -------
        list
            List of candidate host nodes for chain moves.
        """
        self.candidates_to_chain = [h for h in self.T if
                                    not h == self.root_host and not self.out_degree(self.parent(h)) == 1]

        return self.candidates_to_chain

    def get_N_candidates_to_chain(self, recompute=False):
        """
        Get the number of candidate hosts for chain moves, optionally recomputing the list.

        Parameters
        ----------
        recompute : bool, optional
            If True, recompute the list of candidates. Default is False.

        Returns
        -------
        int
            Number of candidate hosts for chain moves.
        """
        if recompute:
            self.N_candidates_to_chain = len(self.get_candidates_to_chain())
        else:
            self.N_candidates_to_chain = len(self.candidates_to_chain)

        return self.N_candidates_to_chain

    def get_root_subtrees(self):
        """
        Retrieve the root subtrees of the transmission tree.

        This method searches for the first sampled siblings of the root host
        in the transmission tree and stores them in the `roots_subtrees` attribute.

        Returns
        -------
        list
            A list of root subtrees.
        """
        self.roots_subtrees = utils.search_firsts_sampled_siblings(self.root_host, self.T)
        return self.roots_subtrees

    def get_unsampled_hosts(self):
        """
        Get the list of unsampled hosts in the transmission tree (excluding the root).

        Returns
        -------
        list
            List of unsampled host nodes.
        """
        self.unsampled_hosts = [h for h in self.T if not h.sampled and h != self.root_host]
        return self.unsampled_hosts

    def get_sampling_model_likelihood(self, hosts=None, T=None, update=False):
        """
        Compute the likelihood of the sampling model.

        Computes the likelihood of the sampling model given a list of hosts.
        If no list is given, the likelihood of the whole transmission tree
        is returned.

        Parameters
        ----------
        hosts: list of host objects

        Returns
        -------
            L: float
                The likelihood of the sampling model given the list of hosts

        """
        if T is None:
            T = self.T
        L = 1
        if hosts is not None:
            if isinstance(hosts,list):
                for h in hosts:
                    if not h.sampled:
                        L*=(1-self.pi)
                    else:
                        L*=self.pi*self.pdf_sampling(h.t_sample-h.t_inf)
            else:
                if not hosts.sampled:
                    L*=(1-self.pi)
                else:
                    L*=self.pi*self.pdf_sampling(hosts.t_sample-hosts.t_inf)
        else:
            for h in T:
                if not h.sampled:
                    L*=(1-self.pi)
                else:
                    L*=self.pi*self.pdf_sampling(h.t_sample-h.t_inf)
            if update:
                self.sampling_likelihood = L
        return L

    def get_sampling_model_log_likelihood(self,hosts=None,T=None, update=False):
        """
        Computes the likelihood of the sampling model given a list of hosts. If no list is given, the likelihood of the
        whole transmission tree is returned.

        Parameters
        ----------
        hosts: list of host objects

        Returns
        -------
            L: float
                The likelihood of the sampling model given the list of hosts

        """
        if T is None:
            T = self.T
        L = 0
        if hosts is not None:
            if isinstance(hosts,list):
                for h in hosts:
                    if not h.sampled:
                        L += np.log((1-self.pi))
                    else:
                        L += np.log(self.pi*self.pdf_sampling(h.t_sample-h.t_inf))
            else:
                if not hosts.sampled:
                    L += np.log((1-self.pi))
                else:
                    L += np.log(self.pi*self.pdf_sampling(hosts.t_sample-hosts.t_inf))
        else:
            for h in T:
                if not h.sampled:
                    L += np.log((1-self.pi))
                else:
                    L += np.log(self.pi*self.pdf_sampling(h.t_sample-h.t_inf))
            if update:
                self.sampling_log_likelihood = L
        return L

    def Delta_log_sampling(self, hosts, T_end, T_ini=None):
        """
        Compute the change in log-likelihood for the sampling model.

        Parameters
        ----------
        hosts : list
            List of host objects.
        T_end : float
            End time.
        T_ini : float, optional
            Initial time. Default is None.

        Returns
        -------
        float
            Change in log-likelihood for the sampling model.

        Notes
        -----
        The function operates as follows:

        1. Computes the log-likelihood for the sampling model at T_end.
        2. If T_ini is provided, subtracts the log-likelihood at T_ini.
        3. Returns the difference.
        """
        L_end = self.get_sampling_model_likelihood(hosts, T_end)

        if T_ini is None:
            T_ini = self.T

        L_ini = self.get_sampling_model_likelihood(hosts, T_ini)
        return np.log(L_end / L_ini)

    def get_offspring_model_likelihood(self,hosts=None,T=None, update=False):
        """
        Computes the likelihood of the offspring model given a list of hosts. If no list is given, the likelihood of the
        whole transmission tree is returned.

        Parameters
        ----------
        hosts: list of host objects

        Returns
        -------
            L: float
                The likelihood of the offspring model given the list of hosts

        """
        if T is None:
            T = self.T
        L = 1
        if hosts is not None:
            if isinstance(hosts,list):
                for h in hosts:
                    L*=self.pmf_offspring(T.out_degree(h))
            else:
                L*=self.pmf_offspring(T.out_degree(hosts))
        else:
            for h in T:
                L*=self.pmf_offspring(T.out_degree(h))
            if update:
                self.offspring_likelihood = L
        return L
    def get_offspring_model_log_likelihood(self,hosts=None,T=None, update=False):
        """
        Computes the likelihood of the offspring model given a list of hosts. If no list is given, the likelihood of the
        whole transmission tree is returned.

        Parameters
        ----------
        hosts: list of host objects

        Returns
        -------
            L: float
                The likelihood of the offspring model given the list of hosts

        """
        if T is None:
            T = self.T
        L = 0
        if hosts is not None:
            if isinstance(hosts,list):
                for h in hosts:
                    L += self.dist_offspring.logpmf(T.out_degree(h))
            else:
                L += self.dist_offspring.logpmf(T.out_degree(hosts))
        else:
            for h in T:
                L += self.dist_offspring.logpmf(T.out_degree(h))
            if update:
                self.offspring_log_likelihood = L
        return L

    def Delta_log_offspring(self, hosts, T_end, T_ini=None):
        """
        Compute the change in log-likelihood for the offspring model.

        Parameters
        ----------
        hosts : list
            List of host objects.
        T_end : float
            End time.
        T_ini : float, optional
            Initial time. Default is None.

        Returns
        -------
        float
            Change in log-likelihood for the offspring model.

        Notes
        -----
        The function operates as follows:

        1. Computes the log-likelihood for the offspring model at T_end.
        2. If T_ini is provided, subtracts the log-likelihood at T_ini.
        3. Returns the difference.
        """
        L_end = self.get_offspring_model_likelihood(hosts, T_end)

        if T_ini is None:
            T_ini = self.T

        L_ini = self.get_offspring_model_likelihood(hosts, T_ini)
        return np.log(L_end / L_ini)


    def get_infection_model_likelihood(self,hosts=None,T=None, update=False):
        """
        Computes the likelihood of the infection model given a list of hosts. If no list is given, the likelihood of the
        whole transmission tree is returned.

        Parameters
        ----------
            hosts: list of host objects

            T: DiGraph object
                Contagious tree which likelihood of the hosts will be computed. If it is None, the network of the model is used.

            update: bool
                If True, the likelihood of the infection model is updated in the model object.

        Returns
        -------
            L: float
                The likelihood of the infection model given the list of hosts


        """
        if T is None:
            T = self.T
        L = 1
        if hosts is not None:
            if isinstance(hosts,list):
                for h in hosts:
                    for j in T.successors(h):
                        L*=self.pdf_infection(j.t_inf-h.t_inf)
                        # print(f"=======>{str(h)=},{str(j)=}\t{j.t_inf-h.t_inf=}\t{self.pdf_infection(j.t_inf-h.t_inf)=}")
            else:
                for j in T.successors(hosts):
                    L*=self.pdf_infection(j.t_inf-hosts.t_inf)
        else:
            for h in T:
                for j in T.successors(h):
                    L*=self.pdf_infection(j.t_inf-h.t_inf)
            if update:
                self.infection_likelihood = L
        return L
    def get_infection_model_log_likelihood(self,hosts=None,T=None, update=False):
        """
        Computes the likelihood of the infection model given a list of hosts. If no list is given, the likelihood of the
        whole transmission tree is returned.

        Parameters
        ----------
            hosts: list of host objects

            T: DiGraph object
                Contagious tree which likelihood of the hosts will be computed. If it is None, the network of the model is used.

            update: bool
                If True, the likelihood of the infection model is updated in the model object.

        Returns
        -------
            L: float
                The likelihood of the infection model given the list of hosts


        """
        if T is None:
            T = self.T
        L = 0
        if hosts is not None:
            if isinstance(hosts,list):
                for h in hosts:
                    for j in T.successors(h):
                        L += np.log(self.pdf_infection(j.t_inf-h.t_inf))
                        # print(f"=======>{str(h)=},{str(j)=}\t{j.t_inf-h.t_inf=}\t{self.pdf_infection(j.t_inf-h.t_inf)=}")
            else:
                for j in T.successors(hosts):
                    L += np.log(self.pdf_infection(j.t_inf-hosts.t_inf))
        else:
            for h in T:
                for j in T.successors(h):
                    L += np.log(self.pdf_infection(j.t_inf-h.t_inf))
            if update:
                self.infection_log_likelihood = L
        return L

    def Delta_log_infection(self, hosts, T_end, T_ini=None):
        """
        Compute the change in log-likelihood for the infection model.

        Parameters
        ----------
        hosts : list
            List of host objects.
        T_end : float
            End time.
        T_ini : float, optional
            Initial time. Default is None.

        Returns
        -------
        float
            Change in log-likelihood for the infection model.

        Notes
        -----
        The function operates as follows:

        1. Computes the log-likelihood for the infection model at T_end.
        2. If T_ini is provided, subtracts the log-likelihood at T_ini.
        3. Returns the difference.
        """
        if T_ini is None:
            T_ini = self.T

        L_end = self.get_infection_model_likelihood(hosts, T_end)


        L_ini = self.get_infection_model_likelihood(hosts, T_ini)

        # Delta2 = (self.k_inf-1)*np.log(h.t_inf/)
        # print()

        # print(f"{Delta=},{np.log(L_end / L_ini)=}")
        return np.log(L_end / L_ini)



    def log_likelihood_host(self, host, T=None):
        """
        Computes the log likelihood of a host given the transmission tree.
        Parameters
        ----------
        host: host object

        T: DiGraph object

        Returns
        -------
            log_likelihood: float
                The log likelihood of the host in the transmission network
        """
        if T is None:
            T = self.T
        L = 0
        # Sampling model
        L += self.get_sampling_model_log_likelihood(host, T)

        # Offspring model
        L += self.get_offspring_model_log_likelihood(host, T)

        # Infection model
        L += self.get_infection_model_log_likelihood(host, T)

        return L

    def Delta_log_likelihood_host(self, hosts, T_end, T_ini=None):
        """
        Compute the change in log-likelihood for a host.

        Parameters
        ----------
        hosts : list
            List of host objects.
        T_end : float
            End time.
        T_ini : float, optional
            Initial time. Default is None.

        Returns
        -------
        float
            Change in log-likelihood for the host.

        Notes
        -----
        The function operates as follows:

        1. Computes the log-likelihood for the host at T_end.
        2. If T_ini is provided, subtracts the log-likelihood at T_ini.
        3. Returns the difference.
        """
        # print(hosts)
        L_end = self.log_likelihood_host(hosts, T_end)

        # Delta = self.Delta_log_infection(hosts, T_end, T_ini) + self.Delta_log_offspring(hosts, T_end, T_ini) + self.Delta_log_sampling(hosts, T_end, T_ini)

        if T_ini is None:
            T_ini = self.T

        L_ini = self.log_likelihood_host(hosts, T_ini)

        # print("----", hosts, L_end - L_ini, Delta)

        return L_end - L_ini

    def log_likelihood_hosts_list(self, hosts, T):
        log_likelihood = 0
        for h in hosts:
            log_likelihood += self.log_likelihood_host(h, T)
        return log_likelihood

    def log_likelihood_transmission_tree(self, T):
        log_likelihood = 0

        for h in T:
            log_likelihood += self.log_likelihood_host(h, T)
        return log_likelihood

    def show_log_likelihoods(self, hosts=None, T=None, verbose=False):
        """
        Print and return the log-likelihoods for the sampling, offspring, and infection models.

        Parameters
        ----------
        hosts : list, optional
            List of host objects to compute log-likelihoods for. If None, computes for all hosts in T.
        T : networkx.DiGraph, optional
            Transmission tree. If None, uses self.T.
        verbose : bool, optional
            If True, prints the log-likelihoods. Default is False.

        Returns
        -------
        tuple
            (LL_sampling, LL_offspring, LL_infection): Log-likelihoods for the sampling, offspring, and infection models.
        """
        if T is None:
            T = self.T

        if hosts is not None:
            LL_sampling = self.get_sampling_model_log_likelihood(hosts, T)
            LL_offspring = self.get_offspring_model_log_likelihood(hosts, T)
            LL_infection = self.get_infection_model_log_likelihood(hosts, T)
        else:
            LL_sampling = 0
            LL_offspring = 0
            LL_infection = 0
            for h in T:
                LL_sampling += self.get_sampling_model_log_likelihood(h, T)
                LL_offspring += self.get_offspring_model_log_likelihood(h, T)
                LL_infection += self.get_infection_model_log_likelihood(h, T)
        if verbose:
            print("Sampling model:", LL_sampling)
            print("Offspring model:", LL_offspring)
            print("Infection model:", LL_infection)
        return LL_sampling, LL_offspring, LL_infection

    def log_likelihood_transmission_tree_old(self, T):
        """
        Compute the log-likelihood of the entire transmission tree using the old method.

        Parameters
        ----------
        T : networkx.DiGraph
            Transmission tree to compute the log-likelihood for.

        Returns
        -------
        float
            The log-likelihood of the transmission tree.
        """
        log_likelihood = 0
        Pi = 1

        for h in T:
            # Sampling model
            if h.sampled:
                sigma = self.pdf_sampling(h.t_sample - h.t_inf)
                Pi = self.pi * (sigma)
            else:
                Pi = (1 - self.pi)

            # Offspring model
            Pi *= self.pmf_offspring(T.out_degree(h))

            # Infection model
            for j in T.successors(h):
                sigma2 = self.pdf_infection(j.t_inf - h.t_inf)
                if sigma2 == 0:
                    self.log_likelihood = -1e30
                    print("Impossible times!!!", int(h), int(j), h.t_inf, j.t_inf, j.t_inf - h.t_inf)
                    return self.log_likelihood
                Pi *= sigma2
            log_likelihood += np.log(Pi)
        return log_likelihood

    def get_log_likelihood_transmission(self):
        self.sampling_likelihood = 1
        self.infection_likelihood = 1
        self.offspring_likelihood = 1
        self.likelihood = 1
        self.log_likelihood = 0
        self.sampling_log_likelihood = 0
        self.infection_log_likelihood = 0
        self.offspring_log_likelihood = 0

        for i,h in enumerate(self.T):
            self.sampling_log_likelihood += self.get_sampling_model_log_likelihood(h)
            self.infection_log_likelihood += self.get_infection_model_log_likelihood(h)
            self.offspring_log_likelihood += self.get_offspring_model_log_likelihood(h)

        self.log_likelihood += self.sampling_log_likelihood + self.infection_log_likelihood + self.offspring_log_likelihood
            # print(i,h,sampling_likelihood,offspring_likelihood,infection_likelihood,self.likelihood)

        self.likelihood = np.exp(self.log_likelihood)

        self.sampling_likelihood = np.exp(self.sampling_log_likelihood)
        self.infection_likelihood = np.exp(self.infection_log_likelihood)
        self.offspring_likelihood = np.exp(self.offspring_log_likelihood)

        return self.log_likelihood


    def add_genetic_prior(self,mu_gen, gen_dist):
        """
        Adds a genetic prior to the model that computes the likelihood that two sampled hosts has a relationship given
        the genetic distance of the virus of the hosts.
        Two nodes are considered that has a relationship if the only hosts that are on they are connected through
        unsampled hosts.

        Parameters
        ----------
        mu_gen: float
            Mutation rate
        gen_dist: np.array
            Genetic distance matrix of the virus of the hosts. The index has to be identical to the index of the hosts.

        """

        self.genetic_prior = genetic_prior_tree(self, mu_gen, gen_dist)
        self.genetic_log_prior = self.genetic_prior.log_prior_T(self.T)


    def add_same_location_prior(self, P_NM, tau, loc_dist):
        """
        Adds a genetic prior to the model that computes the likelihood that two sampled hosts has a relationship given
        the genetic distance of the virus of the hosts.
        Two nodes are considered that has a relationship if the only hosts that are on they are connected through
        unsampled hosts.

        Parameters
        ----------
        log_K: float
            Log probability of two hosts not being in the same location
        gen_dist: np.array
            Genetic distance matrix of the virus of the hosts. The index has to be identical to the index of the hosts.

        """

        self.same_location_prior = same_location_prior_tree(self,P_NM, tau,loc_dist)
        self.same_location_log_prior = self.same_location_prior.log_prior_T(self.T)


    def create_transmision_phylogeny_nets(self, N, mu, P_mut):
        """
            N: Number of hosts
            mu: Mutation rate
            P_mut: Prob of mutation
        """
        n = 1  # Counting hosts for the index
        self.T = nx.DiGraph()
        self.G = nx.DiGraph()

        # Adding first host
        self.T.add_node(self.root_host)
        self.G.add_node(self.root_host.get_genetic_str())
        new_infected = [self.root_host]

        genes = [self.root_host.get_genetic_str()]

        while True:

            last_infected = list(new_infected)
            new_infected = []
            #         print(last_infected)
            for h in last_infected:
                # Generating k infections for host h
                k = nbinom.rvs(self.r, self.p_inf, size=1)[0]
                if k == 0:
                    new_infected = list(last_infected)

                for i in range(k):
                    # sampled with probability Pi
                    rnd = random()
                    if rnd < self.pi:
                        sampled = True
                    else:
                        sampled = False
                        sample_time = None

                    # generate mutations
                    genetic_data = binom_mutation(100, self.p_inf, h.genetic_data)
                    genes.append(genetic_data)

                    # Infections time from gamma distribution
                    t_inf = np.random.gamma(shape=self.k_inf, scale=self.theta_inf, size=1)[0]
                    if sampled:
                        t_u_sampl = np.random.gamma(shape=self.k_samp, scale=self.theta_samp, size=1)[0]
                        sample_time = h.t_inf + t_inf + t_u_sampl
                    #                     print(i+n,h.t_inf+t_inf,t_u_sampl,sample_time,gamma.pdf(t_u_sampl,k,0),gamma.pdf(sample_time,k,h.t_inf+t_inf))
                    # print("host",i+n)
                    new_infected.append(
                        host(str(i + n), i + n, genetic_data, t_inf=h.t_inf + t_inf, t_sample=sample_time))

                    # Adding edges
                    self.T.add_edge(h, new_infected[-1])
                    if h.genetic_data == new_infected[-1].genetic_data: continue

                    ##Adding edges attributes to phylogeny
                    self.G.add_edge(h.get_genetic_str(), new_infected[-1].get_genetic_str(),
                                    hosts=[h, new_infected[-1]])
                n += k
            if last_infected == []: last_infected = [self.root_host]
            # print(k,n,last_infected)
            if n > N: break

        self.host_dict = {int(str(h)): h for h in self.T}

        return self.G, self.T, self.host_dict

    def get_newick(self,lengths=True):
        self.newick = utils.tree_to_newick(self.T, root=self.root_host,lengths=lengths)

        return self.newick

    def save_json(self, filename):
        """
        Save the transmission tree to a JSON file.

        Parameters
        ----------
        filename : str
            Path to the output JSON file.
        """
        utils.tree_to_json(self, filename)

    @classmethod
    def json_to_tree(cls, filename, sampling_params=None, offspring_params=None, infection_params=None):
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
        edge_list = []
        with open(filename, "r") as json_data:
            dict_tree = json.load(json_data)
            json_data.close()
            if sampling_params is None:
                sampling_params = dict_tree["parameters"]["sampling_params"]
            if offspring_params is None:
                offspring_params = dict_tree["parameters"]["offspring_params"]
            if infection_params is None:
                infection_params = dict_tree["parameters"]["infection_params"]

            model = cls(sampling_params, offspring_params, infection_params)
            model.log_likelihood = dict_tree["log_likelihood"]
            model.T = nx.DiGraph()
            model.root_host = utils.get_host_from_dict(dict_tree["tree"])
            edge_list = utils.read_tree_dict(dict_tree["tree"], h1=model.root_host, edge_list=[])
            model.T.add_edges_from(edge_list)

            return model

    def infection_time_from_sampling_step(self, selected_host=None, metHast=True, verbose=False):
        """
        Propose and possibly accept a new infection time for a sampled host using the Metropolis-Hastings algorithm.

        This method samples a new infection time for a selected host (or a random sampled host if not provided),
        computes the acceptance probability, and updates the host's infection time if the proposal is accepted.

        Parameters
        ----------
        selected_host : host, optional
            The host whose infection time will be changed. If None, a random sampled host is selected.
        metHast : bool, optional
            If True, use the Metropolis-Hastings algorithm to accept or reject the proposal. Default is True.
        verbose : bool, optional
            If True, print detailed information about the proposal. Default is False.

        Returns
        -------
        t_inf_new : float
            The proposed new infection time.
        gg : float
            Proposal ratio for the Metropolis-Hastings step.
        pp : float
            Likelihood ratio for the Metropolis-Hastings step.
        P : float
            Acceptance probability for the Metropolis-Hastings step.
        selected_host : host
            The host whose infection time was proposed to change.
        """
        L_old = self.get_log_likelihood_transmission()
        rejects = 0

        ##################################################################
        ##################################################################
        #######                    INFECTION TIME                   ######
        ##################################################################
        ##################################################################
        if selected_host is None: selected_host = sample(list(self.T.nodes()), 1)[0]
        # print(selected_host.t_inf)
        while selected_host.t_sample is None:
            selected_host = sample(list(self.T.nodes()), 1)[0]
            # t_inf_old = selected_host.t_inf
            # print(selected_host.t_sample)

        # print(t_inf_old)
        t_inf_old = selected_host.t_inf
        t_inf_old2 = -selected_host.t_inf + selected_host.t_sample

        # We don't want transmissions happening before the infectors transmission
        acceptable = False
        trys = 0

        # Choosing sampled node
        while not acceptable:
            trys += 1
            # if verbose:print("tries",trys)
            t_inf_new = self.samp_sampling()
            for j in self.T.successors(selected_host):
                if j.t_inf - (selected_host.t_sample - t_inf_new) < 0:
                    acceptable = False
                    break
            else:
                for j in self.T.predecessors(selected_host):
                    if -(j.t_inf - (selected_host.t_sample - t_inf_new)) < 0:
                        acceptable = False
                        break
                else:
                    acceptable = True

        gg = ((t_inf_old2/t_inf_new ) ** (self.k_samp - 1) * np.exp(-(t_inf_old2 - t_inf_new) / self.theta_samp))

        selected_host.t_inf = selected_host.t_sample - t_inf_new

        selected_host.t_inf = t_inf_old
        self.log_likelihood = L_old

        pp = np.exp(L_new - L_old)
        P = gg * pp

        if metHast:
            if P > 1:
                selected_host.t_inf = selected_host.t_sample - t_inf_new
                self.log_likelihood = L_new
            else:
                rnd = random()
                if rnd < P:
                    selected_host.t_inf = selected_host.t_sample - t_inf_new
                    self.log_likelihood = L_new

        return t_inf_new, gg, pp, P, selected_host

    def infection_time_from_infection_model_step(self, selected_host=None, metHast=True, Dt_new=None, verbose=False):
        """
        Method to change the infection time of a host and then accept the change using the Metropolis Hastings algorithm.

        Parameters
        ----------
        selected_host: host object, default=None
            Host whose infection time will be changed. If None, a host is randomly selected.
        metHast: bool, default=True
            If True, the Metropolis Hastings algorithm is used to accept or reject the change.
        Dt_new: float, default=None
            New infection time for the host. If None, a new time is sampled.
        verbose: bool, default=False
            If True, prints the results of the step.
        """
        L_old = self.log_likelihood
        # rejects = 0

        ##################################################################
        ##################################################################
        #######                    INFECTION TIME                   ######
        ##################################################################
        ##################################################################


        if selected_host is None:
            while True:
                selected_host = sample(list(self.T.nodes()), 1)[0]
                if selected_host!=self.root_host:break



        parent = self.parent(selected_host)

        # print(t_inf_old)
        t_inf_old = selected_host.t_inf
        Dt_old = +selected_host.t_inf - parent.t_inf

        # We don't want transmissions happening before the infectors transmission
        t_min = None
        # Choosing sampled node
        if Dt_new is None:
            if self.out_degree(selected_host) == 0 and not selected_host.sampled:
                if selected_host.sampled:
                    t_min = selected_host.t_sample-selected_host.t_inf
                    Dt_new = self.dist_infection.ppf(random() * self.dist_infection.cdf(t_min))
                else:
                    Dt_new = self.samp_infection()
                try:
                    gg = ((Dt_old / Dt_new) ** (self.k_inf - 1) * np.exp(-(Dt_old - Dt_new) / self.theta_inf))
                except RuntimeWarning:
                    print(Dt_old, Dt_new, self.pdf_infection(Dt_new), self.k_inf, self.theta_inf, self.out_degree(selected_host), t_min, self.dist_infection.cdf(t_min))
                    if Dt_old<0:
                        print("NEGATIVE!!!",selected_host.t_inf, parent.t_inf)
                    raise RuntimeWarning

                if verbose:
                    print("An unsampled leaf has been shifted")
                    print(
                        f"\tDt_new: {Dt_new}, Dt_old: {Dt_old}, gg: {gg}, pp: {1/gg}, P: {1}, selected_host: {selected_host}")
                if metHast:

                    DL = np.log(1/gg)
                    pp = np.exp(DL)

                    # DL2 = utils.Delta_log_gamma(Dt_old, Dt_new, self.k_inf, self.theta_inf)


                    # if selected_host.sampled:
                    #     Dt_samp_old = selected_host.t_sample - t_inf_old
                    #     Dt_samp_new = selected_host.t_sample - t_inf_new
                    #     # DL += (self.k_samp-1)*np.log(Dt_samp_new/Dt_samp_old) - ((Dt_samp_new-Dt_samp_old)/self.theta_samp)
                    #     DL2 += utils.Delta_log_gamma(Dt_samp_old, Dt_samp_new, self.k_samp, self.theta_samp)



                    # print(f"----------->{DL=}\t{DL2=}")

                    #Genetic prior
                    if self.genetic_prior is not None:
                        LP_top_old = self.genetic_prior.correction_LL
                        LP_old = self.genetic_log_prior
                        selected_host.t_inf = parent.t_inf + Dt_new
                        LP_new = self.genetic_prior.log_prior_T(self.T,verbose=verbose)
                        selected_host.t_inf = parent.t_inf + Dt_old
                        DL_prior = LP_new-LP_old

                        pp *= np.exp(DL_prior)

                    #location prior
                    if self.same_location_prior is not None:
                        LP_sloc_top_old = self.same_location_prior.correction_LL
                        LP_sloc_old = self.same_location_log_prior
                        selected_host.t_inf = parent.t_inf + Dt_new
                        LP_sloc_new = self.same_location_prior.log_prior_T(self.T)
                        selected_host.t_inf = parent.t_inf + Dt_old
                        DL_prior_same_location = LP_sloc_new - LP_sloc_old

                        pp *= np.exp(DL_prior_same_location)

                    P = gg*pp

                    if P > 1:
                        selected_host.t_inf = parent.t_inf + Dt_new
                        self.log_likelihood += DL
                        accepted = True
                        if self.genetic_prior is not None:
                            self.genetic_log_prior = LP_new

                        if self.same_location_prior is not None:
                            self.same_location_prior.log_prior = LP_sloc_new
                            self.same_location_log_prior = LP_sloc_new

                        if verbose:
                            print("Time shift accepted")
                            print("__" * 50, "\n\n")
                    else:
                        if random() < P:
                            selected_host.t_inf = parent.t_inf + Dt_new
                            self.log_likelihood += DL
                            accepted = True
                            if self.genetic_prior is not None:
                                self.genetic_log_prior = LP_new
                            if self.same_location_prior is not None:
                                self.same_location_prior.log_prior = LP_sloc_new
                                self.same_location_log_prior = LP_sloc_new

                            if verbose:
                                print("Time shift accepted")
                                print("__" * 50, "\n\n")
                        else:
                            accepted = False

                            if self.genetic_prior is not None:
                                self.genetic_prior.correction_LL = LP_top_old
                                self.genetic_prior.log_prior = LP_old
                                self.genetic_log_prior = LP_old
                            if self.same_location_prior is not None:
                                self.same_location_prior.log_prior = LP_sloc_old
                                self.same_location_log_prior = LP_sloc_old
                                self.same_location_prior.correction_LL = LP_sloc_top_old
                            if verbose:
                                print("Time shift rejected")
                                if self.genetic_prior is not None:
                                    print(f"------>\t{DL=}\t{DL_prior=}")
                                else:
                                    print(f"------>\t{DL=}\t{DL_prior=}")

                                print("__" * 50, "\n\n")

                    # self.log_likelihood += DL

                    return Dt_new, gg, pp, P, selected_host, accepted,DL


                # # if selected_host.sampled:
                # DL = (self.k_inf - 1) * np.log(Dt_new / Dt_old) - ((Dt_new - Dt_old) / self.theta_inf)
                # for h in self.successors(selected_host):
                #     Dt_h_old = h.t_inf - t_inf_old
                #     Dt_h_new = h.t_inf - t_inf_new
                #     DL += (self.k_inf - 1) * np.log(Dt_h_new / Dt_h_old) - ((Dt_h_new - Dt_h_old) / self.theta_inf)
                #
                # if selected_host.sampled:
                #     Dt_samp_old = selected_host.t_sample - t_inf_old
                #     Dt_samp_new = selected_host.t_sample - t_inf_new
                #     DL += (self.k_samp - 1) * np.log(Dt_samp_new / Dt_samp_old) - (
                #             (Dt_samp_new - Dt_samp_old) / self.theta_samp)
                #
                #
                # print(
                #     f"A saco {L_new - L_old}, solo cambios {DL}, diff {L_new - L_old - DL}, similar? {np.abs(L_new - L_old - DL) < 1e-9},sampled? {selected_host.sampled}")
                #
                # print(f"\t--> pp={np.exp(L_new - L_old)}, gg={gg}, P={gg * np.exp(L_new - L_old)}")
            else:
                #No leaf
                if self.out_degree(selected_host) > 0:
                    if verbose:
                        print("A no leaf has been shifted")
                    t_min = min(self.T.successors(selected_host), key=lambda j: j.t_inf).t_inf-parent.t_inf
                    #No leaf and sampled
                    if selected_host.sampled:
                        t_min = min(selected_host.t_sample - selected_host.t_inf,t_min)
                else:
                    #Leaf and sampled
                    if verbose:
                        print("A sampled leaf has been shifted")
                    t_min = selected_host.t_sample-selected_host.t_inf
                Dt_new = self.dist_infection.ppf(random() * self.dist_infection.cdf(t_min))



        try:
            gg = ((Dt_old / Dt_new) ** (self.k_inf - 1) * np.exp(-(Dt_old - Dt_new) / self.theta_inf))
        except RuntimeWarning:
            print(Dt_old, Dt_new, self.pdf_infection(Dt_new), self.k_inf, self.theta_inf, selected_host.sampled, self.out_degree(selected_host), t_min, self.dist_infection.cdf(t_min))
            if Dt_old<0:
                print("NEGATIVE!!!",selected_host.t_inf , parent.t_inf)
            raise RuntimeWarning

        # t_inf_old = selected_host.t_inf
        # selected_host.t_inf = parent.t_inf + Dt_new
        t_inf_new = parent.t_inf + Dt_new
        # L_new = self.get_log_likelihood_transmission()

        DL = utils.Delta_log_gamma(Dt_old,Dt_new,self.k_inf,self.theta_inf)
        for h in self.successors(selected_host):
            Dt_h_old = h.t_inf-t_inf_old
            Dt_h_new = h.t_inf-t_inf_new
            # DL += (self.k_inf-1)*np.log(Dt_h_new/Dt_h_old) - ((Dt_h_new-Dt_h_old)/self.theta_inf)
            DL += utils.Delta_log_gamma(Dt_h_old,Dt_h_new,self.k_inf,self.theta_inf)

        if selected_host.sampled:
            Dt_samp_old = selected_host.t_sample-t_inf_old
            Dt_samp_new = selected_host.t_sample-t_inf_new
            # DL += (self.k_samp-1)*np.log(Dt_samp_new/Dt_samp_old) - ((Dt_samp_new-Dt_samp_old)/self.theta_samp)
            DL += utils.Delta_log_gamma(Dt_samp_old,Dt_samp_new,self.k_samp,self.theta_samp)
            # DL += (self.k_inf - 1) * np.log(t_inf_new / Dt_samp_old) - ((t_inf_new - Dt_samp_old) / self.theta_inf)

        selected_host.t_inf = parent.t_inf + Dt_new
        L_new = self.log_likelihood_transmission_tree(self.T)
        selected_host.t_inf = parent.t_inf + Dt_old

        # print(f"A saco {L_new-L_old}, solo cambios {DL}, diff {L_new-L_old-DL}, similar? {np.abs(L_new-L_old-DL)<1e-9},sampled? {selected_host.sampled}")
        pp = np.exp(DL)

        # Genetic prior
        if self.genetic_prior is not None:
            LP_old = self.genetic_log_prior
            LP_top_old = self.genetic_prior.correction_LL
            selected_host.t_inf = parent.t_inf + Dt_new
            LP_new = self.genetic_prior.log_prior_T(self.T)
            selected_host.t_inf = parent.t_inf + Dt_old
            DL_prior = LP_new - LP_old

            pp *= np.exp(DL_prior)

        # location prior
        if self.same_location_prior is not None:
            LP_sloc_top_old = self.same_location_prior.correction_LL
            LP_sloc_old = self.same_location_log_prior
            selected_host.t_inf = parent.t_inf + Dt_new
            LP_sloc_new = self.same_location_prior.log_prior_T(self.T)
            selected_host.t_inf = parent.t_inf + Dt_old
            DL_prior_same_location = LP_sloc_new - LP_sloc_old

            pp *= np.exp(DL_prior_same_location)

        P = gg * pp
        if verbose:
            if self.genetic_prior is not None:
                print(f"\tDt_new: {Dt_new}, Dt_old: {Dt_old}, gg: {gg}, pp: {pp}, P: {P}, selected_host: {selected_host} , L_new: {L_new}, L_old: {L_old}, DL: {DL} vs {L_new-L_old}")
                print(f"\t\tGenetic prior:  LP_new: {LP_new}, LP_old: {LP_old}, DL_prior: {DL_prior}")
            else:
                print(f"\tDt_new: {Dt_new}, Dt_old: {Dt_old}, gg: {gg}, pp: {pp}, P: {P}, selected_host: {selected_host} , L_new: {L_new}, L_old: {L_old}, DL: {DL} vs {L_new-L_old}")
        # pp2 = likelihood_ratio(self,selected_host,t_inf_old,selected_host.t_sample-Dt_new,log=False)
        # L_old = self.log_likelihood_transmission()

        accepted = False

        # Metropolis Hastings
        if metHast:
            if P > 1:
                accepted = True
                selected_host.t_inf = parent.t_inf + Dt_new
                self.log_likelihood += DL

                if verbose:
                    print("Time shift accepted")
                    print("__" * 50, "\n\n")

                if self.genetic_prior is not None:
                    self.genetic_log_prior = LP_new
                if self.same_location_prior is not None:
                    self.same_location_prior.log_prior = LP_sloc_new
                    self.same_location_log_prior = LP_sloc_new
            else:
                rnd = random()
                # print(P,algo)
                if rnd < P:
                    accepted = True
                    selected_host.t_inf = parent.t_inf + Dt_new
                    self.log_likelihood += DL
                    if verbose:
                        print("Time shift accepted")
                        print("__"*50,"\n\n")

                    if self.genetic_prior is not None:
                        self.genetic_log_prior = LP_new
                    if self.same_location_prior is not None:
                        self.same_location_prior.log_prior = LP_sloc_new
                        self.same_location_log_prior = LP_sloc_new
                    # print("rejected",itt)
                else:
                    accepted = False
                    # DL = DL
                    if self.genetic_prior is not None:
                        self.genetic_prior.correction_LL = LP_top_old
                        self.genetic_prior.log_prior = LP_old
                    if self.same_location_prior is not None:
                        self.same_location_prior.log_prior = LP_sloc_old
                        self.same_location_log_prior = LP_sloc_old

                    if verbose:
                        print("Time shift rejected")
                        print("__"*50,"\n\n")

        return Dt_new, gg, pp, P, selected_host, accepted, DL

    def add_unsampled_with_times(self, selected_host=None, P_add=0.5,  P_rewiring=0.5, P_off=0.5, verbose=False,
                                 only_geometrical=False, detailed_probs=False):
        """
        Method to propose the addition of an unsampled host to the transmission tree and get the probability of the proposal.

        Parameters:
        -----------

        selected_host: host object
            Host to which the unsampled host will be added. If None, a host is randomly selected.
        P_add: float
            Probability of proposing to add a new host to the transmission tree.
        P_rewiring: float
            Probability of rewiring the new host to another sibling host.
        P_off: float
            Probability to rewire the new host to be a leaf.
        verbose: bool
            If True, prints the results of the step.
        only_geometrical: bool
            If True, only the proposal of the new topological structure will be considered.
        detailed_probs: bool
            If True, the method will return both probabilities of the proposals, of adding and removing a host.

        Returns:
        --------
        T_new: DiGraph object
            New transmission tree with the proposed changes.
        gg: float
            Ratio of the probabilities of the proposals.
        g_go: float
            Probability of the proposal of adding a host.
        g_ret: float
            Probability of the proposal of removing a host.
        prob_time: float
            Probability of the time of infection of the new host.
        unsampled: host object
            Unsampeld host to be added to the transmission tree.
        added: bool
            If True, the host was added to the transmission tree.

        """
        if selected_host is None:
            selected_host = choice(list(self.T.nodes()))
        #         print(str(selected_host),list(self.T.nodes()))
        k_selected_host = self.out_degree(selected_host)

        # Choosing add in new link or add between two nodes
        if k_selected_host > 0:
            sibling = self.choose_successors(selected_host)[0]
            t_min = sibling.t_inf

            Dt_max = t_min - selected_host.t_inf

            link = [selected_host, sibling]
        else:
            sibling = None
            Dt = self.samp_infection()
            # if only_geometrical:prob_time = 1
            prob_time = self.pdf_infection(Dt)

        s = "U"
        Dt = 0

        try:
            unsampled = host(s, -1, genetic_data=[], t_inf=(selected_host.t_inf + Dt))
        except Exception as e:

            print(sibling, selected_host)
            for h in self.T:
                print(h, type(h))
            raise e
        # N_no_leaves = len([h for h in self.T.nodes() if self.out_degree(h)>0])

        # Rewiring

        if random() < P_rewiring and k_selected_host != 0:
            if verbose:
                print("Trying to rewire too:")
            if random() < P_off or k_selected_host == 1:
                # The unsampled host slices to be with no children (to offspring)
                links = [(selected_host, unsampled)]
                to_remove = []
                if verbose:
                    if k_selected_host > 1:
                        print("\t Unsampled host will no infect no one")
                    else:
                        print("\t Unsampled host CANNOT infect no one")

                if k_selected_host > 1:
                    g_go = (1 / len(self.T)) * P_rewiring * P_off

                    # links.append((unsampled,sibling))
                else:
                    g_go = (1 / len(self.T)) * P_rewiring
                sibling = None
                Dt_max = None
                k_unsampled = 0

                if not only_geometrical:
                    Dt = self.dist_infection.ppf(random())
                    unsampled.t_inf = selected_host.t_inf + Dt
                    prob_time = self.pdf_infection(Dt)
                    g_go *= prob_time

            else:
                k_unsampled = randint(1, k_selected_host - 1)
                links = [(selected_host, unsampled), (unsampled, sibling)]
                to_remove = [(selected_host, sibling)]
                if verbose:
                    print(f"\t Unsampled host will infect {k_unsampled + 1} people.")

                # if k_unsampled>1:
                to_go_down = utils.random_combination(
                    combinations([h for h in self.successors(selected_host) if h != sibling], k_unsampled))[0]
                comb_count = factorial(k_selected_host - 1) / (
                (factorial(k_selected_host - k_unsampled) * factorial(k_unsampled)))

                # Ps = prob_time
                Dt_max = sibling.t_inf - selected_host.t_inf
                for i, h in enumerate(to_go_down):
                    # if k_unsampled ==1:print("-------",h,type(to_go_down))
                    links.append((unsampled, h))
                    to_remove.append((selected_host, h))
                    # print(h,h.t_inf)
                    if not only_geometrical and (h.t_inf - selected_host.t_inf) < Dt_max:
                        Dt_max = h.t_inf - selected_host.t_inf
                # if k_unsampled ==1:

                k_unsampled += 1  # kph=2

                g_go = (1 / len(self.T)) * (1 / (k_selected_host * (k_selected_host - 1)))
                g_go *= k_unsampled * (factorial(k_unsampled - 1) * factorial(k_selected_host - k_unsampled) / (
                    factorial(k_selected_host - 1)))

                g_go *= (1 - P_off) * P_rewiring

                if not only_geometrical:
                    Dt = self.dist_infection.ppf(random() * self.dist_infection.cdf(Dt_max))
                    unsampled.t_inf = selected_host.t_inf + Dt
                    prob_time = self.pdf_infection(Dt) / self.dist_infection.cdf(Dt_max)
                    g_go *= prob_time


        else:
            links = [(selected_host, unsampled)]
            if k_selected_host > 0:
                if verbose: print("No rewire:")
                g_go = (1 / len(self.T)) * (1 - P_rewiring)
                g_go *= (1 / k_selected_host)
                if not only_geometrical:
                    Dt_max = sibling.t_inf - selected_host.t_inf
                    # print("------------",selected_host.t_inf - sibling.t_inf,selected_host.t_inf, sibling.t_inf)
                    Dt = self.dist_infection.ppf(random() * self.dist_infection.cdf(Dt_max))
                    unsampled.t_inf = selected_host.t_inf + Dt
                    prob_time = self.pdf_infection(Dt) / self.dist_infection.cdf(Dt_max)
                    g_go *= prob_time
                k_unsampled = 1
                links.append((unsampled, sibling))
                to_remove = [(selected_host, sibling)]
            else:
                g_go = (1 / len(self.T))
                Dt_max = None
                if not only_geometrical:
                    if verbose: print("No option for rewiring")
                    Dt = self.dist_infection.ppf(random())
                    unsampled.t_inf = selected_host.t_inf + Dt
                    prob_time = self.pdf_infection(Dt)
                    g_go *= prob_time
                k_unsampled = 0
                to_remove = []

        if only_geometrical:
            prob_time = None

        # Ratio of proposals

        g_ret = 1 / (len(self.unsampled_hosts) + 1)

        # If there are no more unsampled hosts, the proposal is divided by P_add
        if len(self.unsampled_hosts) == 0:
            g_ret /= 2

        gg = ((1-P_add)*g_ret) / (P_add*g_go)

        # if abs(gg*prob_time-0.5)<1e-6:
        #     print(gg,g_go,g_ret,len(self.T),(len(self.unsampled_hosts),k_selected_host))

        if verbose:
            if prob_time is not None:
                if sibling is not None:
                    print(
                        f"ADDING:\nadd:\n\t-g(+1){g_go}, Dt {Dt}, Dt_max {Dt_max}, pt {prob_time}, g(+1)/pt {g_go / prob_time}\n\t-g(-1){g_ret}\n\tg_ret/g_go {g_ret / g_go}\nAdding {str(unsampled)} below {str(selected_host)} with degree {k_unsampled} and sibling {sibling} {sibling.t_inf}\n-----------------------------------------------------")
                else:
                    print(
                        f"ADDING:\nadd:\n\t-g(+1){g_go}, Dt {Dt}, Dt_max {Dt_max}, pt {prob_time}, g(+1)/pt {g_go / prob_time}\n\t-g(-1){g_ret}\n\tg_ret/g_go {g_ret / g_go}\nAdding {str(unsampled)} below {str(selected_host)} with degree {k_unsampled} and sibling {sibling}\n-----------------------------------------------------")
            else:
                print(
                    f"ADDING:\nadd:\n\t-g(+1){g_go}, pt {prob_time}, g(+1)/pt {g_go}\n\t-g(-1){g_ret}\n\tg_ret/g_go {g_ret / g_go}\nAdding {str(unsampled)} below {str(selected_host)} with degree {k_unsampled}\n-----------------------------------------------------")



        T_new = nx.DiGraph(self.T)


        T_new.remove_edges_from(to_remove)
        T_new.add_edges_from(links)

        if detailed_probs:
            return T_new, gg, g_go, g_ret, prob_time, unsampled, True
        else:
            return T_new, gg, unsampled, True

    def remove_unsampled_with_times(self, selected_host=None, P_add=0.5, P_rewiring=0.5, P_off=0.5, only_geometrical=False,
                                    detailed_probs=False, verbose=False):
        """
        Method to propose the removal of an unsampled host from the transmission tree and get the probability of the proposal.
        In case that no unsampled hosts are available, a new host is proposed to be added to the transmission tree.

        Parameters:
        -----------
        selected_host: host object
            Unsampled host to be removed from the transmission tree. If None, a host is randomly selected.
        P_add: float
            Probability of proposing to add a new host to the transmission tree.
        P_rewiring: float
            Probability of rewiring the new host to another sibling host.
        P_off: float
            Probability to rewire the new host to be a leaf.
        verbose: bool
            If True, prints the results of the step.
        only_geometrical: bool
            If True, only the proposal of the new topological structure will be considered.
        detailed_probs: bool
            If True, the method will return both probabilities of the proposals, of adding and removing a host.

        Returns:
        --------
        T_new: DiGraph object
            New transmission tree with the proposed changes.
        gg: float
            Ratio of the probabilities of the proposals.
        g_go: float
            Probability of the proposal of adding a host.
        g_ret: float
            Probability of the proposal of removing a host.
        prob_time: float
            Probability of proposing the time of the selected_host.
        added: bool
            If True, the host was added to the transmission tree. Else, the node have been removed
        """
        if selected_host is None:
            if len(self.unsampled_hosts) == 0:
                if detailed_probs:
                    T_new, gg, g_go, g_ret, pt, unsampled, added = self.add_unsampled_with_times(
                                                                                            P_rewiring=P_rewiring,
                                                                                            P_off=P_off,
                                                                                            detailed_probs=detailed_probs,
                                                                                            verbose=verbose)
                    # print("---------",T_new,gg,g_go,g_ret,pt,unsampled,added)
                    return T_new, gg, g_go, g_ret, pt, unsampled, added
                else:
                    T_new, gg, unsampled, added = self.add_unsampled_with_times( P_rewiring=P_rewiring,
                                                                               P_off=P_off,
                                                                               detailed_probs=detailed_probs,
                                                                               verbose=verbose)
                    return T_new, gg, unsampled, added
                # gg /= 2

            selected_host = choice(self.unsampled_hosts)
        try:
            k_selected_host = self.out_degree(selected_host)  # outdefgree selected host
        except NetworkXError as e:
            print(selected_host, type(selected_host))
            raise e

        try:
            parent = self.parent(selected_host)
        except NetworkXError as e:
            # pos = hierarchy_pos_times(self.T,self.root_host)
            pos = nx.spring_layout(self.T)
            tm.plot_transmision_network(self.T, nodes_labels=True, pos=pos, highlighted_nodes=[selected_host])
            raise e

        k_parent = self.out_degree(parent)  # outdefgree parent
        kappa = k_parent + k_selected_host
        children = list(self.successors(selected_host))

        # Ratio of proposals
        ##probability of coming back
        kappa = k_selected_host + k_parent

        ## Probablities of times
        # if only_geometrical:
        ###Probability of having a time with no sibling
        ###Sum of probabilities

        if not only_geometrical:
            Dt = selected_host.t_inf - parent.t_inf
            p_t = self.pdf_infection(Dt)

        g_go = 1 / len(self.unsampled_hosts)

        new_len_prob = 1 / (len(self.T) - 1)
        if k_selected_host == 0:
            if verbose: print("Removing a leave (rewiring to offspring was choosed):")
            if k_parent == 1:
                if verbose: print("\t-It had no rewire")
                g_ret = new_len_prob
            elif k_parent == 2:
                if verbose: print("\t-It had one possible rewire")
                g_ret = new_len_prob * P_rewiring
            elif k_parent > 2:
                if verbose: print("\t-It had a rewire and to become a leave was an option")
                g_ret = new_len_prob * P_rewiring * P_off
            if not only_geometrical:
                if verbose: print(f"\t-p_t {p_t} no Dt_max")
                g_ret *= p_t

                # print("ZERO!!!")
        elif k_selected_host == 1:
            if verbose: print("Removing a host with 1 sibling (no rewiring was choosed):")
            if k_parent == 1:
                g_ret = new_len_prob * (1 - P_rewiring)
            else:
                g_ret = new_len_prob * (1 - P_rewiring) / (kappa - 1)
            if not only_geometrical:
                Dt_max = next(self.successors(selected_host)).t_inf - parent.t_inf
                p_t /= self.dist_infection.cdf(Dt_max)
                if verbose: print(f"\t-p_t {p_t} with Dt {Dt} and Dt_max {Dt_max}")
                g_ret *= p_t

        else:
            # # child = children[0]
            # if k_parent==1:
            if verbose: print("Removing a leave (rewiring to chain other hosts was choosen):")
            g_ret = new_len_prob * P_rewiring * (1 - P_off)
            g_ret *= factorial(k_parent - 1) * factorial(k_selected_host) / (
                    (kappa - 1) * (kappa - 2) * factorial(kappa - 2))
            if not only_geometrical:
                Dt_max = 10000000000
                for h in self.successors(selected_host):
                    Dt_2 = h.t_inf - parent.t_inf
                    if Dt_2 < Dt_max:
                        Dt_max = Dt_2
                p_t /= self.dist_infection.cdf(Dt_max)
                if verbose: print(f"\t-p_t {p_t} with Dt {Dt} and Dt_max {Dt_max}")
                g_ret *= p_t

        if len(self.unsampled_hosts) == 1:
            g_go /= 2

        gg = (P_add * g_ret) / ((1-P_add) * g_go)

        if only_geometrical:
            pt = None

        if verbose:
            if only_geometrical:
                print(
                    f"remove:\n\t-g(+1){g_ret}, g(+1)/pt {g_ret}\n\t-g(-1){g_go}\n\tg_ret/g_go {g_ret / g_go}\n{len(self.unsampled_hosts), len(self.T), k_selected_host, k_parent, new_len_prob}\n-----------------------------------------------------")
            else:
                print(
                    f"remove:\n\t-g(+1){g_ret},Dt {Dt},pt {p_t},g(+1)/pt {g_ret / p_t}\n\t-g(-1){g_go}\n\tg_ret/g_go {g_ret / g_go}\n{len(self.unsampled_hosts), len(self.T), k_selected_host, k_parent, new_len_prob}\n-----------------------------------------------------")
        # rewiring
        T_new = nx.DiGraph(self.T)
        T_new.remove_edge(parent, selected_host)

        for child in children:
            T_new.remove_edge(selected_host, child)
            T_new.add_edge(parent, child)

        T_new.remove_node(selected_host)

        if detailed_probs:
            # print("AQUI!!!")
            return T_new, gg, g_go, g_ret, p_t, selected_host, False
        return T_new, gg, selected_host, False


    def add_remove_step(self, P_add=0.5, P_rewiring=0.5, P_off=0.5, metHast=True, verbose=False):
        """
        Method to propose the addition or removal of an unsampled host to the transmission tree and get the probability of the proposal.

        Parameters:
        -----------
        P_add: float
            Probability of proposing an addition of an unsampled host. Else, an unsampled host is going to be proposed for removal.
        P_rewiring: float
            Probability of rewiring the new host to another sibling host.
        P_off: float
            Probability to rewire the new host to be a leaf.
        metHast: bool
            If True, the Metropolis Hastings algorithm is used to accept or reject the change.
        verbose: bool
            If True, prints the results of the step.

        Returns:
        --------

        """

        if random() < P_add:
            T_new, gg, unsampled, added =  self.add_unsampled_with_times(P_add=P_add, P_rewiring=P_rewiring, P_off=P_off, verbose=verbose)
        else:
            T_new, gg, unsampled, added =  self.remove_unsampled_with_times(P_add=P_add, P_rewiring=P_rewiring, P_off=P_off, verbose=verbose)

        if added:
            affected_hosts = [list(T_new.predecessors(unsampled))[0]] + [h for h in T_new.successors(unsampled)]
            Delta_LL = self.Delta_log_likelihood_host(affected_hosts, T_new) + self.log_likelihood_hosts_list([unsampled], T_new)
        else:
            affected_hosts = [self.parent(unsampled)]+[h for h in self.successors(unsampled)]
            Delta_LL = self.Delta_log_likelihood_host(affected_hosts, T_new) - self.log_likelihood_host(unsampled)

        # L_new = self.log_likelihood_transmission_tree(T_new)

        # nn = dict_replace.get(nn, nn)

        pp = np.exp(Delta_LL)

        if self.genetic_prior is not None:
            LP_top_old = self.genetic_prior.correction_LL
            LP_new = self.genetic_prior.log_prior_T(T_new,verbose=verbose)
            DL_prior = LP_new - self.genetic_log_prior
            pp *= np.exp(DL_prior)

        if self.same_location_prior is not None:
            DL_prior_same_location, LP_sloc_new, LP_sloc_old, LP_sloc_top_old = self.compute_Delta_loc_prior(T_new)
            pp *= np.exp(DL_prior_same_location)
        P = gg * pp

        if metHast:
            if P > 1:
                self.log_likelihood = self.log_likelihood+Delta_LL
                self.T = T_new
                accepted = True

                if not added:
                    if verbose:
                        print(f"Removing host accepted with acceptance probability {P}")
                        print("__" * 50, "\n\n")
                    self.unsampled_hosts.remove(unsampled)
                else:
                    if verbose:
                        print(f"Adding host accepted with acceptance probability {P}")
                        print("__" * 50, "\n\n")
                    self.unsampled_hosts.append(unsampled)
                if self.genetic_prior is not None:
                    self.genetic_log_prior = LP_new
                    self.genetic_prior.log_prior = LP_new
                if self.same_location_prior is not None:
                    self.same_location_prior.log_prior = LP_sloc_new
                    self.same_location_log_prior = LP_sloc_new
            else:
                if random() < P:
                    accepted = True
                    self.log_likelihood = self.log_likelihood+Delta_LL
                    self.T = T_new
                    if self.genetic_prior is not None:
                        self.genetic_log_prior = LP_new
                        self.genetic_prior.log_prior = LP_new
                    if self.same_location_prior is not None:
                        self.same_location_prior.log_prior = LP_sloc_new
                        self.same_location_log_prior = LP_sloc_new

                    if not added:
                        if verbose:
                            print(f"Removing host accepted with acceptance probability {P}")
                            print("__" * 50, "\n\n")
                        self.unsampled_hosts.remove(unsampled)
                    else:
                        if verbose:
                            print(f"Adding host accepted with acceptance probability {P}")
                            print("__" * 50, "\n\n")
                        self.unsampled_hosts.append(unsampled)
                else:
                    if self.genetic_prior is not None:
                        self.genetic_prior.correction_LL = LP_top_old
                        self.genetic_prior.log_prior = self.genetic_log_prior
                    if self.same_location_prior is not None:
                        self.same_location_prior.log_prior = LP_sloc_old
                        self.same_location_log_prior = LP_sloc_old

                    if verbose:
                        if added:
                            print(f"Adding host rejected with acceptance probability {P}")
                            print("__" * 50, "\n\n")
                        else:
                            print(f"Removing host rejected with acceptance probability {P}")
                            print("__" * 50, "\n\n")
                    accepted = False
                    Delta_LL = 0


        return T_new,gg,pp,P,added, accepted, Delta_LL


    def MCMC_step(self, N_steps, verbose=False):

        L_old = self.get_log_likelihood_transmission()

        #Getting sure that the number of candidates to chain is computed
        self.get_N_candidates_to_chain(recompute=True)


        rejects = 0

        for itt in range(N_steps):
            if verbose:
                print(itt, L_old)

            # Choosing type of change
            selected_movement = randint(0, 1)
            if selected_movement == 0:
                self.infection_time_from_sampling_step(verbose=verbose)
            else:
                tree_slicing_step(self)
