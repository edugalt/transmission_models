from itertools import combinations
from math import factorial

from transmission_models import *
import transmission_models.utils as utils
# from transmission_models.utils import tree_to_newick
from transmission_models.models.topology_movements import *

from random import random
from scipy.stats import nbinom, gamma, binom, expon, norm
from scipy.special import gamma as GAMMA

# from ..utils import tree_to_newick


class didelot_unsampled():
    """docstring for didelots_unsampled"""

    def __init__(self, sampling_params, offspring_params, infection_params):
        # super(didelots_unsampled, self).__init__()

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


        self.T = nx.DiGraph()
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
        # def generate_networks(self):

    def samp_t_inf_between(self, h1, h2):
        """
        Samples a time of infection between two hosts, one being the infector and the other the infected.
        It use a rejection sampling method to sample the time of infection of the infected host using the chain model from Didelot et al. 2017.

        Parameters:
        -----------
            h1: host object
                Infector host.
            h2: host object
                Infected host.

        Returns:
        --------
            t: float
                Time of infection of the host infected by h1 and the infector of h2.

        """
        # Dt = abs(h1 - h2)

        Dt = h2.t_inf-h1.t_inf
        return sample_in_between(self,Dt)

    def add_root(self, t_sampl, id="0", genetic_data=[], t_inf=0, t_sample=None):
        self.root_host = host(id, 0, genetic_data, t_inf, t_sampl)
        return self.root_host

    def successors(self, host):
        return self.T.successors(host)

    def parent(self, host):
        return list(self.T.predecessors(host))[0]

    def out_degree(self, host):
        return self.T.out_degree(host)

    def choose_successors(self, host, k=1):
        """
        Chooses k unique successors of a given host.

        Parameters:
            host: host object
                Hosts whose successors will be chosen.
            k: int
                Number of successors to choose.

        Returns:

        """
        return sample(list(self.successors(host)), k)

    def get_candidates_to_chain(self):
        self.candidates_to_chain = [h for h in self.T if
                                    not h == self.root_host and not self.out_degree(self.parent(h)) == 1]

        return self.candidates_to_chain

    def get_N_candidates_to_chain(self, recompute=False):
        if recompute:
            self.N_candidates_to_chain = len(self.get_candidates_to_chain())
        else:
            self.N_candidates_to_chain = len(self.candidates_to_chain)

        return self.N_candidates_to_chain

    def get_unsampled_hosts(self):
        self.unsampled_hosts = [h for h in self.T if not h.sampled or h == self.root_host]
        return self.unsampled_hosts

    def get_sampling_model_likelihood(self,hosts=None,T=None, update=False):
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
        Pi = 1
        # Sampling model
        Pi *= self.get_sampling_model_likelihood(host,T)

        # Offspring model
        Pi *= self.get_offspring_model_likelihood(host,T)

        # Infection model
        Pi *= self.get_infection_model_likelihood(host,T)

        return np.log(Pi)

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

    def log_likelihood_transmission_tree_old(self, T):
        log_likelihood = 0
        Pi = 1

        for h in T:
            # Sampling model
            if h.sampled:
                # sigma = gamma.pdf(h.sample_time,1,h.t_inf,tr)
                sigma = self.pdf_sampling(h.t_sample - h.t_inf)
                # print("--",int(h),sigma,h.t_inf,h.t_sample,h.t_inf-h.t_sample)
                Pi = self.pi * (sigma)
            else:
                Pi = (1 - self.pi)

            # Offspring model
            Pi *= self.pmf_offspring(T.out_degree(h))
            # print("--",int(h),self.pmf_offspring(self.T.out_degree(h)))

            # Infection model
            for j in T.successors(h):
                # print("--",int(h),int(j),j.t_inf-h.t_inf,self.pdf_infection(j.t_inf-h.t_inf))
                sigma2 = self.pdf_infection(j.t_inf - h.t_inf)
                # print("------------------",int(h),int(j),h.t_inf,j.t_inf,j.t_inf-h.t_inf)
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

        for i,h in enumerate(self.T):
        # for h in self.T:
            sampling_likelihood = self.get_sampling_model_likelihood(h)
            infection_likelihood = self.get_infection_model_likelihood(h)
            offspring_likelihood = self.get_offspring_model_likelihood(h)

            self.sampling_likelihood *= sampling_likelihood
            self.infection_likelihood *= infection_likelihood
            self.offspring_likelihood *= offspring_likelihood

            self.log_likelihood += np.log(sampling_likelihood*offspring_likelihood*infection_likelihood)
            # print(i,h,sampling_likelihood,offspring_likelihood,infection_likelihood,self.likelihood)

        self.likelihood = np.exp(self.log_likelihood)

        self.sampling_log_likelihood = np.log(self.sampling_likelihood)
        self.infection_log_likelihood = np.log(self.infection_likelihood)
        self.offspring_log_likelihood = np.log(self.offspring_likelihood)

        return self.log_likelihood


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

    def infection_time_from_sampling_step(self, selected_host=None, metHast=True, verbose=False):
        """
        Method to change the infection time of a host amd then accept the change using the Metropolis Hastings algorithm.

        Parameters
        ----------
        verbose
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
                    # print("kk")
                    # print("...............",intento,int(selected_host),int(j),j.t_inf-(selected_host.t_sample-t_inf_new),j.t_inf,selected_host.t_sample,selected_host.t_sample-t_inf_new)
                    break
            else:
                for j in self.T.predecessors(selected_host):
                    if -(j.t_inf - (selected_host.t_sample - t_inf_new)) < 0:
                        acceptable = False
                        break
                else:
                    acceptable = True

        # gg = self.pdf_sampling(t_inf_new)/self.pdf_sampling(selected_host.t_sample-selected_host.t_inf)
        gg = ((t_inf_old2/t_inf_new ) ** (self.k_samp - 1) * np.exp(-(t_inf_old2 - t_inf_new) / self.theta_samp))
        # pp2 = likelihood_ratio(self,selected_host,t_inf_old,selected_host.t_sample-t_inf_new,log=False)
        # L_old = self.log_likelihood_transmission()

        selected_host.t_inf = selected_host.t_sample - t_inf_new
        L_new = self.get_log_likelihood_transmission()


        selected_host.t_inf = t_inf_old
        self.log_likelihood = L_old

        pp = np.exp(L_new - L_old)
        P = gg * pp
        # P2 = gg*pp2*gg
        # if not P-P2 <1e-9:
        #     for j in self.T.successors(selected_host):
        #         if  j.sampled:print("rejected with non sampled successor",j.sampled)
        # else:
        #     for j in self.T.successors(selected_host):
        #         if  j.sampled:print("accepted with non sampled successor",j.sampled)
        # print(itt,P-P2 <1e-9,selected_host.t_sample-t_inf_new)

        # Metropolis Hastings
        if metHast:
            if P > 1:
                selected_host.t_inf = selected_host.t_sample - t_inf_new
                self.log_likelihood = L_new
            else:
                rnd = random()
                # print(P,algo)
                if rnd < P:
                    selected_host.t_inf = selected_host.t_sample - t_inf_new
                    self.log_likelihood = L_new
                    # print("rejected",itt)\

        return t_inf_new, gg, pp, P, selected_host

    def infection_time_from_infection_model_step(self, selected_host=None, metHast=True, t_inf_new=None, verbose=False):
        """
        Method to change the infection time of a host and then accept the change using the Metropolis Hastings algorithm.

        Parameters
        ----------
        selected_host: host object, default=None
            Host whose infection time will be changed. If None, a host is randomly selected.
        metHast: bool, default=True
            If True, the Metropolis Hastings algorithm is used to accept or reject the change.
        t_inf_new: float, default=None
            New infection time for the host. If None, a new time is sampled.
        verbose: bool, default=False
            If True, prints the results of the step.
        """
        L_old = self.get_log_likelihood_transmission()
        rejects = 0

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
        t_inf_old2 = +selected_host.t_inf - parent.t_inf

        # We don't want transmissions happening before the infectors transmission
        t_min = None
        # Choosing sampled node
        if t_inf_new is None:
            if self.out_degree(selected_host) == 0:
                if selected_host.sampled:
                    t_min = selected_host.t_sample-selected_host.t_inf
                    t_inf_new = self.dist_infection.ppf(random() * self.dist_infection.cdf(t_min))
                else:
                    t_inf_new = self.samp_infection()
                try:
                    gg = ((t_inf_old2 / t_inf_new) ** (self.k_inf - 1) * np.exp(-(t_inf_old2 - t_inf_new) / self.theta_inf))
                except RuntimeWarning:
                    print(t_inf_old2, t_inf_new,self.pdf_infection(t_inf_new), self.k_inf, self.theta_inf,self.out_degree(selected_host),t_min,self.dist_infection.cdf(t_min))
                    if t_inf_old2<0:
                        print("NEGATIVE!!!",selected_host.t_inf , parent.t_inf)
                    raise RuntimeWarning

                if metHast:
                    selected_host.t_inf = parent.t_inf + t_inf_new
                    self.log_likelihood = self.get_log_likelihood_transmission()

                return t_inf_new, gg, 1/gg, 1, selected_host, True
            else:
                t_min = min(self.T.successors(selected_host), key=lambda j: j.t_inf).t_inf-parent.t_inf
                if selected_host.sampled:
                    t_min = min(selected_host.t_sample - selected_host.t_inf,t_min)
                t_inf_new = self.dist_infection.ppf(random() * self.dist_infection.cdf(t_min))



        try:
            gg = ((t_inf_old2/t_inf_new ) ** (self.k_inf - 1) * np.exp(-(t_inf_old2 - t_inf_new) / self.theta_inf))
        except RuntimeWarning:
            print(t_inf_old2,t_inf_new,self.pdf_infection(t_inf_new),self.k_inf,self.theta_inf,selected_host.sampled,self.out_degree(selected_host),t_min,self.dist_infection.cdf(t_min))
            if t_inf_old2<0:
                print("NEGATIVE!!!",selected_host.t_inf , parent.t_inf)
            raise RuntimeWarning


        selected_host.t_inf = parent.t_inf + t_inf_new
        L_new = self.get_log_likelihood_transmission()

        selected_host.t_inf = t_inf_old
        self.log_likelihood = L_old

        pp = np.exp(L_new - L_old)
        P = gg * pp
        if verbose:
            print(f"t_inf_new: {t_inf_new}, t_inf_old: {t_inf_old2}, gg: {gg}, pp: {pp}, P: {P}, selected_host: {selected_host}")
        # pp2 = likelihood_ratio(self,selected_host,t_inf_old,selected_host.t_sample-t_inf_new,log=False)
        # L_old = self.log_likelihood_transmission()

        # Metropolis Hastings
        if metHast:
            if P > 1:
                accepted = True
                selected_host.t_inf = parent.t_inf + t_inf_new
                self.log_likelihood = L_new
            else:
                rnd = random()
                # print(P,algo)
                if rnd < P:
                    accepted = True
                    selected_host.t_inf = parent.t_inf + t_inf_new
                    self.log_likelihood = L_new
                    # print("rejected",itt)
                else:
                    accepted = False
        return t_inf_new, gg, pp, P, selected_host, accepted

    def add_unsampled_with_times(self, selected_host=None, P_rewiring=0.5, P_off=0.5, verbose=False,
                                 only_geometrical=False, detailed_probs=False):
        """
        Method to propose the addition of an unsampled host to the transmission tree and get the probability of the proposal.

        Parameters:
        -----------

        selected_host: host object
            Host to which the unsampled host will be added. If None, a host is randomly selected.
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
                # The unsampled host slices to be with no childs (to offspring)
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

        if len(self.unsampled_hosts) == 0:
            g_ret /= 2

        gg = g_ret / g_go

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

    def remove_unsampled_with_times(self, selected_host=None, P_rewiring=0.5, P_off=0.5, only_geometrical=False,
                                    detailed_probs=False, verbose=False):
        """
        Method to propose the removal of an unsampled host from the transmission tree and get the probability of the proposal.
        In case that no unsampled hosts are available, a new host is proposed to be added to the transmission tree.

        Parameters:
        -----------
        selected_host: host object
            Unsampled host to be removed from the transmission tree. If None, a host is randomly selected.
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

        gg = g_ret / g_go

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
