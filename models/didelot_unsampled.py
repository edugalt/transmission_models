from transmission_models import *
from transmission_models.utils import *
from transmission_models.models.topology_movements import *

from random import random
from scipy.stats import nbinom, gamma, binom, expon, norm
from scipy.special import gamma as GAMMA


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
        self.unsampled_hosts = [h for h in self.T if not h.sampled]
        return self.unsampled_hosts

    def get_sampling_model_likelihood(self,hosts=None):
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
            for h in self.T:
                if not h.sampled:
                    L*=(1-self.pi)
                else:
                    L*=self.pi*self.pdf_sampling(h.t_sample-h.t_inf)
        return L

    def get_offspring_model_likelihood(self,hosts=None):
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
        L = 1
        if hosts is not None:
            if isinstance(hosts,list):
                for h in hosts:
                    L*=self.pmf_offspring(self.T.out_degree(h))
            else:
                L*=self.pmf_offspring(self.T.out_degree(hosts))
        else:
            for h in self.T:
                L*=self.pmf_offspring(self.T.out_degree(h))
        return L

    def get_infection_model_likelihood(self,hosts=None):
        """
        Computes the likelihood of the infection model given a list of hosts. If no list is given, the likelihood of the
        whole transmission tree is returned.

        Parameters
        ----------
        hosts: list of host objects

        Returns
        -------
            L: float
                The likelihood of the infection model given the list of hosts

        """
        L = 1
        if hosts is not None:
            if isinstance(hosts,list):
                for h in hosts:
                    for j in self.T.successors(h):
                        L*=self.pdf_infection(j.t_inf-h.t_inf)
            else:
                for j in self.T.successors(hosts):
                    L*=self.pdf_infection(j.t_inf-hosts.t_inf)
        else:
            for h in self.T:
                for j in self.T.successors(h):
                    L*=self.pdf_infection(j.t_inf-h.t_inf)
        return L

    def log_likelihood_host(self, host, T):
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
        Pi = 1
        # Sampling model
        if host.sampled:
            # sigma = gamma.pdf(h.sample_time,1,h.t_inf,tr)
            sigma = self.pdf_sampling(host.t_sample - host.t_inf)
            # print("--",int(h),sigma,h.t_inf,h.t_sample,h.t_inf-h.t_sample)
            Pi = self.pi * (sigma)
        else:
            Pi = (1 - self.pi)

        # Offspring model
        Pi *= self.pmf_offspring(T.out_degree(host))
        # print("--",int(h),self.pmf_offspring(self.T.out_degree(h)))

        # Infection model
        for j in T.successors(host):
            # print("--",int(h),int(j),j.t_inf-h.t_inf,self.pdf_infection(j.t_inf-h.t_inf))
            sigma2 = self.pdf_infection(j.t_inf - host.t_inf)
            # print("------------------",int(h),int(j),h.t_inf,j.t_inf,j.t_inf-h.t_inf)
            if sigma2 == 0:
                self.log_likelihood = -1e30
                #print("Impossible times!!!", int(host), int(j), host.t_inf, j.t_inf, j.t_inf - host.t_inf)
                return self.log_likelihood
            Pi *= sigma2
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
        self.log_likelihood = self.log_likelihood_transmission_tree(self.T)
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

    def get_newick(self):
        self.newick = tree_to_newick(self.T, root=self.root_host)

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

    def infection_time_from_infection_model_step(self, selected_host=None, metHast=True, verbose=False):
        """
        Method to change the infection time of a host and then accept the change using the Metropolis Hastings algorithm.

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
        if selected_host is None:
            while True:
                selected_host = sample(list(self.T.nodes()), 1)[0]
                if selected_host!=self.root_host:break



            # t_inf_old = selected_host.t_inf
            # print(selected_host.t_sample)

        parent = self.parent(selected_host)

        # print(t_inf_old)
        t_inf_old = selected_host.t_inf
        t_inf_old2 = +selected_host.t_inf - parent.t_inf

        # We don't want transmissions happening before the infectors transmission
        acceptable = False
        trys = 0

        # Choosing sampled node
        while not acceptable:
            trys += 1
            # if verbose:print("tries",trys)
            t_inf_new = self.samp_infection()
            if self.out_degree(selected_host) == 0:
                acceptable= True
            for j in self.T.successors(selected_host):
                if j.t_inf - (parent.t_inf + t_inf_new) < 0:
                    acceptable = False
                    # print("kk")
                    # print("...............",trys,int(selected_host),int(j),j.t_inf-(selected_host.t_sample-t_inf_new),j.t_inf,selected_host.t_sample,selected_host.t_sample-t_inf_new)
                    break
            else:
                acceptable = True

        # gg = self.pdf_sampling(t_inf_new)/self.pdf_sampling(selected_host.t_sample-selected_host.t_inf)
        gg = ((t_inf_old2/t_inf_new ) ** (self.k_inf - 1) * np.exp(-(t_inf_old2 - t_inf_new) / self.theta_inf))
        selected_host.t_inf = parent.t_inf + t_inf_new
        L_new = self.get_log_likelihood_transmission()

        selected_host.t_inf = t_inf_old
        self.log_likelihood = L_old

        pp = np.exp(L_new - L_old)
        P = gg * pp
        # pp2 = likelihood_ratio(self,selected_host,t_inf_old,selected_host.t_sample-t_inf_new,log=False)
        # L_old = self.log_likelihood_transmission()

        # Metropolis Hastings
        if metHast:
            if P > 1:
                selected_host.t_inf = parent.t_inf + t_inf_new
                self.log_likelihood = L_new
            else:
                rnd = random()
                # print(P,algo)
                if rnd < P:
                    selected_host.t_inf = parent.t_inf + t_inf_new
                    self.log_likelihood = L_new
                    # print("rejected",itt)
        return t_inf_new, gg, pp, P, selected_host

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
