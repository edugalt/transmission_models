
from random import choice,randint,random,sample,choices
from scipy.special import gamma as GAMMA
from scipy.stats import nbinom, gamma, binom, expon, poisson
import numpy as np

class genetic_prior_tree():
    def __init__(self, model, mu, distance_matrix):
        self.mu = mu
        self.distance_matrix = distance_matrix

        self.prior_dist = poisson(mu)

        self.model = model

    @staticmethod
    def search_firsts_sampled_siblings(host, T):

        sampled_hosts = []
        for h in T.successors(host):
            if h.sampled:
                sampled_hosts.append(h)
            else:
                sampled_hosts += genetic_prior_tree.search_firsts_sampled_siblings(h, T)

        return sampled_hosts

    @staticmethod
    def search_first_sampleed_parent(host, T, root):

        if host == root:
            return None

        parent = next(T.predecessors(host))

        if not parent.sampled:
            return genetic_prior_tree.search_first_sampleed_parent(parent, T, root)
        else:
            return parent
    @staticmethod
    def get_mut_time_dist(hp, hs):
        return (hs.t_sample + hp.t_sample - 2 * hp.t_inf)

    def prior_host(self, host, T, parent_dist=False):
        log_prior = 0
        for h2 in T[host]:
            if h2.sampled:
                # print(f"{host}-->{h2}")
                Dt = h2.t_sample - host.t_sample
                log_prior += np.log(poisson(self.mu * Dt).pmf(self.distance_matrix[int(host), int(h2)]))
                p = poisson(self.mu * Dt).pmf(self.distance_matrix[int(host), int(h2)])
                # print(int(h),int(h2),Dt,p,np.log(p))
            else:
                siblings = genetic_prior_tree.search_firsts_sampled_siblings(h2, T)
                for hs in siblings:
                    # print(f"{host}-->{hs}")
                    Dt = hs.t_sample - host.t_sample
                    log_prior += np.log(poisson(self.mu * Dt).pmf(self.distance_matrix[int(host), int(hs)]))
        if parent_dist and host != self.model.root_host:
            parent = self.model.parent(host)
            if parent.sampled:
                Dt = host.t_sample - parent.t_sample
                log_prior += np.log(poisson(self.mu * Dt).pmf(self.distance_matrix[int(host), int(parent)]))
            else:
                parent = genetic_prior_tree.search_first_sampleed_parent(host, T, self.model.root_host)
                if parent is not None:
                    # print(f"{parent}-->{host}")
                    Dt = host.t_sample - parent.t_sample
                    log_prior += np.log(poisson(self.mu * Dt).pmf(self.distance_matrix[int(host), int(parent)]))

        return log_prior

    def prior_pair(self, h1, h2):
        log_prior = 0
        if not h1.sampled or not h2.sampled:
            return 0
        Dt = h2.t_sample - h1.t_sample
        return np.log(poisson(self.mu * Dt).pmf(self.distance_matrix[int(h1), int(h2)]))


    def log_prior_T(self, T):
        log_prior = 0
        for h in T:
            if not h.sampled: continue
            for h2 in T[h]:
                if h2.sampled:
                    # print(f"{h}-->{h2}")
                    Dt = self.get_mut_time_dist(h, h2)
                    log_prior += np.log(poisson(self.mu * Dt).pmf(self.distance_matrix[int(h), int(h2)]))
                    # p = poisson(self.mu * Dt).pmf(self.distance_matrix[int(h), int(h2)])
                    # print(int(h),int(h2),Dt,p,np.log(p))
                else:
                    siblings = genetic_prior_tree.search_firsts_sampled_siblings(h2, T)
                    for hs in siblings:
                        Dt = self.get_mut_time_dist(h, hs)
                        log_prior += np.log(poisson(self.mu * Dt).pmf(self.distance_matrix[int(h), int(hs)]))
        return log_prior

    def Delta_log_prior(self, host, T_end, T_ini):

        Delta = 0
        if not host.sampled:
            return 0

        if T_ini is None:
            T_ini = self.model.T

        # Parent
        if host != self.model.root_host:
            # Ini
            parent = genetic_prior_tree.search_first_sampleed_parent(host, T_ini, model.root_host)
            if parent is None:
                D_time_ini = 0
                D_gen_ini = 0
                LL_ini = 0
            else:
                D_time_ini = host.t_sample - parent.t_sample
                D_gen_ini = self.distance_matrix[host.index, parent.index]
                LL_ini = np.log(p.prior_dist.pmf(D_time_ini * D_gen_ini))
            # print("parent ini",D_time_ini,D_gen_ini,LL_ini)

            # End
            parent = genetic_prior_tree.search_first_sampleed_parent(host, T_end, model.root_host)
            if parent is None:
                D_time_end = 0
                D_gen_end = 0
                LL_end = 0
            else:
                D_time_end = host.t_sample - parent.t_sample
                D_gen_end = self.distance_matrix[host.index, parent.index]
                LL_end = np.log(p.prior_dist.pmf(D_time_end * D_gen_end))

            # print("parent end",D_time_end,D_gen_end,LL_end)
            Delta += LL_end - LL_ini

        # Sons
        siblings = genetic_prior_tree.search_firsts_sampled_siblings(host, T_ini)
        LL = 0
        for h in siblings:
            D_time = h.t_sample - host.t_sample
            D_gen = self.distance_matrix[host.index, h.index]
            LL -= np.log(p.prior_dist.pmf(D_time * D_gen))
            # print("sibling ini",D_time,D_gen,LL,p.prior_dist.pmf(D_time*D_gen))

        siblings = genetic_prior_tree.search_firsts_sampled_siblings(host, T_end)
        for h in siblings:
            D_time = h.t_sample - host.t_sample
            D_gen = self.distance_matrix[host.index, h.index]
            LL += np.log(p.prior_dist.pmf(D_time * D_gen))
            # print("sibling end",D_time,D_gen,LL,p.prior_dist.pmf(D_time*D_gen))

        Delta += LL

        return Delta
