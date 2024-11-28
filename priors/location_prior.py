
from random import choice,randint,random,sample,choices
from scipy.special import gamma as GAMMA
from scipy.stats import nbinom, gamma, binom, expon, poisson
import numpy as np
import transmission_models.utils as utils
from itertools import combinations
import networkx as nx

class location_distance_prior_tree():
    def __init__(self, model, mu, distance_matrix):
        self.mu = mu
        self.distance_matrix = distance_matrix

        self.prior_dist = expon(mu)

        self.model = model

        self.correction_LL = 0
        self.log_prior = 0

    @staticmethod
    def search_firsts_sampled_siblings(host, T):

        sampled_hosts = []
        for h in T.successors(host):
            if h.sampled:
                sampled_hosts.append(h)
            else:
                sampled_hosts += location_distance_prior_tree.search_firsts_sampled_siblings(h, T)

        return sampled_hosts

    @staticmethod
    def search_first_sampleed_parent(host, T, root):

        if host == root:
            return None

        parent = next(T.predecessors(host))

        if not parent.sampled:
            return location_distance_prior_tree.search_first_sampleed_parent(parent, T, root)
        else:
            return parent
    @staticmethod
    def get_mut_time_dist(hp, hs):
        return (hs.t_sample + hp.t_sample - 2 * hp.t_inf)

    def get_closest_sampling_siblings(self,T=None):
        if T is None:
            T = self.model.T
            self.model.get_root_subtrees()
            roots_subtrees = self.model.get_root_subtrees()
        else:
            roots_subtrees = utils.search_firsts_sampled_siblings(self.model.root_host, T)
        non_observed = list(roots_subtrees)
        LL_correction = 0
        # print(roots_subtrees[::-1],shuffle(roots_subtrees[::-1]))
        for h in roots_subtrees:
            if h not in non_observed: continue
            N_samp = 0
            parent = h
            relatives = []
            jumped = False
            closest = None
            while N_samp == 0:
                parent = list(T.predecessors(parent))[0]
                closest = None
                # print(h,parent)
                if parent != self.model.root_host:
                    if T.out_degree(parent) == 1:
                        # parent = model.parent(parent)
                        # print(h,parent)
                        jumped = True
                        continue
                    # elif model.out_degree(parent)==2 and not jumped:
                    #     parent = model.parent(parent)
                    #     jumped = True
                    #     continue

                for h2 in T.successors(parent):
                    # print("-"*6,h,parent,h2)
                    if h2.sampled:
                        if h2 == h:
                            continue
                        elif not h2.sampled:#Eliminar el if seria sensato
                            continue
                        else:
                            # if h2 not in non_observed: continue
                            # print("KAKA2",h2,parent)
                            N_samp += 1
                            relatives.append(h2)
                            # non_observed.remove(h2)
                    # else:
                    # print("KAKA",h2,parent)
                if parent == self.model.root_host: break
            # non_observed.remove(h)
            if not relatives: continue
            closest = min(relatives, key=lambda h2: h2.t_sample)
            Dt = h.t_inf#(closest.t_sample + h.t_sample - 2 * parent.t_inf)
            # print(f"\t\t{int(pair[0]),int(pair[1])},{Dt`=}")
            LL_correction += np.log(expon(-self.mu * Dt).pmf(self.distance_matrix[int(h), int(closest)]))

            # print("----->",h,closest,parent)
        return LL_correction


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
                siblings = location_distance_prior_tree.search_firsts_sampled_siblings(h2, T)
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
                parent = location_distance_prior_tree.search_first_sampleed_parent(host, T, self.model.root_host)
                if parent is not None:
                    # print(f"{parent}-->{host}")
                    Dt = host.t_sample - parent.t_sample
                    log_prior += np.log(poisson(self.mu * Dt).pmf(self.distance_matrix[int(host), int(parent)]))

        return log_prior

    def log_prior_T(self, T, update_up=True,verbose=False):
        self.log_prior = 0
        suma = 0
        for h in T:
            if not h.sampled: continue
            for h2 in T[h]:
                if h2.sampled:
                    # print(f"{h}-->{h2}")
                    Dt = h.t_inf#self.get_mut_time_dist(h, h2)

                    log_L = np.log(expon(-self.mu * Dt).pdf(self.distance_matrix[int(h), int(h2)]))
                    if verbose: print(f"{h}-->{h2} {Dt=} {log_L=}")
                    suma += log_L
                    self.log_prior += log_L
                    # p = poisson(self.mu * Dt).pmf(self.distance_matrix[int(h), int(h2)])
                    # print(int(h),int(h2),Dt,p,np.log(p))
                else:
                    siblings = location_distance_prior_tree.search_firsts_sampled_siblings(h2, T)
                    for hs in siblings:
                        Dt = h.t_inf#self.get_mut_time_dist(h, hs)
                        log_L = np.log(expon(-self.mu * Dt).pdf(self.distance_matrix[int(h), int(hs)]))
                        if verbose: print(f"{h}-->{hs} (jumped) {Dt=} {log_L=}")
                        suma += log_L
                        self.log_prior += log_L
        if verbose: print(f"{suma=}")

        if update_up:
            # self.model.get_root_subtrees()
            LL_correction = self.get_closest_sampling_siblings(T)

            self.correction_LL = LL_correction
            self.log_prior += LL_correction
        else:
            self.correction_LL = 0
        # print(f"{self.correction_LL-self.log_prior=},{self.correction_LL=}")
        return self.log_prior


class same_location_prior_tree():
    def __init__(self, model, log_K, distance_matrix):
        self.log_K = log_K
        self.distance_matrix = distance_matrix


        self.model = model

        self.correction_LL = 0
        self.log_prior = 0

    @staticmethod
    def search_firsts_sampled_siblings(host, T):

        sampled_hosts = []
        for h in T.successors(host):
            if h.sampled:
                sampled_hosts.append(h)
            else:
                sampled_hosts += location_distance_prior_tree.search_firsts_sampled_siblings(h, T)

        return sampled_hosts

    @staticmethod
    def search_first_sampleed_parent(host, T, root):

        if host == root:
            return None

        parent = next(T.predecessors(host))

        if not parent.sampled:
            return location_distance_prior_tree.search_first_sampleed_parent(parent, T, root)
        else:
            return parent
    @staticmethod
    def get_mut_time_dist(hp, hs):
        return (hs.t_sample + hp.t_sample - 2 * hp.t_inf)

    def get_closest_sampling_siblings(self,T=None):
        if T is None:
            T = self.model.T
            self.model.get_root_subtrees()
            roots_subtrees = self.model.get_root_subtrees()
        else:
            roots_subtrees = utils.search_firsts_sampled_siblings(self.model.root_host, T)
        non_observed = list(roots_subtrees)
        LL_correction = 0
        # print(roots_subtrees[::-1],shuffle(roots_subtrees[::-1]))
        for h in roots_subtrees:
            if h not in non_observed: continue
            N_samp = 0
            parent = h
            relatives = []
            jumped = False
            closest = None
            while N_samp == 0:
                parent = list(T.predecessors(parent))[0]
                closest = None
                # print(h,parent)
                if parent != self.model.root_host:
                    if T.out_degree(parent) == 1:
                        # parent = model.parent(parent)
                        # print(h,parent)
                        jumped = True
                        continue
                    # elif model.out_degree(parent)==2 and not jumped:
                    #     parent = model.parent(parent)
                    #     jumped = True
                    #     continue

                for h2 in T.successors(parent):
                    # print("-"*6,h,parent,h2)
                    if h2.sampled:
                        if h2 == h:
                            continue
                        elif not h2.sampled:#Eliminar el if seria sensato
                            continue
                        else:
                            # if h2 not in non_observed: continue
                            # print("KAKA2",h2,parent)
                            N_samp += 1
                            relatives.append(h2)
                            # non_observed.remove(h2)
                    # else:
                    # print("KAKA",h2,parent)
                if parent == self.model.root_host: break
            # non_observed.remove(h)
            if not relatives: continue
            closest = min(relatives, key=lambda h2: h2.t_sample)
            Dt = h.t_inf#(closest.t_sample + h.t_sample - 2 * parent.t_inf)
            # print(f"\t\t{int(pair[0]),int(pair[1])},{Dt`=}")
            if self.distance_matrix[int(h), int(closest)]> 0:
                LL_correction += -self.log_K

            # print("----->",h,closest,parent)
        return LL_correction


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
                siblings = location_distance_prior_tree.search_firsts_sampled_siblings(h2, T)
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
                parent = location_distance_prior_tree.search_first_sampleed_parent(host, T, self.model.root_host)
                if parent is not None:
                    # print(f"{parent}-->{host}")
                    Dt = host.t_sample - parent.t_sample
                    log_prior += np.log(poisson(self.mu * Dt).pmf(self.distance_matrix[int(host), int(parent)]))

        return log_prior

    def log_prior_T(self, T, update_up=True,verbose=False):
        self.log_prior = 0
        suma = 0
        for h in T:
            if not h.sampled: continue
            for h2 in T[h]:
                if h2.sampled:
                    log_L = 0
                    # print(f"{h}-->{h2}")
                    Dt = h.t_inf#self.get_mut_time_dist(h, h2)
                    if self.distance_matrix[int(h), int(h2)] > 0:
                        log_L = -self.log_K

                    if verbose: print(f"{h}-->{h2} {Dt=} {log_L=}")
                    suma += log_L
                    self.log_prior += log_L
                    # p = poisson(self.mu * Dt).pmf(self.distance_matrix[int(h), int(h2)])
                    # print(int(h),int(h2),Dt,p,np.log(p))
                else:
                    siblings = location_distance_prior_tree.search_firsts_sampled_siblings(h2, T)
                    for hs in siblings:
                        Dt = h.t_inf#self.get_mut_time_dist(h, hs)
                        log_L = 0
                        if self.distance_matrix[int(h), int(hs)]> 0:
                            log_L = -self.log_K
                        if verbose: print(f"{h}-->{hs} (jumped) {Dt=} {log_L=}")
                        suma += log_L
                        self.log_prior += log_L
        if verbose: print(f"{suma=}")

        if update_up:
            # self.model.get_root_subtrees()
            LL_correction = self.get_closest_sampling_siblings(T)

            self.correction_LL = LL_correction
            self.log_prior += LL_correction
        else:
            self.correction_LL = 0
        # print(f"{self.correction_LL-self.log_prior=},{self.correction_LL=}")
        return self.log_prior

