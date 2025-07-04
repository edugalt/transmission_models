import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from random import choice,randint,random,sample,choices
from scipy.stats import nbinom, gamma, binom, expon, norm
from matplotlib.lines import Line2D
import pandas as pd
from networkx.drawing.nx_pydot import graphviz_layout



try:
    import imageio
except ImportError:
    pass


#CLASSES
class host:
    """
    Represents a host that has been infected with a virus.

    A host object contains information about an infected individual,
    including their genetic data, infection time, sampling time,
    and other attributes.

    Attributes
    ----------
    index : int
        The index of the host.
    sampled : bool
        Indicates whether the host has been sampled or not.
    genetic_data : list
        The genetic data of the host.
    dict_attributes : dict
        A dictionary to store additional attributes.
    t_inf : int
        Time of infection.
    t_sample : int, optional
        The time the host was sampled.
    id : str
        The identifier of the host.

    Methods
    -------
    t_inf : property
        Getter and setter for the time of infection attribute.
    get_genetic_str() : str
        Returns the genetic data as a string.
    __str__() : str
        Returns a string with the id of the host.
    __int__() : int
        Returns the index of the host.

    Examples
    --------
    >>> h = host('host1', 1, ['A', 'T', 'C', 'G'], 10, t_sample=15)
    >>> print(h.t_inf)
    10
    >>> h.t_inf = 20
    >>> print(h.t_inf)
    20
    >>> print(h.get_genetic_str())
    ATCG
    >>> print(h)
    host1

    Notes
    -----
    This class follows the Python naming convention for class names
    (using PascalCase).
    """

    def __init__(self, id, index, genetic_data=[], t_inf=0, t_sample=None):
        """
        Initialize a new instance of the Host class.

        Parameters
        ----------
        id : str
            The id of the host.
        index : int
            The index of the host.
        genetic_data : list, optional
            The genetic data of the host. Defaults to an empty list.
        t_inf : int, optional
            Time of infection. Defaults to 0.
        t_sample : int, optional
            The time the host was sampled. Defaults to None.
        """
        if t_sample is not None:
            self.sampled = True
        else:
            self.sampled = False
        self.t_sample = t_sample
        self.genetic_data = genetic_data
        self._t_inf = t_inf
        self.index = int(index)
        self.id = id
        self.dict_attributes = {}

    @property
    def t_inf(self):
        """
        Getter for the time of infection attribute.

        Returns
        -------
        int
            The time of infection.
        """
        return self._t_inf

    @t_inf.setter
    def t_inf(self, t_inf):
        """
        Setter for the time of infection attribute.

        Parameters
        ----------
        t_inf : int
            The time of infection.
        """
        self._t_inf = t_inf

    def get_genetic_str(self):
        """
        Return the genetic data of the host as a string.

        Returns
        -------
        str
            The genetic data as a string.
        """
        return "".join(self.genetic_data)

    def __str__(self):
        """
        Return a string with the id of the host.

        Returns
        -------
        str
            The id of the host.
        """
        return str(self.id)

    def __int__(self):
        """
        Return the index of the host.

        Returns
        -------
        int
            The index of the host.
        """
        return self.index


#Functions
def create_genome(chain_length):
    """
    Create a random genome sequence of specified length.

    Parameters
    ----------
    chain_length : int
        The length of the genome sequence to create.

    Returns
    -------
    list
        A list of random nucleotides (A, G, C, T) of length chain_length.

    Examples
    --------
    >>> genome = create_genome(10)
    >>> print(genome)
    ['A', 'T', 'C', 'G', 'A', 'T', 'C', 'G', 'A', 'T']
    """
    return [choice("AGCT") for i in range(chain_length)]

def binom_mutation(chain_length, p, genome):
    """
    Perform binomial mutation on a given genome.

    This function generates changes in a genome by randomly selecting 'k' positions
    to mutate, where 'k' follows a binomial distribution with parameters
    'chain_length' and 'p'. The elements at the selected positions are replaced
    with new randomly chosen nucleotides.

    Parameters
    ----------
    chain_length : int
        The length of the genome chain.
    p : float
        The probability of mutation for each element in the chain.
    genome : str or list
        The original genome sequence.

    Returns
    -------
    list
        The mutated genome sequence.

    Notes
    -----
    The function operates as follows:

    1. Calculates the number of positions to mutate, 'k', by sampling from a
       binomial distribution with 'chain_length' trials and success probability 'p'.

    2. Randomly selects 'k' positions from the range [0, chain_length) without replacement.

    3. Creates a new list 'new_genome' from the original genome.

    4. Iterates over the selected positions and replaces the corresponding elements
       in 'new_genome' with randomly chosen nucleotides based on the original
       nucleotide at that position:

       - If the original nucleotide is 'A', it is replaced with a randomly chosen
         nucleotide from 'CTG'.
       - If the original nucleotide is 'C', it is replaced with a randomly chosen
         nucleotide from 'ATG'.
       - If the original nucleotide is 'T', it is replaced with a randomly chosen
         nucleotide from 'ACG'.
       - If the original nucleotide is 'G', it is replaced with a randomly chosen
         nucleotide from 'ACT'.

    5. Returns the mutated genome sequence as 'new_genome'.

    Examples
    --------
    >>> genome = ['A', 'T', 'C', 'G', 'G', 'A', 'T', 'C', 'G', 'A']
    >>> mutated_genome = binom_mutation(len(genome), 0.2, genome)
    >>> print(mutated_genome)
    ['A', 'T', 'C', 'A', 'G', 'A', 'T', 'C', 'G', 'A']

    See Also
    --------
    one_mutation : Perform a single mutation on a genome
    """
    k = np.random.binomial(n=chain_length, p=p, size=1)[0]
    to_change = sample(range(chain_length),k)
    new_genome = list(genome)
    for i in to_change:
        g = new_genome[i]
        if g=="A":
            new_genome[i] = choice("CTG")
        if g=="C":
            new_genome[i] = choice("ATG")
        if g=="T":
            new_genome[i] = choice("ACG")
        if g=="G":
            new_genome[i] = choice("ACT")
    return new_genome


def one_mutation(chain_length, p, genome):
    """
    Perform one mutation on a given genome.

    This function generates a single mutation in a genome by randomly selecting
    one position to mutate. The selected position is replaced with a new
    randomly chosen nucleotide.

    Parameters
    ----------
    chain_length : int
        The length of the genome chain.
    p : float
        The probability of mutation for each element in the chain.
    genome : str or list
        The original genome sequence.

    Returns
    -------
    list
        The mutated genome sequence.

    Notes
    -----
    The function operates as follows:

    1. Randomly selects one position from the range [0, chain_length) to mutate.

    2. Creates a new list 'new_genome' from the original genome.

    3. Checks the original nucleotide at the selected position and replaces it
       with a randomly chosen nucleotide based on the following rules:

       - If the original nucleotide is 'A', it is replaced with a randomly chosen
         nucleotide from 'CTG'.
       - If the original nucleotide is 'C', it is replaced with a randomly chosen
         nucleotide from 'ATG'.
       - If the original nucleotide is 'T', it is replaced with a randomly chosen
         nucleotide from 'ACG'.
       - If the original nucleotide is 'G', it is replaced with a randomly chosen
         nucleotide from 'ACT'.

    4. Returns the mutated genome sequence as 'new_genome'.

    Examples
    --------
    >>> genome = ['A', 'T', 'C', 'G', 'G', 'A', 'T', 'C', 'G', 'A']
    >>> mutated_genome = one_mutation(len(genome), 0.2, genome)
    >>> print(mutated_genome)
    ['A', 'T', 'C', 'A', 'G', 'A', 'T', 'C', 'G', 'T']

    See Also
    --------
    binom_mutation : Perform binomial mutation on a genome
    """
    to_change = sample(range(chain_length),1)
    new_genome = list(genome)
    for i in to_change:
        g = new_genome[i]
        if g=="A":
            new_genome[i] = choice("CTG")
        if g=="C":
            new_genome[i] = choice("ATG")
        if g=="T":
            new_genome[i] = choice("ACG")
        if g=="G":
            new_genome[i] = choice("ACT")
    return new_genome


def average_mutations(mu, P_mut, tau, Dt, host_genetic):
    """
    Generate a list of mutations proportional to a time interval.

    The number of mutations is proportional to a given time interval (Dt)
    where the proportion factor is the mutation rate (mu).

    Parameters
    ----------
    mu : float
        The mutation rate.
    P_mut : float
        The probability of mutation.
    tau : float
        The current time.
    Dt : float
        The time interval.
    host_genetic : list
        The genetic sequence of the host.

    Returns
    -------
    tuple
        A tuple containing:

        - mutations : list
            List of mutated genetic sequences.
        - t_mutations : list
            List of mutation times.

    Notes
    -----
    The function calculates the number of mutations as int(mu * Dt / P_mut)
    and generates that many mutations using the one_mutation function.
    """
    #Mutations
    n_mut = int(mu*(Dt)/P_mut)#number of mutations

    mutations = [one_mutation(len(host_genetic),P_mut,host_genetic)] #First mutation
    tau += t_inf/n_mut
    t_mutations = [tau]#List of mutations times
    if int(np.floor(n_mut)) > 0:
        for l in range(n_mut-1):
            mutations.append(one_mutation(len(host_genetic),P_mut,mutations[-1])) #First mutatiom
            tau += t_inf/n_mut
            t_mutations.append(tau)
    else:
        return [host_genetic],[]
    return mutations,t_mutations


