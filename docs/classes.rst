Classes
=======

Host
----

.. automodule:: transmission_models.classes.host
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: nx,np,plt,choice,randint,random,sample,choices,nbinom,gamma,binom,expon,norm,Line2D,pd,graphviz_layout,imageio

Didelot Unsampled
-----------------

.. autoclass:: transmission_models.classes.didelot_unsampled.didelot_unsampled
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

   **Core Methods**

   .. automethod:: __init__
   .. automethod:: add_root
   .. automethod:: successors
   .. automethod:: parent
   .. automethod:: out_degree
   .. automethod:: choose_successors

   **Tree Structure Methods**

   .. automethod:: get_root_subtrees
   .. automethod:: get_unsampled_hosts
   .. automethod:: get_candidates_to_chain
   .. automethod:: get_N_candidates_to_chain

   **Likelihood Methods**

   .. automethod:: get_sampling_model_likelihood
   .. automethod:: get_sampling_model_log_likelihood
   .. automethod:: get_offspring_model_likelihood
   .. automethod:: get_offspring_model_log_likelihood
   .. automethod:: get_infection_model_likelihood
   .. automethod:: get_infection_model_log_likelihood
   .. automethod:: log_likelihood_host
   .. automethod:: log_likelihood_hosts_list
   .. automethod:: log_likelihood_transmission_tree
   .. automethod:: get_log_likelihood_transmission

   **Delta Methods (for MCMC)**

   .. automethod:: Delta_log_sampling
   .. automethod:: Delta_log_offspring
   .. automethod:: Delta_log_infection
   .. automethod:: Delta_log_likelihood_host

   **MCMC Step Methods**

   .. automethod:: infection_time_from_sampling_step
   .. automethod:: infection_time_from_infection_model_step
   .. automethod:: add_unsampled_with_times
   .. automethod:: remove_unsampled_with_times
   .. automethod:: add_remove_step
   .. automethod:: MCMC_step

   **Prior Methods**

   .. automethod:: add_genetic_prior
   .. automethod:: add_same_location_prior
   .. automethod:: compute_Delta_loc_prior

   **Utility Methods**

   .. automethod:: create_transmision_phylogeny_nets
   .. automethod:: get_newick
   .. automethod:: save_json
   .. automethod:: show_log_likelihoods

MCMC
----

.. automodule:: transmission_models.classes.mcmc.mcmc
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: np,utils,hierarchy_pos,hierarchy_pos_times,plot_transmision_network,tree_to_newick,search_firsts_sampled_siblings,du,tree_slicing_step

Priors
------

.. automodule:: transmission_models.classes.genetic_prior
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: choice,randint,random,sample,choices,GAMMA,nbinom,gamma,binom,expon,poisson,np,utils,combinations,nx

.. automodule:: transmission_models.classes.location_prior
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: choice,randint,random,sample,choices,GAMMA,nbinom,gamma,binom,expon,poisson,np,utils,combinations,nx

Module Documentation
--------------------

Classes Module
~~~~~~~~~~~~~~

.. automodule:: transmission_models.classes
   :members:
   :undoc-members:
   :show-inheritance:

MCMC Module
~~~~~~~~~~~

.. automodule:: transmission_models.classes.mcmc
   :members:
   :undoc-members:
   :show-inheritance: 