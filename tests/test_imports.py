"""
Test that the package can be imported correctly.
"""

def test_import_package():
    """Test that the main package can be imported."""
    import transmission_models
    assert transmission_models is not None

def test_import_modules():
    """Test that all main modules can be imported."""
    from transmission_models.classes import host
    from transmission_models.classes import didelot_unsampled
    from transmission_models.classes import genetic_prior_tree
    from transmission_models import utils
    
    assert host is not None
    assert didelot_unsampled is not None
    assert genetic_prior_tree is not None
    assert utils is not None

def test_import_mcmc():
    """Test that MCMC module can be imported."""
    from transmission_models.classes.mcmc import mcmc
    assert mcmc is not None

def test_import_utils():
    """Test that utility functions can be imported."""
    from transmission_models.utils import tree_to_newick, plot_transmision_network
    from transmission_models.utils import tree_slicing_step
    assert tree_to_newick is not None
    assert plot_transmision_network is not None
    assert tree_slicing_step is not None

if __name__ == "__main__":
    test_import_package()
    test_import_modules()
    test_import_mcmc()
    test_import_utils()
    print("All import tests passed!") 