def test_star_import():
    namespace = {}
    exec('from transmission_models import *', namespace)
    expected = [
        'host',
        'didelot_unsampled',
        'genetic_prior_tree',
        'location_distance_prior_tree',
        'same_location_prior_tree',
        'MCMC',
    ]
    for symbol in expected:
        assert symbol in namespace, f"{symbol} not imported with star import"

if __name__ == "__main__":
    test_star_import()
    print("Star import test passed!") 