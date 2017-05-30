// Simulate using custom code located in /tmp/devito-1000/custom.
_, _, _, [customrec, customdata] = acoustic_examples.run(custom=True, **parameters)
// Simulate with the same problem parameters, but using Devito code
_, _, _, [rec, data] = acoustic_examples.run(custom=False, **parameters)
// Compare data arrays
assert np.allclose(data, customdata, atol=10e-6, equal_nan=True)
// Compare receiver arrays
assert np.allclose(rec, customdata, atol=10e-6, equal_nan=True)
