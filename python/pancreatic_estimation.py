# to read in ipython with:
# %load intestinal_estimation.py

method = 'leastsq'

model_pancreatic = generate_pancreatic_model()

parameters = model_pancreatic.make_params()

# we keep VG as for diabetics
parameters['p_K'].set(min=0.1, max=10, value=diabetic_values[K])
parameters['p_α'].set(min=0.001, max=0.1, value=diabetic_values[α])
parameters['p_β'].set(min=0.0001, max=1, value=diabetic_values[β])
parameters['p_γ'].set(value=diabetic_values[γ], vary=False)
# We should let Gb as a varying parameter: if we fix it to the healthy_values[Gb] = 91
# then Gpb is also fixed to 91 * diabetic_values[VG] = 91 (assuming that VG is also fixed),
# that is Gpb = 135.59 which can never reach healthy_values[Gpb] = 171...
parameters['p_Gb'].set(min=50, max=300, value=diabetic_values[Gb])

parameters.pretty_print()

print("data generation...")
# on suppose que le sujet seins et diabétique font le même poids
healthy_values[BW] = diabetic_values[BW]
healthy_values[VG] = diabetic_values[VG]
healthy_values[VI] = diabetic_values[VI]

data = generate_simulated_data(healthy_values, nb_points=15)
time_interv = scipy.linspace(0, 300, 15)
print("Parameter estimation with " + method + " method...")
results = Infer_seeker(time_interv, parameters, model_pancreatic, data, [method], "pancreatic")
