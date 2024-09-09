# to read in ipython with:
# %load intestinal_estimation.py

method = 'leastsq'

model_gk = generate_glucose_kinetics_model()

parameters = model_gk.make_params()

# we keep VG as for diabetics
parameters['p_VG'].set(value=diabetic_values[VG], vary=False)
parameters['p_k1'].set(min=0.001, max=1, value=diabetic_values[k1])
parameters['p_k2'].set(min=0.001, max=1, value=diabetic_values[k2])
# We should let Gb as a varying parameter: if we fix it to the healthy_values[Gb] = 91
# then Gpb is also fixed to 91 * diabetic_values[VG] = 91 (assuming that VG is also fixed),
# that is Gpb = 135.59 which can never reach healthy_values[Gpb] = 171...
parameters['p_Gb'].set(min=50, max=300, value=healthy_values[Gb])

parameters.pretty_print()

print("data generation...")
data = generate_simulated_data(healthy_values, nb_points=15)
time_interv = scipy.linspace(0, 300, 15)
print("Parameter estimation with " + method + " method...")
results = Infer_seeker(time_interv, parameters, model_gk, data, [method], "glucose_kinetics")
