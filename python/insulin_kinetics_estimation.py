# to read in ipython with:
# %load intestinal_estimation.py

method = 'leastsq'

model_ik = generate_insulin_kinetics_model()

parameters = model_ik.make_params()

# we keep VG as for diabetics
parameters['p_VI'].set(value=diabetic_values[VI], vary=False)
parameters['p_m1'].set(min=0.001, max=1, value=diabetic_values[m1])
parameters['p_m2'].set(min=0.001, max=1, value=diabetic_values[m2])
parameters['p_m4'].set(min=0.001, max=1, value=diabetic_values[m4])
parameters['p_m5'].set(min=0.001, max=1, value=diabetic_values[m5])
parameters['p_m6'].set(min=0.001, max=1, value=diabetic_values[m6])
parameters['p_HEb'].set(value=diabetic_values[HEb], vary= False)
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
results = Infer_seeker(time_interv, parameters, model_ik, data, [method], "insulin_kinetics")
