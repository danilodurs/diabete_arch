# to read in ipython with:
# %load intestinal_estimation.py

method = 'leastsq'

model_intestinal = generate_intestinal_model()

parameters = model_intestinal.make_params()
# parameters['p_kmax'].set(min=0.001, max=1, value=diabetic_values[kmax])
# parameters['p_kmin'].set(min=0.001, max=1, value=diabetic_values[kmin])
# parameters['p_kab'].set(min=0.001, max=1, value=diabetic_values[kab])
# parameters['p_kgri'].set(min=0.001, max=1, value=diabetic_values[kgri])
# parameters['p_f'].set(min=0.01, max=1, value=diabetic_values[f])
# # parameters['p_Gb'].set(min=50, max=300, value=diabetic_values[Gb])
# parameters['p_Gb'].set(value=healthy_values[Gb], vary=False)
# parameters['p_b'].set(min=0.01, max=1, value=diabetic_values[b])
# parameters['p_c'].set(min=0.01, max=1, value=diabetic_values[c])

parameters['p_kmax'].set(min=0.001, max=1, value=diabetic_values[kmax], vary=False)
parameters['p_kmin'].set(min=0.001, max=1, value=diabetic_values[kmin], vary=False)
parameters['p_kab'].set(min=0.001, max=1, value=diabetic_values[kab], vary=False)
parameters['p_kgri'].set(min=0.001, max=1, value=diabetic_values[kgri], vary=False)
parameters['p_f'].set(min=0.01, max=1, value=diabetic_values[f])
# We should let Gb as a varying parameter: if we fix it to the healthy_values[Gb] = 91
# then Gpb is also fixed to 91 * diabetic_values[VG] = 91 (assuming that VG is also fixed),
# that is Gpb = 135.59 which can never reach healthy_values[Gpb] = 171...
parameters['p_Gb'].set(min=50, max=300, value=diabetic_values[Gb])
# parameters['p_Gb'].set(value=healthy_values[Gb])
parameters['p_b'].set(min=0.01, max=1, value=diabetic_values[b], vary=False)
parameters['p_c'].set(min=0.01, max=1, value=diabetic_values[c], vary=False)

parameters.pretty_print()

print("data generation...")
data = generate_simulated_data(healthy_values, nb_points=15)
time_interv = scipy.linspace(0, 300, 15)
print("Parameter estimation with " + method + " method...")
results = Infer_seeker(time_interv, parameters, model_intestinal, data, [method])
