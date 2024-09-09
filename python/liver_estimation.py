# to read in ipython with:
# %load intestinal_estimation.py

method = 'leastsq'

model_liver = generate_liver_model()

parameters = model_liver.make_params()

# we keep VG as for diabetics
parameters['p_kp1'].set(min=2, max=4, value=diabetic_values[kp1])
parameters['p_kp2'].set(min=0.0001, max=0.005, value=diabetic_values[kp2])
parameters['p_kp3'].set(min=0.0001, max=0.005, value=diabetic_values[kp3])
parameters['p_kp4'].set(min=0.001, max=0.5, value=diabetic_values[kp4])
parameters['p_ki'].set(min=0.0001, max=0.05, value=diabetic_values[ki])
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
results = Infer_seeker(time_interv, parameters, model_liver, data, [method], "liver")
