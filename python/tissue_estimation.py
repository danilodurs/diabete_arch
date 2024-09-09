# to read in ipython with:
# %load intestinal_estimation.py

method = 'leastsq'

model_tissue = generate_tissue_model()

parameters = model_tissue.make_params()

# we keep VG as for diabetics
parameters['p_Fcns'].set(value=diabetic_values[Fcns], vary=False)
parameters['p_Vm0'].set(min=1, max=10, value=diabetic_values[Vm0])
parameters['p_Vmx'].set(min=0.01, max=0.1, value=diabetic_values[Vmx])
parameters['p_Km0'].set(min=100, max=1000, value=diabetic_values[Km0])
parameters['p_p2U'].set(min=0.01, max=0.5, value=diabetic_values[p2U])
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
results = Infer_seeker(time_interv, parameters, model_tissue, data, [method], "tissue")
