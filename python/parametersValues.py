from model import *

parameter_values = dict()

# parameter values from Dalla Mann for healthy subject

parameter_values['normal'] = dict()
parameter_values['normal'].update({
    # glucose kinetics
    VG : 1.88,
    k1 : 0.065,
    k2 : 0.079,
    # insulin kinetics
    VI : 0.05,
    m1 : 0.19,
    m2 : 0.484,
    m4 : 0.194,
    m5 : 0.0304,
    m6 : 0.6471,
    HEb: 0.6,
    # rate of appearance
    kmax: 0.0558,
    kmin: 0.008,
    kab: 0.057,
    kgri: 0.0558,
    f: 0.9,
    b: 0.82,
    c: 0.01,
    # Endogenous production
    kp1: 2.7,
    kp2: 0.0021,
    kp3: 0.009,
    kp4: 0.0618,
    ki: 0.0079,
    # Utilization
    Fcns: 1,
    Vm0: 2.5,
    Vmx: 0.047,
    Km0: 225.59,
    Kmx: 0,
    p2U: 0.0331,
    # Secretion
    K: 2.3,
    α: 0.05,
    β: 0.11,
    γ: 0.5,
    # Renal Excretion
    ke1: 0.0005,
    ke2: 339,
    # Body weight (not given in Table 2 but in text, +- 1kg)
    BW: 78,
    # Basal glucose concentration (according to Fig. 5)
    Gb: 91
})

# parameter values from Dalla Mann for diabetic subject

parameter_values['diabetic'] = dict()
parameter_values['diabetic'].update({
    # glucose kinetics
    VG : 1.49,
    k1 : 0.042,
    k2 : 0.071,
    # insulin kinetics
    VI : 0.04,
    m1 : 0.379,
    m2 : 0.673,
    m4 : 0.269,
    m5 : 0.0526,
    m6 : 0.8118,
    HEb: 0.6,
    # rate of appearance
    kmax: 0.0465,
    kmin: 0.0076,
    kab: 0.023,
    kgri: 0.0465,
    f: 0.9,
    b: 0.68,
    c: 0.09,
    # Endogenous production
    kp1: 3.09,
    kp2: 0.0007,
    kp3: 0.005,
    kp4: 0.0786,
    ki: 0.0066,
    # Utilization
    Fcns: 1,
    Vm0: 4.65,
    Vmx: 0.034,
    Km0: 466.21,
    Kmx: 0,
    p2U: 0.0840,
    # Secretion
    K: 0.99,
    α: 0.013,
    β: 0.05,
    γ: 0.5,
    # Renal Excretion
    ke1: 0.0007,
    ke2: 269,
    # Body weight (not given in Table 2 but in text, +- 5kg)
    BW: 91,
    # Basal glucose concentration (according to Fig. 5)
    Gb: 242
})

# for RYGB data, we take the fasting glucose of 125 mg/dl ie. the average clinical data for fat and diabetic subjects

# RYGB_A: we take the parameter values for diabetic subject and just modify the
# fraction of intestinal absorption f from 90% to 50%.

parameter_values['RYGB_A'] = parameter_values['diabetic'].copy()
parameter_values['RYGB_A'].update({f : 0.5, Gb : 125})

# RYGB_B: in addition we take all the parameter for absorption and gastric emptying from healthy subject
parameter_values['RYGB_B'] = parameter_values['diabetic'].copy()
parameter_values['RYGB_B'].update({f : 0.5, Gb : 125, kmax: 0.0558, kmin: 0.008, kgri: 0.0558, b: 0.82, c: 0.01 })

# RYGB_C: same as RYGB except that we increase kmax and keep gastric emptying at kmax
# => kmax = kmin = 0.1
parameter_values['RYGB_C'] = parameter_values['diabetic'].copy()
parameter_values['RYGB_C'].update({f : 0.5, Gb : 125, kmax: 0.1, kmin: 0.1 })
