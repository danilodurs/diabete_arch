import sympy

from sympy.functions.elementary.hyperbolic import tanh

# regular expressions
import re

# if this constant is true then string encoding as latex expression
LATEX = False

def unlatex(exp):
    '''
    Unwrap the latex environment from a string expression.
    Input:  string expression
    Output: string expression
    '''
    patt1 = '(\_\{)|(\\\\[A-Za-z]+\{)|(\})|((?<!DiracDelta)\(t\))|(\_)'
    patt2 = "\\\\alpha"
    patt3 = "\\\\beta"
    patt4 = "\\\\gamma"
    out = re.sub(patt1, "", exp)
    out = re.sub(patt2, "α", out)
    out = re.sub(patt3, "β", out)
    out = re.sub(patt4, "γ", out)
    return(out)

class EmptyObject(): pass

model = EmptyObject()

symbols = EmptyObject()
model.symbols = symbols

def SSymbol(s, **kwargs):
    return sympy.Symbol(s if LATEX else unlatex(s), nonnegative=True, **kwargs)

def SFunction(sf, **kwargs):
    return sympy.Function(sf if LATEX else unlatex(sf), nonnegative=True, **kwargs)

#############################
# Definition of variables
#############################

t      = sympy.Symbol('t')

Gp     = SFunction('G_{\\textit{p}}')(t) # glucose masses in plasma and rapidly equilibrating tissues (mg/kg)
Gt     = SFunction('G_{\\textit{t}}')(t) # glucose masses in plasma and slowly equilibrating tissues (mg/kg)
G      = SFunction('G')(t)   # plasma glucose concentration (mg/dl)
EGP    = SFunction('EGP')(t) # endogenous glucose production (mg/kg/min)
Ra     = SFunction('R_{\\textit{a}}')(t)  # appearance rate of glucose in plasma (mg/kg/min)
E      = SFunction('E')(t)   # renal excretion (mg/kg/min)
Uii    = SFunction('U_{\\textit{ii}}')(t) # insulin-independent glucose utilizations (mg/kg/min)
Uid    = SFunction('U_{\\textit{id}}')(t) # insulin-dependent glucose utilizations (mg/kg/min)
U      = SFunction('U')(t)   # total glucose utilizations (mg/kg/min)

Ip     = SFunction('I_{\\textit{p}}')(t) # Insulin masses in plasma (pmol/kg)
Il     = SFunction('I_{\\textit{l}}')(t) # Insulin masses in liver (pmol/kg)
I      = SFunction('I')(t)   # Insulin concentration (pmol/l)
S      = SFunction('S')(t)   # Insulin secretion (pmol/kg/min)
m3     = SFunction('m_3')(t) # ??
HE     = SFunction('HE')(t)  # Hepatic Extraction of insulin
Ipo    = SFunction('I_{\\textit{po}}')(t) # Amount of insulin in the portal vein (pmol/kg)
Id     = SFunction('I_{\\textit{d}}')(t)  # delayed insulin signal realized with 2 compartments (pmol/l)
I1     = SFunction('I_1')(t) # used by Id (for delayed insulin signal)

Qsto  = SFunction('Q_{\\textit{sto}}')(t)   #  amount of glucose in the stomach (mg) = Qsto1 + Qsto2
Qsto1 = SFunction('Q_{\\textit{sto1}}')(t)  #  amount of solid glucose in the stomach (mg)
Qsto2 = SFunction('Q_{\\textit{sto2}}')(t)  #  amount of liquid glucose in the stomach (mg)
Qgut  = SFunction('Q_{\\textit{gut}}')(t)   #  glucose mass in the intestin (mg)
kempt = SFunction('K_{\\textit{empt}}')(t)  #  rate "constant" of gastric emptying (1/min)

X     = SFunction('X')(t)                   # insulin in the interstitial fluid (pmol/L)
Y     = SFunction('Y')(t)                   # not defined ...
Spo   = SFunction('S_{\\textit{po}}')(t)    # not defined ...

symbols.time = {t}
symbols.variables = [Gp, Gt, G, EGP, Ra, E, Uii, Uid, U, Ip, Il, I, S, m3, HE, Ipo, Id, I1, Qsto, Qsto1,
                     Qsto2, Qgut, kempt, X, Y, Spo]

# Parameters
k1    = SSymbol('k_1')      # rate parameter (1/min) for glucose model
k2    = SSymbol('k_2')      # rate parameter (1/min) for glucose model
VG    = SSymbol('V_G')      # distribution volume of glucose (dl/kg)
VI    = SSymbol('V_I')      # distribution volume of insulin (l/kg)
m1    = SSymbol('m_1')      # rate parameter (1/min) for insulin model
m2    = SSymbol('m_2')      # rate parameter (1/min) for insulin model
m4    = SSymbol('m_4')      # rate parameter (1/min) for insulin model
m5    = SSymbol('m_5')      # rate parameter (1/min) for insulin model
m6    = SSymbol('m_6')      # rate parameter (1/min) for insulin model
kmin  = SSymbol('k_{min}')  # min value of kempt (1/min)
kmax  = SSymbol('k_{max}')  # max value of kempt (1/min)
kab   = SSymbol('k_{ab}')   # rate constant of intestinale absorption (1/min)
kgri  = SSymbol('k_{gri}')  # rate of grinding (1/min)
f     = SSymbol('f')        # fraction of intestinal absorption which actually appears in plasma
b     = SSymbol('b')        # percentage of the dose for which kempt decreases at (kmax-kmin)/2
c     = SSymbol('c')        # percentage of the dose for which kempt is back to (kmax-kmin)/2
kp1   = SSymbol('k_{p1}')   # extrapolated at zero glucose and insulin (mg/kg/min)
kp2   = SSymbol('k_{p2}')   # liver glucose effectiveness (1/min)
kp3   = SSymbol('k_{p3}')   # parameter governing amplitude of insulin action on the liver (mg/kg/min per pmol/l)
kp4   = SSymbol('k_{p4}')   # parameter governing amplitude of portal insulin action on the liver (mg/kg/\min /(pmol/kg))
ki    = SSymbol('k_i')      # rate parameter accounting for delay between insulin signal and insulin action (1/min)
Fcns  = SSymbol('F_{cns}')  # rate of glucose uptake by the brain and erythrocytes (mg/kg/min)
p2U   = SSymbol('p_{\\textit{2U}}') # rate constant of insulin action on the peripheral glucose utilization (1/min)
Vm0   = SSymbol('V_{\\textit{m0}}') # (mg/kg/min)
Vmx   = SSymbol('V_{\\textit{mx}}') # (mg/kg/min per pmol/l)
Km0   = SSymbol('K_{\\textit{m0}}') # (mg/kg)
Kmx   = SSymbol('K_{\\textit{mx}}') # (mg/kg)
γ     = SSymbol('\\gamma')  # transfer rate constant between portal vein and liver (1/min)
K     = SSymbol('K')        # pancreatic responsivity to the glucose rate of change (pmol/kg per mg/dl)
α     = SSymbol('\\alpha')  # delay between glucose signal and insulin secretion (1/min)
β     = SSymbol('\\beta')   # pancreatic responsivity to glucose (pmol/kg/min per mg/dl)
h     = SSymbol('h')        # threshold of glucose above which the cells initiate to produce new insulin (mg/dl)
ke1   = SSymbol('k_{\\textit{e}_1}') # glomerular filtration rate (1/min)
ke2   = SSymbol('k_{\\textit{e}_2}') # renal threshold of glucose (mg/kg)

# memo: h has been set equal to Gb

symbols.parameters = [ k1, k2, VG, VI, m1, m2, m4, m5, m6, kmin, kmax, kab, kgri, f, b, c, kp1, kp2, kp3, kp4,
                      ki, Fcns, p2U, Vm0, Vmx, Km0, Kmx, γ, K, α, β, h, ke1, ke2 ]

# Variables at basal steady state
Gpb  = SSymbol('G_{\\text{pb}}')  # Gp at basal state
Gtb  = SSymbol('G_{\\text{tb}}')  # Gt at basal state
Gb   = SSymbol('G_{\\text{b}}')   # G at basal state
EGPb = SSymbol('EGP_{\\text{b}}') # EGP at basal state
Rab  = SSymbol('R_{\\textit{ab}}')
Eb   = SSymbol('E_{\\text{b}}')   # E at basal state
Uiib = SSymbol('U_{\\text{iib}}') # Uii at basal state
Uidb = SSymbol('U_{\\text{idb}}') # Uid at basal state
Ub   = SSymbol('U_{\\text{b}}')   # Ub at basal state (ie. Uiib+Uid)

Ipb    = SSymbol('I_{\\text{pb}}')   # Ip at basal state
Ilb    = SSymbol('I_{\\text{lb}}')   # Il at basal state
Ib     = SSymbol('I_{\\text{b}}')    # I at basal state
Sb     = SSymbol('S_{\\text{b}}')    # S at basal state
m3b    = SSymbol('m_{\\text{3b}}')   # m3 at basal state
HEb    = SSymbol('HE_{\\text{b}}')   # HE at basal state
Ipob   = SSymbol('I_{\\text{pob}}')  # Ipo at basal state
Idb    = SSymbol('I_{\\text{db}}')   # Id at basal state
I1b    = SSymbol('I_{\\text{1b}}')   # I1 at basal state

Qstob  = SSymbol('Q_{\\textit{stob}}')
Qsto1b = SSymbol('Q_{\\textit{sto1b}}')
Qsto2b = SSymbol('Q_{\\textit{sto2b}}')
Qgutb  = SSymbol('Q_{\\textit{gutb}}')
kemptb = SSymbol('k_{\\textit{emptb}}')

Xb     = SSymbol('X_{\\textit{b}}')
Yb     = SSymbol('Y_{\\textit{b}}')
Spob   = SSymbol('S_{\\textit{pob}}')

symbols.variables_basal = [ Gpb, Gtb, Gb, EGPb, Rab, Eb, Uiib, Uidb, Ub, Ipb, Ilb, Ib, Sb, m3b, HEb,
                            Ipob, Idb, I1b, Qstob, Qsto1b, Qsto2b, Qgutb, kemptb, Xb, Yb, Spob ]

# input variables
D     = SFunction('D')(t) # we set D a function with null derivative instead of a constant parameters to make easier the
                          # updating of its value
# D     = SSymbol('D')        # amount of ingested glucose (mg)
BW    = SSymbol('BW')       # body weight (kg)
symbols.input_variables = { Gb, D, BW }

#############################
# Definition of the equations
#############################

model.equations = dict()

model.equations.update({
    # Equations (1)
    Gp : sympy.Eq(Gp.diff(), EGP + Ra - Uii - E - k1*Gp + k2*Gt),
    Gt : sympy.Eq(Gt.diff(), -Uid + k1*Gp - k2*Gt),
    G  : sympy.Eq(G, Gp / VG),
    # Endogenous Glucose Production (Equation (10))
    EGP  : sympy.Eq(EGP, kp1 - kp2*Gp - kp3*Id - kp4*Ipo),
    # Glucose rate of appearance (Equation (13))
    Qsto  : sympy.Eq(Qsto, Qsto1 + Qsto2),
    Qsto1 : sympy.Eq(Qsto1.diff(), -kgri * Qsto1 + D*sympy.DiracDelta(t)),
    Qsto2 : sympy.Eq(Qsto2.diff(), -kempt*Qsto2 + kgri*Qsto1),
    Qgut  : sympy.Eq(Qgut.diff(), -kab*Qgut + kempt*Qsto2),
    Ra    : sympy.Eq(Ra, f*kab*Qgut/BW),
    # Kempt (Equation (8) from [Dalla Man et al. 2006])
    kempt  : sympy.Eq(kempt, kmin + (kmax-kmin)*(tanh(5*(Qsto-b*D)/(2*D*(1-b)))
                                  - tanh(5*(Qsto-c*D)/(2*D*c))+2)/2),
    # Glucose utilization (Equations (14) (15) (18) and (19))
    Uii : sympy.Eq(Uii, Fcns),
    Uid : sympy.Eq(Uid, (Vm0*Gt + Vmx*X*Gt)/(Km0+Kmx*X+Gt)),
    X   : sympy.Eq(X.diff(), -p2U*X + p2U*(I-Ib)),
    U   : sympy.Eq(U, Uii+Uid),
    # Insulin secretion (Equations (23) (24) (25) and (26)
    S   : sympy.Eq(S,γ*Ipo),
    Ipo : sympy.Eq(Ipo.diff(), -γ*Ipo + Spo),
    Y   : sympy.Eq(Y.diff(),
                   sympy.Piecewise((-α*(Y-β*(G-h)), β*(G-h)>= -Sb),
                                   (-α*Y-α*Sb, True))),
    Spo : sympy.Eq(Spo,
                   sympy.Piecewise((Y+K*G.diff()+Sb, (G.diff()>0)),
                                   (Y+Sb, True))),
    # Glucose renal excretion (Equation (27)
    E : sympy.Eq(E, sympy.Piecewise((ke1*(Gp-ke2), Gp>ke2),
                                    (0, True))),
    # Insulin subsytem (Equations (3), (11))
    Il : sympy.Eq(Il.diff(), -(m1 + m3)*Il + m2*Ip + S),
    Ip : sympy.Eq(Ip.diff(), -(m2 + m4)*Ip + m1*Il),
    I  : sympy.Eq(I,Ip/VI),
    I1 : sympy.Eq(I1.diff(), -ki*(I1-I)),
    Id : sympy.Eq(Id.diff(), -ki*(Id-I1)),
    # Hepatic extraction (Equations (4), (5))
    HE : sympy.Eq(HE, -m5*S + m6),
    m3 : sympy.Eq(m3, HE*m1/(1-HE))
})
