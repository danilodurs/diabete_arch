## intestinal_leastsq_all_params

Name       Value      Min      Max   Stderr     Vary     Expr Brute_Step
p_Gb         242       50      300     None     True     None     None
p_b         0.68     0.01        1     None     True     None     None
p_c         0.09     0.01        1     None     True     None     None
p_f          0.9     0.01        1     None     True     None     None
p_kab      0.023    0.001        1     None     True     None     None
p_kgri    0.0465    0.001        1     None     True     None     None
p_kmax    0.0465    0.001        1     None     True     None     None
p_kmin    0.0076    0.001        1     None     True     None     None

## intestinal_leastsq_all_params2

Name       Value      Min      Max   Stderr     Vary     Expr Brute_Step
p_Gb         242       50      300     None     True     None     None
p_b         0.68     0.01        1     None     True     None     None
p_c         0.09     0.01        1     None     True     None     None
p_f          0.9     0.01        1     None     True     None     None
p_kab      0.023    0.001        1     None     True     None     None
p_kgri    0.0465    0.001        1     None     True     None     None
p_kmax    0.0465    0.001        1     None     True     None     None
p_kmin    0.0076    0.001        1     None     True     None     None

Attention, ici j'essaie de fitter les données d'un patient sein mais avec un poids de diabétique
(BW = 91) et VI et VG (normalement dépendant du poids) de diabétique aussi pour voir si on fit
mieux...
Résultat: pour glucose c'est un peu moins bien, pour l'insuline c'est un peu mieux...

## intestinal_lbfgsb_all_params

Name       Value      Min      Max   Stderr     Vary     Expr Brute_Step
p_Gb         242       50      300     None     True     None     None
p_b         0.68     0.01        1     None     True     None     None
p_c         0.09     0.01        1     None     True     None     None
p_f          0.9     0.01        1     None     True     None     None
p_kab      0.023    0.001        1     None     True     None     None
p_kgri    0.0465    0.001        1     None     True     None     None
p_kmax    0.0465    0.001        1     None     True     None     None
p_kmin    0.0076    0.001        1     None     True     None     None

Le valeurs trouvées sont n'importe quoi (quasiment inchangé):
{'p_kmax': 0.04650000000000003,
 'p_kmin': 0.007600000000000007,
 'p_kab': 0.023000000000000003,
 'p_kgri': 0.04650000000000003,
 'p_f': 0.9,
 'p_b': 0.68,
 'p_c': 0.09000000000000001,
 'p_Gb': 242.0}

Je ne garde pas les résultats...

## intestinal_cobyla_all_params
interruption de l'estimation, trop longue...

## intestinal_dogleg_all_params
"dogleg is an invalid method for this session" ...

## intestinal_shgo_all_params
fait plein de warning division par zéro... J'ai killé l'exécution...

## glucose_kinetics_leastsq_params
à refaire car j'ai mis des bornes pour VG qui l'ont contraint à 1 alors que je voulais le fixer à la
valeur des diabétiques (1.88)
Name     Value      Min      Max   Stderr     Vary     Expr Brute_Step
p_Gb       242       50      300     None     True     None     None
p_VG      1.49     -inf      inf     None    False     None     None
p_k1     0.042    0.001        1     None     True     None     None
p_k2     0.071    0.001        1     None     True     None     None
résultats décevants (Gpb n'est pas du tout fitté, il semble que l'algo ait prévilégié le fitting
de Ipb), à refaire en commençant avec Gb = 91.

## glucose_kinetics_leastsq_params

## liver_leastsq_param
On fit patient sein avec poids, VI et VG de diabetic

Name      Value      Min      Max   Stderr     Vary     Expr Brute_Step
p_Gb        242       50      300     None     True     None     None
p_ki     0.0066   0.0001     0.05     None     True     None     None
p_kp1      3.09        2        4     None     True     None     None
p_kp2    0.0007   0.0001    0.005     None     True     None     None
p_kp3     0.005   0.0001    0.005     None     True     None     None
p_kp4    0.0786    0.001      0.5     None     True     None     None

## insulin_kinetics_leastsq_param
On fit patient sein avec poids, VI et VG de diabetic
Name      Value      Min      Max   Stderr     Vary     Expr Brute_Step
p_Gb        242       50      300     None     True     None     None
p_HEb       0.6     -inf      inf     None    False     None     None
p_VI       0.04     -inf      inf     None    False     None     None
p_m1      0.379    0.001        1     None     True     None     None
p_m2      0.673    0.001        1     None     True     None     None
p_m4      0.269    0.001        1     None     True     None     None
p_m5     0.0526    0.001        1     None     True     None     None
p_m6     0.8118    0.001        1     None     True     None     None

Résultat: complètement à coté de la plaque! il faudrait ressayer cette estimation de paramètres...
!!!!! => il y a avait une erreur dans pancreatic_estimation.py: l'estimation était lancée avec le
modèle model_pacreatic au lieu de model_pancreatic!

(en cours à nouveau sur brassin...)

## pancreatic_leastsq_param
(en cours sur bush...)
Name     Value      Min      Max   Stderr     Vary     Expr Brute_Step
p_Gb       242       50      300     None     True     None     None
p_K       0.99      0.1       10     None     True     None     None
p_α      0.013    0.001      0.1     None     True     None     None
p_β       0.05   0.0001        1     None     True     None     None
p_γ        0.5     -inf      inf     None    False     None     None

## tissue_leastsq_param
Name       Value      Min      Max   Stderr     Vary     Expr Brute_Step
p_Fcns         1     -inf      inf     None    False     None     None
p_Gb         242       50      300     None     True     None     None
p_Km0      466.2      100     1000     None     True     None     None
p_Vm0       4.65        1       10     None     True     None     None
p_Vmx      0.034     0.01      0.1     None     True     None     None
p_p2U      0.084     0.01      0.5     None     True     None     None

## renal_leastsq_param
p_Gb        242       50      300     None     True     None     None
p_ke1    0.0007    1e-05    0.001     None     True     None     None
p_ke2       269      100      500     None     True     None     None
Résultats: à coté de la plaque!
