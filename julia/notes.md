# Répertoires
- `data`: données OBEDIAB, script de curation et d'extraction des données, résultats sérialisés
- `init`: fichiers d'initialisations
- `models`: modèle original et modèle étendus/modifiés
- `myOptims`: encapsulation de diverses fonctions d'optimisation
- `main`: les différents programmes d'estimation de paramètres, analyse d'identifiabilité, etc.
- `results`: les divers résultats
- `misc`: fichiers divers
# Modifications du modèle
## h (19/03/2021)
- ce paramètre doit en fait être égal à Gb (comme écrit dans la sous-section 4. Insulin Secretion).
## Gb (19/03/2021)
- ce paramètre devrait être une propriété du modèle car, intuitivement, la glycémie à jeun est déterminée par les activités de production d'insuline par le pancréas et de consommation par les tissues.
- par Eq (2) et (14), Uidb est fonction de EGPb et Eb (attention, au steady state, Eb=0!)
- par Eq (18), Xb=0
- par Eq (15), Gtb est fonction de Uidb et donc de EGPb et Eb. Mais aussi de Vm0 qui, par (22) dépend de Gtb :( 
- Si on décide plutôt de laisser Vm0 comme paramètre fixe ou estimé, Gtb n'est plus fonction que de EGPb et Eb. 
- il s'ensuit, par Eq (1) que Gb peut être défini comme fonction de EGPb et Eb!
- il faut donc redéfinir Vm0 comme paramètre.
## EGP
- dans [DallaMan 2007], il est dit que "EGP is constrained to be non-negative". => ajout de max(...,0).
  
# Profile Likelihood
## Procedure
- fitted_var: observed variables used to fit the model
- pdev: % of standard deviation used to generate pseudo random data
- st: subject type (normal or diabetic) used for data
- α: seuil de loss 
## Results
1) fitted_var = ["Gp","Ip"], pdev = 0.05, st = normal, α = loss + cquantile(Chisq(1), 0.05)
2) fitted_var = ["Gp","Ip"], pdev = 0.05, st = diabetic, α = loss + cquantile(Chisq(1), 0.05)
3) fitted_var = ["Gp","Ip"], pdev = 0.05, st = normal, α = loss + cquantile(Chisq(1), 0.05)/2
4) fitted_var = ["Gp","Ip"], pdev = 0.05, st = diabetic, α = loss + cquantile(Chisq(1), 0.05)/2
5) fitted_var = ["Gp","Ip","Qgut"], pdev = 0.05, st = normal, α = loss + cquantile(Chisq(1), 0.05)
6) fitted_var = ["Gp","Ip","Qgut"], pdev = 0.05, st = diabetic, α = loss + cquantile(Chisq(1), 0.05)
7) fitted_var = all_vn, pdev = 0.05, st = normal, α = loss + cquantile(Chisq(1), 0.05)
8) fitted_var = all_vn, pdev = 0.05, st = diabetic, α = loss + cquantile(Chisq(1), 0.05)
