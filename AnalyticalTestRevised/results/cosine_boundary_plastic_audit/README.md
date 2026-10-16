# Confronto cosine_boundary plastic: body-fitted, SBM, Gap-SBM

**Aggiornamento successivo:** gli operatori Cut-SBM ora assemblano stress e tangente analitici tramite una macro, su richiesta dell'utente. Vedere [modifica e verifiche](README_cut_analytical.md). La sezione costitutiva sotto descrive la situazione precedente alla macro; il confronto diretto con la legge incrementale resta valido per il ramo alternativo.

Verifica 2026-10-08, sorgenti e librerie **KratosPR**.

## Riferimento e dati prescritti

I tre JSON attuali hanno lo stesso campo di spostamento con fase y+1:

```
ux = -2*t*cos(pi*x/4)*sin(pi*(y+1)/2)
uy =  4*t*sin(pi*x/4)*cos(pi*(y+1)/2)
```

Coincidono esattamente le stringhe delle formule DISPLACEMENT e BODY_FORCE;
GapElements riceve la stessa BODY_FORCE del bulk. I parametri dei solver sono
identici, compresi dt=0.00125, t_final=0.0025, tolleranze e massimo 100 iterazioni.
I materiali risolvono tutti allo stesso materials_plastic.json:
E=67000, nu=0, sigma_y=201, HARDENING_CURVE=3, densità e spessore 1.

La suddivisione del contorno fisico [0,4]×[0,2] coincide:

| Lato | Body-fitted | SBM | Gap-SBM |
|---|---|---|---|
| bottom/top, Neumann | LoadSolidCondition, brep_ids 2,4 | SbmLoadSolidCondition | CutSbmLoadSolidCondition |
| left/right, Dirichlet | SupportSolidCondition, brep_ids 3,5 | SbmSolidCondition | CutSbmSolidCondition |

Il DISPLACEMENT assegnato anche ai gruppi bottom/top nel JSON Gap-SBM
non impone un vincolo Dirichlet: CutSbmLoadSolidCondition assembla soltanto
la trazione prescritta e non legge quella variabile.
Il commento iniziale «Dirichlet esatto su TUTTO il contorno» è obsoleto
in tutti e tre i file attuali.

Le tre condizioni Load nel C++ contengono attualmente il carico analitico
hardcoded cosine_boundary, con fattore cos(pi*y/2), equivalente allo shift y+1
del campo. Questa verifica riguarda questo campo e questi parametri di materiale;
non certifica che i carichi C++ funzionino per altri campi/materiali.

## Verifiche eseguite

`check_cases.py` carica i JSON con il parser Kratos, confronta le impostazioni
e valuta le formule in 5 punti a tre tempi, comprendendo i rami elastico e
plastico. Le espressioni coincidono con analytic_fields.py nei punti verificati.

Inizializza poi i tre modelli con p=2 e mesh breve 8×4, preservando tutte le
altre impostazioni originali. Nessun SolveSolutionStep viene eseguito.
Il RHS nativo di ogni condizione Neumann, allo spostamento numerico nullo,
viene sommato sui coefficienti e diviso per il peso di integrazione;
la partizione dell'unità restituisce la trazione prescritta. Le normali e i punti
sono controllati sul bottom/top. Sono esclusi i punti esattamente ai vertici,
dove la normale vera non è unica.

| Metodo | Condizioni Neumann controllate | Massima differenza assoluta da sigma_analitica·n |
|---|---:|---:|
| body_fitted | 160 | 8.53e-14 |
| sbm | 60 | 1.42e-13 |
| gap_sbm | 72 | 1.42e-13 |

L'inizializzazione dei modelli avviene in directory temporanee; i JSON,
i risultati ordinari e i sorgenti di produzione non vengono modificati.

## Differenze numeriche nei JSON correnti

| Impostazione | Body-fitted | SBM | Gap-SBM |
|---|---|---|---|
| Box della griglia di fondo | [0,4]×[0,2] | [-0.2,4.11]×[-0.3,2.1] | [-0.5,4.2]×[-0.3,2.1] |
| Grado nel JSON | 2,2 | 2,2 | 2,2 |
| Knot span nel JSON | 20,10 | 20,10 | 20,10 |
| Quadratura bulk per direzione | 10 | 20 | 20 |
| Punti iniziali per importazione NURBS | — | 1000 | 10000 |

L'uguaglianza degli span non implica l'uguaglianza di hx/hy su box diversi.
Gli script convergence possono sovrascrivere grado e span con gli argomenti
dell'esecuzione. La quadratura di bulk non imposta automaticamente la quadratura
degli elementi/condizioni del gap, creati da SnakeGapSbmProcess.
Gap-SBM usa gap_approximation_order=1, number_internal_divisions=1
e gap_sbm_type=interpolation: impostazioni proprie del metodo.

## Differenza costitutiva concreta nel KratosPR attuale

SolidElement.cpp definisce ANALYTICAL_TANGENT e sostituisce sia lo stress
sia la tangente assemblati con AnalyticalStress/AnalyticalConstitutiveMatrix.
Questo è il bulk comune ai tre metodi. SupportSolidCondition e le condizioni
SBM usano anch'esse le formule analitiche nel codice corrente.

CutSbmSolidElement e CutSbmSolidCondition usano invece direttamente
CalculateMaterialResponse della legge incrementale
SmallStrainIsotropicPlasticityPlaneStrainVonMisesVonMises.
Il materiale nominale è quindi condiviso, ma la risposta effettivamente
assemblata nel gap e sul suo bordo Dirichlet non segue lo stesso percorso
costitutivo del bulk analitico.

`check_material.py` inizializza quella stessa legge con i parametri del materiale,
impone la deformazione del campo analitico al punto (2,0), e calcola/finalizza
i due passi originali senza risolvere un problema spaziale.
Lo stress analitico è valutato indipendentemente con analytic_fields.py.

| Tempo | sigma_xx, sigma_yy analitiche | sigma_xx, sigma_yy incrementali | Differenza relativa sul vettore (xx,yy,xy) |
|---|---|---|---:|
| 0.00125 | -43.8306, -263.1396 | -62.7085, -224.8788 | 15.9933% |
| 0.0025 | -175.3848, -394.6938 | -170.2304, -221.7432 | 40.0613% |

Questo è un confronto locale alla medesima deformazione, non un errore L2
di una simulazione completa. Dimostra che non si può attribuire ogni differenza
fra Gap-SBM e gli altri metodi alla sola approssimazione geometrica.
Non è stata scelta né modificata una delle due formulazioni costitutive.

Anche il dominio delle norme differisce: L2 include il gap in Gap-SBM,
mentre SBM integra sul surrogato; H1_bulk_relative esclude il gap.

## Riproduzione

```bash
export PYTHONPATH=/home/agorgi/Documents/GECKO/KratosPR/bin/Release
export LD_LIBRARY_PATH=/home/agorgi/Documents/GECKO/KratosPR/bin/Release/libs
export MPLBACKEND=Agg
python3 AnalyticalTestRevised/results/cosine_boundary_plastic_audit/check_cases.py
python3 AnalyticalTestRevised/results/cosine_boundary_plastic_audit/check_material.py
```

`case_audit.json` conserva parametri e controlli dei dati/carichi.
`material_point_check.json` conserva deformazioni, stress e differenze locali.
