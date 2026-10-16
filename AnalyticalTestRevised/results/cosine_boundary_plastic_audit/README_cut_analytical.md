# Macro analitica negli operatori Cut-SBM

KratosPR, modifica del 2026-10-08.

## Scelta della risposta assemblata

La macro è attiva all'inizio di ciascuno dei tre file:

```cpp
#define ANALYTICAL_TANGENT
```

- custom_elements/cut_sbm_solid_element.cpp
- custom_conditions/cut_sbm_solid_condition.cpp
- custom_conditions/cut_sbm_solid_interface_condition.cpp

Con la macro attiva, stress e tangenti assemblati usano le medesime formule
di SolidElement::AnalyticalStress/AnalyticalConstitutiveMatrix. I metodi sono
condivisi nel nuovo custom_utilities/analytical_j2_response_utility.h,
così i tre operatori Cut-SBM usano una sola implementazione.

Per tornare all'assemblaggio con la legge incrementale, commentare la macro
nei tre file e ricompilare/installare KratosIgaCore. I rami `#else` conservano
la risposta costitutiva precedente; i tre file compilano anche con macro spenta.
Questa è una scelta in compilazione, non un parametro dei JSON.

## Punti aggiornati

- CutSbmSolidElement: stress nel RHS, tangente nel LHS, tensioni pubblicate
  in CAUCHY_STRESS_XX/YY/XY durante FinalizeSolutionStep.
- CutSbmSolidCondition: stress nel RHS, tangente nel LHS, operatore aggiunto
  di Nitsche mDBOld calcolato analiticamente e congelato all'inizio del passo;
  tensioni e NORMAL_STRESS pubblicati dopo FinalizeSolutionStep.
- CutSbmSolidInterfaceCondition: stress e tangenti su entrambi i lati
  dell'interfaccia negli operatori LHS/RHS.

CutSbmLoadSolidCondition assembla già il carico analitico prescritto per
cosine_boundary e non richiede questa sostituzione di stress/tangente numerici.
I JSON e study.py non sono stati modificati da questo intervento.

Le leggi incrementali vengono ancora inizializzate, calcolate e finalizzate
per aggiornare il loro stato, come in SolidElement. La macro seleziona la
risposta assemblata e gli stress pubblicati, non elimina il precedente
return mapping. EQUIVALENT_PLASTIC_STRAIN rimane un dato della legge incrementale;
non viene convertito automaticamente nella deformazione plastica del riferimento
analitico a deformazione totale.

Le formule analitiche descrivono la risposta J2 perfettamente plastica
a deformazione totale usata dal benchmark; non introducono hardening.
La lettura dei parametri nella utility è const. YIELD_STRESS assente o nullo
seleziona il ramo elastico, senza inserire la variabile nelle proprietà.

## Verifiche

Build di KratosIgaCore riuscita, libreria installata in KratosPR/bin/Release/libs.
Il precedente binario è conservato in /tmp/libKratosIgaCore.before_cut_analytical.so.
Build ID compilato/installato: 8932c92e553d420f0902da9e73ce07dd60b2cad7.

`check_cut_analytical.py` verifica gli operatori nativi, con campo affine
che produce una deformazione costante e con una proiezione radiale 3D
indipendente come riferimento. Controlla quattro stati per materiale plastico
e quattro per materiale elastico, inclusi taglio e deformazioni oltre snervamento.
Sul gap il lavoro virtuale di campi affini ricostruisce stress e matrice
costitutiva direttamente dal RHS/LHS, senza leggere metodi privati della classe.

- Massima differenza relativa di stress: circa 2.12e-15.
- Massima differenza di tangente costitutiva da differenze finite: 1.73e-10.
- Massima differenza del Jacobiano nativo dal residuo differenziato:
  3.73e-10 sui tre operatori.
- Gli stress pubblicati da elemento e condizione Dirichlet coincidono con
  lo stesso riferimento plastico dopo FinalizeSolutionStep.
- Il caso elastico verifica sia YIELD_STRESS assente, sia il default nullo
  che può essere inserito dagli accessi di altri operatori esistenti.

La verifica è locale: lo schema inizializza gli operatori senza assemblare
il residuo globale. Non costituisce una validazione degli altri elementi
in un problema elastico completo.

Simulazioni brevi del caso reale cosine_boundary plastic, mesh 8×4,
dt=0.00125 e tempo finale 0.0025, materiali e tolleranze originali.
La quadratura del JSON è stata modificata nell'IDE durante la verifica:
sono conservati i parametri effettivi di ciascuna prova.

| Grado | Thread | Quadratura bulk per direzione | Iterazioni nei due passi | L2 relativa |
|---|---:|---:|---|---:|
| p=1 | 4 | 20 | 6,6 | 0.3177844 |
| p=2 | 4 | 30 | 6,5 | 0.0404009 |
| p=2 | 24 | 30 | 6,5 | 0.0404009 |

Le tre prove raggiungono il tempo finale. Il runner di verifica interrompe
l'esecuzione se SolveSolutionStep restituisce false. Nel confronto p=2 i
parametri sono identici fra 4 e 24 thread e la L2 differisce di circa 2.57e-16.
Una sola mesh non dimostra
un ordine di convergenza spaziale.

## Riproduzione

```bash
export PYTHONPATH=/home/agorgi/Documents/GECKO/KratosPR/bin/Release
export LD_LIBRARY_PATH=/home/agorgi/Documents/GECKO/KratosPR/bin/Release/libs
export MPLBACKEND=Agg
python3 AnalyticalTestRevised/results/cosine_boundary_plastic_audit/check_cut_analytical.py
python3 AnalyticalTestRevised/results/cosine_boundary_plastic_audit/check_gap_run.py 2 4
```

Risultati: cut_analytical_checks.json e gap_analytical_*_results.json.
Le prove usano directory temporanee e non sovrascrivono i risultati ordinari.
Patch dei quattro sorgenti: ../../kratos_cut_sbm_analytical_response.patch.
