# Controllo del flusso Neumann e della SbmLoadSolidCondition attuale

Controllo del 18 ottobre 2026: p3, mesh 20×20 e 40×40, caso ciclico misto,
dt=0.00025, penalty-free, lambda_inner=0, lambda_outer=0.5, materiale J2 con
hardening lineare H=10000. Input congelati uguali a quelli delle prove precedenti.

**Non è emerso un errore nel segno, nella normale, nel punto di lettura del
carico o nell'aggiornamento costitutivo che spieghi la mancata convergenza.**
È stato invece verificato un difetto della proiezione geometrica vicino alla
chiusura del cerchio. Correggerlo in una prova isolata non elimina l'arresto.
Questo controllo non dimostra che l'intera formulazione SBM sia stabile o
spazialmente consistente per questa configurazione.

## Versione verificata

È stata usata la libreria installata, con hash registrato in
[build.json](diagnostics/sbm_neumann_current/build.json), senza ricompilarla
o cambiare i sorgenti. Lo
[snapshot della condition](diagnostics/sbm_neumann_current/source_snapshot.cpp)
ha macro `SBM_USE_PROJECTION_POINT_CONSTITUTIVE_LAW` attivo e
`ANALYTICAL_TANGENT` disattivo.

Rispetto alle precedenti librerie diagnostiche, ora la Neumann usa l'ordine
Taylor completo **2p in 2D**, invece di p. La verifica numerica dell'operatore
elastico conferma che la libreria installata esegue effettivamente l'estensione
completa. Non ci si è basati soltanto sul sorgente presente sul disco.

## 1. Dato imposto e punto di lettura

Per ogni condition interna è stato recuperato il suo nodo di proiezione nel
sottomodelpart `skin_model_part.inner.Layer0.projection_data`. Il modeler crea
conditions e nodi in ordine corrispondente; questa associazione è stata poi
verificata mediante ricostruzione dell'intero operatore locale.

La condition legge `FORCE` dal nodo proiettato, non dalla condition. Il processo
del test assegna il dato anche a questi nodi, perché appartengono alla gerarchia
della skin. La formula indipendente usata per il confronto è

```
C(t) = C_FINAL * livello_ciclico(t)
n_true = -(x,y)/r
g_N = -C(t)/r² * e_theta = C(t)*(y,-x)/r³
```

Sono stati controllati sette istanti: 0.0005, 0.001, 0.0015, 0.003, 0.004,
0.005 e 0.006, includendo carico negativo e scarico a trazione nulla.
Il confronto usa la formula polare, indipendente dalla costruzione cartesiana
di `cyclic_analytic.boundary_traction`.

Per tutte le 280 conditions interne della mesh 20 e le 560 della mesh 40:

- `FORCE` coincide con la trazione indipendente entro arrotondamento.
- La normale del punto proiettato coincide con la normale radiale interna
  entro circa 7e-16.
- Il prodotto tra normale surrogate e true è sempre positivo.
- A spostamento nullo, il RHS coincide con `N_i * g_N * (n_s·n_t) * peso`.
- Le prove con coefficienti affini e casuali elastici riproducono anche i
  termini di stress al surrogate e al true.

Questi controlli escludono nel benchmark verificato il precedente problema di
segno della normale interna e un eventuale mancato aggiornamento della FORCE
sui nodi di proiezione. Il campo FORCE assegnato direttamente alle conditions
non viene scambiato per il dato effettivamente letto dalla SBM.

## 2. Gradienti, segni e tangente locale

Il gradiente al true è stato ricostruito **indipendentemente dal Taylor C++**:
si usano B-spline univariate SciPy, convertite in polinomi per ciascuno span,
e il loro prodotto tensoriale. Si valuta il prolungamento polinomiale dello
stesso span del surrogate, anche quando il punto proiettato oltrepassa uno span.
Non viene sostituito con la B-spline di un diverso span.

Prima si verificano valori e gradienti al surrogate; poi si confrontano
gradiente diretto e Taylor al true. Il background del test ha Jacobiano identità,
verificato numericamente su tutte le conditions.

Con D elastica indipendente, il riferimento è:

```
R_i = N_i * [sigma_s*n_s + (g_N - sigma_t*n_t)*(n_s·n_t)] * peso
K_ij = N_i * [-T(n_s)*D*B_s + (n_s·n_t)*T(n_t)*D*B_t]_j * peso
```

| Verifica, massimo errore relativo | Mesh 20 | Mesh 40 |
|---|---:|---:|
| Gradiente true, Taylor completo vs polinomio diretto | 1.05e-14 | 2.04e-14 |
| Matrice locale C++ vs matrice indipendente | 9.11e-15 | 3.25e-14 |
| Gradiente true, vecchio troncamento a p vs diretto | 1.189 | 1.114 |

Il vecchio troncamento a p non ricostruiva completamente il gradiente di una
base tensoriale Qp. La modifica già presente a 2p corregge questo aspetto,
ma non è sufficiente a far convergere il caso attuale.

Il controllo verifica la coerenza dell'assembly con le formule riportate:
non è una dimostrazione che queste formule garantiscano stabilità globale.

## 3. Storia al surrogate e al punto proiettato

Per 12 conditions distribuite sul bordo di ciascuna mesh sono state applicate
otto configurazioni di coefficienti: plasticizzazione, scarico, cambio di
direzione e inversione. Due punti materiali indipendenti, inizialmente vergini,
ricevono rispettivamente `B_s*u` e `B_t*u`, con B_t ottenuta dai polinomi diretti.

Il ciclo di inizializzazione/finalizzazione viene eseguito dallo stesso scheme
Kratos della simulazione. Dopo valutazioni di prova anche plastiche si ripristina
lo stato di spostamento, verificando che il residuo si ripeta identico.
Le storie indipendenti vengono consolidate soltanto a fine passo.

| Massimo scarto, 96 stati per mesh | Mesh 20 | Mesh 40 |
|---|---:|---:|
| Residuo vs due storie indipendenti | 3.59e-14 | 2.64e-14 |
| Tangente vs due storie indipendenti | 6.26e-15 | 1.72e-14 |
| Residuo dopo valutazioni di prova e ripristino | 0 | 0 |
| K·v vs differenza centrale del residuo | 3.03e-9 | 4.96e-9 |

La norma della differenza fra strain true e surrogate arriva a circa 0.00489
sulla mesh 20 e 0.00230 sulla 40: non si tratta quindi di una prova in cui
le due storie coincidono per costruzione. Per questa legge e questo benchmark,
il richiamo alla legge indipendente al true, il suo commit e il suo uso nel
residuo/Jacobiano risultano coerenti.

## 4. Difetto individuato nella proiezione geometrica

L'`IgaModelerSbm` segnala mancata convergenza della proiezione vicino alla
chiusura del cerchio e accetta un punto di ripiego. Rispetto alla proiezione
radiale esatta:

| Mesh | Massimo errore sulla posizione proiettata | Conditions con errore >1e-4 |
|---|---:|---:|
| 20×20 | 0.00508916 | 1 |
| 40×40 | 0.00254460 | 1 |

Sulla 20, per esempio, il punto surrogate `(0.00559813, 1.1)` viene associato
al punto `(0,1)` anziché alla sua proiezione radiale. Il nodo è sul cerchio
e ha normale corretta per la propria posizione, ma non è il closest point.

Il modeler usa una tolleranza `max(1e-6, 1e-3*h)` sia nella ricerca sia nella
selezione tra candidati. La logica di selezione non privilegia automaticamente
un candidato convergente rispetto a un ripiego quasi equivalente in distanza.
È un aspetto da correggere/verificare nella ricerca della proiezione; non
si tratta di un errore di segno nel flusso prescritto.

Per controllarne l'effetto è stata eseguita un'ulteriore prova sulla mesh 20:
prima dell'inizializzazione delle conditions si correggono soltanto coordinate
e normale dei 280 nodi di proiezione interni con la proiezione radiale esatta.
La stessa condition e lo stesso solver restano in uso; nessun sorgente viene
modificato. **Anche questa prova fallisce a t=0.00075.**
Il difetto geometrico non spiega quindi, da solo, l'arresto osservato.

## 5. Solve con la libreria attuale

| Configurazione | Arresto Newton | Ultimo tempo convergente | Ultimo errore L2 |
|---|---:|---:|---:|
| Mesh 20, proiezione del modeler | 0.00075 | 0.0005 | 0.0913041 |
| Mesh 40, proiezione del modeler | 0.00075 | 0.0005 | 0.0258129 |
| Mesh 20, proiezione interna radiale esatta | 0.00075 | 0.0005 | 0.0913041 |

La causa della mancata convergenza globale resta aperta. I controlli restringono
l'indagine alla stabilità/consistenza spaziale dell'operatore misto e alla sua
interazione con la soluzione plastica: non giustificano attribuire l'arresto
a un carico Neumann di segno errato o a un mancato commit della legge al true.
Non dimostrano nemmeno che il problema sia intrinseco al metodo e non a un
ulteriore aspetto d'implementazione non coperto da queste prove.

## Riproduzione e dati

- [Audit completo mesh 20](diagnostics/sbm_neumann_current/audit_p3_n20/summary.json)
- [Audit completo mesh 40](diagnostics/sbm_neumann_current/audit_p3_n40/summary.json)
- [Tutte le conditions, mesh 20](diagnostics/sbm_neumann_current/audit_p3_n20/conditions.json)
- [Tutte le conditions, mesh 40](diagnostics/sbm_neumann_current/audit_p3_n40/conditions.json)
- [Storie indipendenti mesh 20](diagnostics/sbm_neumann_current/audit_p3_n20/history.json)
- [Storie indipendenti mesh 40](diagnostics/sbm_neumann_current/audit_p3_n40/history.json)
- [Solve attuale mesh 20](diagnostics/sbm_neumann_current/new/spans_20x20/result.json)
- [Solve attuale mesh 40](diagnostics/sbm_neumann_current/new/spans_40x40/result.json)
- [Controllo con proiezione esatta](diagnostics/sbm_neumann_current/exact_inner_projection/new/spans_20x20/result.json)

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg python3 AnnulusPlasticity/diagnostics/sbm_neumann_current/audit.py --degree 3 --spans 20
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg python3 AnnulusPlasticity/diagnostics/sbm_neumann_current/audit.py --degree 3 --spans 40
LD_PRELOAD=/home/agorgi/Documents/GECKO/Kratos/bin/Release/libs/libKratosIgaCore.so OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg python3 AnnulusPlasticity/diagnostics/sbm_neumann_current/projection_control.py
```

I test usano cartelle isolate e non riscrivono i JSON di produzione.
Il controllo di proiezione confronta l'hash della libreria con quello registrato
per evitare di mescolare versioni. Lo script di audit usa la libreria installata:
per confrontare nuove compilazioni conservare separatamente questi risultati.
