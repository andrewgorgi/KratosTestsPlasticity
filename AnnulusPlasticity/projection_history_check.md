# Verifica della storia costitutiva al punto proiettato SBM

Verifica del 16 ottobre 2026 su `SbmLoadSolidCondition`, per il caso 2D con
`SmallStrainJ2PlasticityPlaneStrain2DThreeComponentsLaw` e macro
`SBM_USE_PROJECTION_POINT_CONSTITUTIVE_LAW` attiva. `ANALYTICAL_TANGENT` è disattiva.

## Esito

La versione corrente del sorgente e la libreria installata aggiornata superano
il controllo di selezione, valutazione e aggiornamento della legge proiettata.
Non emerge una contaminazione fra storia surrogata e storia proiettata, né un
commit delle iterazioni di Newton/line search rifiutate. Questo non dimostra
la stabilità del problema globale SBM.

Durante la verifica è cambiata la libreria installata: inizialmente aveva
timestamp 14:54:43, successivamente 15:15:36. Il sorgente aveva timestamp
15:02:11. La prima libreria inizializzava `mBasisFunctionsOrder` a **1** anche
per geometrie p2: il disassemblato di `InitializeMemberVariables` contiene
`movq $0x1,0xe0(%rbx)`. Il sorgente corrente ricava invece l'ordine dal numero
di funzioni di base. Sostituendo soltanto `InitializeMemberVariables` con la
versione compilata dal sorgente corrente, il controllo locale p2 passa.
Quindi la prima discrepanza osservata riguarda l'ordine dell'estensione di
Taylor, non un mancato aggiornamento della legge proiettata.

Ho compilato e collegato copie diagnostiche esclusivamente in `/tmp`, senza
modificare o installare sorgenti/librerie Kratos. L'aggiornamento della libreria
installata è avvenuto esternamente durante la verifica. Il controllo finale
passa anche con Python normale, senza `LD_PRELOAD`.

## Percorso nel sorgente

File Kratos: `applications/IgaApplication/custom_conditions/sbm_load_solid_condition.cpp`.

| Fase | Operazione verificata |
|---|---|
| `InitializeMaterial`, righe 49–53 | Due `Clone()` distinti, entrambi inizializzati; nessuna storia condivisa. |
| LHS, righe 294–314 | Calcola `B_sum`, `strain_true = B_sum*u`, chiama `ApplyConstitutiveLaw(..., true)` e usa `D_true*B_sum`. |
| RHS, righe 418–424 | Ricostruisce la stessa deformazione estesa e valuta lo stress con `ApplyConstitutiveLaw(..., true)`. |
| `ApplyConstitutiveLaw`, righe 604–630 | Seleziona `mpConstitutiveLawAtProjectionPoint` per il flag `true`; richiede stress e tangente con deformazione fornita dalla condition. |
| `InitializeSolutionStep`, righe 713–740 | Inizializza le risposte delle due leggi. Per questa legge J2 `InitializeMaterialResponseCauchy` è vuota e non azzera la storia. |
| `UpdateMaterialResponseAtProjectionPoint`, righe 744–767 | Ricostruisce ancora `B_extended*u` con la stessa funzione di Taylor; finalizza soltanto quando `Finalize` è vero. |
| `FinalizeSolutionStep`, righe 694–696 | Finalizza prima la legge surrogata, poi quella proiettata, entrambe con la deformazione finale corrente. |

La legge J2 usa le variabili plastiche memorizzate per costruire la risposta
trial. `CalculateMaterialResponseCauchy` non le committa;
`FinalizeMaterialResponseCauchy` le aggiorna. Il secondo calcolo della risposta
durante la finalizzazione riparte ancora dalla storia precedente: non applica
due incrementi plastici. Le chiamate ripetute LHS/RHS sono quindi ammissibili.

Passare la geometria surrogata a `InitializeMaterial` non introduce un errore
per questa legge: l'override a tre componenti ignora geometria e funzioni di
forma e inizializza soltanto le variabili plastiche. Questa conclusione non
si estende automaticamente a leggi che interpolano campi o usano lunghezze
geometriche. L'estensione di Taylor corrente usa derivate parametriche come
fisiche; il background dell'annulus ha Jacobiano identità, verificato nel test.

La legge al punto proiettato è propria della **Neumann** `SbmLoadSolidCondition`.
La Dirichlet `SbmSolidCondition` mantiene la propria legge al surrogato e
l'operatore aggiunto `mDBOld` congelato all'inizio dello step: non possiede una
seconda storia al true boundary.

## Test riproducibile

```bash
python3 AnnulusPlasticity/check_projection_history.py
```

Il test inizializza le geometrie reali del caso mixed e confronta la condition
con **due istanze indipendenti** della stessa legge materiale, pilotate
rispettivamente da `B*u` e `B_sum*u`. Ricostruisce autonomamente proiezione,
normali, deformazioni, residuo e tangente. Non verifica nuovamente la formula
interna del return mapping: verifica il collegamento e la storia della condition.

Combinazioni: p1/p2, mesh 10/20, H=0/H=10000; otto condition distribuite sul
foro; sei step con carico non proporzionale, scarico, inversione e ricarico.
Totale: **384 confronti di stato**. Fra le valutazioni vengono imposte prove
plastiche perturbate e poi ripristinati i coefficienti, simulando prove
rifiutate. Il commit reale usa lo stesso scheme del solver, senza un solve
Newton, per isolare la gestione della storia.

Massimi su tutte le combinazioni, libreria aggiornata:

| Controllo | Errore massimo |
|---|---:|
| Residuo contro le due storie indipendenti | 2.64e-14 |
| Tangente contro `D_s*B` e `D_true*B_sum` | 2.60e-14 |
| Risposta dopo prove trial rifiutate | 0 |
| Jacobiano del residuo mediante differenze finite, anche dopo commit | 1.25e-8 |

Per p2 le due deformazioni sono effettivamente diverse. Il controllo negativo,
che usa al true la storia committata al surrogato, produce differenze di stress
fino a circa **85.2**: il test distingue davvero le due storie. Per p1
l'operatore di gradiente non ha termini superiori e `B_sum=B`.

Risultati: `diagnostics/projection_history.json`. La prima verifica sulla
libreria precedente è conservata in `diagnostics/projection_history_before_rebuild.json`.

## Risultati globali e limiti

Con la condition corrente compilata in una copia completa della libreria IGA,
sono stati ripetuti i controlli p2, mesh 20, dt=0.001, line search, H e dati
Neumann specifici di ogni riferimento, mantenendo il materiale copiato con
`PENALTY_FACTOR=-1` e lambda inner/outer=1:

| Caso | Esito | L2 relativo |
|---|---|---:|
| Patch affine con hardening | Newton converge | 13.6076 |
| Coseno non affine, perfect plasticity | Newton non converge al primo step | — |
| Annulus con hardening | Newton converge | 67.2305 |

Risultati in `diagnostics/projection_relinked`. Un'ulteriore copia strumentata
conferma `order=2` per tutte le Neumann dell'annulus p2 e restituisce lo stesso
errore globale (`diagnostics/projection_instrumented`). L'ordine di Taylor
corretto e il commit locale corretto non rimuovono quindi il problema globale.
Resta da discriminare stabilità della formulazione, trattamento Dirichlet e
possibili errori di assembly globale. Il fallimento con H=0 non basta a
dimostrare un limite intrinseco di SBM.

Questa verifica riguarda simulazioni nuove 2D con mesh fissa. La serializzazione
`save/load` della condition salva soltanto la classe base, non le due leggi:
non consente di affermare che la storia proiettata venga preservata nei restart.
Non sono stati verificati restart, remeshing e caso 3D.
