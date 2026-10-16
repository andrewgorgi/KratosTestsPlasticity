# Torsione elastoplastica di un anello: SBM e Gap-SBM

Per il problema SBM misto con grado 2 e i controlli in plasticità perfetta,
vedere [p2_mixed_investigation.md](p2_mixed_investigation.md).
Il completamento di Newton non garantisce l'accuratezza: nelle prove p=2
alcuni solve terminano con errori elevati. Non è ancora stabilito se la
causa sia di stabilità della formulazione o di implementazione globale.

Caso separato da `AnalyticalTestRevised`. Dominio fisico: **1 ≤ r ≤ 2**,
centro (0,0), piccole deformazioni, plane strain, J2 con incrudimento isotropo
lineare. Le geometrie vengono lette direttamente da:

- `../AnalyticalTestRevised/geometry/inner_circle.json`, raggio 1;
- `../AnalyticalTestRevised/geometry/outer_circle.json`, raggio 2.

I control point esterni alla circonferenza sono normali per curve NURBS razionali:
il test geometrico verifica i raggi delle curve, non quelli dei control point.

## Dati e riferimento analitico

`analytic.py` definisce E=67000, ν=0, G=33500, σ_y=201, H=10000.
Il carico cresce monotonamente su 0 ≤ t ≤ 0.002:

```
C(t) = C_final * t / 0.002
C_final = 201 * 1.5**2 / sqrt(3) = 261.1066592410083
tau_rtheta = C(t) / r**2
sigma_rr = sigma_thetatheta = sigma_zz = 0
BODY_FORCE = [0, 0, 0]
```

L'equilibrio è `d(tau)/dr + 2*tau/r = 0`, anche nella zona plastica.
Il fronte è `rp(t)=sqrt(sqrt(3)*C(t)/201)`: al carico finale è a r=1.5.
La deformazione plastica equivalente è
`alpha=max((sqrt(3)*C(t)/r**2 - 201)/10000, 0)`.

Con `u_r=0`, `u_theta=r*f(r,t)` e bordo esterno fermo:

```
f_elastic = -C/(2*G) * (1/r**2 - 1/2**2)
f_plastic = f_elastic - 3*C/(2*H)*(1/r**2 - 1/rp**2)
                         + sqrt(3)*201/H*log(rp/r)
```

Si usa il secondo ramo solo per r < rp. In coordinate cartesiane:

```
u = [-f*y, f*x]
sigma_xx = -2*C*x*y/r**4
sigma_yy =  2*C*x*y/r**4
sigma_xy = C*(x*x-y*y)/r**4
```

Al bordo interno, a fine carico: u_theta=-0.010565772032267568,
tau=261.1066592410083, alpha=0.025125. Al bordo esterno: u=0,
tau=65.27666481025207, alpha=0. Il forcing volumetrico è sempre nullo;
la soluzione non è polinomiale. Il riferimento presuppone carico monotono,
materiale inizialmente vergine e deformazioni infinitesime.

## Configurazioni al bordo

**`ProjectParameters_dirichlet.json` (default, entrambi i metodi):**
spostamento analitico sul bordo interno, spostamento nullo su quello esterno.
SBM usa `SbmSolidCondition`; Gap-SBM usa `GapSbmSolidCondition`.
Il runner assegna i valori esatti a ogni passo, compresi quelli ai nodi
di proiezione della skin e alle condizioni Gap-SBM. BODY_FORCE viene impostata
esplicitamente a zero sia nel bulk sia negli elementi del gap.

**`ProjectParameters_mixed.json`:** bordo esterno fermo, trazione interna
`t = sigma*n`, con normale uscente dal materiale `n=-e_r`. Quindi
`t_inner=-(C/r**2)*e_theta`. SBM usa `SbmLoadSolidCondition` sul bordo interno,
che legge FORCE dai nodi proiettati. Questa variante può richiedere passi di
carico più piccoli e la sua convergenza non è garantita sulle mesh grossolane.
Il dato nodale usa la normale fisica analitica, non il NORMAL orientato
della curva: sulla skin interna linearizzata quest'ultimo può puntare
verso +e_r e dare il segno sbagliato alla trazione. La correzione e le
verifiche del caso misto sono in [neumann_data_check.md](neumann_data_check.md).
Audit senza solve: `python3 AnnulusPlasticity/audit_neumann.py`.

**Gap-SBM misto richiede una modifica C++:** la `GapSbmLoadSolidCondition`
locale contiene ancora le tensioni hardcoded del caso cosinusoidale.
`kratos_gap_load_force.patch` prepara la sostituzione con la lettura di
`GetValue(FORCE)`. Il runner assegna la trazione analitica alla condizione.
La patch è fornita come artefatto: **non è applicata né compilata da questi script**.
Solo dopo averla applicata alla versione compatibile di Kratos e ricompilato
IgaApplication, abilitare la variante con `--boundary mixed --gap-load-from-force`.
Il controllo impedisce di eseguire per errore il caso con la trazione cosinusoidale.

Il materiale separato `materials.json` usa **PENALTY_FACTOR=100*E=6700000**.
Questa scelta stabilizza le condizioni di Dirichlet: con PENALTY_FACTOR=-1
le prove iniziali si arrestavano oltre lo snervamento. La scelta della
penalizzazione è parte della configurazione numerica, non della soluzione esatta.
I JSON usano una griglia [-2.2,2.2]², grado 2, mesh 20×20,
4 punti di quadratura per direzione/cella nel bulk. Gap-SBM Dirichlet usa
line search e Δt=0.001. Anche SBM Dirichlet e misto usano le stesse
impostazioni solver e Δt del Gap-SBM Dirichlet; Gap-SBM misto mantiene
le proprie impostazioni nel JSON.
SBM Dirichlet usa ora `lambda_inner=lambda_outer=0.5` per tenere il
surrogate più vicino ai cerchi. La scelta precedente 1/1 falliva sulla
80×80 anche con line search. Cambia il dominio surrogate, non i cerchi
fisici; la variante mista conserva lambda=1. Prove e dettagli:
[sbm_degree1_check.md](sbm_degree1_check.md).

## Esecuzione

Dalla radice del workspace (i runner funzionano anche con percorso assoluto):

```bash
python3 AnnulusPlasticity/sbm/run_and_post.py --degree 2 --spans 20 20
python3 AnnulusPlasticity/gap_sbm/run_and_post.py --degree 2 --spans 20 20

python3 AnnulusPlasticity/sbm/plot.py --degree 2 --spans 20 20 --show
python3 AnnulusPlasticity/gap_sbm/plot.py --degree 2 --spans 20 20 --show

python3 AnnulusPlasticity/sbm/convergence.py --degree 2 --refinements 10 20 40 80
python3 AnnulusPlasticity/gap_sbm/convergence.py --degree 2 --refinements 10 20 40 80

# Variante SBM con carico al bordo interno
python3 AnnulusPlasticity/sbm/run_and_post.py --boundary mixed --degree 2 --spans 20 20 --dt 0.00002
```

Opzioni comuni: `--threads N` (default 1), `--dt`, `--show` e
`--parameters NomeFile.json` (relativo alla directory del metodo).
`--strategy newton_raphson|line_search` permette di confrontare le strategie;
`--echo-level 1` stampa i residui delle iterazioni non lineari.
`--boundary` sceglie `dirichlet` o `mixed`. Il default è esplicito e uguale
per i due metodi: Dirichlet su entrambi i bordi.
`convergence.py` apre sempre il grafico finale, anche senza `--show`, e
mostra le linee di riferimento h^p e h^(p+1), con p uguale al grado scelto.

Gli output sono in `<metodo>/results/<boundary>/p<grado>/`:

- `solution.npz`, `errors.json`, `parameters_used.json`, `conditions.json`;
- contour numerico/analitico/errore di spostamenti, tensioni e deformazione plastica equivalente;
- per la convergenza: `convergence.csv`, `convergence.png` e campioni per ogni mesh.

In caso di errore del solve, `failed_spans_<nx>x<ny>/` conserva
`parameters_used.json` e `failure.json`, con mesh, grado, strategia, Δt,
tempo e passo del fallimento (quando disponibili).

Il foro è mascherato nei contour. I campi sono valutati nelle coordinate
iniziali, senza amplificare o deformare la geometria dei grafici.
`plot.py` inizializza i modeler senza eseguire passi di carico e salva
`geometry.png`, `geometry.npz` e `geometry_parameters.json` nella stessa directory.
Mostra i cerchi esatti, la background mesh, i surrogate interno/esterno e
i centri degli elementi attivi, distinguendo bulk e gap. I centri provengono
da `element.GetGeometry().Center()` nel modello Kratos effettivamente creato.
In caso di Newton non convergente il runner si arresta: non produce norme
che presentino il solve fallito come una soluzione valida.

`L2_relative` confronta lo spostamento nel bulk più il gap; `H1_bulk_relative`
usa solo i punti del bulk, dove Kratos espone il gradiente.
Per SBM il riferimento è esteso analiticamente al dominio surrogate.
Per Gap-SBM si conservano i pesi firmati del gap. Le aree sono salvate:
l'area fisica esatta è 3*pi ≈ 9.42477796. Gli errori sono relativi alla
soluzione continua nel tempo: anche l'errore di integrazione costitutiva
può influire sullo studio spaziale a Δt fissato.

## Verifiche e rigenerazione

```bash
python3 AnnulusPlasticity/test_analytic.py
python3 AnnulusPlasticity/build_parameters.py
```

I test verificano geometrie NURBS, valori al bordo, gradiente, equilibrio
senza forcing e relazione J2 con hardening. Il generatore sovrascrive soltanto
i quattro ProjectParameters e il materiale di questa cartella.
Per modificare i parametri fisici, aggiornare `analytic.py` e rigenerare.
Le geometrie esterne devono rimanere coerenti con RI e RO.

Il runner riusa campionamento e controllo Newton di `AnalyticalTestRevised/study.py`;
il riferimento analitico, i carichi e il postprocessing anulare sono separati.

## Stato delle prove numeriche

Nelle prove iniziali, con lambda=1 per SBM, eseguiti i runner Dirichlet
su mesh 10×10 e i due studi su 10×10 → 20×20,
grado 2, Δt=0.0001, un thread, con la build locale di Kratos. I calcoli
Dirichlet terminano e scrivono i risultati. Errori relativi misurati:

| Metodo | Mesh | L2 | H1 bulk |
|---|---|---:|---:|
| SBM | 10×10 | 0.326935 | 0.645443 |
| SBM | 20×20 | 6.259112 | 7.825722 |
| Gap-SBM | 10×10 | 1.143039 | 2.279714 |
| Gap-SBM | 20×20 | 0.181400 | 0.306665 |

**Il caso non è ancora una validazione della convergenza spaziale di SBM:**
l'errore cresce nella prova 10×10 → 20×20 e richiede una diagnosi della
formulazione/implementazione. Gap-SBM migliora su queste due mesh, ma due
punti grossolani non stabiliscono un ordine asintotico.

Nelle prove iniziali, prima della correzione del segno della trazione,
la variante SBM mista si arrestava per Newton non convergente: a t=0.0013
su 10×10 con Δt=0.0001 e a t=0.0009 su 20×20 con Δt=0.00002.
Questi risultati misti storici non validano il problema con carico corretto.
Dopo la correzione del dato sono state completate le mesh 10/20/40/80/160
con grado 1, line search e Δt=0.001; vedere `neumann_data_check.md`.
Gap-SBM misto non è stato eseguito: richiede la patch C++ descritta sopra.
Questi limiti riguardano le prove numeriche; i test indipendenti del
riferimento analitico e delle due geometrie passano.

### Prove precedenti con Newton standard

Il check con Δt=0.001 completa le mesh 10×10 e 20×20, ma fallisce sulla
40×40 al secondo passo, t=0.002: dopo 100 iterazioni il rapporto del residuo
oscilla circa fra 0.00060 e 0.00064, sopra la tolleranza 1e-9.
Con lo stesso Newton e le stesse tolleranze, Δt=0.0001 completa invece
i 20 passi sulla 40×40; l'ultimo passo richiede 7 iterazioni.
Verificata anche la mesh 80×80 con grado 1, Δt=0.0001 e Newton standard:
completa il caricamento fino a t=0.002.
La riduzione del passo non è però una soluzione robusta per tutte le mesh:
un controllo successivo ha riprodotto il ciclo di Newton sulla 20×20 al
passo 19, t=0.0019, con Δt=0.0001. Questa prova si riproduce con:

```bash
python3 AnnulusPlasticity/gap_sbm/convergence.py --degree 1 --boundary dirichlet --refinements 20 --dt 0.0001 --strategy newton_raphson --echo-level 1
```

Si tratta della convergenza del solve non lineare; non dimostra da sola
la convergenza spaziale verso il riferimento analitico.

Il controllo dei valori assegnati alle condizioni, del contorno numerico e
dei tassi L2 per grado 1 è documentato in [dirichlet_degree1_check.md](dirichlet_degree1_check.md).
Il runner include inoltre un arresto con tolleranza floating point per evitare
un passo spurio oltre il tempo finale; il test è `python3 AnnulusPlasticity/test_time_loop.py`
dalla radice del workspace.

### Configurazione attuale: line search e passo maggiore

`gap_sbm/ProjectParameters_dirichlet.json` e il generatore ora selezionano
**line_search, Δt=0.001**, senza cambiare materiale, penalizzazione,
tolleranze o massimo di 100 iterazioni. La line search ridimensiona
l'aggiornamento di Newton e ha eliminato il ciclo nelle prove eseguite.
Con grado 1, carico finale t=0.002 e un thread, sono state completate
le seguenti mesh con due passi di carico:

| Mesh | Iterazioni primo passo | Iterazioni secondo passo |
|---|---:|---:|
| 10×10 | 4 | 5 |
| 20×20 | 4 | 6 |
| 40×40 | 5 | 8 |
| 80×80 | 5 | 9 |
| 160×160 | 7 | 11 |

```bash
python3 AnnulusPlasticity/gap_sbm/convergence.py --degree 1 --boundary dirichlet
# Per stampare anche i residui:
python3 AnnulusPlasticity/gap_sbm/convergence.py --degree 1 --boundary dirichlet --echo-level 1
# Confronto con un passo più piccolo, mantenendo line search:
python3 AnnulusPlasticity/gap_sbm/convergence.py --degree 1 --boundary dirichlet --dt 0.0001
```

Le mesh sono state verificate anche in invocazioni separate; non è una
garanzia di convergenza per qualsiasi parametro o grado. L'errore temporale
resta distinto dalla convergenza del solve: sulla 80×80, L2 è circa 0.084583
con Δt=0.001, contro 0.083515 nella prova precedente a Δt=0.0001.
