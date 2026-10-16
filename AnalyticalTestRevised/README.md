# Studio analitico rivisto

Tre metodi sul dominio fisico **[0,4] × [0,2]**:

- `body_fitted`: griglia coincidente con il dominio, senza skin separata; `SupportSolidCondition` su tutti i quattro lati.
- `sbm`: stesso rettangolo importato come skin, griglia di fondo estesa a [-0.5,4.2] × [-0.3,2.1], condizioni `SbmCondition` sul surrogate.
- `gap_sbm`: stessa skin e stessa griglia di fondo di SBM; `NurbsGeometryModelerGapSbm`, elementi del gap e condizioni all'interfaccia. Derivato da `../ProjectParametersSolidGap.json`, usando i nomi **GapSbm…** registrati dalla build locale (il file storico usa CutSbm…).

Ogni cartella contiene `run_and_post.py`, `convergence.py`, `plot.py`, `ProjectParameters.json` e otto configurazioni commentate `ProjectParameters_<campo>_<materiale>.json`. I runner condividono `study.py` e `analytic_fields.py`, per evitare formule diverse tra simulazione, grafici e norme di errore.

## Campi e materiale: due scelte separate

| `--field` | Campo di spostamento esatto | Significato |
|---|---|---|
| `linear` | (x t, 0) | Campo spaziale lineare; patch test, BODY_FORCE=0 |
| `polynomial` | (x(x−4)y(y−2)t, 0) | Campo polinomiale degli ultimi JSON storici |
| `cosine_middle` | (−2 cos(πx/4) sin(πy/2)t, 4 sin(πx/4) cos(πy/2)t) | Massimo stress trial nel mezzo, y=1 |
| `cosine_boundary` | Campo precedente con y→y+1 | Massimo stress trial sul bordo, y=0,2 |

`--physics elastic` usa elasticità lineare **plane strain**, E=67000, ν=0, densità=1, spessore=1. `--physics plastic` mantiene la legge J2 plane strain dei file storici, σ_y=201, HARDENING_CURVE=3 (plasticità perfetta).

Il default degli script e di `ProjectParameters.json` è **linear + elastic**, verificato per tutti i metodi. La parola `linear` indica il campo spaziale; non viene usata per nascondere una scelta di materiale.

Ogni file impone lo spostamento esatto su **tutto il contorno** e contiene la forza di volume −div(σ) dello stesso campo. Non riproduce il vecchio mix Dirichlet/Neumann con trazioni nulle incompatibili con il campo scelto. GapElements riceve la stessa BODY_FORCE del bulk. SBM legge lo spostamento non storico dei nodi della skin: il runner lo aggiorna analiticamente prima di ogni passo.

Le formule plastiche riprendono l'ipotesi analitica storica: proiezione radiale J2 del trial elastico in 3D, carico proporzionale monotono t≥0, nessun hardening. **La corrispondenza di questo riferimento con l'integrazione della legge plane strain ridotta deve essere verificata oltre lo snervamento**; le varianti plastiche non sono una validazione già conclusa.

## Esecuzione

Comandi dalla radice del repository; funzionano anche da altre directory usando il percorso completo dello script:

```bash
python3 AnalyticalTestRevised/body_fitted/run_and_post.py --field cosine_middle --physics elastic --show
python3 AnalyticalTestRevised/sbm/convergence.py --field cosine_boundary --physics elastic --degree 2 --refinements 10 20 40 80
python3 AnalyticalTestRevised/gap_sbm/plot.py --field cosine_middle --physics elastic --show
```

Gli stessi argomenti valgono in tutte e tre le cartelle. `--spans 20 10` modifica la singola mesh, `--degree 1` cambia il grado. In `convergence.py` si usano invece gli span quadrati definiti da `--refinements`. `--parameters NomeFile.json` permette di leggere un file modificato a mano: mantenere il `problem_name` coerente con `--field` e `--physics`.

Per **body_fitted, sbm e gap_sbm**, il campo `cosine_boundary` dispone anche del caso
`--physics plastic_hardening`, con J2 plane strain e incrudimento isotropo
lineare **H=10000**, E=67000, ν=0, σ_y iniziale=201. Usa
`ProjectParameters_cosine_boundary_plastic_hardening.json` e il materiale
dedicato `materials_plastic_hardening.json`. Forcing e tensioni di riferimento
sono coerenti con `Utilities/automatic_differentiation_hardening.py` per
carico proporzionale monotono t≥0; non c'è saturazione esponenziale.
Il nuovo caso impone Dirichlet esatto sui quattro lati: non usa le trazioni
manufactured hardcoded nella `LoadSolidCondition` locale.

```bash
python3 AnalyticalTestRevised/body_fitted/run_and_post.py --field cosine_boundary --physics plastic_hardening --degree 2 --threads 1
python3 AnalyticalTestRevised/body_fitted/convergence.py --field cosine_boundary --physics plastic_hardening --degree 2 --refinements 10 20 40 --threads 1
python3 AnalyticalTestRevised/sbm/convergence.py --field cosine_boundary --physics plastic_hardening --degree 2 --refinements 10 20 40 --threads 1
python3 AnalyticalTestRevised/gap_sbm/convergence.py --field cosine_boundary --physics plastic_hardening --degree 2 --refinements 10 20 40 --threads 1
```

I risultati, incluso il contour della deformazione plastica equivalente,
vanno in `<metodo>/results/cosine_boundary_plastic_hardening/`.
Per cambiare H, aggiornare `HARDENING` in `analytic_fields.py` e rigenerare
**solo** queste configurazioni e il loro materiale, senza sovrascrivere gli altri JSON:

```bash
cd AnalyticalTestRevised
python3 - <<'PY'
from build_parameters import write, write_material
write_material('plastic_hardening')
for method in ('body_fitted', 'sbm', 'gap_sbm'):
    write(method, 'cosine_boundary', 'plastic_hardening')
PY
```

Allineare anche `H_VALUE` nello script simbolico se usato separatamente.
SBM aggiorna lo spostamento esatto sulla skin a ogni passo; Gap-SBM riceve
lo stesso forcing con hardening sia nel bulk sia nei GapElements.
Il nuovo caso SBM usa Δt=0.0001: con Δt=0.0002, la prova sulla mesh 10×10
non convergeva a t=0.0014. Body-fitted e Gap-SBM usano Δt=0.0002.

`--threads N` permette di impostare esplicitamente il numero di thread Kratos. Se omesso, si usa il numero di thread definito dall'ambiente.

La quadratura degli elementi si imposta con `number_of_integration_points_per_span` nella voce `SolidElement` di `IgaModelerSbm`: un valore 10 richiede 10 punti per direzione e quindi 100 punti per cella in 2D. `number_of_integration_points` non e letto dal modeler. Nella build KratosPR e stato corretto anche `NurbsSurfaceGeometry::CreateIntegrationPoints`, che ignorava il valore richiesto e usava sempre degree+1; la patch e in `kratos_quadrature_fix.patch`. Le librerie sono state ricompilate/installate e verificato il conteggio effettivo: 100 punti su una cella, 200 su due celle, e il default degree+1 quando il parametro manca.

Correzione del crash parallelo nella build `KratosPR`: `SolidElement` e `SupportSolidCondition` forzavano `ANALYTICAL_TANGENT`, usando formule plastiche anche con la legge elastica. La lettura non const di `YIELD_STRESS`, assente nel materiale elastico, inseriva la variabile nelle proprieta condivise durante i cicli paralleli, causando una race e corruzione della memoria. Sono state rimosse le due definizioni e ricompilata/installata `libKratosIgaCore.so`: tensioni e tangente provengono ora dalla legge costitutiva selezionata. La modifica C++ e conservata in `kratos_parallel_constitutive_fix.patch`. Verificato il caso body-fitted polynomial elastic, p=1, mesh 10/20/40/80, con 24 thread: L2 circa 0.09481/0.02542/0.006535/0.001653, coerente con la prova seriale.

I risultati vengono salvati in `<metodo>/results/<campo>_<materiale>/`: grafici degli spostamenti e delle tensioni numeriche/analitiche/errori, `solution.npz`, `errors.json`, oppure `convergence.csv` e `convergence.png`. `plot.py` salva `geometry.png`. `--show` apre le figure; altrimenti vengono soltanto salvate.

Per i casi plastici, `run_and_post.py` produce anche `equivalent_plastic_strain.png`: la deformazione plastica equivalente numerica all'ultimo passo, letta dalla legge costitutiva nei punti di integrazione di bulk e gap. I valori sono salvati anche in `solution.npz` nella voce `equivalent_plastic_strain`.

Per rigenerare tutti i JSON e i materiali:

```bash
python3 AnalyticalTestRevised/build_parameters.py
```

**La rigenerazione sovrascrive i JSON di questa nuova cartella.** Per modifiche permanenti, cambiare `build_parameters.py` e/o `analytic_fields.py`; le espressioni di carico e gli stress analitici sono generati dalla stessa definizione simbolica.

## Significato delle norme

- `L2_relative`: norma vettoriale dello spostamento sul dominio computazionale, incluso il gap per Gap-SBM.
- `H1_bulk_relative`: norma H1 completa (spostamento + gradiente) sul **bulk**. La build locale non espone NEIGHBOUR_GEOMETRIES a Python; nel gap vengono letti lo spostamento e gli stress già valutati dall'elemento al punto vero, ma non viene inventato un gradiente. Per questo non si dichiara H1 sul dominio completo Gap-SBM.
- `area` e `H1_bulk_area`: aree effettivamente integrate; area esatta del rettangolo=8. SBM integra il surrogate, che cambia con la mesh. Non confondere le sue norme con integrali sul dominio fisico esatto.
- I pesi **firmati** del gap sono conservati: servono anche a sottrarre porzioni del bulk, e non vanno trasformati in valori assoluti.
- `h=max(hx,hy)`; le pendenze sono calcolate tra raffinamenti consecutivi. Il campo lineare è un patch test, con errori vicini alla precisione macchina: per studiare gli ordini usare i campi cosinusoidali.

## Verifiche e limite della build locale

I patch test elastici hanno dato L2 relativo circa 1.4e−16 (body-fitted), 1.5e−16 (SBM), 2.8e−16 (Gap-SBM). I tre script di plot sono stati eseguiti.

Per `cosine_middle`, p=2, due mesh di prova:

| Metodo | Span per direzione | L2 relativo iniziale → finale |
|---|---|---|
| Body-fitted | 4 → 8 | 1.106e−2 → 1.197e−3 |
| SBM | 8 → 16 | 3.511e−3 → 1.776e−4 |
| Gap-SBM | 8 → 16 | 1.001e−2 → 5.262e−4 |

Queste sono verifiche brevi, non uno studio completo degli ordini asintotici. I CSV delle prove sono presenti nelle cartelle dei risultati.

**Plasticità body-fitted e SBM:** la build locale va in segmentation fault durante `FinalizeSolutionStep` delle condizioni `SupportSolidCondition`/`SbmSolidCondition`. Il debugger identifica `GenericSmallStrainIsotropicPlasticity::FinalizeMaterialResponseCauchy`; nel sorgente la condizione passa parametri costitutivi senza i vettori/matrici della risposta. I JSON plastici sono disponibili, ma questi due runner richiedono una correzione/ricompilazione di Kratos prima di usarli. La libreria esterna non è stata modificata.

Gap-SBM plastico, `cosine_middle`, mesh 8×8: il runner termina (L2≈0.200, H1 sul bulk≈0.467). Questi errori su una mesh sola **non dimostrano convergenza**, né risolvono la verifica del riferimento plastico plane strain descritta sopra.
