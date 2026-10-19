# SBM p3: sostituzione della Neumann interna con Dirichlet

Controllo del 18 ottobre 2026 sullo stesso caso ciclico del
[confronto misto p3](sbm_cyclic_old_vs_new_p3_check.md).

**Con penalty-free e dt=0.00025, Dirichlet su entrambi i bordi completa il ciclo
su entrambe le mesh 20×20 e 40×40.** Nel corrispondente caso misto, entrambe
le mesh fallivano a t=0.00075.

## Confronto controllato

È stato usato l'input congelato del caso misto, cambiando soltanto la condition
importata sul cerchio interno da `SbmLoadSolidCondition` a `SbmSolidCondition`.
Sono quindi prescritti gli spostamenti analitici su entrambi i bordi.
Il forcing resta nullo, il materiale resta J2 con hardening lineare H=10000,
e la soluzione di riferimento è la stessa storia ciclica analitica.

Restano p3, quadratura di default, line search, massimo 100 iterazioni,
tolleranze e geometria. In particolare **lambda_inner resta 0**, come nel caso
misto, e lambda_outer=0.5. Non è stato usato direttamente il JSON Dirichlet
di produzione, che ha lambda_inner=0.5 e una quadratura esplicita: ciò avrebbe
modificato più fattori contemporaneamente.

Le conditions create non includono alcuna `SbmLoadSolidCondition`.
I punti e pesi usati per misurare gli errori coincidono con quelli del caso misto.
Le prove usano la stessa libreria diagnostica temporanea dei precedenti controlli,
con hash verificato. Durante le esecuzioni sono cambiati esternamente il sorgente
e la libreria installata di Kratos: il sorgente corrente estende ora l'ordine
della `SbmLoadSolidCondition` a 2p in 2D e 3p in 3D. Questi cambiamenti non sono
inclusi nella libreria congelata usata qui. Nel caso interamente Dirichlet
la `SbmLoadSolidCondition` non viene creata. Non sono state apportate modifiche
alle conditions o alle configurazioni di produzione da questo controllo.

## Risultato principale: penalty-free

Errori normalizzati con la norma analitica al primo picco:

| Mesh | Passi completati | Massimo numero di iterazioni | L2 finale | H1 bulk finale |
|---|---:|---:|---:|---:|
| 20×20 | 24/24 | 8 | 0.404342 | 0.350262 |
| 40×40 | 24/24 | 9 | 0.00822166 | 0.0342641 |

La mesh 20 ha ancora un errore finale elevato; sulla 40 diminuisce sensibilmente.
Due mesh non bastano a stabilire l'ordine asintotico di convergenza.
Il test attraversa effettivamente la plasticità: la massima deformazione plastica
equivalente accumulata finale è circa 0.03143 sulla 20 e 0.04043 sulla 40.

## Sensibilità residua

Con dt dimezzato a 0.000125 e penalty-free:

- La mesh 20 completa 48/48 passi, con errore L2 finale 0.384057.
- La mesh 40 fallisce a t=0.005, dopo 39 passi. L'ultimo passo accettato,
  t=0.004875, ha errore L2 0.00719749.

Quindi Dirichlet elimina il blocco precoce della baseline mista nei due
raffinamenti provati, ma non garantisce convergenza per qualsiasi discretizzazione
temporale con questi criteri e questa strategia.

È stata provata anche penalty positiva 100E, con dt originale: sulla mesh 20
Newton fallisce a t=0.00125. Già al primo picco, t=0.001, l'errore L2 è 6.35166.
Sulla mesh 40 Newton fallisce a t=0.0025, dopo nove passi accettati; l'ultimo
errore L2 a t=0.00225 è 166.189. Anche in questo caso l'avanzamento del solver
non corrisponde a una soluzione accurata.
Questa penalizzazione non costituisce un miglioramento rispetto alla baseline
Dirichlet penalty-free.

## Interpretazione e file

Il confronto mostra che p3 con plasticità e storia ciclica **può completare il
ciclo** nella configurazione corrente. La difficoltà precoce della baseline
è associata al passaggio al bordo interno Neumann. Questo restringe il campo
d'indagine alla formulazione di quel bordo e alla sua interazione con il sistema;
non prova da solo uno specifico errore nella condition, né che una legge
costitutiva al punto proiettato basti a garantire stabilità.

Il macro `SBM_USE_PROJECTION_POINT_CONSTITUTIVE_LAW` riguarda
`SbmLoadSolidCondition`: nel caso interamente Dirichlet questa condition
non viene creata, perciò la distinzione vecchia/nuova della Neumann non interviene.

- [Grafico di confronto](diagnostics/sbm_p3_dirichlet/comparison.png)
- [Tabella completa](diagnostics/sbm_p3_dirichlet/summary.csv)
- [Errori di tutti i passi accettati](diagnostics/sbm_p3_dirichlet/errors_by_time.csv)
- [Verifica dei soli cambiamenti previsti](diagnostics/sbm_p3_dirichlet/checks.json)

Le sottocartelle conservano input congelati, log, esiti, conditions create e
ultimo campo accettato. Per eseguire le prove:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 AnnulusPlasticity/diagnostics/sbm_p3_dirichlet/run.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 AnnulusPlasticity/diagnostics/sbm_p3_dirichlet/run.py --variants free_dt125
OPENBLAS_NUM_THREADS=1 python3 AnnulusPlasticity/diagnostics/sbm_p3_dirichlet/postprocess.py
```

Il driver riutilizza i risultati già presenti nelle cartelle diagnostiche;
per un nuovo confronto indipendente conservare i risultati precedenti in una
cartella separata. Le librerie temporanee precedenti devono essere disponibili
e corrispondere agli hash verificati dal driver. Il controllo iniziale del driver
segnala anche se il sorgente o la libreria installata differiscono dalla baseline,
come accaduto alla fine di queste prove. Il postprocessore consente comunque
di leggere i risultati congelati e registra tale differenza in `checks.json`.
