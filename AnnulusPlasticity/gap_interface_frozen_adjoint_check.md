# Gap-SBM: operatori adjoint congelati sui due lati dell'interfaccia

Correzione applicata e installata in Kratos il 16 ottobre 2026.

In `GapSbmSolidInterfaceCondition` sono stati aggiunti `mDBOldPlus` e
`mDBOldMinus`. A ogni `InitializeSolutionStep` si ricostruiscono i due B estesi
al punto vero e si calcola, con la rispettiva storia costitutiva committata,
`D_alg * B`. I due prodotti rimangono fissi durante le iterazioni dello step.

I termini adjoint usano questi operatori **sia nella LHS sia nella RHS**.
I termini di flusso mantengono stress e tangenti aggiornati alle deformazioni
trial correnti; ciascun lato continua a committare la propria storia soltanto
in `FinalizeSolutionStep`. Il commit non ricalcola gli operatori congelati.

La modifica riguarda la condition d'interfaccia standard usata dal caso
`gap_sbm_type=interpolation`. Non introduce altre varianti e non cambia i
parametri di produzione, il materiale o i criteri di convergenza.

## Build e patch

Modificati soltanto `gap_sbm_solid_interface_condition.h` e `.cpp` nel repository
Kratos. La [patch](diagnostics/gap_interface_frozen_adjoint/change.patch) è
conservata con copie prima/dopo e checksum. La build è stata eseguita con:

```bash
KRATOS_BUILD_TYPE=Release bash build/configure.sh
```

Lo script ufficiale ha ricompilato IgaApplication e installato le librerie.
I test successivi usano la **libreria installata**, senza `LD_PRELOAD`.
Il [log di build](diagnostics/gap_interface_frozen_adjoint/build.log) e
l'[identità della libreria](diagnostics/gap_interface_frozen_adjoint/installed_library.json)
sono salvati con i risultati.

## Verifiche

Si usano gli stessi input congelati del controllo pre-fix: p2, mixed, H=10000,
penalty=-1, quadratura bulk di default e dt=0.00025. Le prove girano in cartelle
isolate per non sovrascrivere i risultati o i file geometrici dell'utente.

| Mesh e strategia | Prima | Dopo |
|---|---|---|
| 10×10, Newton standard | fallisce a t=0.00125 | completa 24 passi, t=0.006 |
| 20×20, line search | fallisce a t=0.00475 nel controllo precedente | completa 24 passi, t=0.006 |

Per la mesh 20, il controllo precedente con line search usava 4 punti/span;
questa prova usa la quadratura default del fixture corrente. La riga 10×10
è il confronto controllato con identico input pre/post.

Prima di finalizzare i passi a t=0.0005, 0.00125, 0.003 e 0.006, si verificano
il Jacobiano globale su tre direzioni, con passi FD 1e-7, 1e-8 e 1e-9, e i
contributi locali per tipo di elemento/condition. Il residuo viene rivalutato
senza committare le perturbazioni della soluzione.

- Errore globale massimo sulle verifiche della mesh 10: **3.56e-9**.
- Errore globale massimo sulla mesh 20: **3.33e-9**.
- Ripetizione delle RHS dopo le perturbazioni: errore **zero**.
- I prodotti delle matrici locali ricostruiscono quello della matrice globale
  entro arrotondamento; i valori puntuali sono nei file `checks.json`.

Il precedente errore globale stabile del 4–9% e quello dell'interfaccia del 48%
non compaiono più in questi controlli. Il primo passo elastico della mesh 10
mantiene esattamente gli stessi errori L2/H1 e le stesse due iterazioni di prima.

Queste verifiche dimostrano la consistenza della linearizzazione negli stati
esaminati e il completamento dei due cicli. Non dimostrano un ordine asintotico:
la mesh 10 resta grossolana e presenta errori spaziali elevati.

## Riproduzione

Dalla radice del workspace, con la build Kratos aggiornata:

```bash
python3 AnnulusPlasticity/check_gap_interface_frozen_adjoint.py --spans 10 --strategy newton_raphson
python3 AnnulusPlasticity/check_gap_interface_frozen_adjoint.py --spans 20 --strategy line_search --output /tmp/gap_frozen_20
```

Lo script conserva i parametri effettivi, gli errori temporali e le verifiche
del Jacobiano. Gli assert controllano la consistenza numerica; lo stato del
solve è riportato in `result.json` e `results/status.json`.

- [Riepilogo](diagnostics/gap_interface_frozen_adjoint/summary.json)
- [Controlli mesh 10](diagnostics/gap_interface_frozen_adjoint/n10_newton/checks.json)
- [Controlli mesh 20](diagnostics/gap_interface_frozen_adjoint/n20_line_search/checks.json)
- [Diagnosi precedente](gap_cyclic_first_unloading_check.md)
