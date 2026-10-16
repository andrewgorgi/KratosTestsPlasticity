# SBM Dirichlet, grado 1: sensibilità al surrogate

## Blocco riprodotto

Con `lambda_inner=lambda_outer=1`, line search standard, Δt=0.001,
PENALTY_FACTOR=6700000 e le tolleranze del JSON:

- le mesh 10×10, 20×20 e 40×40 completano i due passi;
- la 80×80 fallisce al primo passo, t=0.001, dopo 100 iterazioni;
- con Δt=0.0001, la 80×80 fallisce invece a t=0.0006;
- una prova con limiti più restrittivi per la line search non risolve il blocco.

Quindi il solo trasferimento del solver verificato per Gap-SBM non rende
robusto anche questo setup SBM, e ridurre Δt non risolve tutte le mesh.

## Prova controllata sulla scelta del surrogate

Sulla mesh 80×80, h=0.055, misurando la distanza dei punti delle condizioni
surrogate dal cerchio fisico più vicino:

| lambda interno/esterno | Distanza massima / h | Distanza media / h |
|---|---:|---:|
| 1 / 1 | 1.29590 | 0.61170 |
| 0.5 / 0.5 | 0.62666 | 0.24077 |

Mantenendo invariati solver, Δt, materiale, penalizzazione, tolleranze e
quadratura, la seconda scelta completa il calcolo. Questo mostra una
sensibilità alla configurazione del surrogate; non è una dimostrazione
generale della causa matematica dell'instabilità o della correttezza
di tutta l'implementazione non lineare SBM.

## Configurazione salvata e verifiche

`sbm/ProjectParameters_dirichlet.json` e il generatore ora usano
`lambda_inner=lambda_outer=0.5`. Le impostazioni del solver rimangono
uguali a Gap-SBM Dirichlet: line search, Δt=0.001, massimo 100 iterazioni.
Non sono stati cambiati i cerchi fisici, il materiale o la background mesh.
**Cambia il dominio surrogate**, che può attraversare il contorno fisico;
il riferimento analitico è esteso a quel dominio come già previsto dal runner.
I precedenti CSV con lambda=1 non vanno uniti ai nuovi risultati.

Con grado 1, un thread e due passi fino a t=0.002 sono state completate:

| Mesh | Iterazioni primo/secondo passo | L2 relativo | H1 bulk relativo |
|---|---|---:|---:|
| 10×10 | 5 / 6 | 0.299394 | 0.443146 |
| 20×20 | 4 / 7 | 0.187930 | 0.300677 |
| 40×40 | 5 / 9 | 0.0318243 | 0.115224 |
| 80×80 | 5 / 10 | 0.0129349 | 0.0560150 |
| 160×160 | 7 / 11 | 0.00285061 | 0.0262733 |

Le prove sono state eseguite anche in invocazioni separate con gli stessi
parametri. Non garantiscono la convergenza per altri gradi o dati.
La variante SBM mista conserva lambda=1. La successiva correzione del segno
della trazione interna e le sue verifiche sono in `neumann_data_check.md`.

```bash
cd AnnulusPlasticity/sbm
python3 convergence.py --degree 1 --boundary dirichlet
# Per controllare i residui:
python3 convergence.py --degree 1 --boundary dirichlet --echo-level 1
```
