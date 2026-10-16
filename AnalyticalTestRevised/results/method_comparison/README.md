# Confronto cosine_middle plastic: body-fitted / SBM

Verifica del 2026-10-08 con KratosPR/bin/Release, inclusa la correzione della quadratura NURBS. I JSON originali e i sorgenti Kratos non sono stati modificati da questa verifica.

## Problema fisico

I due ProjectParameters_cosine_middle_plastic.json hanno identici solver_settings, processes, output_processes e problem_data salvo problem_name. Identica anche la voce bulk SolidElement: p=2, 30 punti per direzione. Entrambi importano ../materials_plastic.json: E=67000, nu=0, sigma_y=201, penalty=-1. Tempo finale 0.0025, passo 0.00125. Il bordo fisico importato da SBM è il rettangolo [0,4]×[0,2], come nel body-fitted.

AnalyticalStress e AnalyticalConstitutiveMatrix nei tre header SolidElement, SupportSolidCondition e SbmSolidCondition sono identici dopo rimozione di spazi e commenti. Le condizioni di bordo restano operatori diversi: SupportSolidCondition sul bordo vero, SbmSolidCondition con estensione Taylor al nodo skin proiettato. Il runner aggiorna quel nodo con lo stesso campo analitico al tempo corrente.

## Griglia e norme diverse

Inizializzazione con gli span originali 20×10 e p=2; quadratura ridotta a 3 per direzione soltanto per misurare la geometria. DomainSize delle geometrie di quadratura è sommato senza risolvere.

| Quantità | Body-fitted | SBM |
|---|---:|---:|
| Box della griglia | [0,4]×[0,2] | [-0.5,4.2]×[-0.3,2.1] |
| hx | 0.2 | 0.235 |
| hy | 0.2 | 0.24 |
| Knot span attivi | 200 | 112 |
| Area integrata nelle norme | 8 | 6.3168 |

Il rettangolo surrogato di questo SBM è [0.205,3.965]×[0.18,1.86]. Con lambda_outer=1 sono esclusi span tagliati dal bordo fisico. study.errors integra sui punti bulk attivi: il dominio delle norme non coincide fra i metodi, e in SBM cambia con la mesh. I risultati di convergence usano inoltre span n×n, sovrascrivendo 20×10 del JSON.

Sul bordo fisico di cosine_middle il trial von Mises esatto è zero. Sul bordo surrogato della griglia 20×10, valutando lo stesso campo a t=0.0025, 88 dei 230 punti delle condizioni hanno trial von Mises >201, con massimo 336.285. Quindi anche la condizione SBM può campionare una risposta plastica, mentre il bordo body-fitted è elastico. La classificazione qui è del campo esatto, non della soluzione numerica.

## Test breve di isolamento

Tre simulazioni nuove: p=2, span 10×10, quadratura bulk 10×10 per span, 4 thread, stessi materiali, tempi e tolleranze. Per sbm_aligned si cambiano soltanto i quattro estremi xyz/uvw del box, portandoli a [0,4]×[0,2]; il codice delle condizioni resta SbmSolidCondition.

| Variante | L2 relativa | H1 bulk relativa | Area | Iterazioni Newton nei due passi |
|---|---:|---:|---:|---|
| body_fitted | 0.0010514895872 | 0.0043244106355 | 8 | 6, 6 |
| sbm | 0.0037600241517 | 0.0085744272203 | 5.5272 | 7, 7 |
| sbm_aligned | 0.0010549604558 | 0.0043275993132 | 8 | 6, 6 |

Tutti i passi hanno SolveSolutionStep=True; il test interrompe la simulazione in caso contrario. Con il bordo allineato, L2 differisce solo dello 0.33% circa. Questo isola griglia e trattamento del bordo spostato come principali fonti della differenza in questo test; non prova l'assenza di difetti per un bordo spostato e plastico, né stabilisce un ordine asintotico da una sola mesh. La quadratura del fronte elastoplastico resta sensibile alla posizione della griglia.

I CSV di convergenza già presenti non salvano i parametri della singola esecuzione; non è possibile attribuirli con certezza ai JSON e alla build attuali. In particolare il CSV SBM esistente dà L2=0.0334669 a n=10, diverso dal nuovo test con q=10. Non va usato come confronto controllato con quest'ultimo.

Riproduzione, nell'ambiente KratosPR:

```bash
python3 AnalyticalTestRevised/results/method_comparison/run_comparison.py body_fitted
python3 AnalyticalTestRevised/results/method_comparison/run_comparison.py sbm
python3 AnalyticalTestRevised/results/method_comparison/run_comparison.py sbm_aligned
```

I file *_parameters.json, *_solve.json e *.log conservano impostazioni, norme e iterazioni dei test. I file *_geometry.json conservano le misure dell'inizializzazione 20×10.
