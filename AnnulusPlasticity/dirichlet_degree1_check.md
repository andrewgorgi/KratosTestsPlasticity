# Check Gap-SBM, grado 1, solo Dirichlet

## Dati assegnati alle condizioni

Audit del modello effettivamente inizializzato, mesh 20×20:

- 920 `GapSbmSolidCondition`: 320 sul bordo interno, 600 su quello esterno;
- 896 `GapSbmSolidInterfaceCondition`, per l'accoppiamento bulk/gap;
- nessuna condizione di Neumann;
- massimo errore fra `condition.GetValue(DISPLACEMENT)` e il riferimento
  valutato al centro della geometria della condizione: **0**;
- BODY_FORCE su tutti gli elementi: **0**.

Il codice C++ `GapSbmSolidCondition::CalculateRightHandSide` legge proprio
`this->GetValue(DISPLACEMENT)` per il dato `u_D`. Il gradiente e lo spostamento
numerici vengono estesi dal surrogate al punto della condizione mediante
le matrici/funzioni di Taylor. Dirichlet è imposto debolmente con la
penalizzazione/Nitsche della condizione, non fissando i DOF dei control point.

Il contorno Gap-SBM numerico non coincide esattamente con i cerchi NURBS:

| Bordo | Raggio minimo | Raggio massimo |
|---|---:|---:|
| Interno, mesh 20×20 | 0.9951495306 | 0.9998882736 |
| Esterno, mesh 20×20 | 1.9966293897 | 1.9999992383 |

Il runner assegna l'estensione analitica ai punti effettivi di questo contorno.
Di conseguenza, a fine carico il massimo modulo del dato sul bordo esterno
numerico è circa **6.57e-6**, mentre sul cerchio esatto r=2 è zero.
Questo distingue il benchmark discreto sul contorno approssimato dal dato
geometrico ideale; il runner non proietta radialmente il dato su r=1 o r=2.

## Ordine L2: distinguere le mesh grossolane dal regime asintotico

Valori al carico finale (grado 1):

| Mesh | Errore L2 relativo | Ordine rispetto alla precedente |
|---|---:|---:|
| 10×10 | circa 1.080 | — |
| 20×20 | circa 0.835 | 0.37 |
| 40×40 | circa 0.3091 | 1.43 |
| 80×80 | 0.08351453 | 1.89 |

I primi tre valori sono coerenti con la sequenza a Δt=0.0002; il valore
80×80 proviene dal check completato a Δt=0.0001. Il check indipendente
40×40 a Δt=0.0001 dà L2=0.30909874, praticamente uguale al valore a
Δt=0.0002, e conduce ancora a ordine 1.89 rispetto alla 80×80.
Non è quindi una tabella a passo temporale uniformemente identico: serve
un nuovo studio a Δt comune per una misura definitiva dell'ordine.

Ulteriore isolamento: arrestando il caricamento a t=0.0002, con Δt=0.0001,
il materiale resta elastico (alpha massima nulla sulle mesh 10/20/40).

| Mesh | Errore L2 relativo | Ordine |
|---|---:|---:|
| 10×10 | 0.26879489 | — |
| 20×20 | 0.17891730 | 0.59 |
| 40×40 | 0.05290614 | 1.76 |
| 80×80 | 0.01320995 | **2.00** |

Queste prove indicano un regime preasintotico sulle mesh grossolane.
Già sulla 10×10 il bulk occupa solo area 3.872 su circa 9.374 dell'anello
discretizzato; il resto è gap. Sulla 80×80 il bulk occupa area 8.7846 su
circa 9.42478. Questo cambiamento è rilevante per interpretare i primi
raffinamenti; non costituisce da solo una dimostrazione della causa dell'errore.
L'ordine atteso 2 riguarda L2; per H1 l'attesa usuale per grado 1 è 1,
sotto le ipotesi di regolarità e stabilità appropriate.

## Newton e controllo del tempo

Il problema già riprodotto a Δt=0.001 è una mancata convergenza di Newton:
sulla 40×40 il residuo oscilla a t=0.002. Δt=0.0001 ha completato il
caricamento sulle mesh 40×40 e 80×80 senza cambiare le tolleranze.
Questo non garantisce la convergenza per qualunque mesh/passo.

Il check a Δt=0.00005 ha rivelato anche un errore distinto nel runner:
dopo 40 passi, il tempo accumulato era 0.0019999999999999987. Il confronto
stretto con 0.002 faceva partire un passo aggiuntivo, fuori dall'intervallo
del riferimento, generando `ValueError`, non una mancata convergenza NR.
`Analysis.KeepAdvancingSolutionLoop` ora usa una tolleranza floating point;
`test_time_loop.py` verifica i conteggi 40, 20 e 10 passi per i tre Δt.

La sola terminazione di Newton non dimostra l'ordine spaziale. Per uno
studio definitivo occorre mantenere fissi i parametri temporali e verificare
la sensibilità a Δt sulle stesse mesh.

Dopo la correzione del ciclo temporale, anche il check 40×40 a Δt=0.00005
completa esattamente 40 passi. Il confronto temporale sulla stessa mesh è:

| Δt | L2 relativo finale |
|---|---:|
| 0.0002 | 0.3091003819 |
| 0.0001 | 0.3090987405 |
| 0.00005 | 0.3089767619 |

La differenza fra gli ultimi due è circa lo 0.04% dell'errore L2:
in questa prova, dimezzare Δt non rimuove l'errore spaziale osservato.
Non sono stati modificati materiale, penalizzazione o tolleranze durante
questo confronto.
