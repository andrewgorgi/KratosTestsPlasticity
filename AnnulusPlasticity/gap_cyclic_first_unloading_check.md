# Gap-SBM p2: fallimento della prima mesh al primo scarico

**Aggiornamento:** il termine adjoint è stato successivamente congelato sui due
lati dell'interfaccia. Vedere [correzione e verifiche](gap_interface_frozen_adjoint_check.md).
I risultati seguenti descrivono la versione precedente alla correzione.

Controllo del 16 ottobre 2026 con la libreria installata, ricompilata alle 21:30.
Il caso esaminato è mixed, p2, mesh 10×10, dt=0.00025, H=10000,
`PENALTY_FACTOR=-1`, `gap_sbm_type=interpolation` e
`gap_approximation_order=1`.

**I dati al bordo sono corretti nel controllo eseguito, ma il Jacobiano della
`GapSbmSolidInterfaceCondition` non è consistente con il suo residuo plastico.**
Questo impedisce di considerare verificata l'implementazione in generale,
anche se alcune mesh più fini completano il ciclo.

## Il fallimento osservato

Il file del run utente `spans_10x10/status.json` riporta quattro passi completati
e fallimento a t=0.00125: è la prima mesh dello studio che fallisce al primo
scarico, non il primo passo temporale. Il primo picco è a t=0.001.

Il JSON corrente seleziona Newton standard e omette il numero esplicito dei
punti d'integrazione per span. I risultati riusciti già presenti per mesh 40
e 80 sono precedenti e usavano line search e 4 punti per span: non rappresentano
la medesima configurazione. Per isolare questa differenza sono state eseguite
quattro prove con input congelati, fino al primo passo di scarico:

| Strategia | Quadratura bulk | Esito mesh 10×10 |
|---|---|---|
| Newton | default | fallisce a t=0.00125 dopo 4 passi convergenti |
| Newton | 4 punti/span | stesso esito |
| Line search | default | stesso esito |
| Line search | 4 punti/span | stesso esito |

Quindi né cambiare strategia né ripristinare la quadratura risolve da solo
questo fallimento sulla mesh 10.

## Dati di carico

La condition Neumann legge `FORCE`. Nei cinque passi tentati le 120
`GapSbmLoadSolidCondition` ricevono una trazione che concorda con
`C(t)*(y,-x)/r³`, cioè `-C(t)*e_theta/r²` sul cerchio interno.
Errore assoluto massimo 1.14e-13; `BODY_FORCE` è zero su tutti gli elementi.
La legge ciclica e i valori materiali coincidono con il riferimento analitico.
Il macro della `SbmLoadSolidCondition` non seleziona questa condition Gap-SBM.

## Jacobiano globale e localizzazione

Il controllo viene effettuato dopo 100 iterazioni, nello stato fallito a
t=0.00125, senza finalizzare il passo e mantenendo fissa la storia committata.
Si confronta `K*v` con `-[R(u+h*v)-R(u-h*v)]/(2*h)`.

Su tre direzioni casuali, gli errori relativi globali sono 0.09253, 0.04426 e
0.08693. Sono stabili per h=1e-7, 1e-8 e 1e-9: non sono attribuibili al solo
arrotondamento o a un passo FD scelto male.

Sulla prima direzione, il controllo per tipo di entità dà:

| Tipo | Entità | Errore direzionale dei contributi assemblati del tipo |
|---|---:|---:|
| SolidElement | 180 | 1.54e-10 |
| GapSbmSolidElement | 504 | 1.70e-10 |
| GapSbmSolidCondition, Dirichlet esterno | 160 | 7.89e-11 |
| GapSbmLoadSolidCondition, Neumann interno | 120 | 0; carico indipendente da u |
| GapSbmSolidInterfaceCondition | 280 | **0.4847** |

Il peggior errore locale dell'interfaccia è 0.8505, condition 247. La somma
dei prodotti delle matrici locali ricostruisce il prodotto della matrice globale
con errore relativo 4.81e-16. L'errore FD globale ricostruito è 0.09252804156,
uguale al controllo globale indipendente. Ripristinando i coefficienti dopo le
perturbazioni, le RHS locali tornano identiche: errore di ripetizione zero.
Il controllo non rileva quindi una contaminazione della storia da questi trial.

## Termine non linearizzato nell'interfaccia

In `gap_sbm_solid_interface_condition.cpp`, il termine adjoint della RHS usa
`DB_sum_plus/minus = D_alg_plus/minus(u)*B_sum_plus/minus`, rivalutato nello
stato trial corrente, moltiplicato per il salto di spostamento.
La LHS include la derivata del salto, ma tratta quel `DB` come costante.
In forma schematica:

```
R_adj(u) = A(D_alg(u))^T * jump(u)
dR_adj  = A^T * d(jump) + d(A)^T * jump
```

Il secondo termine manca nella linearizzazione implementata. In elasticità
D è costante; in plasticità D_alg può variare con u anche con hardening lineare.
Questo è un difetto concreto della linearizzazione, coerente con la discrepanza
misurata. Non dimostra, da solo, che sia l'unica causa di ogni arresto di Newton
su tutte le mesh.

Per un Newton consistente bisogna includere tale derivata oppure definire
esplicitamente un operatore adjoint congelato, usandolo coerentemente sia nella
RHS sia nella LHS. Nessuna di queste modifiche è stata applicata in questo check.
Il sorgente Kratos, la libreria installata e i JSON di produzione sono invariati.

## Riproduzione e risultati

I quattro casi e i log sono in [diagnostics/gap_cyclic_first_unloading](diagnostics/gap_cyclic_first_unloading).
Lo script diagnostico ferma intenzionalmente il carico a t=0.00125; non è un
runner per il ciclo completo. Usa copie congelate degli input e scrive nella
propria cartella, senza cambiare gli input di produzione:

```bash
python3 AnnulusPlasticity/diagnostics/gap_cyclic_first_unloading/run.py --strategy newton_raphson --quadrature 0 --audit-jacobian
```

- [Identità del sorgente e libreria](diagnostics/gap_cyclic_first_unloading/source_info.json)
- [Jacobiano globale](diagnostics/gap_cyclic_first_unloading/newton_raphson_qdefault/jacobian/audit.json)
- [Controllo per tipo di entità](diagnostics/gap_cyclic_first_unloading/newton_raphson_qdefault/jacobian/entity_audit.json)
- [Dati al bordo](diagnostics/gap_cyclic_first_unloading/newton_raphson_qdefault/boundary_audit.json)
