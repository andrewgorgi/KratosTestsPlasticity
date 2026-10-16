# Ciclo Newton/line search con SBM p2 e proiezione analytical

Verifica del 16 ottobre 2026. Riproduce **gli stessi coefficienti** del log:
0.335052 e 0.374346, alternati dopo prove con alpha=0.1.

## Configurazione riprodotta

SBM mixed, p2, mesh 40×40, proiezione `analytical` su entrambi i bordi,
lambda inner/outer=0.5, materiale J2 con H=10000 e `PENALTY_FACTOR=-1`.
Time step 0.001, massimo 100 iterazioni, criterio AND con le quattro
tolleranze originali. Il fallimento avviene **al primo step, t=0.001**.
Non è quindi necessario un errore nel commit degli step precedenti per
produrre il ciclo.

Le prove sulle mesh 10 e 20 terminano, con errori L2 relativi finali
0.0884387 e 0.0853124. La mesh 40 entra nel ciclo indicato dall'utente.
Risultati: `diagnostics/analytical_cycle`.

Attenzione al confronto con i rapporti precedenti: adesso sono cambiate sia
la proiezione (`analytical`) sia lambda (0.5 anziché 1). L'effetto della sola
proiezione non può essere ricavato confrontando direttamente quelle serie.

## Il residuo non è vicino alla precisione numerica

Nel ciclo il residuo relativo oscilla circa fra 0.00175436 e 0.00209884;
il residuo assoluto normalizzato del criterio Kratos è circa 6.61e-5–7.91e-5,
contro la soglia 1e-10. Anche il criterio di spostamento resta lontano dalla
convergenza. Il residuo globale non normalizzato, ricostruito nel controllo,
è circa 0.13268.

I valori `r1`, `r2`, `rf` stampati dalla line search non sono la norma globale
del residuo: il codice calcola prodotti scalari fra un incremento e il
vettore RHS. Un `rf` dell'ordine di 1e-6 non significa quindi che il criterio
globale sia quasi soddisfatto.

## Tangente globale nello stato oscillante

È stata riassemblata la matrice alla deformazione trial effettiva del ciclo,
senza finalizzare lo step. Il confronto è

`K*d = -[R(u+h*d)-R(u-h*d)]/(2*h)`.

Tre direzioni casuali indipendenti, con storia committata e operatore
Dirichlet congelato invariati:

| h | Massimo errore relativo direzionale |
|---:|---:|
| 1e-7 | 1.27e-11 |
| 1e-8 | 1.18e-10 |
| 1e-9 | 1.16e-9 |

Questi controlli non sono una verifica elemento per elemento di ogni voce
del Jacobiano, ma non mostrano una tangente mancante o incoerente nel ciclo.
Confermano e ampliano il precedente controllo locale delle leggi costitutive.

Test riproducibile, che usa i parametri correnti senza modificarli:

```bash
python3 AnnulusPlasticity/audit_global_cycle.py --spans 40 --degree 2
```

Per un controllo più breve, il ciclo è già presente con
`--max-iterations 30`, usato per i dati in `diagnostics/global_cycle_audit`.
La cartella conserva parametri, esito, audit e matrici sparse.

## Modi negativi concentrati vicino al foro

Fra gli otto autovalori della matrice globale più vicini a zero nello stato
oscillante compaiono quattro valori reali negativi, circa:

`-106.123, -105.775, -105.774, -105.438`.

Fra gli otto autovalori più vicini a zero della matrice valutata a spostamento
nullo, con identica geometria e storia vergine del primo step, si trovano
invece valori reali positivi fra 12.67 e 17.95. Quest'ultimo è un confronto
di tangenti elastiche, non un secondo equilibrio sotto lo stesso carico;
non significa che sia stato verificato l'intero spettro elastico.

I quattro modi negativi hanno **99.4–99.5% della norma quadrata dei coefficienti
su control point con r<1.3**, vicino al bordo interno Neumann. Il raggio medio
pesato è circa 0.866: i coefficienti più attivi appartengono al supporto IGA
che si estende verso l'interno del foro. Avere control point fuori dal dominio
fisico non è di per sé un errore in un metodo embedded.

Questo costituisce un indizio concreto di perdita di stabilità algebrica
nella zona Neumann plastica. La matrice è non simmetrica: il solo spettro
non dimostra una specifica instabilità fisica, né prova che ogni formulazione
SBM con p2 debba fallire. Restano da discriminare formulazione del termine
Neumann, stabilizzazione e implementazione del suo trasferimento.

Matrici, autovalori e localizzazione dettagliata sono in
`diagnostics/analytical_cycle_jacobian`; figura:
`diagnostics/analytical_cycle_jacobian/negative_neumann_mode.png`.

## Prove per eliminare il ciclo

Tutte mantengono le tolleranze originali e la geometria corrente:

| Prova sulla mesh 40 | Esito |
|---|---|
| Newton puro, dt=0.001 | Fallisce a t=0.001 |
| Smorzamento fisso alpha=0.5 | Fallisce a t=0.001 |
| Smorzamento fisso alpha=0.25 | Fallisce a t=0.001 |
| Line search con 12 prove, min alpha=1e-6, max alpha=1, tolleranza interna 0.1 | Fallisce a t=0.001 |
| Line search standard, dt=0.0005 | Fallisce a t=0.001 |
| Line search standard, dt=0.0001 | Fallisce a t=0.001 |
| Penalità positiva 6700000, dt=0.001 | Fallisce a t=0.001 |
| Backtracking diagnostico sulla norma globale RHS | Stagna con norma circa 0.1373; nessuna riduzione accettabile dopo 24 dimezzamenti |

L'ultima prova usa lo stesso assembly e solver lineare, ma una diversa
globalizzazione Newton; non è stata introdotta nel runner di produzione.

Con una stima iniziale costruita campionando lo spostamento analitico sui
control point **dopo** il congelamento dell'operatore Dirichlet, Newton puro
termina i due step. Tuttavia l'errore L2 finale è **0.21013** e il massimo
equivalent plastic strain è 0.08766. Non è una soluzione pratica del problema:
richiede conoscere il riferimento e converge a un'approssimazione poco accurata.
Dimostra che esiste almeno una soluzione algebrica raggiungibile, non che la
branch numerica desiderata sia stabile o accurata.

## Controlli di consistenza aggiuntivi

Con la configurazione analytical corrente, p2 e mesh 40:

- patch affine elastico: L2 relativo **1.38e-15**, H1 **6.83e-15**;
- patch affine con hardening: Newton termina ma L2 relativo **0.01415**,
  H1 **0.08126**, massimo equivalent plastic strain **0.06195**.

Il trasferimento dei dati al bordo passa quindi un controllo affine elastico
molto accurato; la convergenza Newton in plasticità non garantisce da sola
l'accuratezza o la selezione della branch.

Sono state provate, compilando copie della sola libreria IGA in `/tmp`, anche
estensioni più lunghe del gradiente Neumann. Il codice corrente include
derivate della base fino a p, quindi Taylor del gradiente fino a p-1. Sulle
basi tensoriali Qp i termini successivi possono essere non nulli:

- usando derivate fino a p+1, la mesh 40 termina ma con L2 relativo 238.43;
  la mesh 80 fallisce. Variante scartata;
- usando tutte le derivate polinomiali Q2, fino all'ordine 4, l'annulus mesh 40
  termina con L2 0.00871, ma il patch affine plastico dà L2 0.26663 e la mesh
  20 dell'annulus dà L2 0.11283. Non è una correzione validata.

I risultati di queste varianti sono nelle cartelle `gradient_p_plus_1_*` e
`gradient_full_qp*`. Non sono state applicate ai sorgenti o alla libreria
installata. Non basta eliminare il ciclo su una singola mesh per accettare
una modifica dell'operatore.
La prova aggiuntiva della variante Q2 completa sulla mesh 80 è stata interrotta:
il controllo affine plastico aveva già invalidato questa candidata. Non viene
riportato un esito di convergenza per quella prova incompleta.

## Stato del lavoro

Il ciclo è riprodotto e la diagnosi restringe il problema alla stabilità
globale vicino alla Neumann plastica. Non è stata trovata una correzione
verificata che garantisca insieme convergenza e accuratezza. Non sono state
allentate le tolleranze, cambiati i parametri di produzione o modificati i
sorgenti Kratos. Il prossimo intervento deve verificare i termini Neumann
e la stabilizzazione dei modi di supporto vicino al foro, usando sia il patch
affine plastico sia il test non affine come controlli.
