# Correzione del dato di Neumann SBM sull'anello

## Errore trovato nel runner

Con `projection_type="linealized"`, lo snake discretizza le curve e salva
ai nodi la normale derivata dall'orientamento NURBS. Entrambi i file dei
cerchi hanno orientamento antiorario: il loro `NORMAL` nodale è +e_r.
Sul bordo interno questa non è la normale uscente dal materiale, che è -e_r.

Il runner leggeva quel `NORMAL` per costruire FORCE anche sui nodi interni.
Prima della correzione, al carico finale imponeva quindi **+261.106659 e_theta**
invece di **-261.106659 e_theta**. L'errore massimo vettoriale era circa
522.213318. Il campo FORCE sulle condizioni era già orientato correttamente,
ma `SbmLoadSolidCondition` legge FORCE dal nodo della skin scelto come
proiezione: correggere solo il dato sulla condizione non sarebbe bastato.

Nel ramo linearizzato, la condition ricostruisce inoltre la propria
`mTrueNormal` dai segmenti e dall'identificatore inner/outer. Non si può
quindi assumere che coincida con il `NORMAL` memorizzato al nodo.
Il ramo con proiezione NURBS analitica è diverso: il modeler orienta già
la normale del nodo proiettato in base a inner/outer.

## Dato corretto

`analytic.boundary_traction` ora costruisce la normale fisica dell'anello
direttamente dalla posizione, indipendentemente dal campo nodale NORMAL:

```
n_inner = -(x,y)/r
n_outer = +(x,y)/r
g_N = sigma_exact(x,y,t) * n
```

Il runner usa questa funzione sia per la FORCE nodale letta da SBM,
sia per quella sulla condizione letta da Gap-SBM. I campi NORMAL geometrici
di Kratos non vengono modificati. Per esempio, a fine carico:

- nel punto (1,0): g_N = (0,-261.1066592410083);
- nel punto (0,1): g_N = (261.1066592410083,0).

La classificazione interno/esterno usa il cerchio più vicino e vale per
questo benchmark anulare; non è una regola generale per altre geometrie.

## Verifica riproducibile

```bash
python3 AnnulusPlasticity/audit_neumann.py --degree 1 --spans 20 20
python3 AnnulusPlasticity/test_analytic.py
```

L'audit inizializza il modello effettivo, assegna i dati al tempo richiesto
e controlla la FORCE nodale con una formula indipendente, senza risolvere
il problema meccanico. Salva `sbm/results/mixed/p1/boundary_audit.json`.
Sul modello 20×20 verificato:

- 96 SbmLoadSolidCondition interne e 216 SbmSolidCondition esterne;
- 999 nodi della skin interna con normale memorizzata circa +e_r;
- FORCE tangenziale interna corretta: -261.1066592410083;
- errore massimo del dato interno: circa 2.3e-13;
- errore massimo del dato esterno: circa 5.8e-14.

L'audit riporta anche la risultante e il momento del RHS delle sole
condizioni Neumann a spostamento nullo. Sulla 20×20 il momento è circa
-1314.57, mentre il momento fisico sul cerchio esatto è -1640.58.
Questo confronto non è validato come identità: il RHS SBM usa test e
pesi sul surrogate e la formulazione contiene ulteriori termini di flusso
spostato. La verifica puntuale del dato non valida da sola conservazione,
consistenza o tasso spaziale dell'intero assembly.

## Prove del solve dopo la correzione

Completate le mesh 10,20,40,80,160 con grado 1, line search, Δt=0.001,
lambda interno/esterno=1 e le tolleranze originali. Nessun cambio del
surrogate è stato usato per questa verifica. Le prove sono state eseguite
anche in invocazioni separate con gli stessi parametri.
Questo non garantisce la convergenza per altri gradi, mesh o parametri.

| Mesh | L2 relativo finale |
|---|---:|
| 10×10 | 0.676592 |
| 20×20 | 0.593100 |
| 40×40 | 0.367900 |
| 80×80 | 0.229598 |
| 160×160 | 0.126394 |

Il tasso L2 dell'ultimo raffinamento è circa 0.86: il completamento di
Newton dopo la correzione del carico non dimostra l'ordine spaziale atteso.
Questo limite della verifica del caso misto rimane aperto.

I risultati misti SBM prodotti prima della correzione avevano il segno
del carico interno sbagliato e devono essere rigenerati. La correzione
non modifica i dati di spostamento dei casi Dirichlet.
