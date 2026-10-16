# SBM misto p=2: controlli di implementazione e plasticità perfetta

Verifica successiva del ciclo costitutivo al punto proiettato:
[projection_history_check.md](projection_history_check.md). Il confronto su
384 stati con due storie indipendenti passa nella libreria aggiornata;
i problemi globali p2 restano riproducibili. Il rapporto distingue anche una
precedente build con ordine Taylor Neumann forzato a 1 dal sorgente corrente.

## Configurazione effettivamente letta

Le prove leggono il JSON misto salvato, non lo rigenerano:

- `PENALTY_FACTOR=-1` nel materiale attuale;
- `lambda_inner=lambda_outer=1` nel JSON misto;
- proiezione `linealized`, grado 1 o 2, quadratura bulk 4×4;
- line search, massimo 100 iterazioni, tolleranze originali;
- Δt=0.001, salvo il controllo iniziale elastico che usa Δt=0.0001
  e si arresta a t=0.0002.

Questi parametri non sono quelli di alcune vecchie prove con penalizzazione
positiva. Inoltre il generatore attuale contiene lambda=0.5 per il misto,
mentre il JSON salvato ha ancora 1: rigenerarlo cambia il confronto.
I controlli non modificano materiale/JSON di produzione o sorgenti Kratos.
Usano copie isolate, salvando i parametri di ogni prova.

## Un H=0 diretto nell'anello non è un riferimento valido

La torsione anulare corrente impone tau=C/r². Al carico finale:

```
q(r=1) = sqrt(3)*C_final = 452.25 > sigma_y = 201
```

Con H=0, questo stress viola il criterio di plasticità perfetta. Nel campo
puramente tangenziale il carico ammissibile è C ≤ 201/sqrt(3), circa
116.05; al limite snerva il bordo interno, senza la fascia plastica
1 < r < 1.5 del caso con hardening. Le formule dello spostamento con 1/H
non si possono riutilizzare ponendo H=0.

Per il confronto sono stati quindi costruiti due riferimenti validi:

- `patch_perfect`: u=(10*t*x,0), stato uniforme, b=0, stress J2 perfettamente
  plastico da carico proporzionale; solo test di implementazione, non uno
  studio di convergenza non polinomiale;
- `cosine_perfect`: spostamento manufactured cosinusoidale `cosine_middle`
  di `AnalyticalTestRevised/analytic_fields.py`, con stress J2 H=0,
  **forcing b=-div(sigma) e Neumann sigma*n coerenti** sull'anello.
  Questo secondo controllo è non affine, ma non ha forcing nullo.

In entrambi i casi la normale fisica è -e_r all'interno e +e_r all'esterno.
La deformazione è proporzionale monotona nel tempo, come richiesto dal
riferimento radiale. I test delle funzioni analitiche controllano anche
la divergenza numerica delle tensioni del caso cosinusoidale.

## Tangente locale e residuo del patch esatto

`check_local_tangent.py` confronta K*d con -dR/du*d mediante differenze
centrali, h=1e-7 e 1e-8, su tre condizioni per tipo e per grado, in uno
stato trial affine plastico. Non finalizza la storia durante le perturbazioni.

| Grado | SbmSolidCondition: massimo errore relativo | SbmLoadSolidCondition |
|---|---:|---:|
| 1 | 7.4e-10 | 2.2e-9 |
| 2 | 1.6e-9 | 1.7e-9 |

`check_patch_residual.py` assegna i coefficienti della soluzione affine
esatta e assembla il RHS completo, senza risolvere Newton. Nelle prove
elastiche, con hardening e con plasticità perfetta, per p=1 e p=2 su
10×10 e 20×20 il residuo relativo resta circa 1.8e-6–2.5e-6.
La normalizzazione usa la norma del vettore RHS delle condizioni di carico.
Il contorno linearizzato e le normali ricostruite non coincidono esattamente
con i cerchi, perciò non si pretende uno zero a precisione macchina.

Questi controlli non mostrano un errore macroscopico locale nella derivata
del residuo o nella riproduzione affine a p=2. Non verificano però ogni
stato, l'intera storia al punto proiettato o la stabilità globale.

## Confronti numerici sulla mesh 20×20

| Controllo | p=1 | p=2 |
|---|---|---|
| Anello prima dello snervamento | L2=0.1761, alpha=0 | L2=0.02528, alpha=0 |
| Patch affine elastico | completato | completato |
| Patch affine con H=10000 | L2≈6.0e-7 | Newton termina, L2≈13.61 |
| Patch affine perfettamente plastico | L2≈7.1e-7 | Newton termina, L2≈2180 |
| Coseno non affine perfettamente plastico | L2≈0.09730 | fallisce al primo passo |
| Anello con H=10000 | L2≈0.5931 | Newton termina, L2≈67.23 |

Il caso anulare p=2 è stato anche riprodotto sulla 80×80: fallisce al primo
passo t=0.001. Una piccola norma del residuo non garantisce che la soluzione
raggiunta sia vicina a quella fisica: i patch p=2 mostrano proprio questo.
Il patch completamente perfettamente plastico può inoltre avere una
tangente degenere e meccanismi non unici; da solo non è una prova decisiva.

Due ulteriori variazioni diagnostiche, senza cambiare i file di produzione:

- PENALTY_FACTOR=6700000, p=2, 20×20: il patch con hardening fallisce
  al secondo passo e l'anello al primo;
- proiezione NURBS `analytical`, p=2, 20×20, penalizzazione -1: il patch
  con hardening termina con L2≈4.45; il coseno H=0 fallisce ancora;
  l'anello termina con L2≈91.95.

Quindi né il solo cambio di penalizzazione né la sola proiezione esatta
rimuovono il problema nelle prove eseguite.

## Conclusione e limiti

Le prove confermano un problema p=2 nella risposta plastica globale del
setup/implementazione attuale, anche su un riferimento diverso dall'anello
con hardening. **Non confermano un difetto esclusivo dell'hardening**:
anche il controllo non affine H=0 fallisce.

Non è dimostrato che sia un limite intrinseco del metodo SBM. Restano
da distinguere stabilità della formulazione discreta, aggiornamento della
storia durante percorsi non proporzionali e comportamento globale del solver.
Il superamento delle verifiche locali non esclude un errore globale o di
storia. Nessuna correzione C++ è stata applicata sulla base di queste prove.

## Riproduzione

Dalla radice del workspace:

```bash
python3 AnnulusPlasticity/check_local_tangent.py
python3 AnnulusPlasticity/check_patch_residual.py

# Confronto perfettamente plastico non affine, Neumann e forcing coerenti
python3 AnnulusPlasticity/diagnose_mixed.py --cases cosine_perfect --degrees 1 2 --meshes 10 20

# Tutti i controlli con i parametri attualmente salvati
python3 AnnulusPlasticity/diagnose_mixed.py

# Variazioni isolate, senza alterare materials.json
python3 AnnulusPlasticity/diagnose_mixed.py --cases patch_hardening annulus_hardening --degrees 2 --penalty-factor 6700000 --output AnnulusPlasticity/diagnostics/positive_penalty
python3 AnnulusPlasticity/diagnose_mixed.py --cases patch_hardening cosine_perfect annulus_hardening --degrees 2 --meshes 20 --projection-type analytical --output AnnulusPlasticity/diagnostics/analytical_projection
```

Ogni controllo salva `parameters.json`, `materials.json`, `result.json` e,
solo se il solve termina, i campioni in `solution.npz`. In caso di errore
salva il traceback e continua con gli altri controlli: consultare `status`
nel riepilogo, non soltanto il codice di uscita dello script.
Le prove riportate sono state eseguite anche in invocazioni separate;
i riepiloghi si trovano nelle sottocartelle di `diagnostics/`.
