# Pulizia chimica del grafene trasferito con PMMA/PPC

Ricerca sulla rimozione di residui dalla faccia esposta di grafene CVD monostrato già trasferito su Si/SiO₂. Il processo riferito usa PMMA AR-P 672.045 e PPC; il materiale viene poi impiegato per la raccolta a secco e la costruzione di strutture sovrapposte.

**La pulizia del campione reale non è ancora dimostrata.** Questo repository conserva teoria, calcoli eseguiti, controlli, risultati negativi e proposte sperimentali. Le energie di piccoli aggregati nel vuoto, il distacco imposto da una forza e i controlli numerici non costituiscono una soluzione sperimentalmente validata.

Il [paper, revisione 7](outputs/pulizia-grafene-pmma-ppc.pdf) contiene 21 pagine e 26 riferimenti. Il [sorgente LaTeX](outputs/pulizia-grafene-pmma-ppc.tex) conserva la tipografia AMS/Computer Modern del documento di riferimento. La compilazione è stata verificata, insieme all'impaginazione delle pagine.

La proposta operativa principale confronta la sequenza già provata con DMF seguito direttamente da THF, separando l'effetto del risciacquo intermedio in IPA. HFIP rimane un candidato motivato, da confrontare sul materiale. Nessuna probabilità di successo è stata ricavata dai calcoli.

| Contenuto | Percorso |
| --- | --- |
| Paper, figure e archivio completo | [outputs](outputs) |
| Codice molecolare, parametri, geometrie, dati e registri | [work/molecular](work/molecular) |
| Modelli teorici e verifiche | [work](work) |
| Note di ricerca e letture critiche | [notes](notes) |
| Versioni precedenti degli elaborati | [history](history) |
| Impronte dei file nell'archivio corrente | [MANIFEST-ARCHIVIO.json](MANIFEST-ARCHIVIO.json) |

I principali risultati comprendono calcoli quantistici con controlli della base e della griglia; simulazioni di PMMA4 su grafene in cinque solventi; 27 finestre guidate da 500 ps; identificazione di un lento cambiamento di orientazione e contatto della catena; controlli dei modelli appresi; correzione locale della derivata della pressione; campionamenti molecolari di HFIP e controlli dell'autoassociazione dei solventi. Nessun profilo di energia libera del distacco o confronto quantitativo definitivo fra pulenti è stato accettato.

Le proprietà del liquido sono ancora da qualificare. Rimangono da completare i confronti fra partenze e dimensioni della cella, il trattamento di PPC e catene più lunghe, la selettività rispetto all'interfaccia grafene–SiO₂ e la verifica sul campione reale. I dati interrotti o rifiutati conservano la relativa etichetta; i registri sono istantanee datate. Le traiettorie ancora attive al momento dell'archiviazione non sono presentate come complete.

Gli archivi ZIP sono gestiti con Git LFS. Dopo il recupero del repository, `git lfs pull` recupera gli oggetti disponibili sul server. I file estratti dell'archivio corrente sono già consultabili nel repository.

L'ultima serie molecolare e il successivo controllo energetico sono terminati dopo la chiusura dell'archivio corrente. I file sono aggiunti separatamente e descritti in [MANIFEST-DATI-SUCCESSIVI.json](MANIFEST-DATI-SUCCESSIVI.json), con impronte verificabili. Le proposte più lunghe riducono la differenza di energia interna fra le due configurazioni confrontate, da circa 20 a 6 kJ/mol per molecola, mentre le densità restano circa 1,265 e 0,756 g/cm³. Sono confronti di configurazioni singole: nessuna proprietà del liquido all'equilibrio ne è stata dedotta.

La [nota sul campionamento dell'associazione](outputs/campionamento-associazione.txt) deriva una proposta che sposta e orienta le molecole verso possibili partner, con il rapporto di Metropolis-Hastings completo. Quattro controlli su distribuzioni note sono stati eseguiti, compresi i controlli negativi. La proposta geometrica non è ancora stata applicata a HFIP o all'interfaccia del campione.

Per controllare la corrispondenza fra archivio, impronte e file estratti:

```sh
python3 scripts/verifica_archivio.py
```

Le istruzioni scientifiche sono in [work/molecular/LEGGIMI.txt](work/molecular/LEGGIMI.txt). Sono conservate le versioni dei pacchetti e dei motori di calcolo. Gli ambienti installati, le cache e i pesi dei modelli scaricati sono dipendenze esterne: le fonti e le impronte necessarie per recuperarli sono in [work/molecular/sources](work/molecular/sources). Le modifiche locali devono essere applicate a un ambiente isolato; i pesi dei modelli rimangono invariati.

La ricerca è stata redatta ed eseguita con Codex su richiesta di Mattia Apicella. I componenti esterni conservano riferimenti e licenze; il codice ByteFF-Pol incluso è accompagnato dalla sua [licenza](work/molecular/vendor/byteff-pol/LICENSE).
