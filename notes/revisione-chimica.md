# Revisione indipendente della proposta chimica

Ricerca primaria verificata il 6 ottobre 2026. Oggetto: residui presumibilmente PMMA AR-P 672.045 su grafene CVD trasferito da Cu a SiO2/Si. Questa nota distingue risultati pubblicati, trasferibilità ancora da verificare e scelte progettuali originali. Non certifica un campione mai misurato.

## Aggiornamento decisivo: processo reale con PMMA e PPC

Dati forniti dall'utente dopo la prima ricerca: PMMA trattato a 90°C per 2 min; distacco elettrochimico dal rame in un bagno descritto come acqua, elettrolita ancora ignoto; ulteriore strato sottile di PPC in anisolo, inizialmente non a contatto diretto con il grafene. Il processo reale va quindi discusso come trasferimento mediante delaminazione elettrochimica e supporto PMMA/PPC, non come attacco umido FeCl3 o APS. Non introdurre quei reagenti nella storia del campione.

Il candidato HCl sotto riportato è un'estrapolazione da una diversa formulazione/processo; manca verifica sul sistema PMMA950K/PPC. Il dato 90°C/2 min non prova carbonizzazione né reticolazione. PPC significa qui poli(carbonato di propilene), da confermare sulla scheda del prodotto, e non il PC aromatico bisfenolo-A usato in altri studi di trasferimento.

Proposta diagnostica prioritaria:

- Preparare campioni di riferimento con lo stesso PMMA, lo stesso PPC e la miscela/sovrapposizione, con uguali solventi, trattamento termico e rimozione. La posizione iniziale del PPC sopra il PMMA riduce il contatto diretto atteso, ma non esclude contaminazione durante rimozione, ai bordi o attraverso discontinuità: questa è un'ipotesi, non un fatto osservato.
- Confrontare l'impronta chimica completa con questi riferimenti: XPS C1s e O1s insieme, più ToF-SIMS basata su molteplici frammenti o IR su residui sufficientemente abbondanti/tecnica locale adatta. Il solo picco carbonilico non identifica PMMA perché PMMA e PPC contengono entrambi carbonio ossigenato e carbonile. In IR, gruppi estere e carbonato e l'intero insieme di bande possono aiutare, ma non promettere discriminazione automatica di un residuo subnanometrico e non usare una sola frequenza di tabella come verdetto.
- Se il risultato resta ambiguo, marcatura isotopica di uno dei supporti in un nuovo trasferimento è un controllo più discriminante, ma costituisce una prova nuova e impegnativa. La fattibilità della marcatura PMMA e ToF-SIMS è pubblicata da Wang et al., *Direct Observation of Poly(Methyl Methacrylate) Removal from a Graphene Surface*, [Chemistry of Materials, DOI 10.1021/acs.chemmater.6b03875](https://doi.org/10.1021/acs.chemmater.6b03875). Non dire che questo lavoro abbia distinto PPC da PMMA nel campione attuale.
- L'assenza di contatto iniziale col PPC non autorizza a escluderlo, così come l'uso del PMMA non autorizza a identificare ogni particella come PMMA. Finché manca l'identità, il test HCl va descritto come prova di recupero empirica con esiti interpretabili, non reazione selettiva già fondata sul contaminante certo.

## Identità del materiale

La [scheda ufficiale Allresist](https://www.allresist.de/wp-content/uploads/2020/03/AR-P630-670_Deutsch_Allresist_Produktinformation.pdf), pp. 1 e 6, identifica AR-P 672.045 come PMMA della famiglia 950K, 4,5% di solidi in anisolo. A 4000 rpm riporta 0,23 µm. Il codice fornito dall'utente AR672.045 è coerente con questo prodotto, ma la bottiglia e il lotto restano da verificare. Non è un copolimero MMA/MAA. La stessa scheda descrive la scissione da elettroni e la reticolazione alle dosi molto elevate: non trasferire questa spiegazione a un campione senza sapere se è stato irradiato. La scheda di un resist per litografia non garantisce pulizia molecolare del grafene.

## Candidato principale: HCl acquoso, con confini molto precisi

Fonte primaria: Z. Xiao, Q. Wan, C. Durkan, *Cleaning Transferred Graphene for Optimization of Device Performance*, Advanced Materials Interfaces 6 (2019), 1801794, [DOI 10.1002/admi.201801794](https://doi.org/10.1002/admi.201801794). [Manoscritto accettato nel deposito Cambridge](https://www.repository.cam.ac.uk/items/ac1b3655-d7c2-4eb1-8dad-f36233a7a676); [testo integrale istituzionale](https://api.repository.cam.ac.uk/server/api/core/bitstreams/310c554f-60e4-4413-b4ce-24d8e0f0e11b/content).

Condizioni: grafene CVD/Cu → SiO2 300 nm/Si; PMMA 495K A4 in anisolo, rimosso inizialmente con acetone 3 h e IPA. Successivo HCl 1 M per 48 h, risciacquo DI e asciugatura N2. I solventi di confronto comprendevano acetone, cloroformio, DMSO e Remover PG. Ra passa da 2,4 a 0,4 nm con HCl; solventi circa 1 nm. Nei dispositivi trattati prima della metallizzazione: Rc 0,54→0,26 kΩ, mobilità 3060→3478 cm²/Vs, VCNP 8→10 V. Sui residui invecchiati tre mesi e riscaldati 200°C/3 h, la copertura cala soltanto 75→34%. L'affermazione aggiuntiva HCl concentrato/6 h non specifica pienamente condizioni e metrologia. L'idrolisi degli esteri è un'interpretazione, non dimostrazione. Nessun test sul particolare AR-P 672.045; nessuna prova di pulizia atomica. L'interpretazione proposta per NaCl/clorometano non va ripetuta come chimica stabilita. Temperatura del bagno non esplicitata nel testo verificato: se si sceglie temperatura ambiente, dichiararla scelta progettuale.

**Giudizio:** miglior candidato trovato per il recupero chimico di campioni già trasferiti, entro il vincolo acidi e senza fluoruri o basi. È un esperimento motivato da una pubblicazione, non una soluzione nuova già dimostrata. Se Armando ha già provato esattamente HCl 1 M/48 h sul suo materiale, questo candidato è già falsificato per quelle condizioni; non nascondere la circostanza cambiando nome alla proposta.

## Alternative: prova pubblicata e limite

- Acido acetico: M. Her, R. Beams, L. Novotny, *Graphene transfer with reduced residue*, Physics Letters A 377 (2013), 1455–1458, [DOI 10.1016/j.physleta.2013.04.015](https://doi.org/10.1016/j.physleta.2013.04.015), [preprint autori](https://arxiv.org/abs/1301.4106). Dimostra un miglioramento rispetto all'acetone, non garantisce residuo nullo. Il lavoro [Frontiers 2023, DOI 10.3389/fmats.2023.1279939](https://www.frontiersin.org/journals/materials/articles/10.3389/fmats.2023.1279939/full) usa 3 h a temperatura ambiente contro acetone 50°C/3 h e ammette PMMA residuo dopo entrambi. Il coordinatore ha identificato inoltre un'incoerenza fra tabella XPS e interpretazione testuale: non adottare il testo come conferma chimica incontestata.
- Formammide: il lavoro di Ruoff *Reducing Extrinsic Performance-Limiting Factors in Graphene Grown by Chemical Vapor Deposition*, Nano Letters (2013), [manoscritto autori](https://utw10193.utweb.utexas.edu/Archive/RuoffsPDFs/340.pdf), afferma che la formammide solvata nel residuo compensa il drogaggio p. È esattamente il motivo per cui miglioramento elettrico non equivale a rimozione del polimero. Non proporla come pulizia provata soltanto perché il punto di neutralità migliora.
- Trattamento NaNO2 acidificato: K. Lee et al., Polymers 17 (2025), 689, [DOI 10.3390/polym17050689](https://doi.org/10.3390/polym17050689), [testo completo](https://pmc.ncbi.nlm.nih.gov/articles/PMC11902524/). Il risciacquo pH 3,5/10 min è effettuato subito dopo l'attacco FeCl3 e **prima** del trasferimento sul substrato, seguito da cloroformio 1 h, monoclorobenzene 30 min, cloroformio 30 min. Non dimostra recupero del campione già trasferito. Le attribuzioni a NO, ossidazione PMMA e neutralizzazione di Cl richiedono conferme indipendenti: non ereditarle acriticamente. Il protocollo introduce sodio e un'altra variabile contaminante, e non è la priorità di questa proposta.
- Trattamento assistito da protoni: G. Yuan et al., Nature Communications 14 (2023), 5457, [DOI 10.1038/s41467-023-41296-5](https://doi.org/10.1038/s41467-023-41296-5), [testo](https://pmc.ncbi.nlm.nih.gov/articles/PMC10482836/). Il nome non indica un semplice bagno acido: è un trattamento con plasma di H2; quindi fuori dal quesito chimico liquido.
- Miscela IPA/acetone/MIBK: H. J. Jeong et al., Carbon 66 (2014), 612–618, [DOI 10.1016/j.carbon.2013.09.050](https://doi.org/10.1016/j.carbon.2013.09.050). Il risultato è associato a irraggiamento UV precedente. Non separare la miscela da quel passaggio per attribuirle l'effetto completo.
- Acidi di Lewis/scissione selettiva: la ricerca mirata non ha reperito una dimostrazione primaria abbastanza specifica su PMMA950K/grafene/SiO2 da giustificare BBr3 o simili come ricetta. Una conversione dell'estere laterale non è una depolimerizzazione della catena carboniosa; potrebbe lasciare un altro polimero adsorbito. Nessuna miscela nuova va presentata come funzionante senza esperimenti.

## Perché le particelle non identificano da sole PMMA

[Lupina et al., ACS Nano 9 (2015), DOI 10.1021/acsnano.5b01261](https://arxiv.org/abs/1505.00889) misurano metalli residui Cu/Fe oltre 10^13 atomi/cm² anche dopo trasferimenti umidi o delaminazione elettrochimica. Questo stabilisce una categoria plausibile di contaminanti, non l'identità delle particelle di Armando.

[Li et al., Langmuir 34 (2018), 1827–1833, DOI 10.1021/acs.langmuir.7b03117](https://pubmed.ncbi.nlm.nih.gov/29303580/) ricavano, dopo acetone/IPA, un modello di strato interno 17 Å e strato esterno diffuso 31 Å, mediante riflettometria neutronica e AFM. Una superficie liscia non dimostra assenza di un film sottile.

Ipotesi da distinguere: polimero sulla faccia esposta; residui metallici/ossidi/sali; contaminante sotto il grafene; bolle/pieghe della membrana; carbonio amorfo o PMMA chimicamente modificato. XPS C1s deve essere interpretata insieme ai segnali elementali e ad altri controlli: sp3 e carbonio ossigenato non sono marcatori esclusivi del PMMA. AFM altezza da sola non identifica la composizione né quale lato del grafene ospita il contaminante. Un bagno dall'alto non ha necessariamente accesso a un'interfaccia sepolta.

## Disegno di prova: contributo proposto, non dati

1. Registrare formulazione, solvente, lotto, spessore del supporto, ogni trattamento termico, esposizione a elettroni/UV, etchant del rame, tempi trascorsi e bagni già tentati. Il codice PMMA non sostituisce la storia del campione.
2. Suddividere lo stesso lotto in campioni indipendenti randomizzati: nessun nuovo trattamento; acqua DI per lo stesso tempo; solvente di riferimento concordato; HCl 1 M/48 h. Prevedere più campioni per braccio e più campi per campione; i campi non sono repliche indipendenti. Ulteriori tempi sono esplorativi, non il parametro pubblicato. Conservare aree gemelle non prescansionate.
3. Verificare sul medesimo insieme di aree la copertura residua, volume e distribuzione delle altezze, non solo Ra. Scansioni AFM delicate e aree non prescansionate evitano scambiare la pulizia meccanica della punta per quella chimica. Mappare anche campi più ampi e bordi per intercettare semplice spostamento/aggregazione.
4. Prima/dopo: segnali chimici compatibili con PMMA (XPS o ToF-SIMS su provini gemelli), Cu/Fe/Cl e altri elementi effettivamente introdotti. Analizzare il risciacquo con metodo adatto se disponibile: la scomparsa topografica non basta. Controllo SiO2 senza grafene sottoposto agli stessi bagni.
5. Escludere asportazione del grafene: continuità spaziale Raman/ottica/AFM, area coperta, comparsa di rotture; la diminuzione del segnale polimerico con perdita della membrana è fallimento. Misurare difetti Raman, drogaggio/isteresi e prestazioni, senza confondere migliore resistenza con superficie più pulita.
6. Stabilire **prima** dell'esperimento soglie compatibili con l'uso e con la ripetibilità strumentale. Qualsiasi riduzione del 90%, tolleranza 10%, o soglia AFM espressa in nm deve essere etichettata come criterio progettuale proposto, non standard della letteratura né dato ottenuto.
7. Ripetere su un secondo lotto se il primo supera tutti i criteri. Se composizione o posizione smentiscono PMMA esposto, interrompere l'inseguimento del solvente e modificare l'ipotesi. Se resta un film uniforme, non proclamare successo per la scomparsa delle isole.

## Formulazione sostenibile per la conclusione del paper

«Identifichiamo un trattamento chimico documentato che resta da validare sul sistema PMMA AR-P 672.045/PPC ottenuto per delaminazione elettrochimica e formuliamo un confronto capace di distinguerne l'efficacia da redistribuzione, drogaggio e perdita di grafene. Non riportiamo nuovi dati sperimentali né rivendichiamo pulizia universale. La prova sul campione richiesto è ancora necessaria.»

# Revisione della bozza completa (6 ottobre 2026)

File letto: outputs/pulizia-grafene-pmma-ppc.tex. Nessuna modifica apportata al sorgente. La struttura epistemica è corretta: non ho trovato una falsa dichiarazione di successo, di depolimerizzazione dimostrata o di prova diretta sul PPC. Le correzioni sotto rendono effettivi alcuni controlli ancora incompleti.

## 1. La scomparsa della firma iniziale non prova rimozione del polimero

Punto: identità/firme chimiche, poi criterio «diminuzione coerente della firma chimica», righe 278–291 e 325–328 della versione letta.

L'articolo stesso contempla trasformazioni esteriche. Se HCl converte una frazione in materiale più ricco di COOH, il marcatore iniziale PMMA può diminuire mentre resta materiale carbonioso, eventualmente appiattito e invisibile come particella. Questa combinazione potrebbe superare impropriamente i criteri morfologici e chimici.

Correzione minima: aggiungere riferimenti PMMA/PPC sottoposti agli stessi trattamenti e precisare che la perdita del solo marcatore esterico non vale come rimozione. Ricercare segnali dei prodotti trasformati e residuo organico complessivo mediante misure complementari. Confrontare quantità, limiti e normalizzazioni dichiarati, non soltanto le percentuali relative delle componenti XPS.

## 2. Controllo anti-asportazione: il 99% globale non basta

Punto: riduzione90% e conservazione99%, righe 325–328.

È possibile che l'1% di area di grafene perso contenga il90% del residuo. La soglia globale non impedisce di contare la rimozione del grafene sporco come pulizia.

Correzione minima: calcolare riduzione del residuo sull'intersezione registrata delle aree in cui il grafene è presente sia prima sia dopo; escludere le aree asportate dalla categoria pulizia e riportarle separatamente come danno. Questa è una definizione della metrica, non una nuova soglia arbitraria.

## 3. La sequenza W/S/A/AS non è ancora interamente fissata

Punto: righe 252–270.

«Stesso tempo di manipolazione» non specifica dove stazionano W e A durante i30min di anisolo degli altri campioni, né quale sia il liquido finale. Attesa a secco, acqua e IPA non sono equivalenti. Lo scambio IPA/aniso/IPA può incidere su precipitazione, rigonfiamento, riorganizzazione e asciugatura. Un bianco pulito è necessario, ma non esclude rideposizione selettiva sul grafene quando il bianco non ha lo stesso supporto e le stesse condizioni.

Correzione concreta: fase iniziale e finale IPA e asciugatura comuni; nella fase centrale S/AS ricevono anisolo fresco2×15min, W/A IPA fresco2×15min, con ricambi e volumi abbinati. Dichiarare quindi che si misura l'effetto della sostituzione dell'IPA centrale con anisolo nell'intera sequenza. L'effetto di una sequenza non coincide con solubilità intrinseca o scissione. In alternativa, fissare un altro percorso abbinato ma rimuovere l'espressione «solo contributo estrattivo».

Se si desidera la riproduzione più fedele di Xiao, mantenerla come prima prova HCl/acqua +DI/N2; il fattoriale con IPA è un esperimento distinto e più elaborato. Questo è già abbastanza chiaro nella struttura, basta non confondere i due risultati.

## 4. Misure chimiche distruttive e gemelli

Punto: ToF-SIMS insieme a generiche misure iniziali/finali.

ToF-SIMS bombarda il materiale; l'area analizzata non è un campione intatto per il successivo bagno di pulizia. Per XPS va almeno documentata la dose/esposizione se si usa il medesimo provino.

Correzione minima: misure distruttive su provini gemelli dedicati; le aree bombardate non entrano nella prova prima/dopo. AFM/Raman restano registrabili nelle stesse aree con controlli dose e interazione.

## 5. Decisione causale: confronto con il controllo, non solo prima/dopo

Punto: «Se A e AS migliorano mentre S fallisce, il bagno acido ha un contributo misurabile».

Correzione minima: «Se A migliora rispetto a W e AS rispetto a S, nelle condizioni e nell'incertezza misurate, si sostiene un contributo dell'acido». Miglioramenti di A e AS rispetto al proprio valore iniziale non bastano: anche W potrebbe migliorare ugualmente. Nessun bisogno di promettere sinergia; se la si volesse studiare servirebbe una differenza delle differenze con incertezza e scala della metrica dichiarate.

## 6. Refuso

Riga 291 della versione letta: «Quando disponibile,+una» → «Quando disponibile, una».

## Valutazione complessiva

Il precedente HCl è riportato fedelmente, con la differenza PMMA495K/950K e la temperatura del bagno non ricostruita. La cautela su PPC, residui sepolti, pirolisi90°C, elettricità e XPS è adeguata. L'anisolo è ammissibile come fattore empirico, ma non ha evidenza di sinergia post-HCl: la bozza lo dice correttamente. Le correzioni1–5 sono sufficienti; non serve aggiungere altri solventi o una reazione speculativa.
