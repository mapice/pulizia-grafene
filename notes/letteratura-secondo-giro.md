# Secondo giro: pulizia chimica PMMA/PPC su grafene — fonti e revisione critica

Ricerca effettuata il 6 ottobre 2026. Nessuna sperimentazione sul campione. Fonti nuove rispetto alla prima nota; le ipotesi di processo in fondo sono deduzioni da verificare, non risultati degli articoli.

## Dati del caso che cambiano le priorità

Armando riferisce di avere già provato HCl per 48 h, acido acetico glaciale, acetato di etile, cloroformio, DMF, THF, MPK (sigla da mantenere tale), acetone, IPA/acqua 75/25, AR600.7 e AR300-76 a temperatura ambiente e a 80 °C. L'unico effetto scarso è associato alla sequenza DMF 30 min → IPA 5 min → IPA 3 min → THF 30 min → IPA 5 min → IPA 3 min → N2, con agitazione a 500 rpm. Concentrazioni, quantità e preparazione esatta dei bagni non vanno inventate.

**Conseguenza:** ritirare HCl 48 h come nuova proposta per questo caso; non riproporre gli stessi solventi come se fossero inediti. L'agitazione rende inappropriato descrivere i bagni come statici. Non misura però direttamente la velocità locale vicino alla superficie o la concentrazione polimerica del liquido trascinato.

## 1. Ricambio del solvente: risultato pertinente con PMMA 950k

Olubunmi O. Ayodele, Sajedeh Pourianejad, Anthony Trofe, Aleksandrs Prokofjevs, Tetyana Ignatova, **Application of Soxhlet Extractor for Ultra-clean Graphene Transfer**, ACS Omega 7 (2022), 7297–7303. DOI: https://doi.org/10.1021/acsomega.1c07113
Testo integrale: https://pmc.ncbi.nlm.nih.gov/articles/PMC8892648/

PMMA 950k al 4%, distacco elettrochimico del grafene da Cu, deposizione su SiO2/Si. Acetone distillato riciclato tramite Soxhlet per 4 h, poi N2, contro immersione in acetone a 25 °C per 2 h e risciacquo IPA. Rq finale dichiarato 1,26 nm; microscopia e Raman migliori. Non usa XPS o prova atomica di assenza di residui. Il confronto confonde rinnovo, temperatura, durata, manipolazione e risciacquo. Storia termica molto più intensa del caso Armando: anche 150 °C per 2 h e 135 °C per 30 min. Sostiene il ricambio come variabile sperimentale; non dimostra pulizia universale né recupero di residui già resistenti a molti solventi.

## 2. Controevidenza nuova: ricambio e omissione IPA non bastano alla pulizia chimica

Zian Tang et al., **A Closed-Loop Solvent Recycling Device for Polymer Removal in Graphene Transfer Process**, Separations 12 (2025), 295, pubblicato 26 ottobre 2025. DOI: https://doi.org/10.3390/separations12110295
Testo integrale indicizzato: https://www.mdpi.com/2297-8739/12/11/295

Confronta cinque cicli di pulizia con distillazione/ricambio e immersioni convenzionali. PMMA 495k, acetone circa 38 °C nella camera campione, immersione continua fra cicli, successiva asciugatura sotto vuoto a 180 °C. Miglioramento ottico dei residui grandi, ma XPS quasi identico: rapporti C=C:C–O:O–C=O 100:28:8 (immersione) e 100:30:7 (ciclico). Entrambi lasciano residui. Nel controllo c'è IPA finale, nel ciclico no: è una controevidenza alla promessa “eliminare IPA risolve tutto”. La caratterizzazione non coincide con la storia PMMA/PPC di Armando. Anche le prove di purezza del solvente sono svolte in cicli senza campioni, quindi non dimostrano da sole assenza di ogni contaminante durante il processo reale. Il paper contiene incongruenze redazionali/quantitative nelle stime energetiche; qui non si utilizzano tali stime.

## 3. Miscela solvente/non solvente: base sperimentale, non diagramma del nostro sistema

J. Manjkow, J. S. Papanu, D. S. Soong, D. W. Hess, A. T. Bell, **An in situ study of dissolution and swelling behavior of poly-(methyl methacrylate) thin films in solvent/nonsolvent binary mixtures**, J. Appl. Phys. 62 (1987), 682–688. DOI verificato mediante Crossref: https://doi.org/10.1063/1.339742
Abstract originale riprodotto: https://cir.nii.ac.jp/crid/1360861290933469952
PDF editore identificato: https://pubs.aip.org/aip/jap/article-pdf/62/2/682/18611722/682_1_online.pdf

Ellissometria in situ su PMMA 1,2 micrometri: transizione stretta fra dissoluzione e insolubilità nelle miscele MEK/IPA e MIBK/metanolo vicino a 50:50. In regime insolubile il film può gonfiarsi fino a tre volte. A 50:50 MEK/IPA, raffreddare da 24,8 a 18,4 °C rende il film soltanto parzialmente solubile. Pertinenza: la traiettoria di composizione può contare tanto quanto il nome del solvente. **Limite:** non sono DMF/IPA, THF/IPA, acetone/IPA né PPC; non esportare soglie numeriche o diagrammi di fase al nostro caso. Letto abstract primario, non testo integrale.

## 4. IPA non è universalmente e assolutamente un non solvente

M. J. Rooks et al., **Low stress development of poly(methylmethacrylate) for high aspect ratio structures**, JVST B 20 (2002), 2937–2941. DOI: https://doi.org/10.1116/1.1524971
PDF primario letto: https://nano.yale.edu/sites/default/files/files/pmma_develop_ipa_water.pdf

IPA è un solvente molto debole per PMMA a basso peso molecolare; la miscela IPA/acqua può sviluppare PMMA esposto, benché i componenti singoli siano poco efficaci. Studio su resist irradiato e film spessi, non pulizia di PMMA intatto adsorbito su grafene. Serve a evitare l'affermazione assoluta che qualsiasi residuo precipiti inevitabilmente durante un risciacquo IPA. Per AR-P 672.045/PPC la precipitazione deve essere verificata con l'eluato reale e i polimeri di riferimento.

## 5. Ottimizzazione 2024: più calore non garantisce meno residui

Ahmed F. Abdelaal et al., **Polymethylmethacrylate (PMMA) deposition and removal optimization in CVD-grown graphene transfer: A Taguchi technique study**, Diamond and Related Materials 149 (2024), 111660. DOI verificato: https://doi.org/10.1016/j.diamond.2024.111660
Pagina autori: https://pure.kfupm.edu.sa/en/publications/polymethylmethacrylate-pmma-deposition-and-removal-optimization-i/
Pagina editore: https://www.sciencedirect.com/science/article/abs/pii/S0925963524008732

Condizioni selezionate per bassa copertura residua: PMMA 4,5%, 3000 rpm, acetone a 40 °C, 60 min; residuo ancora 3,42% dell'area. Un altro campione con più residui aveva proprietà elettriche migliori. I passaggi accessibili dell'editore descrivono maggiori particelle non disciolte nei campioni a 60 °C. Non è una prova causale isolata dell'effetto termico perché variano altri fattori; è comunque contro un'estrapolazione monotona ingenua. Non è disponibile qui il testo integrale: non inventare tipo/massa molecolare PMMA o tutti i dettagli del substrato.

## 6. Formammide: esempio preciso di miglioramento elettrico senza rimozione

Ji Won Suk et al., **Enhancement of the Electrical Properties of Graphene Grown by Chemical Vapor Deposition via Controlling the Effects of Polymer Residue**, Nano Letters 13 (2013), 1462–1467. DOI: https://doi.org/10.1021/nl304420b
PDF autori letto: https://utw10193.utweb.utexas.edu/Archive/RuoffsPDFs/340.pdf

CVD su SiO2/Si; esposizione notturna a formammide migliora mobilità e sposta Dirac verso zero. Gli autori attribuiscono l'effetto a donazione elettronica della formammide assorbita nel residuo. La morfologia AFM del residuo rimane simile. Questo lavoro **non** è prova che un altro solvente rimuova PMMA; mostra perché mobilità/Dirac da soli non misurano la pulizia. Non proporre formammide come soluzione sulla base di quei numeri.

## 7. Vapore di acetone: esiste, ma evidenza distante dal recupero del caso

Sakib Ishraq, **Cleaning and Characterization of Chemical Vapor Deposited Graphene for Nanoelectronic Device Development**, tesi MS, West Virginia University (2024). DOI: https://doi.org/10.33915/etd.12515
https://researchrepository.wvu.edu/etd/12515/

Abstract riferisce rimozione con vapore di acetone e caratterizzazione AFM/Raman su grafene trasferito a vetro. Collegato protocollo JoVE degli autori: https://www.jove.com/v/63393/development-and-functionalization-of-electrolyte-gated-graphene-field-effect-transistor-for-biomarker-detection . Il protocollo indica esposizione a vapore per 4 min, seguita da acetone per 5 min. Il valore 70 °C riportato è una regolazione del riscaldamento, non prova della temperatura effettiva di vapore o campione. Non ho letto la tesi integrale. Nessuna prova comparativa recuperata su residui PMMA/PPC già trattati con DMF/THF e remover caldi. Da non promuovere a ricetta risolutiva.

## 8. Solubilità del PMMA massivo non basta a predire pulizia su grafene

Qingling Hang, Davide A. Hill, Gary H. Bernstein, **Efficient removers for poly(methylmethacrylate)**, JVST B 21 (2003), 91–97. DOI: https://doi.org/10.1116/1.1532734
Abstract primario disponibile sulla pagina della pubblicazione autore: https://www.researchgate.net/publication/260309520_Efficient_removers_for_polymethylmethacrylate

Confronto di solventi su PMMA/SiO2: 1,2-dicloroetano riporta rugosità vicina al substrato di partenza, meglio di acetone. Impiega interazione Flory–Huggins e energia dell'interfaccia solvente/superficie come descrittori. **Non è grafene.** Non giustifica estrapolare un parametro di solubilità a un residuo ignoto o promettere DCE dopo il fallimento di solventi già efficaci per PMMA massivo. Utile soltanto come fondamento per separare dissoluzione del polimero e distacco dall'interfaccia.

## 9. Particelle a base di silicio possono coesistere con PMMA su CVD-grafene

**CF4/H2 Plasma Cleaning of Graphene Regenerates Electronic Properties of the Pristine Material**, ACS Applied Nano Materials (2019). DOI: https://doi.org/10.1021/acsanm.8b02249
https://pubs.acs.org/doi/10.1021/acsanm.8b02249

Abstract primario descrive nanoparticelle a base di Si insieme a PMMA; H2 plasma da solo frammenta il grafene e attacca il SiO2 esposto. La rilevanza è diagnostica, non un suggerimento di usare plasma: non tutte le particelle osservate dopo trasferimento sono necessariamente PMMA. La sola presenza di Si2p su un substrato SiO2 non distingue nanoparticelle, supporto e PDMS. Serve discriminazione spaziale/chimica con riferimento appropriato.

## 10. Identificazione chimica locale dei residui di stampo

Jeffrey J. Schwartz et al., **Chemical Identification of Interlayer Contaminants within van der Waals Heterostructures**, ACS Applied Materials & Interfaces 11 (2019), 25578–25585. DOI: https://doi.org/10.1021/acsami.9b06594
Testo integrale: https://pmc.ncbi.nlm.nih.gov/articles/PMC6903401/
Autori NIST: https://www.nist.gov/publications/chemical-identification-interlayer-contaminants-within-van-der-waals-heterostructures

PTIR/AFM-IR identifica chimicamente PDMS e policarbonato in aggregati intrappolati dentro eterostrutture WS2/WSe2/hBN. Non si basa soltanto sulla forma AFM. Il policarbonato di questo paper è PC, **non PPC**. La pulizia preventiva degli stampi PDMS riduce contaminazione; non dimostra rimozione postuma dal campione di Armando. Supporta la scelta di confrontare spettri locali con materiali reali PMMA/PPC/PDMS e controlli del processo. La presenza della cornice PDMS nel metodo di letteratura non prova una contaminazione PDMS nel campione reale.

## 11. Residui PDMS e difficoltà dei solventi: prova su altro materiale 2D

Achint Jain et al., **Minimizing residues and strain in 2D materials transferred from PDMS**, Nanotechnology 29 (2018), 265203. DOI: https://doi.org/10.1088/1361-6528/aabd90
Preprint autore: https://arxiv.org/abs/1801.02971
PDF editore in repository: https://www.research-collection.ethz.ch/bitstreams/b7accb1c-4375-4be6-93f7-13c230a5bdfe/download

Dimostra residui significativi su MoS2 trasferito con PDMS, anche in zone apparentemente pulite a bassa risoluzione. Discute limiti di solventi come toluene e diclorometano e della sola ricottura. Metodi risolutivi esplorati includono UV/ozono preventivo e calore, fuori dal vincolo del caso. Pertinente per diagnosi differenziale e limiti della topografia; analogia, non replica sul nostro grafene.

## Deduzioni operative originali, condizionali

1. La debole risposta alla sequenza DMF/IPA/THF non identifica il componente attivo: DMF, THF, tempo complessivo, agitazione, passaggi IPA e asciugatura variano insieme. È necessario distinguere fra rimozione e redistribuzione/rigonfiamento.
2. La prova immediata di maggior valore è un confronto fattoriale che conserva tutto il procedimento già efficace in misura scarsa e cambia separatamente il passaggio IPA intermedio e quello finale. Nel ramo senza IPA intermedio, sostituire con un solvente realmente capace di mantenere in soluzione entrambi i polimeri, verificato con materiali reali; non inventare un rapporto DMF/THF ottimale.
3. Prima di consumare molti campioni, usare eluato reale e controlli PMMA, PPC e miscela preparati con la storia termica del processo. Aggiungere IPA fuori dal campione e misurare torbidità/particelle e loro reversibilità al solvente iniziale. Precipitazione positiva supporta un meccanismo; un esito negativo non esclude fenomeni localizzati nell'ultimo film liquido.
4. Confrontare chimicamente residui del campione e residui raccolti, preservando il campione dagli effetti della misura: SEM/elettroni o scansioni AFM aggressive possono modificare polimeri o spostarli. Non attribuire tutto a reticolazione senza sapere esposizioni e dose.
5. Se il materiale rimosso ha firma PMMA/PPC e il residuo persistente ha altra firma, insistere con solventi PMMA è strategicamente sbagliato. Se anche il controllo polimerico sul supporto inerte non si dissolve, cercare alterazioni del materiale, lotti, contaminanti e storia del processo; non invocare subito legame eccezionale col grafene.
6. Tenere HCl nel registro delle prove fallite. Né la letteratura HCl né un calcolo di solubilità possono annullare quel dato. La ricerca non ha identificato una soluzione chimica già dimostrata sul caso PMMA/PPC di Armando.

## Lacune esplicite

Non recuperata una fonte primaria diretta che dimostri precipitazione specifica di PPC su grafene durante IPA nella sequenza esatta DMF→IPA→THF. La purificazione PPC con antisolventi è citata da review e brevetti, insufficienti per trasformarla in causa certa del caso. Non recuperata una replica indipendente positiva o negativa esattamente di Xiao/HCl 1 M/48 h su AR-P 672.045/PPC. Non trovata nuova chimica 2024–2026 validata che superi tutti i vincoli e tutti i fallimenti riportati. Queste lacune vanno preservate nella conclusione.
