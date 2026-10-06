# Terzo giro — posizione e origine della contaminazione

Ricerca del 6 ottobre 2026. Scopo: capire se gli oggetti resistenti ai solventi siano sopra il grafene, all'interfaccia con SiO2 o già presenti dopo crescita/distacco. Non si assume contaminazione PDMS. I trattamenti falliti riferiti da Armando, compreso HCl per 48 h, rimangono dati e non vengono riproposti come nuove soluzioni.

## Risultato principale

La letteratura primaria documenta almeno tre motivi per cui una protrusione AFM dopo trasferimento non identifica PMMA accessibile dall'alto: sacche di acqua sotto il foglio; contaminanti carboniosi già presenti sul grafene/Cu; metalli residui dopo delaminazione elettrochimica. Esistono metodi per distinguere questi casi, ma non è stata trovata una nuova pulizia con soli solventi/acidi che garantisca l'estrazione di particelle intrappolate sotto un foglio intatto e aderente a SiO2.

## Fonti primarie nuove e quanto permettono di dire

### 1. Magnozzi et al.: rilievi AFM che erano acqua intrappolata

M. Magnozzi, N. Haghighian, V. Mišeikis, O. Cavalleri, C. Coletti, F. Bisio, M. Canepa, **Fast detection of water nanopockets underneath wet-transferred graphene**, Carbon 118 (2017), 208–214. DOI: https://doi.org/10.1016/j.carbon.2017.03.022
PDF completo letto: https://arxiv.org/pdf/1805.01175

CVD/Cu, PMMA, FeCl3, SiO2 circa 290 nm. Una parte dei campioni mostra due livelli AFM distanti circa 2–2,5 nm. XPS ed ellissometria sostengono l'interpretazione con acqua interposta. Dopo trattamento termico il segnale graphitico C1s resta sostanzialmente costante, mentre quello del substrato aumenta: argomento per materiale rimosso sotto il grafene. **Limiti:** etching diverso dal caso; il trattamento chiamato “mild” arriva a 350 °C in circa 10^-7 mbar, non 90 °C. Lo spessore ottenuto da ellissometria (~5 nm) differisce dall'AFM (~2 nm): modello non perfettamente identificato. Non suggerire la ricottura come soluzione nel vincolo chimico.

### 2. Lee, Ahn, Ryu: accesso ai bordi e movimento dell'acqua

D. Lee, G. Ahn, S. Ryu, **Two-Dimensional Water Diffusion at a Graphene–Silica Interface**, JACS 136 (2014), 6634–6642. DOI: https://doi.org/10.1021/ja4121988
PDF autore letto: https://arxiv.org/pdf/1404.6733

Raman in tempo reale segue acqua che entra dal bordo e avanza sotto il grafene; AFM osserva strati di circa 3,5 Å. Cinetica dipendente dallo stato della silice. Rilevanza: l'agitazione nel bagno non equivale a ricambio del liquido confinato. **Limite:** acqua, grafene e substrati con trattamento specifico; non fornisce coefficienti per DMF/THF o tempi di estrazione dei residui polimerici del caso. Variazioni del Raman possono misurare intercalazione/drogaggio, non rimozione del polimero.

### 3. Lupina et al.: anche il distacco elettrochimico lascia metalli

G. Lupina et al., **Residual Metallic Contamination of Transferred Chemical Vapor Deposited Graphene**, ACS Nano (2015). DOI: https://doi.org/10.1021/acsnano.5b01261
PDF completo con supplemento letto: https://arxiv.org/pdf/1505.00889

Confronta trasferimenti chimici ed elettrochimici su Si/SiO2. ToF-SIMS/TXRF rilevano Cu e Fe, anche quando XPS non vede metalli. Nessuno dei procedimenti studiati porta il Cu sotto 10^13 atomi/cm². Alcuni campioni elettrochimici hanno aggregazioni di Cu molto maggiori; queste aggregazioni scompaiono dopo trattamenti con HCl, ma il fondo residuo persiste. Le zone vicino ai bordi possono essere contaminate dalla lavorazione e non sono rappresentative. **Limite:** metallo in tracce non implica che tutte le particelle otticamente visibili siano metalliche; HCl fallito non discrimina da solo metallo, accessibilità e composizione.

### 4. Amontree et al. 2024: contaminazione già durante la crescita

J. Amontree et al., **Reproducible graphene synthesis by oxygen-free chemical vapour deposition**, Nature 630 (2024), 636–642. DOI: https://doi.org/10.1038/s41586-024-07454-5
Fonte editore: https://www.nature.com/articles/s41586-024-07454-5
PDF autore/NIST: https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=936656

Esperimenti con quantità controllate di ossigeno mostrano variazioni di crescita, contaminazione superficiale, Raman e trasporto. AFM/XPS confrontano grafene ancora sul Cu: la contaminazione non richiede contatto con PMMA. La qualità migliore deriva dalla crescita senza ossigeno, non dalla pulizia successiva. **Deduzione per il caso:** conservare una porzione dello stesso grafene/Cu prima di depositare polimeri è un controllo ad alto valore. Non attribuire però l'impurità di Armando a ossigeno o crescita senza caratterizzarla.

### 5. Lin et al.: isotopi per distinguere contaminazione di crescita e polimero

**Towards super-clean graphene**, Nature Communications 10 (2019), 1912. DOI: https://doi.org/10.1038/s41467-019-09565-4
https://pmc.ncbi.nlm.nih.gov/articles/PMC6478734/

Lo studio usa marcatura 12C/13C per indagare l'origine del carbonio amorfo e PMMA deuterato con ToF-SIMS per rilevare residui di trasferimento. Grafene inizialmente più pulito lascia meno PMMA dopo un trasferimento identico. Sono quindi possibili contributi accoppiati: sporco precedente può favorire residui successivi. **Limite:** cambiare la crescita o usare PMMA marcato non recupera retroattivamente il campione già trattato. Una contaminazione carboniosa senza firma PMMA non va automaticamente chiamata “PMMA carbonizzato”.

### 6. Zemek et al.: profondità senza asportare materiale

J. Zemek, J. Houdková, P. Jiříček, T. Ižák, M. Kalbáč, **Non-destructive depth profile reconstruction of single-layer graphene using angle-resolved X-ray photoelectron spectroscopy**, Applied Surface Science 491 (2019), 16–23. DOI: https://doi.org/10.1016/j.apsusc.2019.06.083
Pagina primaria letta: https://www.sciencedirect.com/science/article/pii/S0169433219317921

ARXPS di C1s, O1s e Cu3p ricostruisce contaminazione carboniosa/ossigenata sopra il grafene e ossigeno associato al rame all'interfaccia inferiore. Il sistema è grafene/Cu, non PMMA/PPC/grafene/SiO2. L'applicazione proposta al caso richiede confrontare modelli di strato superiore, inferiore e misto, vincolati da copertura/rugosità AFM. Un semplice rapporto misurato a due angoli non è una localizzazione certa: isole e superfici non uniformi possono produrre tendenze ambigue. Non sostituire la ricostruzione con una sola deconvoluzione C1s.

### 7. Holroyd et al.: Raman “pulito” e residuo polimerico

C. Holroyd, A. Horn, C. Casiraghi, S. Koehler, **Vibrational fingerprints of residual polymer on transferred CVD-graphene**, Carbon 117 (2017), 473–475. DOI: https://doi.org/10.1016/j.carbon.2017.03.008
PDF completo letto: https://serwiss.bib.hs-hannover.de/files/2260/holroyd_etal2017-fingerprints_cvd-graphene.pdf

Grafene CVD trasferito su oro. RAIRS e spettroscopia di somma di frequenze identificano firme PMMA che il Raman ordinario non rivela. Non dimostra che il Raman non possa mai trovare PMMA: dimostra che l'assenza dei suoi segnali non basta a escluderlo in un residuo sottile. Serve uno spettro completo e confronto con il materiale reale; una singola banda C–H non è specifica per PMMA né distingue automaticamente PPC. Trasferibilità strumentale a SiO2 e sensibilità locale da qualificare.

### 8. Wang et al.: tracciare il polimero, non soltanto la forma

X. Wang et al., **Direct Observation of Poly(Methyl Methacrylate) Removal from a Graphene Surface**, Chemistry of Materials 29 (2017), 2033–2039. DOI: https://doi.org/10.1021/acs.chemmater.6b03875
https://pubs.acs.org/doi/10.1021/acs.chemmater.6b03875

La marcatura isotopica PMMA con ToF-SIMS identifica, localizza e quantifica il residuo; la ricottura sotto vuoto può trasformarlo in carbonio amorfo senza eliminarlo. Letto abstract primario e descrizione supplemento, non testo principale completo. Il messaggio pertinente è analitico: scomparsa dell'esterico non equivale a scomparsa di tutto il carbonio estraneo. Non trasferire quella trasformazione termica alla breve cottura di Armando a 90 °C.

### 9. Matruglio et al.: “non ha toccato inizialmente il polimero” non significa immune

A. Matruglio et al., **Contamination-free suspended graphene structures by a Ti-based transfer method**, Carbon 103 (2016), 305–310. DOI: https://doi.org/10.1016/j.carbon.2016.03.023
https://www.sciencedirect.com/science/article/abs/pii/S0008622316302123

Il lavoro segnala contaminazione della faccia posteriore del grafene sospeso da parte di polimero già dissolto: proteggere soltanto la faccia che inizialmente tocca PMMA non protegge l'altra durante i bagni. **Limite:** membrane sospese e altra strategia di trasferimento; non prova che PMMA/PPC sia penetrato sotto il campione aderente di Armando. Smentisce soltanto l'inferenza automatica “PPC era sopra PMMA, quindi non può contaminare grafene”. Non proporre il Ti come pulizia di recupero.

### 10. Jang et al.: rimuovere l'interfaccia invece di attaccare il substrato

D. J. Jang et al., **A Modified Wet Transfer Method for Eliminating Interfacial Impurities in Graphene**, Nanomaterials 13 (2023), 1494. DOI: https://doi.org/10.3390/nano13091494
https://pmc.ncbi.nlm.nih.gov/articles/PMC10179892/

PMMA/grafene è appoggiato a SiO2 temporaneo idrofilo; acqua lo separa nuovamente e parte delle impurità rimane sul supporto. Dimostra che il trasferimento può separare contaminanti dall'interfaccia senza dissolvere l'ossido. **Ma** il supporto è pretrattato UV e il procedimento è preventivo, prima della rimozione PMMA: non soddisfa automaticamente il vincolo del caso né dimostra recupero del campione ormai processato. Modificarlo con preparazione chimica dell'ossido sarebbe una nuova proposta da provare, non una ricetta già validata. Non si raccomanda qui questo procedimento come soluzione.

## Un esperimento discriminante concreto, senza nuova famiglia di solventi

### Confronto per stadio con faccia inferiore esposta

Usare tre piccole porzioni adiacenti dello stesso grafene/Cu, conservando identiche condizioni fino al punto in cui si dividono:

1. **Prima del polimero:** conservare una porzione sul Cu per verificare se la classe di contaminazione esiste già dopo crescita. Registrare zone e metodi; non confrontare soltanto rugosità assoluta Cu contro SiO2.
2. **Dopo il distacco e i risciacqui, prima della deposizione normale:** depositare una porzione capovolta su supporto diagnostico, lasciando la faccia precedentemente sul Cu esposta e il PMMA/PPC fra grafene e supporto. Lo scopo è osservare la faccia inferiore, non dissolvere da questa geometria il supporto polimerico. Applicare solo le misure che non confondano inevitabilmente il PMMA sottostante col contaminante superiore: localizzazione delle particelle, loro spettro locale e controlli del supporto.
3. **Campione finale standard:** completare la terza porzione con la procedura di Armando, incluso il risciacquo realmente usato. Confrontare classe morfologica e firma chimica, non assumere che particelle simili siano lo stesso materiale.

Il capovolgimento è già documentato come idea da Liu 2025 (citato nella prima nota); il valore nuovo di questo disegno è combinarlo con un controllo prima dei polimeri, il distacco elettrochimico effettivo e la sequenza reale fallita. Non usa la pulizia con acidi/perossido di quel lavoro. La sua esecuzione deve controllare introduzione di impurità dal supporto e dalla manipolazione.

**Interpretazioni ammissibili:**

- Contaminazione con firma corrispondente già nella porzione sul Cu: origine almeno parzialmente precedente al trasferimento; non dimostra che tutta la contaminazione finale provenga da lì.
- Particelle già sulla faccia ex-Cu dopo distacco ma non nel controllo iniziale osservabile: restringe l'origine a faccia inferiore, distacco/rinsciacquo o zona del Cu inizialmente inaccessibile; non assegna ancora una chimica.
- Solo il campione finale le mostra: diventano prioritarie dissoluzione, risciacqui, substrato e asciugatura; non è prova automatica di precipitazione IPA.
- La faccia inferiore nuova è pulita: non esclude acqua o polimeri introdotti dopo la deposizione, né dimostra che i rilievi del vecchio campione siano sopra.

### Sul campione già preparato

La posizione richiede una prova aggiuntiva: ARXPS con modelli concorrenti e vincoli AFM, oppure una separazione diagnostica su campione sacrificabile con registrazione di ciò che rimane sul substrato e sulla membrana. Non chiamare una semplice immagine AFM o un segnale Si2p “prova del lato inferiore”. Anche una bolla interfaciale può spostarsi durante scansioni AFM; il solo movimento sotto la punta non identifica residuo superficiale.

## Perché il lato cambia il problema chimico

Un solvente che dissolve il polimero esposto può non estrarlo da un'interfaccia chiusa: occorrono accesso, solubilizzazione e una via di uscita dei prodotti. Il grafene ideale è una barriera molecolare (Bunch et al., Nano Lett. 8 (2008), 2458–2462, DOI https://doi.org/10.1021/nl801457b), ma i campioni CVD possono avere bordi, difetti e rotture. Non assumere né impermeabilità assoluta del campione né passaggio uniforme dei bagni attraverso il foglio.

L'eventuale andamento della risposta alla pulizia con distanza da bordi/aperture è un indizio di trasporto laterale, soltanto dopo aver misurato lo stato iniziale: Lupina mostra che i bordi possono partire più contaminati. Non estrarre una costante di diffusione da tempi arbitrari né esportare quella dell'acqua a PMMA/PPC.

## Decisione sull'alternativa chimica

- Se il residuo è sopra e ha firma polimerica, resta sensato il confronto della traiettoria dei solventi con IPA separato descritto nel secondo giro.
- Se il residuo è sotto, ripetere gli stessi solventi senza dimostrare accesso ha poco valore. La scelta concreta è esporre/prevenire la contaminazione prima della chiusura dell'interfaccia, oppure ammettere che il recupero del campione resta irrisolto entro i vincoli.
- Se la firma è minerale, metallica o carboniosa di crescita, nessun risultato qui autorizza a chiamarla PMMA e nessun trattamento già fallito va rimesso in lista senza una nuova discriminante.

Non è stata individuata una ricetta nuova con soli acidi/solventi organici che rimuova in modo dimostrato tutti questi possibili contaminanti preservando SiO2 e grafene nel campione descritto. Il risultato del terzo giro è rendere distinguibili scenari prima confusi, non dichiararne vero uno.
