# Verifica fonti sulla pulizia chimica del grafene — 6 ottobre 2026

Appunti bibliografici, non risultati sperimentali. Fonti primarie lette in questa sessione.

## Aggiornamento del caso fornito dall'utente

Dati nuovi: riscaldamento a 90 °C per 2 min; distacco elettrochimico del grafene da Cu in bagno acquoso con elettrodi, elettrolita non noto; ulteriore strato di PPC in anisolo sopra PMMA, senza contatto diretto intenzionale con grafene. Questo è coerente con il ramo principale semi-secco di Tyagi, NON autorizza ad attribuire al caso un etching con FeCl3 o APS. Neppure dimostra NaOH o esatta polarità del processo.

SI Tyagi letto tramite testo indicizzato: https://www.rsc.org/suppdata/d1/nr/d1nr05904a/d1nr05904a1.pdf
Riporta PMMA 100 nm e 90 °C/2 min; PPC 1,5 micrometri e altro trattamento 90 °C/2 min; cornice PDMS; delaminazione in NaOH 1 M, circa 2,4 V con controelettrodo Pt e circa 3 mA, due risciacqui DI, asciugatura, laminazione sul substrato a 90 °C. Il SI chiama Cu/SLG anodo ma parla anche di limitare bolle di H2: non estrapolare polarità e meccanismo del caso reale da questo solo testo. Nel paper PPC è 15 volte più spesso del PMMA; il dato utente "molto fine" non permette assumere questi spessori.

Implicazioni da segnalare come inferenze: PPC introduce un'altra possibile origine dei residui dopo dissoluzione, anche se inizialmente non tocca grafene. I 90 °C/2 min, da soli, non documentano carbonizzazione del PMMA. Servono riferimenti analitici del PMMA, del PPC e di entrambi sottoposti alla stessa storia; XPS C1s da solo può non separarli agevolmente dato che entrambi contengono carbonio e ossigeno.

## Prodotto indicato da Armando

AR-P 672.045 Allresist: PMMA nominale 950K, solidi 4,5%, solvente anisolo. Non confondere con il copolimero PMMA/MA AR-P 617. Il nome del prodotto non dimostra che le particelle residue siano PMMA non modificato.

Fonti del produttore:
- https://www.allresist.de/portfolio-item/e-beam-resist-ar-p-672-serie/
- https://www.allresist.com/wp-content/uploads/sites/2/2024/12/Allresist_Product-information-E-Beamresist-AR-P-630-670-English-web.pdf
- https://www.allresist.com/wp-content/uploads/sites/2/2021/05/Allresist_Product-information-E-beamresist-English-web.pdf
- SDS: https://www.allresist.com/wp-content/uploads/sites/2/2021/12/AR-P672-series_SDB_GB_1.1_1496.pdf

## Paper probabilmente richiamato da Armando

Ayush Tyagi et al., “Ultra-clean high-mobility graphene on technologically relevant substrates”, Nanoscale 14 (2022), 2167–2176. DOI https://doi.org/10.1039/D1NR05904A
PDF primario: https://pubs.rsc.org/en/content/articlepdf/2022/nr/d1nr05904a

Autori comprendono Vaidotas Mišeikis e Camilla Coletti. Campioni principali: monocristalli 200–250 micrometri, CVD su Cu 25 micrometri, trasferimento semi-secco elettrochimico su SiO2 285 nm/Si; prova supplementare anche grafene policristallino trasferito umido. 1SC: acetone 2h, IPA 5min, N2. 2SC aggiunge AR600-71 3min, acqua DI10s,N2. Remover dichiarato 70%1,3-diossolano +30%1-metossi-2-propanolo. Non attribuire automaticamente queste percentuali al prodotto commerciale attuale senza SDS attuale.

AFM: particelle da810 a34 su100 micrometri quadrati (>95% riduzione), non zero. Per campione umido Rq2,2→0,7nm. L'autore non dimostra rimozione assoluta di ogni molecola. Il corpo principale è semisecco, ma non è corretto dire che manchino del tutto dati sul trasferimento umido.

## Candidato post-trasferimento veramente pertinente

Zhuocong Xiao, Qifang Wan, Colm Durkan, “Cleaning Transferred Graphene for Optimization of Device Performance”, Advanced Materials Interfaces 6 (2019), 1801794. DOI https://doi.org/10.1002/admi.201801794
Repository autore: https://www.repository.cam.ac.uk/handle/1810/291256
Pagina completa autore usata per lettura: https://www.researchgate.net/publication/332784821_Cleaning_Transferred_Graphene_for_Optimization_of_Device_Performance

PMMA495KA4(4%in anisolo), CVD/Cu→SiO2 300nm/Si. Dopo acetone3h+IPA, bagnoHCl1M2giorni, risciacquoDI,N2. Temperatura del bagno non esplicitata nella sezione sperimentale letta: non attribuire25°C come dato pubblicato. AFM Ra2,4nm iniziale, ~1nm con solventi organici,0,4nm conHCl/NaCl. Raman nessun incrementoD, ma doping aumenta. Valori elettrici provengono da pulizia dopo sviluppo e prima metallizzazione, non dallo stesso confronto AFM posttrasferimento: non mescolare i due esperimenti.

Limite più importante: sul lotto invecchiato3mesi e ricotto200°C3h, copertura residui75%→34% conHCl1M2giorni (DI70%,NaCl40%). Gli autori riferiscono rimozione completa dopoHCl concentrato6h senza specificare una molarità riproducibile in quel passaggio: NON usarlo come ricetta alternativa pronta. Trasferibilità da495K a950K e dalla loro storia termica a quella di Armando non dimostrata. Idrolisi acida delle catene laterali è interpretazione, non provata con analisi dei prodotti. Ipotesi degli autori di clorurazione daNaCl non sufficientemente fondata e da non ripetere come certezza. NaCl può aggiungere contaminazione mobile/doping suSiO2.

## Fonte 2025: prevenzione distinta dal recupero

Hao Liu, Yingzhi Li, Kun Yang, Lei Guo, Zebing Zeng, Yifan Yao, “Revealing key surface contaminants via stack-flipping strategy: investigating rinsing protocols for clean graphene transfer and enhanced electrical performance”, Journal of Materials Chemistry C13 (2025),19606–19614. DOI https://doi.org/10.1039/D5TC02259B
Testo integrale indicizzato: https://pubs.rsc.org/be/content/articlehtml/2025/tc/d5tc02259b?page=search

Studia risciacqui del film galleggiante PMMA/grafene DOPO dissoluzioneCu e PRIMA deposizioneSiO2. Le etichette abbreviateHCl,HNO3,H2SO4 indicano misceleACIDO/H2O2/H2O, non acidi soli. Tempo2h; concentrazioni e dettagli nelSI non recuperato in questa sessione. Non usare come protocollo riproducibile senzaSI. Testo menziona anche ammoniaca residua dall'ultimo risciacquo: non promettere processo integralmente privo di basi senza verificareSI.

Con capovolgimento caratterizza lato inferiore esposto. Residui carboniosi probabilmente amorfi, non soltantoPMMA. Cu basso simile eFe non rilevato anche col controlloDI dopo2h, quindi differenze non spiegate da soli metalli. FTIRPMMA simile fra trattamenti: integrità macroscopica del polimero mantenuta. Dati sostengono prevenzione contaminazione e miglior successiva rimozionePMMA, non recupero di campione già suSiO2, nonAR-P672.045 specifico.

## Limiti delle fonti e della conclusione

La fonte più direttamente utilizzabile èXiao2019 come candidato controllato da adattare, non risultato nuovo. Non emerge una garanzia universale di pulizia completa mediante solvente/acido per campione ignoto. La conclusione corretta richiede caratterizzazione chimica delle particelle, verifica che siano sul lato accessibile e confronto controllato su campioni gemelli con stessa storia.
