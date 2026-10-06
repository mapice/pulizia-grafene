# Cambio di obiettivo: recupero per successivo prelievo e impilamento

Dati umani: particelle presenti subito dopo trasferimento, prima della litografia. Il “90% PMMA” è una stima di identità, non copertura misurata. La membrana PC e la cupola PDMS sono previste per il passaggio futuro: non possono spiegare retroattivamente i residui già presenti. PC, PPC e PMMA sono tre materiali distinti.

## Due fonti primarie pertinenti, con limiti

**Pizzocchero et al., The hot pick-up technique for batch assembly of van der Waals heterostructures**, Nature Communications 7 (2016), 11894. DOI https://doi.org/10.1038/ncomms11894 ; PDF completo https://arxiv.org/pdf/1605.02334 . Su grafene ESFOLIATO da SiO2, una struttura PPC/PDMS trasporta hBN che preleva il grafene. A 110 °C viene riferita resa di impilamento prossima al 100%; i contatti hanno resa distinta, 88%, in 22 dispositivi. Il contatto lento (<1 micrometro/s) influenza l'intrappolamento. Gli autori dichiarano che PPC da solo non preleva il grafene come fa hBN: la superficie che effettivamente tocca il cristallo è cruciale. Questi numeri non validano cristalli CVD trasferiti umidi, sporchi, né un prelievo diretto con PC. Non esportare temperatura o velocità come prescrizione universale.

**Purdie et al., Cleaning interfaces in layered materials heterostructures**, Nature Communications 9 (2018), 5387. DOI https://doi.org/10.1038/s41467-018-07558-3 ; PDF https://www-g.eng.cam.ac.uk/nms/publications/pdf/Purdie2018.pdf . Stampo PC/PDMS e hBN per grafene ESFOLIATO. Lo studio include contaminazione intenzionale con PMMA 495k all'8% e successiva rimozione in acetone/IPA prima dell'incapsulamento. È possibile ottenere elevate prestazioni finali anche dopo questa storia; tuttavia le bolle dei campioni contaminati restano immobili a 180 °C e richiedono 250 °C nella pulizia tramite avanzamento del contatto. Questo è un trattamento fisico/termico dell'interfaccia, non pulizia chimica e non dimostrazione di resa elevata su CVD/PMMA/PPC. Contaminazione iniziale e impossibilità di ottenere uno stack utile non sono equivalenti.

## Metrica coerente con lo scopo di Armando — proposta

Definire prima del confronto un insieme di cristalli idonei per dimensione, integrità e geometria. Non selezionare soltanto quelli riusciti dopo la pulizia. L'esito primario è:

Y_utile = numero di stack che soddisfano i requisiti finali / numero di cristalli idonei assegnati al trattamento.

Registrare separatamente: sopravvivenza alla pulizia; cristalli prelevati integri; area recuperata rispetto all'area iniziale; rilascio e allineamento riusciti; area finale priva di bolle/particelle nei limiti della misura; eventuale prestazione elettrica richiesta. I requisiti vanno fissati dall'applicazione, senza inventare una soglia percentuale universale.

Il confronto chimico deve mantenere uguali membrana PC, geometria della cupola, superficie effettivamente a contatto con grafene, velocità, temperatura e tempi del prelievo. La pulizia vincente deve aumentare Y_utile, non soltanto abbassare rugosità o conteggio delle particelle prima del prelievo. Un trattamento che rende la superficie più bella ma aumenta adesione al SiO2 o fragilità può peggiorare l'obiettivo. È un rischio da misurare, non un effetto già accertato nel caso.

Non usare questi due lavori per rispondere al vincolo chimico con “basta fare dry transfer”. Servono a scegliere la metrica finale e a evitare di scartare un campione soltanto perché non atomisticamente pulito prima del prelievo.
