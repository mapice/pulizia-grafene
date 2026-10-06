"""Append actual molecular/bulk results to the existing standalone article."""
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[2]
TARGET=ROOT/'outputs/pulizia-grafene-pmma-ppc.tex'
SOURCE=ROOT/'work/versione-2/pulizia-grafene-pmma-ppc.tex'
HERE=Path(__file__).resolve().parent
s=SOURCE.read_text()
s=s.replace('revisione 2','revisione 3').replace('Revisione 2: teoria cinetica, simulazioni condizionali e prove per la sequenza DMF-IPA-THF','Revisione 3: calcoli quantistici, modello molecolare appreso e simulazione del liquido HFIP')
s=s.replace(r'\title[Grafene: teoria della pulizia a solventi sequenziali]{Pulizia chimica del grafene con residui di PMMA/PPC:\\ teoria della sequenza dei solventi e prove discriminanti}',r'\title[Grafene: teoria e calcoli molecolari della pulizia]{Pulizia chimica del grafene gi\`a trasferito con PMMA/PPC:\\ teoria, calcoli molecolari e verifica dei modelli}')
s=s.replace(r'L\textquotesingle{}efficacia sul campione reale rimane una domanda sperimentale.','') if False else s
s=s.replace("L'efficacia sul campione reale rimane una domanda sperimentale.","Aggiungiamo calcoli molecolari quantistici e con un potenziale appreso, oltre a una simulazione atomistica del liquido HFIP. Il confronto con un acido gi\\`a fallito mostra che una forte interazione locale con l'estere non basta a scegliere il pulente. Il controllo del liquido impedisce di promuovere un modello dalla sola densit\\`a. L'efficacia sul campione reale rimane da verificare.")
intro=r'''
\paragraph{Aggiornamento molecolare.}
HFIP \`e la nuova candidata chimica discussa nella Sezione~\ref{sec:molecular}.
Esiste evidenza di dissoluzione del PMMA in HFIP a temperatura ambiente,
ma la sua efficacia sui residui aderenti del caso non \`e dimostrata.
I calcoli eseguiti hanno anzi scartato il criterio semplicistico
``interazione pi\`u forte con il carbonile, quindi pulizia migliore''.
La sequenza DMF--THF rimane il controllo operativo; una classifica
quantitativa fra DMF, THF e HFIP richiede ancora un modello interfaciale
qualificato, oltre alla verifica sul materiale.

'''
s=s.replace(r'\section{Il caso reale e i risultati negativi da conservare}',intro+r'\section{Il caso reale e i risultati negativi da conservare}')
s=s.replace("chimica delle particelle non sono ancora noti.",r'''chimica delle particelle non sono ancora noti.
Le particelle sono riferite subito dopo il trasferimento, prima della
litografia. La membrana PC e il supporto PDMS sono impiegati in un
successivo prelievo a secco: questo uso futuro non identifica l'origine
dei depositi gi\`a presenti. La stima ``90\% PMMA'' comunicata dall'esperto
\`e una valutazione di probabilit\`a, non una misura della composizione.''')

data=[]
for pair,label in [('HFIP__MethylAcetate','HFIP'),('AceticAcid__MethylAcetate','Acido acetico'),('IPA__MethylAcetate','IPA')]:
    for basis in ['def2-svp','def2-tzvp']:
        r=json.loads((HERE/f'results/dft/{pair}_pbe0_{basis}_interaction.json').read_text())
        data.append((label,basis,r['interaction_counterpoise_kJ_mol'],r['BSSE_correction_kJ_mol']))
table='\n'.join(f'{label} & {"SVP" if basis=="def2-svp" else "TZVP"} & {e:.2f} & {bsse:.2f} \\\\' for label,basis,e,bsse in data)
bulk=json.loads((HERE/'results/gromacs/bulk_scale0.5_seed20261006/summary.json').read_text())
epsilon=bulk['dielectric_GROMACS'];rho=bulk['density_g_cm3']
chapter=r'''
\section{Calcoli molecolari eseguiti e ipotesi scartate}\label{sec:molecular}

\subsection{Il bersaglio locale non determina il distacco}

Il PMMA presenta gruppi estere; PPC presenta gruppi carbonato. Abbiamo
calcolato interazioni di HFIP con acetato di metile, pivalato di metile
e carbonato di dimetile come analoghi locali, oltre a coppie di solventi.
Queste molecole non rappresentano la catena 950K, le sue molteplici
ancore o la storia del deposito. Sono controlli di un'ipotesi locale.

Sono state ottimizzate 12 molecole e 13 coppie, con tre geometrie iniziali
per coppia, usando GFN2-xTB. Tutte le 39 geometrie delle coppie sono state
raffinate fino a convergenza con forza massima sotto
0,008\,eV/\AA. Non rivendichiamo un minimo globale.
Le energie di interazione a geometria fissata sono
\begin{equation}
 E_{\rm int}=E_{AB}(\mathbf R_A,\mathbf R_B)
             -E_A(\mathbf R_A)-E_B(\mathbf R_B),
\end{equation}
separate dall'energia di legame che comprende la deformazione dei frammenti.

Abbiamo poi eseguito 50 calcoli autoconsistenti PBE0-D3(BJ), tutti
convergenti, per sette coppie con base def2-SVP e tre controlli con
def2-TZVP. La Tabella~\ref{tab:quantum} corregge l'errore dovuto
all'uso della base del partner tramite il metodo \emph{counterpoise}.
I frammenti rimangono nella geometria del complesso ottenuta con GFN2.
Sono energie elettroniche nel vuoto, non energie libere nel liquido.

\begin{table}[htbp]\centering
\caption{Interazione con acetato di metile a geometria fissata,
PBE0-D3(BJ); kJ/mol. La correzione di base \`e positiva e riduce
l'apparente attrazione. Le cifre esprimono il calcolo numerico,
non un'accuratezza fisica al centesimo.}\label{tab:quantum}
\small
\begin{tabular}{lcrr}\toprule
Partner & Base & $E_{\rm int}$ corretto & Correzione di base\\\midrule
TABLE_QUANTUM
\bottomrule\end{tabular}
\end{table}

\textbf{Il confronto produce una confutazione utile: anche l'acido
acetico interagisce almeno altrettanto fortemente con questo analogo,
pur essendo gi\`a fallito sul materiale riferito.} Per queste geometrie
la differenza fra SVP e TZVP dopo correzione \`e circa 1\,kJ/mol;
non dimostra per\`o convergenza rispetto a geometrie, funzionale,
entropia, liquido o superficie. Il solo legame con l'estere non
discrimina un trattamento efficace.

La letteratura documenta PMMA di massa 99.200\,g/mol disciolto in HFIP
a temperatura ambiente \cite{WangHFIP2020}. Dimostra la solvatazione
di quel materiale libero. La selezione di HFIP perch\'e dissolva anche
metilcellulosa non implica che gli altri solventi non dissolvano PMMA.
Non dimostra rimozione del PMMA 950K aderente al grafene o del PPC.

Un eventuale primo confronto operativo usa HFIP puro in due bagni
freschi di 10 minuti a temperatura ambiente, con il secondo destinato
a esportare il materiale estratto, inizialmente evitando la mescolanza
con IPA. \textbf{Quei tempi sono una proposta pilota, non condizioni
validate o derivate dai calcoli.} Occorrono un controllo del materiale
libero dello stesso lotto, del substrato e della contaminazione
fluorurata dopo il trattamento, insieme alle misure della Sezione~9.
La sequenza abituale conserva la funzione di controllo operativo.

\subsection{Perch\'e serve il liquido e una selettivit\`a interfaciale}

Siano $P$ il polimero, $G$ il grafene, $O$ l'ossido e $L$ il liquido.
Nel limite reversibile e non reattivo i lavori di separazione sono
\begin{align}
 W_{PG}^{L}&=\gamma_{PL}+\gamma_{GL}-\gamma_{PG},\\
 W_{GO}^{L}&=\gamma_{GL}+\gamma_{OL}-\gamma_{GO}.
\end{align}
La differenza cancella $\gamma_{GL}$: una maggiore affinit\`a del
solvente per il grafene non garantisce selettivit\`a fra rimozione
del polimero e distacco del foglio dall'ossido. Questi lavori non
sono, da soli, barriere cinetiche o energie di frattura del processo
di prelievo, che dipende anche da geometria, temperatura e velocit\`a.

Il medesimo problema appare a scala molecolare. Se il carbonile \`e
solvatato nello stesso modo sia nel polimero aderente sia nel polimero
liberato, quel contributo si cancella nel costo del distacco.
Bisogna calcolare la differenza fra i due stati, includendo liquido,
contatti, infiltrazione, conformazioni e riaggancio. Moltiplicare
un'energia di coppia per circa 9.500 unit\`a non fornisce una barriera
corretta: la catena pu\`o staccarsi progressivamente.

\subsection{Un controllo del liquido che limita le conclusioni}

Abbiamo trascritto i parametri HFIP pubblicati nelle Tabelle I--V
del supplemento di Marchelli et al. \cite{Marchelli2022}, convertendo
esplicitamente unit\`a e prefattori armonici. Il fattore 1--4, non
specificato nel testo consultato, \`e assunto pari a 0,5 e deve essere
oggetto di sensibilit\`a. Sono state simulate 125 molecole:
50\,ps a volume fissato e 1\,ns a pressione fissata a 298,15\,K.
L'analisi usa 200--1000\,ps; una sola realizzazione e questo intervallo
non bastano a certificare convergenza della risposta dielettrica.

\begin{table}[htbp]\centering
\caption{Controllo preliminare del modello HFIP. I riferimenti
sperimentali sono riportati da Casoria et al. \cite{Casoria2024}.
La deviazione della densit\`a \`e uno scarto temporale, non un
intervallo di confidenza di repliche indipendenti.}
\begin{tabular}{lcc}\toprule
Quantit\`a & Simulazione pilota & Riferimento sperimentale\\\midrule
Densit\`a, g/cm$^3$ & RHO_BULK ($s=0,022$) & 1,607\\
Permittivit\`a relativa & EPS_BULK & 16,7\\\bottomrule
\end{tabular}
\end{table}

L'espressione delle fluttuazioni del dipolo \`e stata verificata
indipendentemente in unit\`a SI:
\begin{equation}
 \epsilon_r=1+\frac{\langle\mathbf M^2\rangle
                         -\langle\mathbf M\rangle^2}
                       {3\epsilon_0 k_B T\langle V\rangle},
\end{equation}
con condizioni elettrostatiche conduttrici al contorno e cariche fisse.
Il dipolo dell'ultimo fotogramma, ricostruito separatamente dalle
coordinate e dai parametri pubblicati, concorda entro
0,0005\,Debye con il programma. Il contributo elettronico non \`e
incluso da un modello a cariche fisse; non basta a giustificare
automaticamente il divario osservato.

\textbf{La vicinanza della densit\`a non qualifica il modello per
una classifica quantitativa del distacco.} Vanno controllate durata,
repliche, convenzioni 1--4, conformazione del solvente e polarizzazione.
Una conversione errata del volume nelle prime esecuzioni \`e stata
rilevata mediante controllo fra motori e derivazione indipendente
in SI: quei tentativi sono archiviati come invalidi ed esclusi.

\subsection{Potenziali appresi e catene realistiche: stato della ricerca}

MACE-OFF23-small \cite{MACE2025} \`e stato effettivamente eseguito
sulle 13 coppie e confrontato con i calcoli precedenti. Usa un diverso
livello quantistico di addestramento, quindi le discrepanze non sono
una misura universale dell'errore rispetto a PBE0. Il modello non
comprende silicio: non descrive l'intero supporto Si/\SiO.
La sua esecuzione Metal a precisione singola \`e stata confrontata
con CPU a precisione doppia; la differenza quadratica media delle
forze sul complesso di controllo \`e circa
$1,9\times10^{-6}$\,eV/\AA. Una modifica locale, registrata,
riguarda soltanto il tipo numerico delle energie atomiche diagnostiche
su Metal, senza modificare l'espressione di energia totale e forze.

Sono inoltre state costruite topologie di PMMA a 4 e 8 unit\`a con
ogni termine di legame, angolo e torsione vincolato alla Tabella S2
di Behbahani e Harmandaris \cite{Behbahani2021}. Mancanze avrebbero
prodotto un errore esplicito. La neutralizzazione di 0,0051 elettroni
in eccesso ai due estremi \`e una scelta di chiusura documentata,
ancora da verificare. Quelle catene non sono ancora una simulazione
del distacco e non sostituiscono la massa 950K.

\textbf{Stato: calcoli elettronici e controllo del liquido completati;
energia libera del distacco interfaciale non ancora qualificata.}
\`E stato eseguito anche MACE-POLAR-1-small \cite{MACEPolar2026},
con interazioni elettrostatiche e 83 elementi. Per gli stessi complessi
con acetato di metile otteniamo $-46{,}19$\,kJ/mol con HFIP,
$-44{,}83$ con acido acetico e $-28{,}27$ con IPA. Il livello di
addestramento \`e diverso da PBE0: la variazione dell'ordine fra
HFIP e acido, di pochi kJ/mol, non produce una previsione di pulizia.
La versione della libreria elettrostatica \`e stata fissata a 0.4.0
per corrispondere all'interfaccia del modello rilasciato.

Un controllo ulteriore applica un campo statico alla stessa molecola
HFIP, a nuclei fissi. PBE0/def2-TZVP produce una polarizzabilit\`a
lungo $z$ di circa $0{,}4242\,e\,\text{\AA}^2/\mathrm V$, coerente fra
ampiezze 0,002 e 0,01\,V/\AA. Nel percorso di esecuzione del
modello appreso otteniamo invece circa $-0{,}0335\,e\,\text{\AA}^2/\mathrm V$
a 0,01\,V/\AA. La stima usa
\begin{equation}
 \alpha_{zz}\simeq -\frac{E(F)+E(-F)-2E(0)}{F^2}.
\end{equation}
Per uno stato fondamentale elettronico stabile, l'energia ottenuta
minimizzando un funzionale con accoppiamento lineare al campo \`e
concava in $F$. Il segno negativo ottenuto dal modello appreso
richiede chiarire convenzioni, implementazione e dominio di validit\`a;
la risposta al campo di questa configurazione non viene promossa a
previsione fisica del liquido. Una buona energia di coppia non
sostituisce questo controllo.

\textbf{Restano aperti i controlli di sensibilit\`a 1--4 del liquido,
l'energia libera polimero--grafene, la selettivit\`a sull'ossido e
la misura sul campione.} Una classifica di efficacia o una soluzione
non viene dedotta dalla semplice disponibilit\`a dei programmi.

'''
chapter=chapter.replace('TABLE_QUANTUM',table).replace('RHO_BULK',f'{rho:.3f}'.replace('.',',')).replace('EPS_BULK',f'{epsilon:.3f}'.replace('.',','))
s=s.replace(r'\section{Decisione scientifica e prossimo risultato da cercare}',chapter+r'\section{Decisione scientifica e prossimo risultato da cercare}')
s=s.replace(r'\end{thebibliography}',r'''
\bibitem{WangHFIP2020}
Y. Wang, T. Duo, X. Xu et al., \emph{Eco-Friendly High-Performance
Poly(methyl methacrylate) Film Reinforced with Methylcellulose},
ACS Omega \textbf{5} (2020), 24256--24261.
\href{https://doi.org/10.1021/acsomega.0c02249}{doi:10.1021/acsomega.0c02249}.

\bibitem{Marchelli2022}
G. Marchelli, J. Ingenmey, O. Holl\'oczki, A. Chaumont e B. Kirchner,
\emph{Hydrogen Bonding and Vaporization Thermodynamics in
Hexafluoroisopropanol-Acetone and -Methanol Mixtures.
A Joined Cluster Analysis and Molecular Dynamic Study},
ChemPhysChem \textbf{23} (2022), e202100620.
\href{https://doi.org/10.1002/cphc.202100620}{doi:10.1002/cphc.202100620}.
Parametri dal supplemento consultato integralmente.

\bibitem{Casoria2024}
M. Casoria, M. Macchiagodena, P. Rovero et al.,
\emph{Upgrading of the general AMBER force field 2 for fluorinated
alcohol biosolvents: A validation for water solutions and melittin solvation},
Journal of Peptide Science \textbf{30} (2024), e3543.
\href{https://doi.org/10.1002/psc.3543}{doi:10.1002/psc.3543}.

\bibitem{Behbahani2021}
A. Foroozani Behbahani e V. Harmandaris,
\emph{Gradient of Segmental Dynamics in Stereoregular Poly(methyl
methacrylate) Melts Confined between Pristine or Oxidized Graphene Sheets},
Polymers \textbf{13} (2021), 830.
\href{https://doi.org/10.3390/polym13050830}{doi:10.3390/polym13050830}.
Consultati Tabelle S1--S2 e correzione
\href{https://doi.org/10.3390/polym14010106}{doi:10.3390/polym14010106}.

\bibitem{MACE2025}
D. P. Kov\'acs, J. H. Moore, N. J. Browning et al.,
\emph{MACE-OFF: Transferable Short Range Machine Learning Force Fields
for Organic Molecules}, arXiv:2312.15211, versione 5 (2025).
\href{https://arxiv.org/abs/2312.15211v5}{Manoscritto};
\href{https://github.com/ACEsuit/mace-off}{Modello ufficiale}.

\bibitem{MACEPolar2026}
I. Batatia, W. J. Baldwin, D. Kuryla et al.,
\emph{MACE-POLAR-1: A Polarisable Electrostatic Foundation Model
for Molecular Chemistry}, arXiv:2602.19411 (2026).
\href{https://arxiv.org/abs/2602.19411}{Manoscritto};
\href{https://github.com/ACEsuit/mace-foundations/releases/tag/mace_polar_1}{Modelli ufficiali}.
\end{thebibliography}''')
TARGET.write_text(s)
print(TARGET)
