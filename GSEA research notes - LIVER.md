# GSEA research notes - LIVER

10/06/2024

Created sampleTables and DGE objects for liver and head-kidney.

Wrote vaccine constrasts and extracted results tables for 10wpi and 4wpc.

*There were 53 pathways for which P-values were not calculated properly due to unbalanced (positive and negative) gene-level statistic values. For such pathways pval, padj, NES, log2err are set to NA.*

This kept happening, so I added an *nPermSimple = 100000* flag to *gseGO* and *nPermSimple = 100000* to *gsePathway*.

## 10WPI

### DNA vaccine

Signaling by interleukins shows upregulated.

| Description                                                  | NES       | setSize | Count |
| ------------------------------------------------------------ | --------- | ------- | ----- |
| rRNA modification in the nucleus and cytosol                 | -2.772010 | 54      | 39    |
| rRNA processing                                              | -2.453593 | 160     | 65    |
| rRNA processing in the nucleus and cytosol                   | -2.407239 | 154     | 61    |
| Major pathway of rRNA processing in the nucleolus and cytosol | -2.239763 | 145     | 53    |
| tRNA processing                                              | -2.218107 | 64      | 29    |
| tRNA modification in the nucleus and cytosol                 | -2.091645 | 21      | 15    |
| PERK regulates gene expression                               | -2.048805 | 24      | 11    |
| ATF4 activates genes in response to endoplasmic reticulum  stress | -1.987790 | 19      | 12    |
| Formation of tubulin folding intermediates by CCT/TriC       | -1.887469 | 10      | 4     |
| Transcriptional regulation by small RNAs                     | -1.805579 | 40      | 20    |
| RNA Polymerase I Promoter Escape                             | -1.793334 | 23      | 15    |
| RNA Pol II CTD phosphorylation and interaction with CE during HIV infection | -1.781392 | 23      | 13    |
| RNA Pol II CTD phosphorylation and interaction with CE       | -1.781392 | 23      | 13    |
| mRNA Capping                                                 | -1.750957 | 25      | 14    |
| Unfolded Protein Response (UPR)                              | -1.743634 | 70      | 28    |
| TAK1-dependent IKK and NF-kappa-B activation                 | -1.741067 | 18      | 7     |
| Gene Silencing by RNA                                        | -1.733720 | 52      | 24    |
| snRNP Assembly                                               | -1.727473 | 36      | 18    |
| Metabolism of non-coding RNA                                 | -1.727473 | 36      | 18    |
| KSRP (KHSRP) binds and destabilizes mRNA                     | -1.727124 | 11      | 8     |
| Formation of TC-NER Pre-Incision Complex                     | -1.720563 | 41      | 22    |
| HIV Life Cycle                                               | -1.711677 | 107     | 47    |
| Late Phase of HIV Life Cycle                                 | -1.710021 | 101     | 45    |
| mRNA decay by 3' to 5' exoribonuclease                       | -1.707305 | 13      | 9     |
| Nucleotide-binding domain, leucine rich repeat containing receptor (NLR) signaling pathways | -1.669651 | 30      | 9     |
| tRNA processing in the nucleus                               | -1.667210 | 38      | 16    |
| RNA Polymerase III Chain Elongation                          | -1.663220 | 12      | 7     |
| Viral Messenger RNA Synthesis                                | -1.653304 | 34      | 17    |
| Inflammasomes                                                | -1.653136 | 12      | 4     |
| SUMOylation of RNA binding proteins                          | -1.653020 | 29      | 15    |
| Mitochondrial translation elongation                         | -1.648479 | 70      | 33    |
| Nuclear import of Rev protein                                | -1.645935 | 24      | 13    |
| RNA polymerase II transcribes snRNA genes                    | -1.624966 | 58      | 27    |
| Mitochondrial translation                                    | -1.612842 | 75      | 34    |
| Interleukin-4 and Interleukin-13 signaling                   | -1.607890 | 39      | 19    |
| Formation of the Early Elongation Complex                    | -1.596029 | 28      | 14    |
| Formation of the HIV-1 Early Elongation Complex              | -1.596029 | 28      | 14    |
| tRNA Aminoacylation                                          | -1.587923 | 33      | 18    |
| Mitochondrial translation initiation                         | -1.573332 | 70      | 32    |
| Mitochondrial translation termination                        | -1.563183 | 69      | 31    |
| HCMV Late Events                                             | -1.551581 | 40      | 19    |
| HIV Infection                                                | -1.361084 | 170     | 60    |
| Processing of Capped Intron-Containing Pre-mRNA              | -1.346116 | 222     | 69    |
| Cell Cycle                                                   | 1.314339  | 363     | 127   |
| Diseases of signal transduction by growth factor receptors and second messengers | 1.316558  | 266     | 77    |
| Axon guidance                                                | 1.344675  | 289     | 118   |
| SARS-CoV Infections                                          | 1.349494  | 247     | 74    |
| Signaling by Receptor Tyrosine Kinases                       | 1.349672  | 257     | 71    |
| Nervous system development                                   | 1.390382  | 298     | 122   |
| Cell Cycle, Mitotic                                          | 1.405877  | 306     | 110   |
| Signaling by Nuclear Receptors                               | 1.426774  | 132     | 33    |
| Hemostasis                                                   | 1.432707  | 270     | 101   |
| MAPK family signaling cascades                               | 1.442332  | 159     | 41    |
| Signaling by WNT                                             | 1.476225  | 150     | 54    |
| Transmission across Chemical Synapses                        | 1.478456  | 70      | 24    |
| Mitotic Metaphase and Anaphase                               | 1.501343  | 146     | 58    |
| Mitotic Anaphase                                             | 1.501343  | 146     | 58    |
| Cell Cycle Checkpoints                                       | 1.507114  | 156     | 61    |
| GTP hydrolysis and joining of the 60S ribosomal subunit      | 1.514329  | 91      | 47    |
| L13a-mediated translational silencing of Ceruloplasmin expression | 1.517851  | 90      | 44    |
| Signaling by GPCR                                            | 1.521175  | 154     | 54    |
| GPCR downstream signalling                                   | 1.529661  | 147     | 52    |
| Beta-catenin independent WNT signaling                       | 1.541891  | 80      | 33    |
| Muscle contraction                                           | 1.559529  | 51      | 15    |
| Response of EIF2AK4 (GCN2) to amino acid deficiency          | 1.562418  | 80      | 42    |
| Formation of a pool of free 40S subunits                     | 1.583754  | 82      | 41    |
| Integration of energy metabolism                             | 1.584293  | 38      | 14    |
| RHOG GTPase cycle                                            | 1.586506  | 51      | 16    |
| Platelet activation, signaling and aggregation               | 1.590911  | 125     | 39    |
| RAC2 GTPase cycle                                            | 1.594233  | 59      | 19    |
| RND1 GTPase cycle                                            | 1.594742  | 27      | 9     |
| Parasite infection                                           | 1.607800  | 42      | 22    |
| Leishmania phagocytosis                                      | 1.607800  | 42      | 22    |
| FCGR3A-mediated phagocytosis                                 | 1.607800  | 42      | 22    |
| ECM proteoglycans                                            | 1.610333  | 30      | 10    |
| Neurotransmitter receptors and postsynaptic signal transmission | 1.624353  | 53      | 20    |
| Nonsense-Mediated Decay (NMD)                                | 1.626122  | 87      | 48    |
| Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC) | 1.626122  | 87      | 48    |
| Semaphorin interactions                                      | 1.627524  | 35      | 15    |
| GPVI-mediated activation cascade                             | 1.628625  | 20      | 9     |
| PLC beta mediated events                                     | 1.629346  | 16      | 7     |
| RAF activation                                               | 1.631989  | 21      | 10    |
| Myogenesis                                                   | 1.640567  | 13      | 7     |
| Eukaryotic Translation Termination                           | 1.642670  | 72      | 40    |
| Separation of Sister Chromatids                              | 1.648325  | 115     | 48    |
| MAP2K and MAPK activation                                    | 1.650966  | 23      | 11    |
| Meiotic synapsis                                             | 1.654659  | 16      | 7     |
| Oncogenic MAPK signaling                                     | 1.662560  | 44      | 18    |
| Potential therapeutics for SARS                              | 1.665109  | 63      | 27    |
| Gene and protein expression by JAK-STAT signaling after Interleukin-12 stimulation | 1.671597  | 26      | 12    |
| CDC42 GTPase cycle                                           | 1.674033  | 86      | 22    |
| Sensory Perception                                           | 1.690622  | 60      | 28    |
| Collagen biosynthesis and modifying enzymes                  | 1.693640  | 23      | 7     |
| Degradation of cysteine and homocysteine                     | 1.694713  | 10      | 4     |
| RHO GTPase cycle                                             | 1.708945  | 264     | 84    |
| TRAF6 mediated IRF7 activation                               | 1.709489  | 11      | 6     |
| NR1H2 and NR1H3-mediated signaling                           | 1.711277  | 26      | 7     |
| Viral mRNA Translation                                       | 1.717562  | 69      | 38    |
| Mitotic Prometaphase                                         | 1.720485  | 107     | 42    |
| Extracellular matrix organization                            | 1.730960  | 117     | 37    |
| SRP-dependent cotranslational protein targeting to membrane  | 1.731785  | 88      | 48    |
| Signaling by high-kinase activity BRAF mutants               | 1.732551  | 20      | 11    |
| Notch-HLH transcription pathway                              | 1.733203  | 18      | 9     |
| Interleukin-12 signaling                                     | 1.735706  | 32      | 14    |
| Cholesterol biosynthesis                                     | 1.735860  | 23      | 14    |
| Mitotic Spindle Checkpoint                                   | 1.748624  | 62      | 28    |
| RHO GTPases Activate Formins                                 | 1.748758  | 68      | 32    |
| EPHB-mediated forward signaling                              | 1.754624  | 23      | 13    |
| G alpha (12/13) signalling events                            | 1.760218  | 36      | 20    |
| Interleukin-12 family signaling                              | 1.761629  | 35      | 15    |
| G-protein beta:gamma signalling                              | 1.770287  | 14      | 8     |
| EML4 and NUDC in mitotic spindle formation                   | 1.788346  | 61      | 29    |
| L1CAM interactions                                           | 1.793282  | 53      | 18    |
| Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC) | 1.799826  | 74      | 43    |
| RHOD GTPase cycle                                            | 1.802760  | 39      | 19    |
| Signaling by Rho GTPases, Miro GTPases and RHOBTB3           | 1.818813  | 382     | 120   |
| RHOA GTPase cycle                                            | 1.821890  | 87      | 26    |
| Peptide chain elongation                                     | 1.824288  | 69      | 38    |
| Signaling by Rho GTPases                                     | 1.832108  | 372     | 117   |
| Selenocysteine synthesis                                     | 1.833169  | 72      | 39    |
| RHOC GTPase cycle                                            | 1.847258  | 50      | 17    |
| RND3 GTPase cycle                                            | 1.859686  | 30      | 15    |
| Signaling by RAF1 mutants                                    | 1.860056  | 23      | 13    |
| Eukaryotic Translation Elongation                            | 1.861776  | 71      | 40    |
| Amplification of signal from the kinetochores                | 1.872992  | 54      | 26    |
| Amplification  of signal from unattached  kinetochores via a MAD2  inhibitory signal | 1.872992  | 54      | 26    |
| Cardiac conduction                                           | 1.921692  | 26      | 8     |
| Signaling by moderate kinase activity BRAF mutants           | 1.930290  | 25      | 14    |
| Signaling by RAS mutants                                     | 1.930290  | 25      | 14    |
| Paradoxical activation of RAF signaling by kinase inactive BRAF | 1.930290  | 25      | 14    |
| Signaling downstream of RAS mutants                          | 1.930290  | 25      | 14    |
| RHO GTPase Effectors                                         | 1.934879  | 146     | 72    |
| RHOQ GTPase cycle                                            | 1.973296  | 38      | 14    |
| Resolution of Sister Chromatid Cohesion                      | 1.978613  | 69      | 36    |
| RHOB GTPase cycle                                            | 2.086523  | 42      | 16    |

### EOMES

Signaling by interleukins shows upregulated.

| Description                                                  | NES       | setSize | Count |
| ------------------------------------------------------------ | --------- | ------- | ----- |
| rRNA modification in the nucleus and cytosol                 | -2.629134 | 57      | 38    |
| Mitochondrial translation                                    | -2.510041 | 82      | 52    |
| Mitochondrial translation elongation                         | -2.498719 | 76      | 49    |
| Mitochondrial translation termination                        | -2.465083 | 76      | 48    |
| Mitochondrial translation initiation                         | -2.452399 | 76      | 47    |
| Cobalamin (Cbl, vitamin B12) transport and metabolism        | -2.190745 | 15      | 4     |
| Activation of Matrix Metalloproteinases                      | -2.169753 | 15      | 2     |
| Mitochondrial tRNA aminoacylation                            | -2.086415 | 20      | 14    |
| tRNA Aminoacylation                                          | -2.004145 | 40      | 28    |
| rRNA processing                                              | -1.981503 | 174     | 60    |
| rRNA processing in the mitochondrion                         | -1.974662 | 10      | 8     |
| Removal of the Flap Intermediate                             | -1.917535 | 13      | 8     |
| Processive synthesis on the lagging strand                   | -1.902146 | 14      | 8     |
| Kidney development                                           | -1.896871 | 19      | 7     |
| Extension of Telomeres                                       | -1.887146 | 45      | 20    |
| rRNA processing in the nucleus and cytosol                   | -1.882236 | 164     | 54    |
| Translation                                                  | -1.855510 | 250     | 87    |
| Chromosome Maintenance                                       | -1.846192 | 80      | 32    |
| Formation of ATP by chemiosmotic coupling                    | -1.827672 | 12      | 10    |
| tRNA processing                                              | -1.795627 | 89      | 35    |
| RNA Polymerase III Chain Elongation                          | -1.783728 | 16      | 10    |
| Telomere Maintenance                                         | -1.765739 | 62      | 30    |
| Major pathway of rRNA processing in the nucleolus and cytosol | -1.724958 | 154     | 46    |
| E2F mediated regulation of DNA replication                   | -1.701744 | 16      | 8     |
| DNA strand elongation                                        | -1.700657 | 30      | 17    |
| Mitochondrial protein import                                 | -1.698303 | 54      | 24    |
| Cytosolic tRNA aminoacylation                                | -1.680270 | 23      | 13    |
| Cristae formation                                            | -1.672193 | 24      | 14    |
| Hemostasis                                                   | 1.246532  | 399     | 152   |
| Signaling by Interleukins                                    | 1.278670  | 277     | 102   |
| Diseases of signal transduction by growth factor receptors and second messengers | 1.299682  | 349     | 123   |
| MAPK family signaling cascades                               | 1.305292  | 225     | 78    |
| Nervous system development                                   | 1.321692  | 403     | 130   |
| Axon guidance                                                | 1.324522  | 388     | 91    |
| Intracellular signaling by second messengers                 | 1.348471  | 230     | 72    |
| MAPK1/MAPK3 signaling                                        | 1.354531  | 193     | 68    |
| RAF/MAP kinase cascade                                       | 1.362480  | 190     | 67    |
| Signaling by Nuclear Receptors                               | 1.368205  | 181     | 61    |
| ESR-mediated signaling                                       | 1.369432  | 129     | 43    |
| Signaling by Receptor Tyrosine Kinases                       | 1.370449  | 383     | 184   |
| Chromatin modifying enzymes                                  | 1.398007  | 174     | 77    |
| Chromatin organization                                       | 1.398007  | 174     | 77    |
| MITF-M-regulated melanocyte development                      | 1.419168  | 85      | 34    |
| Sensory Perception                                           | 1.426483  | 112     | 56    |
| Epigenetic regulation of gene expression                     | 1.428471  | 116     | 37    |
| Adipogenesis                                                 | 1.431171  | 88      | 26    |
| RHO GTPase cycle                                             | 1.437374  | 374     | 144   |
| Signaling by NTRKs                                           | 1.447580  | 99      | 42    |
| Platelet activation, signaling and aggregation               | 1.450350  | 186     | 81    |
| Cell-Cell communication                                      | 1.457617  | 96      | 50    |
| Activation of NMDA receptors and postsynaptic events         | 1.462013  | 55      | 22    |
| Leishmania infection                                         | 1.466123  | 117     | 48    |
| Parasitic Infection Pathways                                 | 1.466123  | 117     | 48    |
| Sphingolipid metabolism                                      | 1.469900  | 77      | 23    |
| RAC1 GTPase cycle                                            | 1.470595  | 157     | 80    |
| Signaling by TGF-beta Receptor Complex                       | 1.474899  | 77      | 25    |
| Sensory processing of sound by inner hair cells of the cochlea | 1.475068  | 45      | 27    |
| Sensory processing of sound                                  | 1.475908  | 48      | 28    |
| Death Receptor Signaling                                     | 1.479910  | 118     | 42    |
| Inositol phosphate metabolism                                | 1.484833  | 40      | 20    |
| Regulation of insulin secretion                              | 1.497956  | 44      | 17    |
| RHOC GTPase cycle                                            | 1.503305  | 64      | 27    |
| Signaling by TGFB family members                             | 1.509704  | 97      | 32    |
| NRAGE signals death through JNK                              | 1.510095  | 45      | 21    |
| Signaling by NTRK1 (TRKA)                                    | 1.515012  | 85      | 37    |
| CDC42 GTPase cycle                                           | 1.520473  | 130     | 73    |
| Ca2+ pathway                                                 | 1.523515  | 42      | 25    |
| Platelet Aggregation (Plug Formation)                        | 1.524694  | 32      | 18    |
| Constitutive Signaling by Aberrant PI3K in Cancer            | 1.526696  | 42      | 22    |
| Pre-NOTCH Transcription and Translation                      | 1.527756  | 30      | 12    |
| G alpha (s) signalling events                                | 1.529277  | 69      | 36    |
| TGF-beta receptor signaling activates SMADs                  | 1.532312  | 38      | 15    |
| Potential therapeutics for SARS                              | 1.537806  | 77      | 42    |
| RUNX1 regulates genes involved in megakaryocyte differentiation and platelet function | 1.538548  | 33      | 16    |
| Nephrin family interactions                                  | 1.539485  | 19      | 9     |
| NOTCH1 Intracellular Domain Regulates Transcription          | 1.542501  | 38      | 13    |
| Sema4D induced cell migration and growth-cone collapse       | 1.543528  | 17      | 6     |
| Downregulation of SMAD2/3:SMAD4 transcriptional activity     | 1.544392  | 25      | 8     |
| Transcriptional regulation by RUNX3                          | 1.544888  | 79      | 24    |
| Collagen degradation                                         | 1.545320  | 35      | 18    |
| Post NMDA receptor activation events                         | 1.552244  | 46      | 19    |
| Neurotransmitter release cycle                               | 1.552616  | 31      | 15    |
| MAP2K and MAPK activation                                    | 1.556544  | 28      | 16    |
| Signaling by GPCR                                            | 1.556555  | 302     | 133   |
| MET promotes cell motility                                   | 1.557673  | 27      | 17    |
| PI3K events in ERBB2 signaling                               | 1.560929  | 11      | 6     |
| HDMs demethylate histones                                    | 1.562580  | 17      | 8     |
| NR1H2 and NR1H3-mediated signaling                           | 1.567284  | 36      | 17    |
| Regulation of CDH11 Expression and Function                  | 1.570134  | 21      | 13    |
| Adrenaline,noradrenaline inhibits insulin secretion          | 1.571358  | 14      | 6     |
| G alpha (q) signalling events                                | 1.575934  | 100     | 41    |
| Adenylate cyclase inhibitory pathway                         | 1.576767  | 11      | 6     |
| G alpha (12/13) signalling events                            | 1.579723  | 55      | 25    |
| Oncogenic MAPK signaling                                     | 1.583526  | 60      | 28    |
| Formation of WDR5-containing histone-modifying complexes     | 1.584590  | 35      | 17    |
| MITF-M-dependent gene expression                             | 1.585433  | 63      | 27    |
| Synthesis of IP3 and IP4 in the cytosol                      | 1.587776  | 21      | 12    |
| L1CAM interactions                                           | 1.590207  | 72      | 39    |
| Uptake and actions of bacterial toxins                       | 1.590709  | 18      | 9     |
| Interleukin-4 and Interleukin-13 signaling                   | 1.593440  | 58      | 27    |
| SUMOylation of intracellular receptors                       | 1.595570  | 25      | 14    |
| Elastic fibre formation                                      | 1.597056  | 35      | 13    |
| Regulation of beta-cell development                          | 1.600646  | 24      | 10    |
| G alpha (i) signalling events                                | 1.602068  | 115     | 51    |
| GPVI-mediated activation cascade                             | 1.608847  | 27      | 13    |
| Polo-like kinase mediated events                             | 1.609342  | 14      | 6     |
| RHOA GTPase cycle                                            | 1.615091  | 127     | 51    |
| Glucagon-like Peptide-1 (GLP1) regulates insulin secretion   | 1.617192  | 21      | 12    |
| Non-integrin membrane-ECM interactions                       | 1.622412  | 42      | 17    |
| Interferon gamma signaling                                   | 1.626992  | 40      | 24    |
| FOXO-mediated transcription of cell cycle genes              | 1.638358  | 13      | 8     |
| Sphingolipid catabolism                                      | 1.638525  | 10      | 5     |
| GPCR downstream signalling                                   | 1.639542  | 275     | 125   |
| Signaling by NODAL                                           | 1.641624  | 11      | 8     |
| Neurotransmitter receptors and postsynaptic signal transmission | 1.642256  | 107     | 52    |
| Signaling by high-kinase activity BRAF mutants               | 1.643647  | 25      | 16    |
| Integrin signaling                                           | 1.645587  | 24      | 17    |
| ECM proteoglycans                                            | 1.645736  | 46      | 22    |
| Bacterial Infection Pathways                                 | 1.646065  | 50      | 18    |
| HSF1-dependent transactivation                               | 1.648852  | 17      | 8     |
| Anti-inflammatory response favouring Leishmania parasite infection | 1.655132  | 49      | 21    |
| Leishmania parasite growth and survival                      | 1.655132  | 49      | 21    |
| Interaction between L1 and Ankyrins                          | 1.657753  | 13      | 7     |
| Integration of energy metabolism                             | 1.660360  | 71      | 29    |
| RHOQ GTPase cycle                                            | 1.660770  | 50      | 27    |
| NR1H3 & NR1H2 regulate gene expression linked to cholesterol transport and efflux | 1.665302  | 28      | 14    |
| Aquaporin-mediated transport                                 | 1.665910  | 29      | 14    |
| RHO GTPases activate CIT                                     | 1.666394  | 15      | 8     |
| G alpha (z) signalling events                                | 1.667774  | 29      | 16    |
| FCGR3A-mediated IL10 synthesis                               | 1.671625  | 29      | 13    |
| Long-term potentiation                                       | 1.675547  | 10      | 7     |
| Molecules associated with elastic fibres                     | 1.678294  | 28      | 13    |
| Signaling by BMP                                             | 1.682198  | 22      | 10    |
| ADORA2B mediated anti-inflammatory cytokines production      | 1.682606  | 24      | 13    |
| RUNX3 regulates NOTCH signaling                              | 1.688974  | 11      | 7     |
| RHOJ GTPase cycle                                            | 1.690497  | 47      | 26    |
| Neuronal System                                              | 1.697738  | 206     | 98    |
| Regulation of MITF-M-dependent genes involved in cell cycle and proliferation | 1.704705  | 13      | 7     |
| Voltage gated Potassium channels                             | 1.706622  | 13      | 9     |
| NCAM signaling for neurite out-growth                        | 1.708120  | 40      | 17    |
| Notch-HLH transcription pathway                              | 1.712344  | 22      | 9     |
| PKA activation in glucagon signalling                        | 1.717508  | 13      | 9     |
| Signaling by RAF1 mutants                                    | 1.722014  | 29      | 17    |
| Transmission across Chemical Synapses                        | 1.727897  | 144     | 68    |
| Glucagon signaling in metabolic regulation                   | 1.729771  | 20      | 12    |
| TP53 Regulates Transcription of Genes Involved in G1 Cell Cycle Arrest | 1.731729  | 12      | 5     |
| Regulation of gene expression in late stage (branching morphogenesis) pancreatic bud precursor cells | 1.741200  | 11      | 7     |
| PKA-mediated phosphorylation of CREB                         | 1.744077  | 16      | 8     |
| Myogenesis                                                   | 1.750783  | 20      | 13    |
| GPER1 signaling                                              | 1.764746  | 27      | 15    |
| Signaling by moderate kinase activity BRAF mutants           | 1.768421  | 32      | 19    |
| Signaling by RAS mutants                                     | 1.768421  | 32      | 19    |
| Paradoxical activation of RAF signaling by kinase inactive BRAF | 1.768421  | 32      | 19    |
| Signaling downstream of RAS mutants                          | 1.768421  | 32      | 19    |
| Vasopressin regulates renal water homeostasis via Aquaporins | 1.769831  | 24      | 13    |
| PKA activation                                               | 1.771608  | 14      | 8     |
| Phase 0 - rapid depolarisation                               | 1.777203  | 12      | 7     |
| G-protein mediated events                                    | 1.881901  | 39      | 24    |
| Opioid Signalling                                            | 1.900240  | 55      | 30    |
| DAG and IP3 signaling                                        | 1.907237  | 30      | 18    |
| PLC beta mediated events                                     | 1.910717  | 36      | 23    |
| Calmodulin induced events                                    | 1.970152  | 24      | 16    |
| CaM pathway                                                  | 1.970152  | 24      | 16    |
| Ca-dependent events                                          | 2.019479  | 26      | 18    |

### GATA3

| Description                                                  | NES       | setSize | Count |
| ------------------------------------------------------------ | --------- | ------- | ----- |
| rRNA processing in the mitochondrion                         | -2.379128 | 10      | 9     |
| Metabolism of Angiotensinogen to Angiotensins                | -2.372706 | 10      | 3     |
| rRNA modification in the nucleus and cytosol                 | -2.293831 | 57      | 34    |
| Mitochondrial tRNA aminoacylation                            | -2.253875 | 20      | 9     |
| Mitochondrial translation                                    | -1.936104 | 82      | 41    |
| Mitochondrial translation initiation                         | -1.887474 | 76      | 38    |
| Mitochondrial translation elongation                         | -1.875701 | 76      | 39    |
| Mitochondrial translation termination                        | -1.836147 | 76      | 38    |
| rRNA processing                                              | -1.637892 | 174     | 62    |
| rRNA processing in the nucleus and cytosol                   | -1.436040 | 164     | 55    |
| Signaling by Receptor Tyrosine Kinases                       | 1.237291  | 397     | 172   |
| RHO GTPase cycle                                             | 1.256677  | 382     | 199   |
| GPCR downstream signalling                                   | 1.330100  | 309     | 135   |
| Extracellular matrix organization                            | 1.330252  | 196     | 93    |
| Transmission across Chemical Synapses                        | 1.330932  | 159     | 75    |
| Signaling by GPCR                                            | 1.346442  | 339     | 135   |
| RAC1 GTPase cycle                                            | 1.362416  | 159     | 104   |
| Neurotransmitter receptors and postsynaptic signal transmission | 1.406526  | 116     | 62    |
| PLC beta mediated events                                     | 1.549178  | 36      | 29    |
| ECM proteoglycans                                            | 1.578233  | 50      | 18    |
| Non-integrin membrane-ECM interactions                       | 1.581483  | 44      | 17    |
| Myogenesis                                                   | 1.661284  | 23      | 13    |

### IV-HD

### IV-LD

## 4WPC

### DNA vaccine

| Description                                                  | NES       | setSize | Count |
| ------------------------------------------------------------ | --------- | ------- | ----- |
| SRP-dependent cotranslational protein targeting to membrane  | 2.746160  | 88      | 49    |
| Viral mRNA Translation                                       | 2.562936  | 69      | 43    |
| Eukaryotic Translation Elongation                            | 2.560727  | 71      | 44    |
| Peptide chain elongation                                     | 2.526258  | 69      | 42    |
| Selenocysteine synthesis                                     | 2.465764  | 73      | 44    |
| Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC) | 2.459142  | 74      | 43    |
| Formation of a pool of free 40S subunits                     | 2.438645  | 82      | 43    |
| L13a-mediated translational silencing of Ceruloplasmin expression | 2.412538  | 90      | 52    |
| GTP hydrolysis and joining of the 60S ribosomal subunit      | 2.389440  | 91      | 52    |
| Eukaryotic Translation Initiation                            | 2.260711  | 97      | 47    |
| Cap-dependent Translation Initiation                         | 2.260711  | 97      | 47    |
| Eukaryotic Translation Termination                           | 2.230120  | 73      | 39    |
| Response of EIF2AK4 (GCN2) to amino acid deficiency          | 2.215357  | 80      | 42    |
| Formation of the ternary complex, and subsequently, the 43S complex | 2.175515  | 44      | 23    |
| SARS-CoV-1 modulates host translation machinery              | 2.158552  | 29      | 15    |
| Nonsense-Mediated Decay (NMD)                                | 2.128311  | 89      | 43    |
| Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC) | 2.128311  | 89      | 43    |
| Formation of ATP by chemiosmotic coupling                    | 2.069048  | 12      | 7     |
| Ribosomal scanning and start codon recognition               | 2.001157  | 50      | 23    |
| Translation initiation complex formation                     | 1.997164  | 49      | 23    |
| Activation of the mRNA upon binding of the cap-binding complex and eIFs, and subsequent binding to 43S | 1.985995  | 50      | 23    |
| SARS-CoV-2 modulates host translation machinery              | 1.942443  | 38      | 16    |
| Selenoamino acid metabolism                                  | 1.923872  | 93      | 45    |
| Mitochondrial Fatty Acid Beta-Oxidation                      | 1.897481  | 27      | 10    |
| Cellular response to starvation                              | 1.895605  | 117     | 51    |
| Defects in vitamin and cofactor metabolism                   | 1.882877  | 17      | 11    |
| Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins. | 1.826944  | 96      | 45    |
| Regulation of expression of SLITs and ROBOs                  | 1.720593  | 132     | 49    |
| Fatty acid metabolism                                        | 1.692700  | 108     | 36    |
| Respiratory electron transport                               | 1.650109  | 80      | 36    |
| Signaling by ROBO receptors                                  | 1.603262  | 161     | 54    |
| Translation                                                  | 1.598214  | 237     | 84    |
| Aerobic respiration and respiratory electron transport       | 1.559828  | 151     | 61    |
| Metabolism of amino acids and derivatives                    | 1.373951  | 264     | 70    |
| Cytokine Signaling in Immune system                          | -1.333362 | 352     | 103   |
| RHO GTPase cycle                                             | -1.341345 | 305     | 81    |
| Signaling by Interleukins                                    | -1.450596 | 224     | 63    |
| Hemostasis                                                   | -1.510786 | 305     | 103   |
| RAC1 GTPase cycle                                            | -1.542343 | 128     | 36    |
| Platelet activation, signaling and aggregation               | -1.544752 | 145     | 46    |
| SLC-mediated transmembrane transport                         | -1.547072 | 104     | 36    |
| Chromosome Maintenance                                       | -1.593103 | 66      | 29    |
| VEGFA-VEGFR2 Pathway                                         | -1.608699 | 70      | 26    |
| RAC3 GTPase cycle                                            | -1.653358 | 64      | 20    |
| Mitochondrial tRNA aminoacylation                            | -1.655229 | 16      | 11    |
| Transport of inorganic cations/anions and amino acids/oligopeptides | -1.657739 | 45      | 15    |
| RHO GTPases Activate NADPH Oxidases                          | -1.741948 | 17      | 8     |
| Cyclin D associated events in G1                             | -1.765330 | 28      | 14    |
| G1 Phase                                                     | -1.765330 | 28      | 14    |
| Cell surface interactions at the vascular wall               | -1.768807 | 66      | 22    |
| RAF-independent MAPK1/3 activation                           | -1.791937 | 14      | 6     |
| Amino acid transport across the plasma membrane              | -1.807066 | 18      | 8     |
| Interleukin-4 and Interleukin-13 signaling                   | -1.848273 | 43      | 24    |

### EOMES

There is a marked difference between the number of enriched GO terms and enriched pathways from ReactomePA.

While 257 GO terms were enriched with *gseGO*, 400 pathways were tagged in *ReactomePA*.

In this case, I arranged them in ascending NES, as the downregulated pathways were more interesting.

| Description                                                  | NES       | setSize | Count |
| ------------------------------------------------------------ | --------- | ------- | ----- |
| Interleukin-4 and Interleukin-13 signaling                   | -2.262873 | 38      | 22    |
| GPER1 signaling                                              | -2.096039 | 15      | 7     |
| Fcgamma receptor (FCGR) dependent phagocytosis               | -2.089789 | 52      | 19    |
| FCGR3A-mediated IL10 synthesis                               | -2.083127 | 16      | 9     |
| Integrin signaling                                           | -2.057657 | 19      | 9     |
| G alpha (s) signalling events                                | -2.040387 | 24      | 10    |
| Platelet Aggregation (Plug Formation)                        | -2.038182 | 22      | 11    |
| RHO GTPase Effectors                                         | -2.018924 | 135     | 49    |
| Cell surface interactions at the vascular wall               | -2.016042 | 53      | 20    |
| GPVI-mediated activation cascade                             | -1.987619 | 19      | 8     |
| Hemostasis                                                   | -1.979233 | 261     | 87    |
| Interferon gamma signaling                                   | -1.945199 | 29      | 18    |
| Integrin cell surface interactions                           | -1.937112 | 36      | 15    |
| Leishmania infection                                         | -1.929806 | 79      | 32    |
| Parasitic Infection Pathways                                 | -1.929806 | 79      | 32    |
| Signaling by CSF1 (M-CSF) in myeloid cells                   | -1.928144 | 18      | 10    |
| Signal transduction by L1                                    | -1.918460 | 15      | 8     |
| Anti-inflammatory response favouring Leishmania parasite infection | -1.912667 | 29      | 17    |
| Leishmania parasite growth and survival                      | -1.912667 | 29      | 17    |
| Parasite infection                                           | -1.911090 | 39      | 13    |
| Leishmania phagocytosis                                      | -1.911090 | 39      | 13    |
| FCGR3A-mediated phagocytosis                                 | -1.911090 | 39      | 13    |
| RHO GTPases Activate Formins                                 | -1.883883 | 61      | 24    |
| RAC1 GTPase cycle                                            | -1.871119 | 102     | 38    |
| GPCR downstream signalling                                   | -1.864692 | 128     | 54    |
| Signaling by GPCR                                            | -1.854200 | 134     | 57    |
| Semaphorin interactions                                      | -1.840478 | 35      | 17    |
| RAC2 GTPase cycle                                            | -1.835605 | 54      | 20    |
| EML4 and NUDC in mitotic spindle formation                   | -1.834083 | 56      | 23    |
| FLT3 Signaling                                               | -1.832020 | 22      | 11    |
| Processive synthesis on the C-strand of the telomere         | -1.823964 | 12      | 7     |
| Regulation of actin dynamics for phagocytic cup formation    | -1.822127 | 39      | 12    |
| MAPK targets/ Nuclear events mediated by MAP kinases         | -1.804437 | 19      | 12    |
| Metabolism of nitric oxide: NOS3 activation and regulation   | -1.803939 | 11      | 4     |
| Platelet activation, signaling and aggregation               | -1.802865 | 121     | 37    |
| Extra-nuclear estrogen signaling                             | -1.801008 | 29      | 8     |
| VEGFA-VEGFR2 Pathway                                         | -1.793506 | 57      | 28    |
| DNA strand elongation                                        | -1.792021 | 24      | 15    |
| Collagen degradation                                         | -1.783196 | 28      | 9     |
| Amplification of signal from the kinetochores                | -1.782581 | 49      | 19    |
| Amplification  of signal from unattached  kinetochores via a MAD2  inhibitory signal | -1.782581 | 49      | 19    |
| Interleukin-3, Interleukin-5 and GM-CSF signaling            | -1.781255 | 26      | 13    |
| Mitotic Prometaphase                                         | -1.780742 | 96      | 43    |
| RAC3 GTPase cycle                                            | -1.776373 | 53      | 18    |
| Signaling by VEGF                                            | -1.771072 | 61      | 30    |
| Toll Like Receptor 7/8 (TLR7/8) Cascade                      | -1.770673 | 46      | 25    |
| MyD88 dependent cascade initiated on endosome                | -1.770673 | 46      | 25    |
| G1/S-Specific Transcription                                  | -1.764554 | 14      | 9     |
| Toll Like Receptor 9 (TLR9) Cascade                          | -1.764530 | 49      | 26    |
| Muscle contraction                                           | -1.761417 | 48      | 20    |
| Extracellular matrix organization                            | -1.760416 | 108     | 32    |
| Signaling by Interleukins                                    | -1.759556 | 197     | 71    |
| Cytokine Signaling in Immune system                          | -1.757916 | 308     | 105   |
| Antigen activates B Cell Receptor (BCR) leading to generation of second messengers | -1.756051 | 16      | 9     |
| EPHA-mediated growth cone collapse                           | -1.753683 | 10      | 8     |
| MAP kinase activation                                        | -1.752945 | 36      | 19    |
| Processive synthesis on the lagging strand                   | -1.749877 | 12      | 8     |
| Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell | -1.744003 | 19      | 11    |
| TRAF6 mediated induction of NFkB and MAP kinases upon TLR7/8 or 9 activation | -1.736590 | 45      | 24    |
| Mitotic Spindle Checkpoint                                   | -1.735430 | 57      | 21    |
| SMAD2/SMAD3:SMAD4 heterotrimer regulates transcription       | -1.728680 | 19      | 9     |
| Cell-extracellular matrix interactions                       | -1.727461 | 10      | 5     |
| MyD88-independent TLR4 cascade                               | -1.725437 | 54      | 26    |
| TRIF (TICAM1)-mediated TLR4 signaling                        | -1.725437 | 54      | 26    |
| Lagging Strand Synthesis                                     | -1.721041 | 16      | 11    |
| Telomere C-strand (Lagging Strand) Synthesis                 | -1.719281 | 19      | 11    |
| G-protein beta:gamma signalling                              | -1.709217 | 12      | 5     |
| Fanconi Anemia Pathway                                       | -1.699367 | 17      | 11    |
| FLT3 signaling in disease                                    | -1.698379 | 17      | 8     |
| G alpha (12/13) signalling events                            | -1.691303 | 32      | 20    |
| Toll Like Receptor 10 (TLR10) Cascade                        | -1.687917 | 45      | 23    |
| Toll Like Receptor 5 (TLR5) Cascade                          | -1.687917 | 45      | 23    |
| MyD88 cascade initiated on plasma membrane                   | -1.687917 | 45      | 23    |
| RHO GTPases Activate NADPH Oxidases                          | -1.686677 | 12      | 7     |
| Removal of the Flap Intermediate from the C-strand           | -1.684042 | 11      | 6     |
| Regulation of CDH11 Expression and Function                  | -1.682777 | 14      | 8     |
| Regulation of Homotypic Cell-Cell Adhesion                   | -1.682777 | 14      | 8     |
| Regulation of Expression and Function of Type II Classical Cadherins | -1.682777 | 14      | 8     |
| Translesion synthesis by Y family DNA polymerases bypasses lesions on DNA template | -1.681754 | 23      | 12    |
| Toll Like Receptor 4 (TLR4) Cascade                          | -1.681057 | 64      | 29    |
| Smooth Muscle Contraction                                    | -1.679834 | 20      | 10    |
| ERK/MAPK targets                                             | -1.679352 | 13      | 8     |
| Toll Like Receptor 3 (TLR3) Cascade                          | -1.676690 | 52      | 25    |
| Resolution of Sister Chromatid Cohesion                      | -1.674593 | 63      | 24    |
| Signaling by TGFB family members                             | -1.672626 | 63      | 23    |
| Sema4D induced cell migration and growth-cone collapse       | -1.670284 | 11      | 7     |
| Signaling by RAF1 mutants                                    | -1.668716 | 22      | 8     |
| Signaling by high-kinase activity BRAF mutants               | -1.668466 | 19      | 8     |
| Nuclear Envelope Breakdown                                   | -1.667025 | 39      | 11    |
| Regulation of signaling by CBL                               | -1.666557 | 15      | 9     |
| DAP12 interactions                                           | -1.664445 | 14      | 5     |
| DAP12 signaling                                              | -1.664445 | 14      | 5     |
| Removal of the Flap Intermediate                             | -1.664370 | 11      | 8     |
| Signaling by TGF-beta Receptor Complex                       | -1.663848 | 53      | 21    |
| Uptake and actions of bacterial toxins                       | -1.659376 | 12      | 2     |
| Signaling by ERBB2                                           | -1.658931 | 24      | 15    |
| Interleukin-17 signaling                                     | -1.646827 | 37      | 19    |
| ECM proteoglycans                                            | -1.645292 | 29      | 9     |
| Signaling by moderate kinase activity BRAF mutants           | -1.635207 | 24      | 8     |
| Signaling by RAS mutants                                     | -1.635207 | 24      | 8     |
| Paradoxical activation of RAF signaling by kinase inactive BRAF | -1.635207 | 24      | 8     |
| Signaling downstream of RAS mutants                          | -1.635207 | 24      | 8     |
| Depolymerization of the Nuclear Lamina                       | -1.629480 | 12      | 5     |
| Activation of ATR in response to replication stress          | -1.626161 | 23      | 17    |
| Signaling by Rho GTPases                                     | -1.625659 | 346     | 101   |
| Syndecan interactions                                        | -1.622241 | 17      | 7     |
| Elastic fibre formation                                      | -1.617773 | 17      | 9     |
| Role of phospholipids in phagocytosis                        | -1.613648 | 10      | 5     |
| FCERI mediated MAPK activation                               | -1.610961 | 18      | 12    |
| Cell Cycle, Mitotic                                          | -1.604524 | 282     | 106   |
| Activation of Matrix Metalloproteinases                      | -1.604274 | 12      | 4     |
| FCERI mediated Ca+2 mobilization                             | -1.601557 | 13      | 9     |
| Constitutive Signaling by Ligand-Responsive EGFR Cancer Variants | -1.597643 | 13      | 6     |
| Signaling by EGFR in Cancer                                  | -1.597643 | 13      | 6     |
| Signaling by Ligand-Responsive EGFR Variants in Cancer       | -1.597643 | 13      | 6     |
| Basigin interactions                                         | -1.588490 | 13      | 5     |
| Potential therapeutics for SARS                              | -1.587549 | 60      | 20    |
| Termination of translesion DNA synthesis                     | -1.586378 | 18      | 10    |
| Signaling by Rho GTPases, Miro GTPases and RHOBTB3           | -1.582732 | 354     | 101   |
| Signaling by FGFR1                                           | -1.581086 | 17      | 8     |
| HDR through Homologous Recombination (HRR)                   | -1.579300 | 27      | 13    |
| ADORA2B mediated anti-inflammatory cytokines production      | -1.578305 | 12      | 6     |
| Notch-HLH transcription pathway                              | -1.577985 | 17      | 7     |
| Miscellaneous transport and binding events                   | -1.573612 | 12      | 6     |
| Interleukin-2 family signaling                               | -1.572867 | 16      | 7     |
| Signaling by Receptor Tyrosine Kinases                       | -1.571384 | 243     | 89    |
| Cell junction organization                                   | -1.568458 | 43      | 17    |
| Gap-filling DNA repair synthesis and ligation in GG-NER      | -1.566696 | 16      | 9     |
| MAP2K and MAPK activation                                    | -1.565525 | 22      | 8     |
| G0 and Early G1                                              | -1.564954 | 15      | 10    |
| Mitotic Prophase                                             | -1.562224 | 56      | 13    |
| DDX58/IFIH1-mediated induction of interferon-alpha/beta      | -1.560684 | 36      | 13    |
| Signaling by CSF3 (G-CSF)                                    | -1.556781 | 20      | 8     |
| Translesion Synthesis by POLH                                | -1.555846 | 14      | 8     |
| L1CAM interactions                                           | -1.555558 | 50      | 17    |
| Signaling by FGFR3                                           | -1.551163 | 13      | 7     |
| Signaling by FGFR4                                           | -1.551163 | 13      | 7     |
| Death Receptor Signaling                                     | -1.548836 | 77      | 30    |
| MyD88:MAL(TIRAP) cascade initiated on plasma membrane        | -1.545719 | 48      | 24    |
| Toll Like Receptor TLR1:TLR2 Cascade                         | -1.545719 | 48      | 24    |
| Toll Like Receptor TLR6:TLR2 Cascade                         | -1.545719 | 48      | 24    |
| Toll Like Receptor 2 (TLR2) Cascade                          | -1.545719 | 48      | 24    |
| DARPP-32 events                                              | -1.545263 | 13      | 5     |
| TNFR1-induced NF-kappa-B signaling pathway                   | -1.544483 | 16      | 8     |
| Constitutive Signaling by EGFRvIII                           | -1.543304 | 11      | 5     |
| Signaling by EGFRvIII in Cancer                              | -1.543304 | 11      | 5     |
| G-protein mediated events                                    | -1.542183 | 17      | 11    |
| Cell-Cell communication                                      | -1.540384 | 57      | 20    |
| Signaling by NTRK2 (TRKB)                                    | -1.537126 | 11      | 7     |
| Downregulation of ERBB2 signaling                            | -1.536510 | 12      | 7     |
| Amino acid transport across the plasma membrane              | -1.536036 | 15      | 6     |
| Signaling by FGFR                                            | -1.533343 | 39      | 16    |
| RHOB GTPase cycle                                            | -1.529351 | 40      | 22    |
| NRAGE signals death through JNK                              | -1.527677 | 27      | 17    |
| Signaling by FLT3 fusion proteins                            | -1.527151 | 12      | 8     |
| Chromosome Maintenance                                       | -1.526339 | 55      | 19    |
| Extension of Telomeres                                       | -1.525929 | 31      | 19    |
| VEGFR2 mediated vascular permeability                        | -1.525305 | 18      | 12    |
| Resolution of Abasic Sites (AP sites)                        | -1.520702 | 23      | 11    |
| Signaling by ALK in cancer                                   | -1.514907 | 61      | 26    |
| Signaling by ALK fusions and activated point mutants         | -1.514907 | 61      | 26    |
| Homology Directed Repair                                     | -1.507547 | 47      | 19    |
| Nuclear events stimulated by ALK signaling in cancer         | -1.506139 | 20      | 11    |
| Inactivation of CSF3 (G-CSF) signaling                       | -1.506066 | 16      | 7     |
| Degradation of the extracellular matrix                      | -1.505435 | 54      | 13    |
| Platelet degranulation                                       | -1.503440 | 65      | 14    |
| Cell Cycle                                                   | -1.495341 | 334     | 121   |
| p75 NTR receptor-mediated signalling                         | -1.492332 | 47      | 24    |
| Toll-like Receptor Cascades                                  | -1.490329 | 76      | 31    |
| Assembly of collagen fibrils and other multimeric structures | -1.488964 | 22      | 7     |
| Interferon alpha/beta signaling                              | -1.487745 | 26      | 14    |
| Response to elevated platelet cytosolic Ca2+                 | -1.486235 | 69      | 14    |
| Deactivation of the beta-catenin transactivating complex     | -1.485505 | 19      | 9     |
| DNA Damage Bypass                                            | -1.485166 | 29      | 13    |
| Factors involved in megakaryocyte development and platelet production | -1.484846 | 71      | 22    |
| Opioid Signalling                                            | -1.479913 | 28      | 15    |
| Signaling by SCF-KIT                                         | -1.476386 | 25      | 11    |
| DNA Repair                                                   | -1.465392 | 150     | 40    |
| Transcriptional activity of SMAD2/SMAD3:SMAD4 heterotrimer   | -1.464608 | 29      | 10    |
| Negative regulation of the PI3K/AKT network                  | -1.464134 | 39      | 23    |
| Signaling by BRAF and RAF1 fusions                           | -1.461940 | 32      | 9     |
| AURKA Activation by TPX2                                     | -1.459428 | 35      | 10    |
| Cell death signalling via NRAGE, NRIF and NADE               | -1.457175 | 38      | 20    |
| TGF-beta receptor signaling activates SMADs                  | -1.456626 | 25      | 12    |
| RHO GTPase cycle                                             | -1.450703 | 245     | 64    |
| Pre-NOTCH Transcription and Translation                      | -1.448507 | 22      | 15    |
| Base Excision Repair                                         | -1.445906 | 27      | 11    |
| Interferon Signaling                                         | -1.445318 | 110     | 39    |
| Negative regulation of MAPK pathway                          | -1.444211 | 22      | 13    |
| NOTCH1 Intracellular Domain Regulates Transcription          | -1.439979 | 27      | 10    |
| SUMOylation of DNA damage response and repair proteins       | -1.439858 | 43      | 9     |
| Non-integrin membrane-ECM interactions                       | -1.438765 | 27      | 11    |
| Rab regulation of trafficking                                | -1.436890 | 67      | 17    |
| Signaling by FGFR2                                           | -1.436754 | 35      | 14    |
| G alpha (q) signalling events                                | -1.435915 | 48      | 21    |
| PI5P, PP2A and IER3 Regulate PI3K/AKT Signaling              | -1.424922 | 37      | 22    |
| Signaling by MET                                             | -1.414098 | 39      | 20    |
| SARS-CoV-2 activates/modulates innate and adaptive immune responses | -1.406048 | 60      | 12    |
| DNA Double-Strand Break Repair                               | -1.399285 | 63      | 22    |
| Generic Transcription Pathway                                | -1.396730 | 473     | 164   |
| CDC42 GTPase cycle                                           | -1.393908 | 80      | 33    |
| Diseases of signal transduction by growth factor receptors and second messengers | -1.391798 | 256     | 71    |
| M Phase                                                      | -1.377932 | 201     | 66    |
| SUMO E3 ligases SUMOylate target proteins                    | -1.372799 | 92      | 23    |
| Cilium Assembly                                              | -1.357340 | 71      | 15    |
| ESR-mediated signaling                                       | -1.340824 | 91      | 38    |
| SUMOylation                                                  | -1.336488 | 95      | 23    |
| Adaptive Immune System                                       | -1.335264 | 341     | 112   |
| Mitotic G1 phase and G1/S transition                         | -1.324119 | 96      | 40    |
| Transcriptional Regulation by TP53                           | -1.301197 | 208     | 78    |
| Innate Immune System                                         | -1.296696 | 480     | 125   |
| Diseases of metabolism                                       | 1.316776  | 115     | 44    |
| Respiratory electron transport                               | 1.353684  | 74      | 34    |
| Aerobic respiration and respiratory electron transport       | 1.380650  | 140     | 52    |
| Mitochondrial protein degradation                            | 1.384518  | 73      | 28    |
| Autodegradation of the E3 ubiquitin ligase COP1              | 1.412268  | 43      | 24    |
| Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs) | 1.419517  | 62      | 31    |
| Stabilization of p53                                         | 1.443104  | 44      | 25    |
| Nonsense-Mediated Decay (NMD)                                | 1.459617  | 86      | 45    |
| Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC) | 1.459617  | 86      | 45    |
| Metabolism of lipids                                         | 1.463174  | 381     | 95    |
| Signaling by ROBO receptors                                  | 1.470582  | 150     | 76    |
| Golgi-to-ER retrograde transport                             | 1.476496  | 79      | 39    |
| SARS-CoV-1 modulates host translation machinery              | 1.485222  | 29      | 15    |
| ABC transporter disorders                                    | 1.489243  | 59      | 30    |
| COPI-dependent Golgi-to-ER retrograde traffic                | 1.490854  | 52      | 17    |
| Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins. | 1.495957  | 90      | 41    |
| Protein localization                                         | 1.498337  | 106     | 33    |
| Disorders of transmembrane transporters                      | 1.502778  | 104     | 48    |
| Cellular response to starvation                              | 1.505566  | 112     | 54    |
| Hh mutants are degraded by ERAD                              | 1.510001  | 46      | 25    |
| Biosynthesis of the N-glycan precursor (dolichol lipid-linked oligosaccharide, LLO) and transfer to a nascent protein | 1.511547  | 39      | 18    |
| Complement cascade                                           | 1.512167  | 26      | 5     |
| Plasma lipoprotein assembly, remodeling, and clearance       | 1.524190  | 43      | 18    |
| rRNA processing in the nucleus and cytosol                   | 1.534583  | 151     | 75    |
| rRNA processing                                              | 1.535307  | 157     | 61    |
| Defective CFTR causes cystic fibrosis                        | 1.541313  | 50      | 27    |
| Fatty acid metabolism                                        | 1.545336  | 92      | 18    |
| ABC-family proteins mediated transport                       | 1.549515  | 66      | 32    |
| Post-translational protein phosphorylation                   | 1.554096  | 54      | 27    |
| COPII-mediated vesicle transport                             | 1.562566  | 38      | 16    |
| Lysine catabolism                                            | 1.565716  | 12      | 7     |
| Regulation of Glucokinase by Glucokinase Regulatory Protein  | 1.567169  | 22      | 2     |
| Defective TPR may confer susceptibility towards thyroid papillary carcinoma (TPC) | 1.567169  | 22      | 2     |
| Cytochrome P450 - arranged by substrate type                 | 1.575572  | 21      | 12    |
| Metabolism of vitamins and cofactors                         | 1.578185  | 103     | 37    |
| Hh mutants abrogate ligand secretion                         | 1.578914  | 47      | 26    |
| Major pathway of rRNA processing in the nucleolus and cytosol | 1.581543  | 142     | 72    |
| Intra-Golgi and retrograde Golgi-to-ER traffic               | 1.592137  | 124     | 57    |
| Activation of the mRNA upon binding of the cap-binding complex and eIFs, and subsequent binding to 43S | 1.603380  | 49      | 27    |
| Translation initiation complex formation                     | 1.618414  | 48      | 27    |
| Maturation of spike protein                                  | 1.626331  | 19      | 7     |
| Glutathione conjugation                                      | 1.627009  | 13      | 6     |
| Retrograde transport at the Trans-Golgi-Network              | 1.627279  | 35      | 14    |
| Hedgehog ligand biogenesis                                   | 1.632777  | 49      | 28    |
| Sensory Perception                                           | 1.633976  | 57      | 17    |
| Ribosomal scanning and start codon recognition               | 1.645786  | 49      | 28    |
| Visual phototransduction                                     | 1.679929  | 34      | 12    |
| Regulation of expression of SLITs and ROBOs                  | 1.681093  | 128     | 71    |
| Fatty acyl-CoA biosynthesis                                  | 1.689937  | 23      | 7     |
| Response of EIF2AK4 (GCN2) to amino acid deficiency          | 1.703857  | 78      | 46    |
| Phase II - Conjugation of compounds                          | 1.721793  | 38      | 14    |
| Metabolism of water-soluble vitamins and cofactors           | 1.722364  | 63      | 28    |
| Formation of ATP by chemiosmotic coupling                    | 1.740552  | 12      | 11    |
| Mitochondrial translation                                    | 1.757240  | 71      | 33    |
| Metabolism of amino acids and derivatives                    | 1.792006  | 246     | 120   |
| Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC) | 1.796156  | 73      | 42    |
| Transport to the Golgi and subsequent modification           | 1.823875  | 101     | 44    |
| Mitochondrial translation elongation                         | 1.828171  | 67      | 33    |
| Mitochondrial translation initiation                         | 1.830089  | 66      | 19    |
| Mitochondrial translation termination                        | 1.839537  | 66      | 19    |
| Nucleotide biosynthesis                                      | 1.841813  | 13      | 7     |
| Cytosolic tRNA aminoacylation                                | 1.842301  | 21      | 11    |
| Phase I - Functionalization of compounds                     | 1.842301  | 34      | 19    |
| Asparagine N-linked glycosylation                            | 1.864055  | 172     | 74    |
| Biological oxidations                                        | 1.886675  | 77      | 32    |
| Eukaryotic Translation Initiation                            | 1.913272  | 95      | 57    |
| Cap-dependent Translation Initiation                         | 1.913272  | 95      | 57    |
| Formation of the ternary complex, and subsequently, the 43S complex | 1.917351  | 43      | 21    |
| Selenoamino acid metabolism                                  | 1.919174  | 90      | 49    |
| Defects in vitamin and cofactor metabolism                   | 1.919440  | 14      | 11    |
| Eukaryotic Translation Elongation                            | 1.920686  | 70      | 42    |
| Cristae formation                                            | 1.928859  | 24      | 13    |
| COPI-mediated anterograde transport                          | 1.929901  | 60      | 31    |
| Peptide chain elongation                                     | 1.931277  | 68      | 40    |
| Eukaryotic Translation Termination                           | 1.957626  | 71      | 47    |
| The canonical retinoid cycle in rods (twilight vision)       | 1.978377  | 10      | 8     |
| L13a-mediated translational silencing of Ceruloplasmin expression | 2.003959  | 88      | 53    |
| ER to Golgi Anterograde Transport                            | 2.008324  | 88      | 42    |
| GTP hydrolysis and joining of the 60S ribosomal subunit      | 2.015188  | 89      | 54    |
| Selenocysteine synthesis                                     | 2.036777  | 70      | 41    |
| Purine ribonucleoside monophosphate biosynthesis             | 2.048663  | 10      | 7     |
| Viral mRNA Translation                                       | 2.074703  | 68      | 41    |
| Formation of a pool of free 40S subunits                     | 2.149929  | 80      | 49    |
| Translation                                                  | 2.201753  | 223     | 121   |
| Regulation of cholesterol biosynthesis by SREBP (SREBF)      | 2.294680  | 45      | 17    |
| Metabolism of steroids                                       | 2.298662  | 86      | 32    |
| SRP-dependent cotranslational protein targeting to membrane  | 2.401320  | 87      | 58    |
| Activation of gene expression by SREBF (SREBP)               | 2.409997  | 35      | 17    |
| Cholesterol biosynthesis                                     | 2.573984  | 23      | 18    |

### GATA3

| Description                                                  | NES       | setSize | Count |
| ------------------------------------------------------------ | --------- | ------- | ----- |
| Transport of inorganic cations/anions and amino acids/oligopeptides | -1.929944 | 44      | 16    |
| RHO GTPases Activate NADPH Oxidases                          | -1.885335 | 15      | 8     |
| Cell surface interactions at the vascular wall               | -1.877235 | 65      | 30    |
| Amino acid transport across the plasma membrane              | -1.855301 | 18      | 10    |
| Platelet Aggregation (Plug Formation)                        | -1.799747 | 25      | 12    |
| Interleukin-4 and Interleukin-13 signaling                   | -1.785542 | 42      | 25    |
| Signaling by FGFR1                                           | -1.783707 | 24      | 11    |
| Binding and Uptake of Ligands by Scavenger Receptors         | -1.728694 | 23      | 5     |
| Negative regulation of FGFR1 signaling                       | -1.724222 | 13      | 5     |
| Basigin interactions                                         | -1.718560 | 15      | 7     |
| Signaling by GPCR                                            | -1.712427 | 167     | 56    |
| GPCR downstream signalling                                   | -1.708012 | 157     | 48    |
| Integrin signaling                                           | -1.704380 | 21      | 10    |
| Signaling by CSF1 (M-CSF) in myeloid cells                   | -1.702810 | 21      | 11    |
| RAC1 GTPase cycle                                            | -1.702467 | 120     | 46    |
| Cyclin D associated events in G1                             | -1.698183 | 27      | 12    |
| G1 Phase                                                     | -1.698183 | 27      | 12    |
| Signaling by FLT3 fusion proteins                            | -1.687607 | 13      | 9     |
| Signaling by ALK in cancer                                   | -1.674926 | 65      | 37    |
| Signaling by ALK fusions and activated point mutants         | -1.674926 | 65      | 37    |
| PRC2 methylates histones and DNA                             | -1.667802 | 16      | 7     |
| SLC-mediated transmembrane transport                         | -1.667667 | 100     | 26    |
| Class A/1 (Rhodopsin-like receptors)                         | -1.664491 | 47      | 13    |
| Role of LAT2/NTAL/LAB on calcium mobilization                | -1.661828 | 10      | 8     |
| RAC3 GTPase cycle                                            | -1.660552 | 60      | 23    |
| HSF1 activation                                              | -1.657144 | 11      | 2     |
| FLT3 signaling in disease                                    | -1.651469 | 18      | 10    |
| Signaling by FGFR3                                           | -1.638074 | 18      | 7     |
| Signaling by FGFR4                                           | -1.638074 | 18      | 7     |
| Transcriptional activity of SMAD2/SMAD3:SMAD4 heterotrimer   | -1.629110 | 37      | 13    |
| RHOB GTPase cycle                                            | -1.623881 | 45      | 22    |
| Downstream signaling of activated FGFR2                      | -1.620808 | 11      | 7     |
| Signaling by FGFR1 in disease                                | -1.619984 | 18      | 13    |
| Hemostasis                                                   | -1.618392 | 289     | 96    |
| RAF-independent MAPK1/3 activation                           | -1.615505 | 13      | 10    |
| Signaling by ERBB2 ECD mutants                               | -1.615229 | 13      | 8     |
| Syndecan interactions                                        | -1.613136 | 17      | 6     |
| Signaling by ERBB2 in Cancer                                 | -1.607487 | 14      | 8     |
| Signaling by ERBB2 KD Mutants                                | -1.607487 | 14      | 8     |
| G alpha (q) signalling events                                | -1.603283 | 58      | 23    |
| Platelet activation, signaling and aggregation               | -1.603022 | 134     | 47    |
| Signaling by ALK                                             | -1.599328 | 17      | 11    |
| GPVI-mediated activation cascade                             | -1.595708 | 22      | 11    |
| ROS and RNS production in phagocytes                         | -1.590936 | 18      | 9     |
| Signaling by FGFR                                            | -1.581584 | 46      | 17    |
| Integrin cell surface interactions                           | -1.576134 | 39      | 14    |
| Downstream signaling of activated FGFR3                      | -1.575990 | 10      | 6     |
| Downstream signaling of activated FGFR4                      | -1.575990 | 10      | 6     |
| Signaling by TGF-beta Receptor Complex                       | -1.571323 | 63      | 25    |
| Epigenetic regulation of gene expression                     | -1.565754 | 104     | 38    |
| Signaling by Type 1 Insulin-like Growth Factor 1 Receptor (IGF1R) | -1.565482 | 20      | 10    |
| IGF1R signaling cascade                                      | -1.565482 | 20      | 10    |
| Interferon gamma signaling                                   | -1.549480 | 30      | 15    |
| SMAD2/SMAD3:SMAD4 heterotrimer regulates transcription       | -1.540638 | 25      | 12    |
| Interleukin-37 signaling                                     | -1.539840 | 10      | 8     |
| G alpha (12/13) signalling events                            | -1.535679 | 39      | 17    |
| Signaling by CSF3 (G-CSF)                                    | -1.529843 | 22      | 12    |
| GPCR ligand binding                                          | -1.529233 | 61      | 14    |
| G alpha (i) signalling events                                | -1.520189 | 61      | 13    |
| rRNA modification in the nucleus and cytosol                 | -1.512065 | 54      | 26    |
| Fcgamma receptor (FCGR) dependent phagocytosis               | -1.502643 | 59      | 18    |
| Signaling by Interleukins                                    | -1.501554 | 216     | 86    |
| Diseases of signal transduction by growth factor receptors and second messengers | -1.492363 | 279     | 82    |
| Extracellular matrix organization                            | -1.475897 | 124     | 30    |
| Signaling by Receptor Tyrosine Kinases                       | -1.437949 | 277     | 76    |
| RHO GTPase cycle                                             | -1.426424 | 289     | 78    |
| CDC42 GTPase cycle                                           | -1.417491 | 97      | 34    |
| Cytokine Signaling in Immune system                          | -1.333779 | 338     | 112   |
| Metabolism of lipids                                         | 1.362471  | 436     | 114   |
| SLC transporter disorders                                    | 1.500487  | 53      | 13    |
| Diseases of metabolism                                       | 1.515116  | 131     | 36    |
| Complex I biogenesis                                         | 1.525126  | 42      | 17    |
| Metabolism of amino acids and derivatives                    | 1.530165  | 260     | 99    |
| Ubiquitin Mediated Degradation of Phosphorylated Cdc25A      | 1.530984  | 42      | 25    |
| p53-Independent DNA Damage Response                          | 1.530984  | 42      | 25    |
| p53-Independent G1/S DNA damage checkpoint                   | 1.530984  | 42      | 25    |
| AUF1 (hnRNP D0) binds and destabilizes mRNA                  | 1.554179  | 44      | 27    |
| Regulation of activated PAK-2p34 by proteasome mediated degradation | 1.559095  | 42      | 25    |
| Biological oxidations                                        | 1.559514  | 85      | 43    |
| Somitogenesis                                                | 1.577357  | 40      | 24    |
| Autodegradation of the E3 ubiquitin ligase COP1              | 1.579467  | 44      | 27    |
| Disorders of transmembrane transporters                      | 1.580149  | 114     | 22    |
| Regulation of Apoptosis                                      | 1.587880  | 43      | 26    |
| Mitochondrial protein degradation                            | 1.594150  | 75      | 30    |
| Sulfur amino acid metabolism                                 | 1.619725  | 22      | 10    |
| ABC transporter disorders                                    | 1.625230  | 61      | 31    |
| Regulation of beta-cell development                          | 1.625786  | 19      | 4     |
| Respiratory electron transport                               | 1.626783  | 78      | 33    |
| Metabolism of nucleotides                                    | 1.651396  | 69      | 18    |
| Asymmetric localization of PCP proteins                      | 1.651439  | 47      | 23    |
| Aquaporin-mediated transport                                 | 1.658361  | 17      | 5     |
| Pyruvate metabolism                                          | 1.667465  | 31      | 14    |
| RHO GTPases activate CIT                                     | 1.667608  | 10      | 3     |
| Degradation of DVL                                           | 1.667783  | 45      | 27    |
| Diseases of carbohydrate metabolism                          | 1.681234  | 22      | 8     |
| Drug ADME                                                    | 1.682511  | 36      | 10    |
| Degradation of AXIN                                          | 1.686945  | 45      | 27    |
| Regulation of Glucokinase by Glucokinase Regulatory Protein  | 1.693556  | 24      | 6     |
| Defective TPR may confer susceptibility towards thyroid papillary carcinoma (TPC) | 1.693556  | 24      | 6     |
| Nucleotide catabolism                                        | 1.696699  | 22      | 8     |
| RHO GTPases activate PAKs                                    | 1.697307  | 13      | 4     |
| Triglyceride catabolism                                      | 1.699835  | 10      | 4     |
| Lysine catabolism                                            | 1.713094  | 12      | 8     |
| Hedgehog ligand biogenesis                                   | 1.729542  | 49      | 30    |
| Regulation of pyruvate dehydrogenase (PDH) complex           | 1.729713  | 10      | 5     |
| Cristae formation                                            | 1.742941  | 24      | 13    |
| Visual phototransduction                                     | 1.748230  | 35      | 14    |
| Hh mutants abrogate ligand secretion                         | 1.749674  | 47      | 29    |
| Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins. | 1.751881  | 94      | 41    |
| Aerobic respiration and respiratory electron transport       | 1.753255  | 147     | 58    |
| Hh mutants are degraded by ERAD                              | 1.758956  | 46      | 29    |
| Peroxisomal protein import                                   | 1.759622  | 44      | 28    |
| Formation of ATP by chemiosmotic coupling                    | 1.788182  | 12      | 7     |
| Branched-chain amino acid catabolism                         | 1.846272  | 17      | 9     |
| Defective CFTR causes cystic fibrosis                        | 1.891005  | 51      | 32    |
| Degradation of cysteine and homocysteine                     | 1.894221  | 11      | 6     |
| Triglyceride metabolism                                      | 1.898246  | 18      | 7     |
| Initial triggering of complement                             | 1.912913  | 12      | 3     |
| Complement cascade                                           | 1.960510  | 28      | 5     |
| The canonical retinoid cycle in rods (twilight vision)       | 2.335860  | 10      | 8     |

### IV-HD

| Description                                                  | NES       | setSize | Count |
| ------------------------------------------------------------ | --------- | ------- | ----- |
| tRNA modification in the nucleus and cytosol                 | -2.177561 | 21      | 14    |
| tRNA processing                                              | -2.114216 | 64      | 33    |
| Activation of the pre-replicative complex                    | -2.092571 | 23      | 13    |
| Interleukin-4 and Interleukin-13 signaling                   | -2.024606 | 39      | 13    |
| DNA strand elongation                                        | -1.991399 | 25      | 13    |
| rRNA modification in the nucleus and cytosol                 | -1.951206 | 54      | 23    |
| Activation of ATR in response to replication stress          | -1.881421 | 26      | 13    |
| Synthesis of DNA                                             | -1.771895 | 89      | 30    |
| Regulation of actin dynamics for phagocytic cup formation    | -1.769697 | 41      | 11    |
| Fcgamma receptor (FCGR) dependent phagocytosis               | -1.733032 | 57      | 15    |
| Amplification of signal from the kinetochores                | -1.721040 | 54      | 22    |
| Amplification  of signal from unattached  kinetochores via a MAD2  inhibitory signal | -1.721040 | 54      | 22    |
| Recruitment of mitotic centrosome proteins and complexes     | -1.718997 | 41      | 13    |
| Centrosome maturation                                        | -1.718997 | 41      | 13    |
| Cell Cycle, Mitotic                                          | -1.711076 | 306     | 100   |
| RHO GTPases Activate Formins                                 | -1.688922 | 68      | 29    |
| DNA Replication                                              | -1.674303 | 95      | 30    |
| S Phase                                                      | -1.672746 | 116     | 37    |
| Mitotic Prometaphase                                         | -1.671646 | 107     | 38    |
| Response to elevated platelet cytosolic Ca2+                 | -1.664125 | 70      | 17    |
| Cell Cycle                                                   | -1.606754 | 363     | 110   |
| Platelet degranulation                                       | -1.604579 | 66      | 15    |
| G1/S Transition                                              | -1.598515 | 94      | 30    |
| Mitotic G1 phase and G1/S transition                         | -1.573187 | 103     | 32    |
| Cell Cycle Checkpoints                                       | -1.571586 | 156     | 52    |
| M Phase                                                      | -1.538726 | 216     | 67    |
| Sphingolipid metabolism                                      | 1.648678  | 51      | 20    |
| NCAM signaling for neurite out-growth                        | 1.803092  | 21      | 10    |
| Regulation of beta-cell development                          | 1.899808  | 17      | 5     |
| FOXO-mediated transcription of oxidative stress, metabolic and neuronal genes | 1.952205  | 14      | 5     |
| FOXO-mediated transcription                                  | 2.146193  | 33      | 14    |

### IV-LD

| Description                                                  | NES       | setSize | Count |
| ------------------------------------------------------------ | --------- | ------- | ----- |
| Antigen activates B Cell Receptor (BCR) leading to generation of second messengers | -1.963323 | 21      | 10    |
| Downregulation of SMAD2/3:SMAD4 transcriptional activity     | -1.922826 | 21      | 6     |
| RHO GTPases Activate NADPH Oxidases                          | -1.921551 | 17      | 8     |
| Transcriptional activity of SMAD2/SMAD3:SMAD4 heterotrimer   | -1.917286 | 37      | 9     |
| Fcgamma receptor (FCGR) dependent phagocytosis               | -1.880702 | 64      | 28    |
| Amino acid transport across the plasma membrane              | -1.879409 | 18      | 7     |
| Transport of inorganic cations/anions and amino acids/oligopeptides | -1.856819 | 45      | 15    |
| Signaling by ERBB2 in Cancer                                 | -1.838310 | 14      | 6     |
| Signaling by ERBB2 KD Mutants                                | -1.838310 | 14      | 6     |
| FCGR3A-mediated IL10 synthesis                               | -1.836638 | 21      | 9     |
| Platelet activation, signaling and aggregation               | -1.835522 | 145     | 52    |
| Cell surface interactions at the vascular wall               | -1.833250 | 66      | 26    |
| Signaling by ERBB2 ECD mutants                               | -1.804882 | 13      | 6     |
| Regulation of insulin secretion                              | -1.802926 | 26      | 9     |
| FCERI mediated Ca+2 mobilization                             | -1.802784 | 19      | 10    |
| ROS and RNS production in phagocytes                         | -1.792430 | 19      | 7     |
| Signaling by ERBB2 TMD/JMD mutants                           | -1.765167 | 11      | 5     |
| VEGFA-VEGFR2 Pathway                                         | -1.754364 | 70      | 18    |
| RHO GTPases Activate ROCKs                                   | -1.741422 | 12      | 9     |
| PI3K/AKT Signaling in Cancer                                 | -1.737552 | 46      | 26    |
| G alpha (q) signalling events                                | -1.729742 | 62      | 19    |
| Netrin-1 signaling                                           | -1.724029 | 17      | 6     |
| Signaling by VEGF                                            | -1.718570 | 74      | 18    |
| SMAD2/SMAD3:SMAD4 heterotrimer regulates transcription       | -1.714377 | 25      | 7     |
| Semaphorin interactions                                      | -1.709944 | 40      | 20    |
| G alpha (i) signalling events                                | -1.706816 | 65      | 19    |
| Platelet calcium homeostasis                                 | -1.702921 | 11      | 5     |
| Effects of PIP2 hydrolysis                                   | -1.700834 | 11      | 6     |
| Cholesterol biosynthesis                                     | -1.695521 | 23      | 11    |
| EPHA-mediated growth cone collapse                           | -1.688049 | 14      | 11    |
| Interleukin-4 and Interleukin-13 signaling                   | -1.686119 | 43      | 25    |
| Hemostasis                                                   | -1.685229 | 305     | 119   |
| RAF-independent MAPK1/3 activation                           | -1.677323 | 14      | 10    |
| G-protein mediated events                                    | -1.672555 | 22      | 6     |
| Sema4D induced cell migration and growth-cone collapse       | -1.668424 | 13      | 9     |
| Signaling by FGFR1                                           | -1.667585 | 25      | 14    |
| Role of phospholipids in phagocytosis                        | -1.666463 | 16      | 11    |
| Nuclear Receptor transcription pathway                       | -1.663007 | 32      | 11    |
| Sema4D in semaphorin signaling                               | -1.658031 | 17      | 10    |
| Smooth Muscle Contraction                                    | -1.646729 | 24      | 14    |
| RAC3 GTPase cycle                                            | -1.644313 | 64      | 24    |
| Signaling by ERBB2                                           | -1.644090 | 29      | 12    |
| Signaling by CSF1 (M-CSF) in myeloid cells                   | -1.643930 | 22      | 8     |
| PLC beta mediated events                                     | -1.640126 | 19      | 5     |
| GPCR downstream signalling                                   | -1.636954 | 170     | 49    |
| Constitutive Signaling by Aberrant PI3K in Cancer            | -1.634688 | 28      | 14    |
| Signaling by GPCR                                            | -1.634416 | 184     | 53    |
| Muscle contraction                                           | -1.633358 | 64      | 24    |
| Negative regulation of the PI3K/AKT network                  | -1.617022 | 49      | 25    |
| PI5P, PP2A and IER3 Regulate PI3K/AKT Signaling              | -1.615472 | 46      | 23    |
| SLC-mediated transmembrane transport                         | -1.615031 | 104     | 29    |
| Signaling by ALK in cancer                                   | -1.614893 | 66      | 26    |
| Signaling by ALK fusions and activated point mutants         | -1.614893 | 66      | 26    |
| Activation of gene expression by SREBF (SREBP)               | -1.612002 | 38      | 17    |
| Class A/1 (Rhodopsin-like receptors)                         | -1.609446 | 51      | 14    |
| Platelet degranulation                                       | -1.603596 | 74      | 22    |
| HSP90 chaperone cycle for steroid hormone receptors (SHR) in the presence of ligand | -1.599137 | 31      | 5     |
| Negative regulation of MAPK pathway                          | -1.593495 | 27      | 14    |
| Parasite infection                                           | -1.591010 | 46      | 18    |
| Leishmania phagocytosis                                      | -1.591010 | 46      | 18    |
| FCGR3A-mediated phagocytosis                                 | -1.591010 | 46      | 18    |
| Fc epsilon receptor (FCERI) signaling                        | -1.573421 | 93      | 24    |
| Anti-inflammatory response favouring Leishmania parasite infection | -1.571348 | 36      | 17    |
| Leishmania parasite growth and survival                      | -1.571348 | 36      | 17    |
| Response to elevated platelet cytosolic Ca2+                 | -1.569407 | 78      | 31    |
| Ca2+ pathway                                                 | -1.567243 | 27      | 10    |
| Signaling by Interleukins                                    | -1.561478 | 224     | 82    |
| Costimulation by the CD28 family                             | -1.547630 | 35      | 16    |
| RAC2 GTPase cycle                                            | -1.542744 | 66      | 16    |
| Signaling by Receptor Tyrosine Kinases                       | -1.523483 | 296     | 95    |
| Transcriptional regulation by RUNX1                          | -1.500291 | 122     | 41    |
| RHO GTPase Effectors                                         | -1.491957 | 170     | 58    |
| RAC1 GTPase cycle                                            | -1.474237 | 128     | 49    |
| Diseases of signal transduction by growth factor receptors and second messengers | -1.456128 | 289     | 105   |
| Cytokine Signaling in Immune system                          | -1.443442 | 352     | 125   |
| Adaptive Immune System                                       | -1.396534 | 409     | 120   |
| RHO GTPase cycle                                             | -1.353014 | 305     | 101   |
| Signaling by Rho GTPases                                     | -1.336077 | 429     | 134   |
| Signaling by Rho GTPases, Miro GTPases and RHOBTB3           | -1.332690 | 439     | 139   |
| rRNA processing in the nucleus and cytosol                   | 1.276087  | 157     | 75    |
| Major pathway of rRNA processing in the nucleolus and cytosol | 1.422418  | 147     | 77    |
| Metabolism of amino acids and derivatives                    | 1.426549  | 264     | 122   |
| Metabolism of nucleotides                                    | 1.469857  | 72      | 22    |
| Translation                                                  | 1.475773  | 237     | 117   |
| Cellular response to starvation                              | 1.561049  | 117     | 45    |
| Aerobic respiration and respiratory electron transport       | 1.664542  | 151     | 77    |
| Regulation of expression of SLITs and ROBOs                  | 1.668917  | 132     | 50    |
| Selenoamino acid metabolism                                  | 1.730473  | 93      | 56    |
| Retrograde neurotrophin signalling                           | 1.807523  | 11      | 6     |
| Response of EIF2AK4 (GCN2) to amino acid deficiency          | 1.812629  | 80      | 54    |
| Respiratory electron transport                               | 1.818833  | 80      | 44    |
| The canonical retinoid cycle in rods (twilight vision)       | 1.819757  | 10      | 8     |
| Formation of ATP by chemiosmotic coupling                    | 1.848615  | 12      | 9     |
| SARS-CoV-2 modulates host translation machinery              | 1.861811  | 38      | 25    |
| SARS-CoV-1 modulates host translation machinery              | 1.885818  | 29      | 21    |
| Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins. | 1.982287  | 96      | 54    |
| Activation of the mRNA upon binding of the cap-binding complex and eIFs, and subsequent binding to 43S | 2.004224  | 50      | 31    |
| Ribosomal scanning and start codon recognition               | 2.013422  | 50      | 31    |
| Nonsense-Mediated Decay (NMD)                                | 2.013986  | 89      | 57    |
| Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC) | 2.013986  | 89      | 57    |
| Cristae formation                                            | 2.016761  | 24      | 16    |
| Translation initiation complex formation                     | 2.024032  | 49      | 31    |
| Eukaryotic Translation Termination                           | 2.154457  | 73      | 51    |
| Formation of the ternary complex, and subsequently, the 43S complex | 2.157697  | 44      | 30    |
| SRP-dependent cotranslational protein targeting to membrane  | 2.237501  | 88      | 58    |
| Eukaryotic Translation Initiation                            | 2.271027  | 97      | 63    |
| Cap-dependent Translation Initiation                         | 2.271027  | 97      | 63    |
| Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC) | 2.283325  | 74      | 52    |
| Selenocysteine synthesis                                     | 2.335805  | 73      | 52    |
| L13a-mediated translational silencing of Ceruloplasmin expression | 2.347708  | 90      | 61    |
| GTP hydrolysis and joining of the 60S ribosomal subunit      | 2.348593  | 91      | 61    |
| Peptide chain elongation                                     | 2.375799  | 69      | 51    |
| Eukaryotic Translation Elongation                            | 2.390348  | 71      | 53    |
| Formation of a pool of free 40S subunits                     | 2.400952  | 82      | 58    |
| Viral mRNA Translation                                       | 2.415355  | 69      | 51    |

