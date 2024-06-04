There are no explicitly immune related pathways highly ranked in GSEA at 1WPC. However, they show up prominently at 4WPC.

Can this be related to the fact that the heart is not an 'immune organ'? What is the case in the liver, head-kidney, and spleen?

However, in EOMES, more general GO terms like *negative regulation of response to external stimulus*, *positive regulation of response to external stimulus*, and *regulation of response to external stimulus* are present, with relatively large gene sets, generally upregulated. Is this evidence of an early antiviral response?

Responses to external stimuli show up in the DNA vaccine, EOMES, and GATA3, but not in IV-HD? Evidence of viral infection recognition?



## 1WPC

### DNA vaccine

Running GSEA results through ReactomePA, I get a set of downregulated pathways for DNA vaccine:

| Description | NES |
| --- | --- |
| Interferon alpha/beta signaling | -2.511935 |
| Potassium Channels | -2.462609 |
| p130Cas linkage to MAPK signaling for integrins | -2.066685 |
| Formation of Fibrin Clot (Clotting Cascade) | -2.022246 |
| ATF4 activates genes in response to endoplasmic reticulum  stress | -1.990491 |
| TNFR1-induced NF-kappa-B signaling pathway | -1.947870 |
| Activation of IRF3, IRF7 mediated by TBK1, IKKε (IKBKE) | -1.925654 |
| PERK regulates gene expression | -1.912257 |
| Inactivation of CSF3 (G-CSF) signaling | -1.868618 |
| Membrane binding and targetting of GAG proteins | -1.867574 |
| Synthesis And Processing Of GAG, GAGPOL Polyproteins | -1.867574 |
| Nicotinamide salvaging | -1.867505 |
| Nicotinate metabolism | -1.866540 |
| Assembly Of The HIV Virion | -1.859421 |
| Synthesis of active ubiquitin: roles of E1 and E2 enzymes | -1.812187 |
| GRB2:SOS provides linkage to MAPK signaling for Integrins | -1.793079 |
| IKK complex recruitment mediated by RIP1 | -1.783806 |
| Interleukin-6 family signaling | -1.776561 |
| Regulation of NF-kappa B signaling | -1.772008 |
| Chaperone Mediated Autophagy | -1.762597 |
| TAK1-dependent IKK and NF-kappa-B activation | -1.753485 |
| Nucleotide catabolism | -1.750191 |
| Response of EIF2AK1 (HRI) to heme deficiency | -1.730220 |
| Complement cascade | -1.727122 |
| Nucleotide salvage | -1.724230 |
| Negative regulators of DDX58/IFIH1 signaling | -1.719055 |
| TICAM1, RIP1-mediated IKK complex recruitment | -1.697703 |
| Regulated Necrosis | -1.686290 |
| Signaling by CSF3 (G-CSF) | -1.671001 |
| Negative regulation of FLT3 | -1.659995 |
| Negative regulation of MAPK pathway | -1.631777 |
| Regulation of TNFR1 signaling | -1.594833 |
| rRNA modification in the nucleus and cytosol | -1.574979 |
| DDX58/IFIH1-mediated induction of interferon-alpha/beta | -1.552165 |
| ER-Phagosome pathway | -1.521858 |
| Interleukin-4 and Interleukin-13 signaling | -1.463792 |



Relevant pathways for immune response and inflammation, cell stress and homeostasis, and viral infection and response:

### Immune Response and Inflammation

1. **Interferon Alpha/Beta Signaling**:
   - Critical for antiviral defense, activating immune cells, and upregulating antigen presentation.
2. **TNFR1-induced NF-kappa-B Signaling Pathway**:
   - Involved in inflammation, immune response, and cell survival.
3. **Activation of IRF3, IRF7 Mediated by TBK1, IKKε (IKBKE)**:
   - Plays a role in antiviral responses by activating interferon-stimulated genes.
4. **Regulation of NF-kappa B Signaling**:
   - Central in controlling immune responses, inflammation, and cell survival.
5. **Interleukin-6 Family Signaling**:
   - Important for immune responses, inflammation, and hematopoiesis.
6. **Interleukin-4 and Interleukin-13 Signaling**:
   - Involved in the regulation of immune responses, particularly in allergic reactions.
7. **Complement Cascade**:
   - Part of the innate immune system, enhancing the ability to clear microbes and damaged cells.



### Cell Stress and Homeostasis

1. **ATF4 Activates Genes in Response to Endoplasmic Reticulum Stress**:
   - Manages stress responses within the endoplasmic reticulum, crucial for cell survival under stress.
2. **PERK Regulates Gene Expression**:
   - Involved in the unfolded protein response, helping cells cope with endoplasmic reticulum stress.
3. **Chaperone Mediated Autophagy**:
   - Selective degradation of damaged or misfolded proteins, crucial for cellular homeostasis.
4. **Response of EIF2AK1 (HRI) to Heme Deficiency**:
   - Regulates protein synthesis in response to heme availability, affecting red blood cell production.



### Viral Infection and Response

1. **Assembly of the HIV Virion**:
   - Pathway detailing the steps in the assembly of HIV, a key process in the viral life cycle.
2. **Membrane Binding and Targeting of GAG Proteins**:
   - Involved in the assembly and budding of retroviruses like HIV.
3. **Synthesis and Processing of GAG, GAGPOL Polyproteins**:
   - Crucial for the production of viral proteins necessary for viral replication and assembly.



### Summary

These pathways reflect a wide range of cellular functions, including immune responses, stress responses, metabolic processes, cell signaling, viral infection mechanisms, and apoptosis. This combination suggests a coordinated effort to maintain cellular homeostasis, respond to external and internal stressors, and regulate immune and inflammatory responses. Understanding these pathways can provide insights into various physiological conditions and potential therapeutic targets for diseases such as infections, cancer, immune disorders, and metabolic syndromes.



### EOMES

Running GSEA results through ReactomePA, I get a set of downregulated pathways:

| Description | NES |
| --- | --- |
| Potassium Channels | -2.584500 |
| Inwardly rectifying K+ channels | -2.124522 |
| Activation of IRF3, IRF7 mediated by TBK1, IKKε (IKBKE) | -2.007038 |
| Nicotinamide salvaging | -1.994544 |
| Neuronal System | -1.967783 |
| Nicotinate metabolism | -1.955000 |
| GABA B receptor activation | -1.936191 |
| Activation of GABAB receptors | -1.936191 |
| Response of EIF2AK1 (HRI) to heme deficiency | -1.902384 |
| ATF4 activates genes in response to endoplasmic reticulum  stress | -1.870037 |
| GABA receptor activation | -1.861655 |
| mRNA decay by 3' to 5' exoribonuclease | -1.853049 |
| Response of Mtb to phagocytosis | -1.784642 |
| Mitochondrial tRNA aminoacylation | -1.739713 |



This combination of pathways suggests an integration of diverse cellular functions:

- **Neuronal Function and Neurotransmission**: Pathways involving potassium channels, GABA receptor activation, and the neuronal system highlight critical mechanisms in maintaining neuronal excitability and neurotransmission.
- **Metabolic Processes**: Nicotinamide salvaging and nicotinate metabolism are crucial for NAD+ production, affecting cellular energy metabolism.
- **Immune Response and Stress**: Pathways involving IRF3/IRF7 activation, ATF4 in ER stress, and response to phagocytosis are indicative of cellular responses to stress and infection.
- **RNA Processing and Mitochondrial Function**: mRNA decay mechanisms and mitochondrial tRNA aminoacylation point to the importance of RNA regulation and mitochondrial protein synthesis.

Overall, this combination of pathways illustrates a complex network of interactions that are vital for cellular homeostasis, neuronal function, metabolic processes, and immune responses. Understanding these pathways provides insights into how cells respond to various physiological and pathological conditions.



### GATA3

| Description | NES |
| --- | --- |
| Keratinization | -2.378933 |
| Formation of the cornified envelope | -2.378933 |
| Interferon alpha/beta signaling | -2.337436 |
| Potassium Channels | -2.330262 |
| Nicotinate metabolism | -2.132847 |
| Nicotinamide salvaging | -2.116462 |
| Evasion by RSV of host interferon responses | -2.077602 |
| Activation of IRF3, IRF7 mediated by TBK1, IKKε (IKBKE) | -1.972133 |
| TRAF3-dependent IRF activation pathway | -1.945352 |
| Translation of Structural Proteins | -1.886217 |
| Inactivation of CSF3 (G-CSF) signaling | -1.866082 |
| Maturation of nucleoprotein | -1.817987 |
| TRAF6 mediated IRF7 activation | -1.786926 |
| Negative regulators of DDX58/IFIH1 signaling | -1.749586 |
| RIPK1-mediated regulated necrosis | -1.740030 |
| Regulation of necroptotic cell death | -1.740030 |
| Maturation of nucleoprotein | -1.727373 |
| Synthesis of PE | -1.713567 |
| Assembly Of The HIV Virion | -1.680388 |
| Interleukin-6 family signaling | -1.679486 |
| Chondroitin sulfate biosynthesis | -1.662545 |
| Membrane binding and targetting of GAG proteins | -1.653289 |
| Synthesis And Processing Of GAG, GAGPOL Polyproteins | -1.653289 |
| Signaling by CSF3 (G-CSF) | -1.652082 |



### IV-HD

| Description | NES |
| --- | --- |
| Interferon alpha/beta signaling | -2.493979 |
| Keratinization | -2.373031 |
| Formation of the cornified envelope | -2.373031 |
| RHO GTPase cycle |  1.357184 |
| M Phase |  1.377808 |
| Signaling by GPCR |  1.404715 |
| G2/M Checkpoints |  1.419820 |
| Cellular Senescence |  1.423657 |
| DNA Double-Strand Break Repair |  1.442063 |
| GPCR downstream signalling |  1.462783 |
| Signaling by Rho GTPases, Miro GTPases and RHOBTB3 |  1.463620 |
| Signaling by Rho GTPases |  1.487763 |
| RHOJ GTPase cycle |  1.514199 |
| trans-Golgi Network Vesicle Budding |  1.520059 |
| COPI-dependent Golgi-to-ER retrograde traffic |  1.535652 |
| RAC1 GTPase cycle |  1.540232 |
| Assembly of collagen fibrils and other multimeric structures |  1.547597 |
| G2/M DNA damage checkpoint |  1.553502 |
| Fcgamma receptor (FCGR) dependent phagocytosis |  1.557622 |
| RHO GTPase Effectors |  1.561024 |
| Kinesins |  1.563174 |
| Collagen formation |  1.587237 |
| FCGR3A-mediated IL10 synthesis |  1.591024 |
| SUMOylation of DNA damage response and repair proteins |  1.593563 |
| SUMOylation of DNA methylation proteins |  1.599738 |
| Diseases of mitotic cell cycle |  1.602410 |
| Transcription of E2F targets under negative control by DREAM complex |  1.606372 |
| RHOC GTPase cycle |  1.607240 |
| Sensory Perception |  1.611937 |
| Synthesis of bile acids and bile salts |  1.613874 |
| E2F mediated regulation of DNA replication |  1.615655 |
| Aquaporin-mediated transport |  1.620156 |
| Resolution of D-loop Structures through Holliday Junction Intermediates |  1.621492 |
| RND2 GTPase cycle |  1.621869 |
| Presynaptic phase of homologous DNA pairing and strand exchange |  1.627912 |
| G alpha (i) signalling events |  1.630012 |
| Arachidonic acid metabolism |  1.630762 |
| Leishmania infection |  1.631708 |
| Parasitic Infection Pathways |  1.631708 |
| VEGFR2 mediated cell proliferation |  1.633267 |
| RHOQ GTPase cycle |  1.635536 |
| SUMOylation |  1.639222 |
| Vasopressin regulates renal water homeostasis via Aquaporins |  1.642338 |
| Ca-dependent events |  1.644998 |
| RND3 GTPase cycle |  1.648496 |
| SUMO E3 ligases SUMOylate target proteins |  1.652404 |
| RHOU GTPase cycle |  1.653308 |
| Removal of the Flap Intermediate |  1.653557 |
| Aberrant regulation of mitotic cell cycle due to RB1 defects |  1.655203 |
| Defective pyroptosis |  1.659783 |
| Processive synthesis on the lagging strand |  1.664029 |
| G1/S Transition |  1.666620 |
| Cell Cycle |  1.670658 |
| PCNA-Dependent Long Patch Base Excision Repair |  1.672704 |
| PLC beta mediated events |  1.673461 |
| Resolution of D-Loop Structures |  1.682975 |
| Ion homeostasis |  1.684722 |
| Mitotic Prometaphase |  1.687451 |
| Chromosome Maintenance |  1.688812 |
| Cell Cycle Checkpoints |  1.690703 |
| Calmodulin induced events |  1.693557 |
| CaM pathway |  1.693557 |
| HDR through Homologous Recombination (HRR) |  1.694333 |
| EML4 and NUDC in mitotic spindle formation |  1.696535 |
| Telomere Maintenance |  1.697580 |
| HDR through Single Strand Annealing (SSA) |  1.699025 |
| RHO GTPases Activate Formins |  1.699854 |
| Homology Directed Repair |  1.706192 |
| Impaired BRCA2 binding to RAD51 |  1.708318 |
| Mitotic G1 phase and G1/S transition |  1.711530 |
| Cholesterol biosynthesis |  1.711996 |
| DNA Replication |  1.712553 |
| Opioid Signalling |  1.714581 |
| Glucagon-like Peptide-1 (GLP1) regulates insulin secretion |  1.717834 |
| HDR through Homologous Recombination (HRR) or Single Strand Annealing (SSA) |  1.719656 |
| Homologous DNA Pairing and Strand Exchange |  1.734755 |
| G0 and Early G1 |  1.738882 |
| Synthesis of DNA |  1.742812 |
| Diseases of DNA repair |  1.744355 |
| Anti-inflammatory response favouring Leishmania parasite infection |  1.753147 |
| Leishmania parasite growth and survival |  1.753147 |
| Resolution of AP sites via the multiple-nucleotide patch replacement pathway |  1.755019 |
| G-protein mediated events |  1.757053 |
| Mitotic Spindle Checkpoint |  1.764308 |
| Extension of Telomeres |  1.765192 |
| Regulation of insulin secretion |  1.772276 |
| Cell Cycle, Mitotic |  1.781071 |
| Polymerase switching |  1.784076 |
| Leading Strand Synthesis |  1.784076 |
| Resolution of Abasic Sites (AP sites) |  1.790265 |
| Amplification of signal from the kinetochores |  1.792241 |
| Amplification  of signal from unattached  kinetochores via a MAD2  inhibitory signal |  1.792241 |
| Diseases of DNA Double-Strand Break Repair |  1.795437 |
| Defective homologous recombination repair (HRR) due to BRCA2 loss of function |  1.795437 |
| G1/S-Specific Transcription |  1.797571 |
| Base Excision Repair |  1.809282 |
| Aberrant regulation of mitotic G1/S transition in cancer due to RB1 defects |  1.813502 |
| Defective binding of RB1 mutants to E2F1,(E2F2, E2F3) |  1.813502 |
| Resolution of Sister Chromatid Cohesion |  1.820693 |
| S Phase |  1.820888 |
| Lagging Strand Synthesis |  1.831390 |
| Activation of the pre-replicative complex |  1.833038 |
| DAG and IP3 signaling |  1.833079 |
| CDC42 GTPase cycle |  1.841011 |
| Activation of ATR in response to replication stress |  1.869308 |
| Unwinding of DNA |  1.886716 |
| Polymerase switching on the C-strand of the telomere |  1.888445 |
| Telomere C-strand (Lagging Strand) Synthesis |  1.891149 |
| DNA strand elongation |  2.134875 |



The pathways cover a wide range of biological processes, including cell  cycle regulation ("M Phase", "G2/M Checkpoints", "G1/S Transition"), DNA repair ("DNA Double-Strand Break Repair", "Base Excision Repair"),  signaling pathways ("Signaling by GPCR", "Interferon alpha/beta  signaling"), and metabolism ("Synthesis of bile acids and bile salts",  "Arachidonic acid metabolism").





### IV-LD

The only *Reactome* pathway is 'Neuronal System'.



## 4WPC

### DNA vaccine

Running GSEA with the whole gene list, without pre-selecting by p-value, yields the following downregulated pathways on *ReactomePA*.

| Description | NES |
| --- | --- |
| Interferon alpha/beta signaling | -2.245783 |
| Signaling by CSF1 (M-CSF) in myeloid cells | -2.183906 |
| GPVI-mediated activation cascade | -2.156946 |
| Interleukin-3, Interleukin-5 and GM-CSF signaling | -2.118311 |
| ROS and RNS production in phagocytes | -2.073070 |
| Antigen activates B Cell Receptor (BCR) leading to generation of second messengers | -2.063355 |
| RHO GTPases Activate NADPH Oxidases | -1.964348 |
| Activation of IRF3, IRF7 mediated by TBK1, IKKε (IKBKE) | -1.956229 |
| Inactivation of CSF3 (G-CSF) signaling | -1.950036 |
| Regulation of signaling by CBL | -1.915487 |
| Cytosolic tRNA aminoacylation | -1.915362 |
| Signaling by CSF3 (G-CSF) | -1.897145 |
| FCERI mediated MAPK activation | -1.881534 |
| Parasite infection | -1.788290 |
| Leishmania phagocytosis | -1.788290 |
| FCGR3A-mediated phagocytosis | -1.788290 |
| Generation of second messenger molecules | -1.780454 |
| Antigen processing-Cross presentation | -1.775111 |
| Signaling by the B Cell Receptor (BCR) | -1.749496 |
| Costimulation by the CD28 family | -1.702488 |
| Antigen Presentation: Folding, assembly and peptide loading of class I MHC | -1.690269 |
| Hh mutants are degraded by ERAD | -1.689402 |
| Hh mutants abrogate ligand secretion | -1.689402 |
| DAP12 interactions | -1.671403 |
| DAP12 signaling | -1.671403 |
| FCERI mediated NF-kB activation | -1.665623 |
| Defective CFTR causes cystic fibrosis | -1.665085 |
| Regulation of actin dynamics for phagocytic cup formation | -1.648707 |
| SCF-beta-TrCP mediated degradation of Emi1 | -1.615164 |
| Fcgamma receptor (FCGR) dependent phagocytosis | -1.603441 |
| Hedgehog ligand biogenesis | -1.602925 |
| TCR signaling | -1.595014 |
| Degradation of GLI2 by the proteasome | -1.585871 |
| Activation of NF-kappaB in B cells | -1.577746 |
| Fc epsilon receptor (FCERI) signaling | -1.544485 |
| Interleukin-1 family signaling | -1.528250 |
| Signaling by Interleukins | -1.498810 |
| Cytokine Signaling in Immune system | -1.396404 |
| Adaptive Immune System | -1.241338 |



Overall, these pathways highlight key processes in the immune response,  cell signaling, and regulation, emphasizing their roles in maintaining  immune system function and responding to pathogens. Negative NES values  indicate downregulation, suggesting a potential suppression of these  pathways under certain conditions.

On the opposite side of the spectrum, the upregulated pathways:

| Description | NES |
| --- | --- |
| Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins. | 2.166720 |
| Aerobic respiration and respiratory electron transport | 2.127555 |
| Respiratory electron transport | 2.114959 |
| Post-translational protein phosphorylation | 2.086596 |
| Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs) | 2.085496 |
| Formation of Fibrin Clot (Clotting Cascade) | 2.065645 |
| Intrinsic Pathway of Fibrin Clot Formation | 2.043202 |
| Complex I biogenesis | 1.923962 |
| Heme signaling | 1.898826 |
| Mitochondrial biogenesis | 1.863493 |
| Branched-chain amino acid catabolism | 1.830416 |
| Muscle contraction | 1.819388 |
| Potassium Channels | 1.815271 |
| Platelet degranulation | 1.802439 |
| Cristae formation | 1.799254 |
| Voltage gated Potassium channels | 1.788866 |
| Mitochondrial protein degradation | 1.783941 |
| Response to elevated platelet cytosolic Ca2+ | 1.782834 |
| Formation of ATP by chemiosmotic coupling | 1.775668 |
| Citric acid cycle (TCA cycle) | 1.772693 |
| Nuclear Receptor transcription pathway | 1.759821 |
| FOXO-mediated transcription of oxidative stress, metabolic and neuronal genes | 1.731914 |
| GABA receptor activation | 1.727551 |
| Striated Muscle Contraction | 1.712729 |
| Transcriptional activation of mitochondrial biogenesis | 1.683622 |
| FOXO-mediated transcription | 1.682785 |
| Collagen chain trimerization | 1.670550 |
| Metabolism of nitric oxide: NOS3 activation and regulation | 1.669765 |
| Gluconeogenesis | 1.665098 |
| Opioid Signalling | 1.652053 |
| Sensory Perception | 1.651697 |
| Collagen biosynthesis and modifying enzymes | 1.622587 |
| Neuronal System | 1.596390 |
| Cardiac conduction | 1.584142 |
| Organelle biogenesis and maintenance | 1.496958 |
| Hemostasis | 1.421507 |



Overall, these pathways collectively impact:

- Cellular energy production and metabolism.
- Cell signaling and regulatory mechanisms.
- Blood clotting and cardiovascular health.
- Muscle function and movement.
- Structural maintenance and tissue integrity.
- Protein synthesis, modification, and degradation.
- Growth, development, and response to environmental stimuli.



If one starts from a list of genes pre-selected by p-value (*p < 0.05*), and runs GSEA against a genome wide annotation for human (org.Hs.eg.db), the overall result is the same, with immune-related pathways downregulated and energy production upregulated.



| Description | NES |
| --- | --- |
| aerobic respiration |  3.157926 |
| oxidative phosphorylation |  2.937657 |
| proton motive force-driven mitochondrial ATP synthesis |  2.624725 |
| dicarboxylic acid metabolic process |  2.361743 |
| mitochondrial respiratory chain complex assembly |  2.300209 |
| muscle system process |  2.252134 |
| regulation of muscle system process |  2.094244 |
| rhythmic process |  2.085570 |
| nucleoside triphosphate biosynthetic process |  1.995967 |
| ribonucleoside triphosphate biosynthetic process |  1.995967 |
| amino acid catabolic process |  1.936152 |
| nucleoside triphosphate metabolic process |  1.881641 |
| ribonucleoside triphosphate metabolic process |  1.881641 |
| ATP metabolic process |  1.810749 |
| small molecule catabolic process |  1.722197 |
| purine-containing compound metabolic process |  1.656750 |
| cellular response to oxygen-containing compound |  1.647138 |
| small molecule metabolic process |  1.549061 |
| multicellular organismal process | -1.425700 |
| cell communication | -1.510842 |
| signaling | -1.521286 |
| positive regulation of biological process | -1.575174 |
| negative regulation of macromolecule biosynthetic process | -1.634718 |
| negative regulation of biosynthetic process | -1.663773 |
| negative regulation of cellular biosynthetic process | -1.663773 |
| import into cell | -1.684992 |
| lymphocyte apoptotic process | -1.774158 |
| regulation of leukocyte mediated cytotoxicity | -1.787597 |
| positive regulation of cell migration | -1.792857 |
| positive regulation of locomotion | -1.792857 |
| positive regulation of cell motility | -1.792857 |
| locomotion | -1.801708 |
| positive regulation of cell development | -1.812092 |
| regulation of protein phosphorylation | -1.822433 |
| negative regulation of T cell activation | -1.848189 |
| collagen metabolic process | -1.854287 |
| leukocyte homeostasis | -1.856544 |
| response to cytokine | -1.858390 |
| regulation of cell killing | -1.860695 |
| regulated exocytosis | -1.866868 |
| positive regulation of catalytic activity | -1.872830 |
| regulation of chemotaxis | -1.875409 |
| negative regulation of protein transport | -1.894241 |
| negative regulation of establishment of protein localization | -1.894241 |
| peptidyl-tyrosine modification | -1.894578 |
| import across plasma membrane | -1.895777 |
| regulation of leukocyte migration | -1.897036 |
| response to external stimulus | -1.933134 |
| positive regulation of phosphorylation | -1.949875 |
| phagocytosis, engulfment | -1.953162 |
| membrane invagination | -2.007528 |
| plasma membrane invagination | -2.007528 |
| regulation of binding | -2.009462 |
| leukocyte mediated cytotoxicity | -2.029608 |
| epithelial to mesenchymal transition | -2.038571 |
| cytokine-mediated signaling pathway | -2.057145 |
| steroid biosynthetic process | -2.079742 |
| skeletal system development | -2.103461 |
| production of molecular mediator of immune response | -2.113634 |
| regulation of production of molecular mediator of immune response | -2.113634 |
| T cell mediated immunity | -2.115681 |
| positive regulation of response to stimulus | -2.123837 |
| positive regulation of multicellular organismal process | -2.132066 |
| regulation of ERK1 and ERK2 cascade | -2.137176 |
| regulation of innate immune response | -2.137273 |
| cell killing | -2.172378 |
| leukocyte activation involved in immune response | -2.183926 |
| biological process involved in interspecies interaction between organisms | -2.186225 |
| regulation of lymphocyte differentiation | -2.199522 |
| phagocytosis | -2.208789 |
| regulation of defense response | -2.223763 |
| chemotaxis | -2.247164 |
| taxis | -2.247164 |
| response to biotic stimulus | -2.281739 |
| leukocyte migration | -2.302119 |
| cell chemotaxis | -2.305048 |
| positive regulation of T cell activation | -2.309675 |
| regulation of cell adhesion | -2.366362 |
| leukocyte cell-cell adhesion | -2.372385 |
| cytokine production | -2.409554 |
| regulation of cytokine production | -2.409554 |
| regulation of mononuclear cell proliferation | -2.450445 |
| defense response to other organism | -2.457111 |
| lymphocyte proliferation | -2.480956 |
| mononuclear cell differentiation | -2.565588 |
| leukocyte proliferation | -2.596309 |
| inflammatory response | -2.599778 |
| immune system process | -2.732992 |
| cell activation | -2.764861 |
| lymphocyte activation | -2.795845 |
| regulation of immune system process | -2.824275 |
| positive regulation of immune system process | -2.826171 |
| regulation of leukocyte activation | -2.841657 |
| T cell activation | -2.872483 |
| regulation of immune response | -2.915157 |
| leukocyte activation | -2.962128 |
| immune response | -3.085458 |

### EOMES




From ChatGPT when prompted to explain a negative *Normalised Enrichment Score* (NES):

*In summary, a negative NES in GSEA indicates that the gene set is  significantly downregulated in the condition of interest compared to the control, reflecting a negative correlation with the phenotype being  studied. This information can provide valuable insights into the  biological processes that are less active or suppressed in the  experimental condition.*



### Interpretation of GSEA Results:

- **Positive NES:**
  - Indicates that the gene set is enriched at the top of the ranked list, suggesting that genes in the set are upregulated in the condition of interest.
- **Negative NES:**
  - Indicates that the gene set is enriched at the bottom of the ranked list, suggesting that genes in the set are downregulated in the condition of interest.







Most studies related to viral infection use qPCR to 'measure' gene expression of target genes in relation to a *housekeeping gene*. The approach of using RNAseq as a tool to study immune responses is inherently different, in the sense that it does not compare gene expression against a housekeeping gene. Instead, it requires using a *reference*, such as a control group. This control group provides the baseline gene expression against which all other groups are compared to.

Because with RNAseq one looks at thousands of genes simultaneously, it makes more sense to have a holistic view of gene expression. Any emerging patterns in gene expression tend to stand out during pathway analysis, since pathways, by definition, incorporate hundreds, or even thousands, of differentially regulated genes.



WHAT GENES ARE EXPECTED TO BE DIFFERENTIALLY REGULATED BY EOMES, GATA3, AND DNA VACCINE?

CREATE A TIMELINE PLOT OF GENE EXPRESSION, I.E. 1WPC VS 10WPI, 4WPC VS 1WPC, 6WPC VS 4WPC, 10WPC VS 6WPC FOR CONTROL



Immune-related genes

Housekeeping: β-actin (ACTB): AF012125, ENSSSAG00000116649. A lot of reads everywhere.

SasaIFN-α: AY216594, no Ensembl

Μx: U66475, ENSSSAG00070007970 (mxa). No hits in read files.

ISG15: AY926456, ENSSSAG00000105254 (UBIL). Mapped but mostly no reads.

IFN-γ: AY795563, ENSSSAG00000105299. Mapped but no reads.

ΤΝF-α (1&2): NM_001123589 & NM_001123590, ENSSSAG00000065312 (mapped but no reads) & ENSSSAG00000053783 (no hits)

IL12-β: BT049114, ENSSSAG00000009655. Mapped but with few reads.

IL-10: EF165029, ENSSSAG00000107544. Mapped but with no reads.

IL-8: NM_001140710, ENSSSAG00000006498. Mapped but with no reads.

CD3-ε: NM_001123622, ENSSSAG00000076824. A lot of reads in s and hk. Fewer reads in h.

CD4: NM_001124539, rainbow trout gene

CD8-α: NM_001123583, ENSSSAG00000065860. Mapped but with few reads.

TCR-α: BT050114, no Ensembl. ENSSSAG00070030723 for T-cell receptor alpha/delta variable 31.0. No hits.

MHC I: AF508864, no Ensembl. ENSSSAG00000077407 for a novel gene. A lot of hits in h, hk, and s. 

MHC II: BT049430, no Ensembl. ENSSSAG00000004635 for a novel gene. A LOT of hits everywhere.



How did you look at specific genes in the RNAseq data (genes in Table 1, for example)? Are these manually curated, and you get their respective raw counts from gene count files?

So these were manually curated and then extracted into a table where we included the FC values from Deseq2 and then used total counts and the house keeping gene EF1a to make a normalised expression. It's maybe not the best way of doing it, but was a useful proxy of the genes expression. (Shahmir Naseer)



About normalising expression

Expression values were obtained by normalising gene count against total raw counts and then the ratio to ELF1a multiplied by 1000, trout (ENSOMYG00000038328) and salmon (ENSSSAG00000077892) counts,

blank cells = no counts for the gene.
