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



![dnavaccine_4wpc](/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_4wpc/gsea_definitive_plots/dnavaccine_4wpc.png)



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



If one starts from a list of genes pre-selected by p-value (*p < 0.05*), and runs GSEA against a genome wide annotation for human (not ReactomePA, org.Hs.eg.db), the overall result is the same, with immune-related pathways downregulated and energy production upregulated.

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

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| Aerobic respiration and respiratory electron transport |  1.949627 | 164 | 106 |
| Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins. |  1.895979 |  98 |  68 |
| Formation of definitive endoderm |  1.894218 |  10 |   4 |
| Respiratory electron transport |  1.884847 |  83 |  57 |
| Complex I biogenesis |  1.784491 |  44 |  32 |
| Mitochondrial biogenesis |  1.754441 |  79 |  41 |
| Muscle contraction |  1.745952 | 117 |  60 |
| NCAM1 interactions |  1.719882 |  16 |   7 |
| Adaptive Immune System | -1.363424 | 491 | 125 |
| Cell Cycle | -1.394166 | 496 | 155 |
| Cytokine Signaling in Immune system | -1.425435 | 424 | 143 |
| Neutrophil degranulation | -1.470707 | 298 |  84 |
| Cell Cycle Checkpoints | -1.518823 | 211 | 107 |
| Cholesterol biosynthesis | -1.777279 |  20 |   9 |
| ABC transporter disorders | -1.902475 |  62 |  33 |
| Antigen processing-Cross presentation | -1.915112 |  75 |  32 |
| Diseases of Immune System | -1.920536 |  15 |   6 |
| Diseases associated with the TLR signaling cascade | -1.920536 |  15 |   6 |
| Inwardly rectifying K+ channels | -1.944344 |  14 |   4 |
| Chromosome Maintenance | -1.961067 |  77 |  29 |
| RHO GTPases Activate NADPH Oxidases | -1.973060 |  19 |   7 |
| Regulation of Complement cascade | -2.085535 |  17 |   7 |
| Interferon alpha/beta signaling | -2.096606 |  34 |  18 |
| Potassium Channels | -2.124389 |  34 |   6 |



#### Key regulated pathways

##### Upregulated

- Aerobic respiration and respiratory electron transport:
  - This pathway is highly upregulated. It involves processes crucial for energy production via aerobic respiration, indicating an increase in metabolic activity.
- Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins
  - This pathway is also significantly upregulated. It points to active ATP synthesis and possible thermogenesis through uncoupling proteins, suggesting heightened mitochondrial activity.
- Respiratory electron transport
  - Similar to the first pathway, this one focuses on electron transport, further emphasizing the increased mitochondrial function and energy production.
- Complex I biogenesis
  - Upregulation of this pathway indicates an increase in the formation of *Complex I* of the mitochondrial electron transport chain, supporting enhanced respiration and ATP production.
- Mitochondrial biogenesis
  - Increased biogenesis of mitochondria, which correlates with the upregulation of respiratory pathways and overall cellular energy production.



##### Downregulated

- Adaptive immune system
  - Downregulation suggests a suppressed adaptive immune response, potentially indicating an immune evasion mechanism or reduced immune activity.
- Cell cycle
  - Suggests a reduction in cell proliferation activities, possibly indicating a shift from cell division to other cellular activities.
- Cytokine signaling in immune system
  - Downregulation here indicates reduced cytoking-mediated signaling, which could be part of broader suppression of immune responses.
- Neutrophil degranulation
  - Indicative of reduced neutrophil activity, which play a key role in innate immunity.
- Antigen processing-cross presentation
  - Reduced antigen presentation capability, further suggesting immune suppression.

- Diseases of Immune System and Diseases associated with the TLR signaling cascade
- RHO GTPases activate NADPH oxidases
  - Reduced activity in oxidative stress responses mediated by RHO GTPases.

- Regulation of Complement cascade
  - Deregulation of complement regulation could impair immune complex clearance and inflammatory response.

- Interferon alpha/beta signaling
  - Reduced interferon signaling, affecting antiviral responses and immune modulation.


### Overall Interpretation

- **Upregulation**: The upregulated pathways largely involve mitochondrial activity, energy production, and muscle contraction. This suggests a state of increased metabolic activity and energy demand.
- **Downregulation**: The downregulated pathways are predominantly associated with immune responses, cell cycle processes, and certain biosynthetic pathways like cholesterol synthesis. This indicates a suppression of immune activity and cell proliferation, possibly in favor of energy conservation or other cellular priorities.

These changes might reflect a specific physiological or pathological state where the organism is prioritizing energy production and mitochondrial function over immune responses and cell division. This could be relevant in contexts such as metabolic adaptation, stress responses, or certain disease states.





### GATA3

| Description | NES | setSize | Count |
| :-- | --- | --- | --- |
| Aerobic respiration and respiratory electron transport |  2.204298 | 164 | 104 |
| Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins. |  2.082230 |  98 |  77 |
| Collagen chain trimerization |  2.024129 |  18 |   5 |
| Respiratory electron transport |  2.002613 |  83 |  65 |
| Collagen biosynthesis and modifying enzymes |  1.906599 |  35 |   6 |
| Formation of definitive endoderm |  1.899582 |  10 |   1 |
| Muscle contraction |  1.857070 | 115 |  39 |
| Complex I biogenesis |  1.837222 |  44 |  36 |
| Mitochondrial biogenesis |  1.795692 |  79 |  37 |
| Pyruvate metabolism |  1.793286 |  40 |  19 |
| Regulation of pyruvate metabolism |  1.775754 |  29 |  12 |
| Cristae formation |  1.757520 |  24 |  18 |
| Eukaryotic Translation Elongation |  1.738635 |  73 |  53 |
| Estrogen-dependent nuclear events downstream of ESR-membrane signaling |  1.738402 |  17 |  11 |
| Citric acid cycle (TCA cycle) |  1.737412 |  30 |  17 |
| Peptide chain elongation |  1.730130 |  71 |  52 |
| Formation of ATP by chemiosmotic coupling |  1.729062 |  12 |  10 |
| Mitochondrial protein degradation |  1.712164 |  80 |  49 |
| FOXO-mediated transcription of oxidative stress, metabolic and neuronal genes |  1.710878 |  17 |   9 |
| Striated Muscle Contraction |  1.701422 |  25 |   9 |
| FOXO-mediated transcription |  1.697493 |  46 |  21 |
| TP53 Regulates Transcription of Cell Death Genes |  1.663520 |  27 |   6 |
| Formation of a pool of free 40S subunits |  1.651648 |  83 |  56 |
| Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC) |  1.648478 |  75 |  52 |
| Aberrant regulation of mitotic cell cycle due to RB1 defects |  1.647869 |  30 |  15 |
| Viral mRNA Translation |  1.634192 |  70 |  50 |
| Platelet Adhesion to exposed collagen |  1.631415 |  10 |   4 |
| Nonsense-Mediated Decay (NMD) |  1.617930 |  92 |  60 |
| Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC) |  1.617930 |  92 |  60 |
| Eukaryotic Translation Initiation |  1.597933 |  98 |  61 |
| Cap-dependent Translation Initiation |  1.597933 |  98 |  61 |
| L13a-mediated translational silencing of Ceruloplasmin expression |  1.577314 |  91 |  58 |
| GTP hydrolysis and joining of the 60S ribosomal subunit |  1.565948 |  92 |  58 |
| Cardiac conduction |  1.557649 |  65 |  24 |
| SRP-dependent cotranslational protein targeting to membrane |  1.556804 |  89 |  55 |
| Selenocysteine synthesis |  1.554998 |  75 |  54 |
| Eukaryotic Translation Termination |  1.527636 |  74 |  52 |
| Organelle biogenesis and maintenance |  1.518678 | 218 |  61 |
| Factors involved in megakaryocyte development and platelet production |  1.497255 | 102 |  20 |
| Adaptive Immune System | -1.272833 | 486 |  96 |
| GPCR downstream signalling | -1.338223 | 244 |  44 |
| Hemostasis | -1.343668 | 375 |  75 |
| Signaling by GPCR | -1.375192 | 271 |  47 |
| Signaling by Interleukins | -1.383535 | 264 |  84 |
| Metabolism of steroids | -1.431832 |  95 |  26 |
| Interferon Signaling | -1.444722 | 147 |  41 |
| Signaling by the B Cell Receptor (BCR) | -1.469307 |  92 |  34 |
| SARS-CoV-2 activates/modulates innate and adaptive immune responses | -1.469872 |  74 |  18 |
| Cytokine Signaling in Immune system | -1.474240 | 417 | 120 |
| G alpha (q) signalling events | -1.511764 |  89 |  16 |
| Interleukin-1 family signaling | -1.515449 | 104 |  34 |
| Neutrophil degranulation | -1.522535 | 296 |  83 |
| Disorders of transmembrane transporters | -1.565124 | 120 |  35 |
| Plasma lipoprotein assembly, remodeling, and clearance | -1.577641 |  43 |  12 |
| TCR signaling | -1.592106 |  83 |  31 |
| Toll Like Receptor 4 (TLR4) Cascade | -1.601942 |  93 |  22 |
| GPCR ligand binding | -1.628377 | 112 |  26 |
| ABC transporter disorders | -1.679921 |  62 |  21 |
| Toll-like Receptor Cascades | -1.706463 | 109 |  27 |
| Class A/1 (Rhodopsin-like receptors) | -1.706959 |  80 |  21 |
| tRNA Aminoacylation | -1.714003 |  40 |  14 |
| ER-Phagosome pathway | -1.738654 |  63 |  24 |
| Post-translational protein phosphorylation | -1.754338 |  61 |  18 |
| Signaling by FGFR1 in disease | -1.762609 |  24 |  15 |
| Plasma lipoprotein remodeling | -1.773702 |  18 |   6 |
| Binding and Uptake of Ligands by Scavenger Receptors | -1.788881 |  24 |   6 |
| Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs) | -1.794771 |  69 |  23 |
| Antigen activates B Cell Receptor (BCR) leading to generation of second messengers | -1.820756 |  26 |  13 |
| FCERI mediated MAPK activation | -1.833991 |  25 |  11 |
| Platelet activation, signaling and aggregation | -1.844287 | 179 |  39 |
| Scavenging by Class A Receptors | -1.849698 |  13 |   3 |
| Signaling by CSF1 (M-CSF) in myeloid cells | -1.864712 |  26 |  13 |
| Response to elevated platelet cytosolic Ca2+ | -1.989799 |  90 |  23 |
| Cytosolic tRNA aminoacylation | -2.003514 |  23 |  12 |
| Antigen processing-Cross presentation | -2.005494 |  75 |  33 |
| Platelet degranulation | -2.005702 |  86 |  22 |
| GPVI-mediated activation cascade | -2.009896 |  26 |  11 |
| ROS and RNS production in phagocytes | -2.063678 |  22 |  11 |
| Peptide ligand-binding receptors | -2.074588 |  39 |  11 |
| RHO GTPases Activate NADPH Oxidases | -2.092266 |  19 |   8 |
| Diseases of Immune System | -2.198914 |  14 |   8 |
| Diseases associated with the TLR signaling cascade | -2.198914 |  14 |   8 |
| Complement cascade | -2.263055 |  20 |  10 |
| Regulation of Complement cascade | -2.307570 |  16 |   7 |
| Interferon alpha/beta signaling | -2.342834 |  33 |  17 |















### IV-HD

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins. |  2.373994 |  98 |  68 |
| Aerobic respiration and respiratory electron transport |  2.360927 | 164 | 104 |
| Respiratory electron transport |  2.283299 |  83 |  56 |
| Muscle contraction |  2.067722 | 112 |  49 |
| Collagen chain trimerization |  2.041581 |  18 |   9 |
| Collagen biosynthesis and modifying enzymes |  2.029762 |  35 |  11 |
| Complex I biogenesis |  2.018103 |  44 |  31 |
| Mitochondrial biogenesis |  1.996827 |  77 |  51 |
| Mitochondrial protein degradation |  1.921150 |  79 |  44 |
| Cristae formation |  1.903527 |  24 |  19 |
| Formation of ATP by chemiosmotic coupling |  1.872484 |  12 |  11 |
| Cardiac conduction |  1.842621 |  62 |  25 |
| Nuclear Receptor transcription pathway |  1.828501 |  34 |  16 |
| Citric acid cycle (TCA cycle) |  1.817385 |  30 |  18 |
| Striated Muscle Contraction |  1.809174 |  25 |  12 |
| Glyoxylate metabolism and glycine degradation |  1.784874 |  14 |   9 |
| Signaling by Retinoic Acid |  1.776750 |  26 |  15 |
| Branched-chain amino acid catabolism |  1.775862 |  19 |  18 |
| Pyruvate metabolism |  1.740393 |  40 |  20 |
| Maturation of TCA enzymes and regulation of TCA cycle |  1.733603 |  17 |  11 |
| Signaling by BMP |  1.728433 |  20 |   9 |
| Collagen formation |  1.722355 |  50 |  14 |
| Lysine catabolism |  1.688240 |  12 |   7 |
| Defective B3GALTL causes PpS |  1.687668 |  18 |  11 |
| FOXO-mediated transcription |  1.684639 |  46 |  22 |
| Processing of SMDT1 |  1.680166 |  16 |   3 |
| Regulation of pyruvate dehydrogenase (PDH) complex |  1.679770 |  13 |  10 |
| Transcriptional activation of mitochondrial biogenesis |  1.673736 |  46 |  21 |
| Mitochondrial protein import |  1.656725 |  50 |  35 |
| FOXO-mediated transcription of oxidative stress, metabolic and neuronal genes |  1.654292 |  17 |   9 |
| Myogenesis |  1.650532 |  21 |  16 |
| Triglyceride metabolism |  1.639665 |  19 |  11 |
| Assembly of collagen fibrils and other multimeric structures |  1.632868 |  32 |  11 |
| Ion homeostasis |  1.628526 |  33 |  13 |
| O-glycosylation of TSR domain-containing proteins |  1.618838 |  19 |  11 |
| Regulation of pyruvate metabolism |  1.605760 |  29 |  15 |
| ADORA2B mediated anti-inflammatory cytokines production |  1.593179 |  25 |  12 |
| Diseases associated with O-glycosylation of proteins |  1.590013 |  28 |  15 |
| Peroxisomal lipid metabolism |  1.588763 |  21 |  13 |
| Regulation of TP53 Activity through Association with Co-factors |  1.588492 |  10 |   7 |
| Receptor-type tyrosine-protein phosphatases |  1.577917 |  14 |   5 |
| Factors involved in megakaryocyte development and platelet production |  1.575224 | 102 |  28 |
| Extra-nuclear estrogen signaling |  1.564253 |  47 |  29 |
| Estrogen-dependent nuclear events downstream of ESR-membrane signaling |  1.558617 |  16 |  13 |
| Phase 0 - rapid depolarisation |  1.554122 |  14 |   9 |
| Metabolism of nitric oxide: NOS3 activation and regulation |  1.544521 |  11 |   7 |
| Transcriptional Regulation by MECP2 |  1.539376 |  40 |  15 |
| Collagen degradation |  1.534987 |  33 |  10 |
| RUNX1 regulates genes involved in megakaryocyte differentiation and platelet function |  1.532652 |  34 |  11 |
| Transmission across Chemical Synapses |  1.530830 | 125 |  44 |
| G alpha (s) signalling events |  1.528299 |  55 |  16 |
| Energy dependent regulation of mTOR by LKB1-AMPK |  1.527262 |  24 |  17 |
| Heme biosynthesis |  1.506504 |  12 |   6 |
| Protein localization |  1.496549 | 126 |  70 |
| Signaling by Nuclear Receptors |  1.485743 | 165 |  60 |
| ESR-mediated signaling |  1.403965 | 123 |  50 |
| Adaptive Immune System | -1.266639 | 476 | 118 |
| Mitotic Metaphase and Anaphase | -1.312167 | 184 |  61 |
| Cytokine Signaling in Immune system | -1.345557 | 406 | 130 |
| SARS-CoV-2-host interactions | -1.351659 | 136 |  40 |
| Platelet activation, signaling and aggregation | -1.392250 | 177 |  46 |
| Downstream signaling events of B Cell Receptor (BCR) | -1.398520 |  67 |  25 |
| Interferon Signaling | -1.409395 | 145 |  52 |
| APC/C:Cdc20 mediated degradation of Securin | -1.434706 |  58 |  24 |
| TNFR2 non-canonical NF-kB pathway | -1.436748 |  57 |  24 |
| Neutrophil degranulation | -1.441234 | 293 | 111 |
| Regulation of PTEN stability and activity | -1.447049 |  60 |  23 |
| Respiratory Syncytial Virus Infection Pathway | -1.447476 |  75 |  21 |
| Dual incision in TC-NER | -1.450489 |  53 |  25 |
| Fcgamma receptor (FCGR) dependent phagocytosis | -1.464650 |  73 |  21 |
| APC/C:Cdh1 mediated degradation of Cdc20 and other APC/C:Cdh1 targeted proteins in late mitosis/early G1 | -1.471462 |  64 |  27 |
| Presynaptic phase of homologous DNA pairing and strand exchange | -1.482260 |  29 |  21 |
| Fc epsilon receptor (FCERI) signaling | -1.483509 | 105 |  39 |
| Interleukin-1 family signaling | -1.489433 | 102 |  40 |
| SARS-CoV-2 modulates host translation machinery | -1.494387 |  39 |  24 |
| Regulation of RAS by GAPs | -1.498015 |  56 |  22 |
| CLEC7A (Dectin-1) signaling | -1.503632 |  78 |  29 |
| Platelet degranulation | -1.510700 |  84 |  26 |
| Toll Like Receptor 4 (TLR4) Cascade | -1.512454 |  90 |  22 |
| Cdc20:Phospho-APC/C mediated degradation of Cyclin A | -1.513434 |  63 |  26 |
| Dectin-1 mediated noncanonical NF-kB signaling | -1.518192 |  52 |  22 |
| RSV-host interactions | -1.521409 |  57 |  15 |
| Switching of origins to a post-replicative state | -1.529257 |  79 |  32 |
| Chromosome Maintenance | -1.539670 |  74 |  25 |
| FCERI mediated Ca+2 mobilization | -1.540980 |  24 |  10 |
| TNF signaling | -1.541100 |  45 |  12 |
| Interferon gamma signaling | -1.542065 |  38 |   9 |
| DNA Replication Pre-Initiation | -1.545875 |  89 |  35 |
| APC:Cdc20 mediated degradation of cell cycle proteins prior to satisfation of the cell cycle checkpoint | -1.547625 |  64 |  27 |
| Cytosolic tRNA aminoacylation | -1.547764 |  23 |  12 |
| Signaling by FGFR1 in disease | -1.552353 |  23 |  13 |
| Ubiquitin-dependent degradation of Cyclin D | -1.557170 |  46 |  20 |
| Costimulation by the CD28 family | -1.562173 |  42 |  14 |
| Regulation of ornithine decarboxylase (ODC) | -1.576657 |  46 |  21 |
| UCH proteinases | -1.576768 |  74 |  26 |
| Response of Mtb to phagocytosis | -1.577899 |  16 |   3 |
| APC/C:Cdc20 mediated degradation of mitotic proteins | -1.580759 |  65 |  28 |
| Disorders of transmembrane transporters | -1.581266 | 117 |  44 |
| HDR through Single Strand Annealing (SSA) | -1.585122 |  29 |  19 |
| Binding and Uptake of Ligands by Scavenger Receptors | -1.591402 |  24 |   6 |
| Cell Cycle Checkpoints | -1.591653 | 209 |  77 |
| AUF1 (hnRNP D0) binds and destabilizes mRNA | -1.598520 |  47 |  20 |
| Regulation of signaling by CBL | -1.598687 |  20 |   8 |
| MAP2K and MAPK activation | -1.605994 |  29 |   9 |
| Ovarian tumor domain proteases | -1.610943 |  28 |   9 |
| Signaling by FGFR3 | -1.611010 |  23 |   5 |
| Signaling by FGFR4 | -1.611010 |  23 |   5 |
| Dual Incision in GG-NER | -1.615932 |  33 |  17 |
| Processing of DNA double-strand break ends | -1.616790 |  50 |  16 |
| Downstream TCR signaling | -1.616901 |  67 |  26 |
| Regulation of actin dynamics for phagocytic cup formation | -1.617697 |  50 |  16 |
| Regulation of Apoptosis | -1.619726 |  46 |  20 |
| Activation of APC/C and APC/C:Cdc20 mediated degradation of mitotic proteins | -1.620192 |  66 |  29 |
| Stabilization of p53 | -1.626019 |  48 |  20 |
| FCERI mediated NF-kB activation | -1.628735 |  64 |  24 |
| Removal of the Flap Intermediate | -1.633338 |  13 |  12 |
| Regulation of APC/C activators between G1/S and early anaphase | -1.634746 |  70 |  30 |
| Translesion synthesis by Y family DNA polymerases bypasses lesions on DNA template | -1.634863 |  31 |  16 |
| G2/M DNA damage checkpoint | -1.636487 |  49 |  26 |
| Activation of NF-kappaB in B cells | -1.641091 |  56 |  24 |
| Fanconi Anemia Pathway | -1.646247 |  32 |  17 |
| Glycosphingolipid metabolism | -1.648253 |  35 |  18 |
| Ubiquitin Mediated Degradation of Phosphorylated Cdc25A | -1.652383 |  45 |  20 |
| p53-Independent DNA Damage Response | -1.652383 |  45 |  20 |
| p53-Independent G1/S DNA damage checkpoint | -1.652383 |  45 |  20 |
| Host Interactions of HIV factors | -1.661886 | 105 |  40 |
| Signaling by high-kinase activity BRAF mutants | -1.663672 |  26 |   9 |
| SARS-CoV-2 activates/modulates innate and adaptive immune responses | -1.664612 |  74 |  25 |
| APC/C-mediated degradation of cell cycle proteins | -1.666051 |  76 |  33 |
| Regulation of mitotic cell cycle | -1.666051 |  76 |  33 |
| FCERI mediated MAPK activation | -1.666453 |  25 |   9 |
| Lagging Strand Synthesis | -1.667606 |  18 |  15 |
| Processive synthesis on the lagging strand | -1.670080 |  14 |  13 |
| DNA Replication | -1.670569 | 113 |  44 |
| NRIF signals cell death from the nucleus | -1.670941 |  11 |   4 |
| Termination of translesion DNA synthesis | -1.672201 |  25 |  19 |
| NIK-->noncanonical NF-kB signaling | -1.672782 |  51 |  22 |
| Cross-presentation of soluble exogenous antigens (endosomes) | -1.673384 |  42 |  20 |
| Degradation of AXIN | -1.676835 |  48 |  20 |
| Toll-like Receptor Cascades | -1.678746 | 106 |  30 |
| TCR signaling | -1.692458 |  83 |  34 |
| Class A/1 (Rhodopsin-like receptors) | -1.702797 |  73 |  26 |
| Parasite infection | -1.718058 |  51 |  17 |
| Leishmania phagocytosis | -1.718058 |  51 |  17 |
| FCGR3A-mediated phagocytosis | -1.718058 |  51 |  17 |
| MAP3K8 (TPL2)-dependent MAPK1/3 activation | -1.718704 |  14 |   4 |
| FBXL7 down-regulates AURKA during mitotic entry and in early mitosis | -1.727324 |  48 |  22 |
| Regulation of activated PAK-2p34 by proteasome mediated degradation | -1.728387 |  44 |  20 |
| Unwinding of DNA | -1.729156 |  11 |  10 |
| Synthesis of DNA | -1.729487 | 104 |  42 |
| Negative regulators of DDX58/IFIH1 signaling | -1.732429 |  25 |   9 |
| Degradation of DVL | -1.733938 |  50 |  21 |
| Negative regulation of FGFR3 signaling | -1.734380 |  14 |   3 |
| Negative regulation of FGFR4 signaling | -1.734380 |  14 |   3 |
| GSK3B and BTRC:CUL1-mediated-degradation of NFE2L2 | -1.737255 |  46 |  20 |
| Autodegradation of the E3 ubiquitin ligase COP1 | -1.739149 |  46 |  20 |
| p130Cas linkage to MAPK signaling for integrins | -1.749088 |  12 |   5 |
| Impaired BRCA2 binding to RAD51 | -1.749773 |  26 |  20 |
| Degradation of GLI1 by the proteasome | -1.757204 |  50 |  22 |
| Plasma lipoprotein remodeling | -1.764388 |  17 |   9 |
| Suppression of phagosomal maturation | -1.764965 |  10 |   3 |
| Signaling by FGFR3 in disease | -1.768984 |  10 |   4 |
| Alpha-protein kinase 1 signaling pathway | -1.773696 |  10 |   3 |
| GLI3 is processed to GLI3R by the proteasome | -1.776562 |  51 |  22 |
| Diseases associated with glycosylation precursor biosynthesis | -1.777302 |  13 |  11 |
| STING mediated induction of host immune responses | -1.785335 |  10 |   5 |
| Signaling by the B Cell Receptor (BCR) | -1.790076 |  92 |  40 |
| Defective CFTR causes cystic fibrosis | -1.790669 |  52 |  23 |
| Vif-mediated degradation of APOBEC3G | -1.794264 |  47 |  21 |
| Vpu mediated degradation of CD4 | -1.795407 |  45 |  21 |
| Evasion by RSV of host interferon responses | -1.795968 |  17 |   8 |
| GRB2:SOS provides linkage to MAPK signaling for Integrins | -1.801383 |  12 |   5 |
| Activation of ATR in response to replication stress | -1.801752 |  32 |  22 |
| Regulation of RUNX3 expression and activity | -1.806489 |  47 |  22 |
| Degradation of GLI2 by the proteasome | -1.812854 |  50 |  22 |
| Amyloid fiber formation | -1.813502 |  35 |  14 |
| Translesion Synthesis by POLH | -1.818860 |  16 |  12 |
| TRAF3-dependent IRF activation pathway | -1.831612 |  10 |   7 |
| Integrin signaling | -1.840796 |  24 |   5 |
| RHO GTPases Activate NADPH Oxidases | -1.851174 |  18 |   6 |
| GPVI-mediated activation cascade | -1.851624 |  26 |  13 |
| Peptide ligand-binding receptors | -1.858229 |  34 |  14 |
| Diseases of Immune System | -1.860082 |  12 |   7 |
| Diseases associated with the TLR signaling cascade | -1.860082 |  12 |   7 |
| DAP12 interactions | -1.861983 |  18 |   8 |
| DAP12 signaling | -1.861983 |  18 |   8 |
| Hedgehog ligand biogenesis | -1.873010 |  51 |  24 |
| Formation of Fibrin Clot (Clotting Cascade) | -1.876662 |  17 |   8 |
| ER-Phagosome pathway | -1.879640 |  61 |  25 |
| Activation of Matrix Metalloproteinases | -1.882756 |  12 |   8 |
| Hh mutants are degraded by ERAD | -1.883445 |  48 |  23 |
| Hh mutants abrogate ligand secretion | -1.883445 |  48 |  23 |
| SCF-beta-TrCP mediated degradation of Emi1 | -1.888642 |  48 |  23 |
| Orc1 removal from chromatin | -1.895336 |  62 |  29 |
| G2/M Checkpoints | -1.908061 | 112 |  60 |
| Inactivation of CSF3 (G-CSF) signaling | -1.909885 |  20 |  12 |
| ABC transporter disorders | -1.934190 |  62 |  28 |
| Antigen activates B Cell Receptor (BCR) leading to generation of second messengers | -1.935602 |  26 |  13 |
| Signaling by CSF3 (G-CSF) | -1.948166 |  25 |  12 |
| Interleukin-3, Interleukin-5 and GM-CSF signaling | -1.948460 |  33 |  15 |
| Activation of IRF3, IRF7 mediated by TBK1, IKKε (IKBKE) | -1.975521 |  10 |   5 |
| Signaling by CSF1 (M-CSF) in myeloid cells | -2.003615 |  26 |   9 |
| DNA strand elongation | -2.010772 |  29 |  24 |
| Activation of the pre-replicative complex | -2.033192 |  27 |  22 |
| Glycosphingolipid catabolism | -2.090416 |  26 |  16 |
| Antigen processing-Cross presentation | -2.098850 |  73 |  33 |
| Potassium Channels | -2.126779 |  28 |   4 |
| Regulation of Complement cascade | -2.143450 |  16 |  10 |
| Complement cascade | -2.147672 |  19 |  11 |
| Cholesterol biosynthesis | -2.205904 |  20 |  12 |
| ROS and RNS production in phagocytes | -2.257629 |  22 |  11 |
| Interferon alpha/beta signaling | -2.287254 |  33 |  22 |



### IV-LD

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| Interferon alpha/beta signaling | -2.558410 |  33 |  19 |
| Antigen activates B Cell Receptor (BCR) leading to generation of second messengers | -2.256919 |  26 |  14 |
| Signaling by CSF1 (M-CSF) in myeloid cells | -2.247666 |  25 |   9 |
| Signaling by CSF3 (G-CSF) | -2.180512 |  25 |   8 |
| Evasion by RSV of host interferon responses | -2.177209 |  17 |   9 |
| Cytosolic tRNA aminoacylation | -2.161461 |  22 |  14 |
| Regulation of RUNX3 expression and activity | -2.148208 |  47 |  29 |
| Hh mutants are degraded by ERAD | -2.132418 |  48 |  31 |
| Hh mutants abrogate ligand secretion | -2.132418 |  48 |  31 |
| Glycosphingolipid catabolism | -2.129369 |  20 |  14 |
| ROS and RNS production in phagocytes | -2.114074 |  21 |  15 |
| GPVI-mediated activation cascade | -2.104005 |  26 |  13 |
| FCERI mediated Ca+2 mobilization | -2.089354 |  23 |   9 |
| Inactivation of CSF3 (G-CSF) signaling | -2.080420 |  20 |   8 |
| Signaling by the B Cell Receptor (BCR) | -2.078370 |  92 |  47 |
| Antigen processing-Cross presentation | -2.055188 |  73 |  40 |
| Hedgehog ligand biogenesis | -2.049071 |  51 |  32 |
| Vif-mediated degradation of APOBEC3G | -2.047086 |  47 |  28 |
| TCR signaling | -2.021049 |  81 |  42 |
| RHO GTPases Activate NADPH Oxidases | -2.020000 |  17 |   7 |
| Defective CFTR causes cystic fibrosis | -2.005344 |  52 |  31 |
| Interleukin-3, Interleukin-5 and GM-CSF signaling | -1.991710 |  33 |  12 |
| SARS-CoV-2 activates/modulates innate and adaptive immune responses | -1.988786 |  71 |  23 |
| Vpu mediated degradation of CD4 | -1.988303 |  45 |  27 |
| TRAF6 mediated IRF7 activation | -1.984811 |  12 |   6 |
| Activation of NF-kappaB in B cells | -1.980958 |  56 |  35 |
| GSK3B and BTRC:CUL1-mediated-degradation of NFE2L2 | -1.966287 |  46 |  25 |
| Diseases associated with glycosylation precursor biosynthesis | -1.953961 |  13 |  10 |
| Regulation of signaling by CBL | -1.938136 |  20 |   9 |
| Degradation of GLI2 by the proteasome | -1.932326 |  50 |  27 |
| FCERI mediated MAPK activation | -1.931653 |  25 |  10 |
| Degradation of AXIN | -1.920396 |  48 |  28 |
| DAP12 interactions | -1.906521 |  18 |   8 |
| DAP12 signaling | -1.906521 |  18 |   8 |
| TRAF3-dependent IRF activation pathway | -1.905810 |  10 |   6 |
| FBXL7 down-regulates AURKA during mitotic entry and in early mitosis | -1.905041 |  48 |  27 |
| Autodegradation of the E3 ubiquitin ligase COP1 | -1.896513 |  46 |  27 |
| Stabilization of p53 | -1.885473 |  48 |  27 |
| Regulation of activated PAK-2p34 by proteasome mediated degradation | -1.882216 |  44 |  26 |
| Generation of second messenger molecules | -1.880302 |  14 |   7 |
| NRIF signals cell death from the nucleus | -1.871491 |  10 |   5 |
| FCERI mediated NF-kB activation | -1.865545 |  63 |  34 |
| DDX58/IFIH1-mediated induction of interferon-alpha/beta | -1.859041 |  45 |  18 |
| Cross-presentation of soluble exogenous antigens (endosomes) | -1.850929 |  42 |  25 |
| Parasite infection | -1.845931 |  49 |  15 |
| Leishmania phagocytosis | -1.845931 |  49 |  15 |
| FCGR3A-mediated phagocytosis | -1.845931 |  49 |  15 |
| SCF-beta-TrCP mediated degradation of Emi1 | -1.845172 |  47 |  27 |
| SARS-CoV-1 modulates host translation machinery | -1.837814 |  30 |  16 |
| AUF1 (hnRNP D0) binds and destabilizes mRNA | -1.835794 |  47 |  28 |
| Formation of the ternary complex, and subsequently, the 43S complex | -1.832124 |  45 |  27 |
| Regulation of Apoptosis | -1.827798 |  46 |  26 |
| Ubiquitin Mediated Degradation of Phosphorylated Cdc25A | -1.827420 |  45 |  26 |
| p53-Independent DNA Damage Response | -1.827420 |  45 |  26 |
| p53-Independent G1/S DNA damage checkpoint | -1.827420 |  45 |  26 |
| Fc epsilon receptor (FCERI) signaling | -1.826947 | 103 |  46 |
| Interferon gamma signaling | -1.826233 |  37 |  16 |
| Negative regulators of DDX58/IFIH1 signaling | -1.824630 |  24 |  10 |
| Activation of IRF3, IRF7 mediated by TBK1, IKKε (IKBKE) | -1.813470 |  10 |   6 |
| Degradation of GLI1 by the proteasome | -1.801153 |  50 |  27 |
| Regulation of KIT signaling | -1.794670 |  13 |   6 |
| Orc1 removal from chromatin | -1.793490 |  62 |  29 |
| ER-Phagosome pathway | -1.793228 |  61 |  32 |
| Cholesterol biosynthesis | -1.789364 |  18 |  10 |
| SARS-CoV-2 modulates host translation machinery | -1.775869 |  39 |  21 |
| Interferon Signaling | -1.771713 | 140 |  49 |
| Negative regulation of NOTCH4 signaling | -1.766764 |  47 |  25 |
| Cytokine Signaling in Immune system | -1.759874 | 387 | 120 |
| Degradation of DVL | -1.759245 |  49 |  27 |
| Interleukin-1 signaling | -1.739892 |  81 |  41 |
| TNFR2 non-canonical NF-kB pathway | -1.730040 |  56 |  31 |
| NIK-->noncanonical NF-kB signaling | -1.728636 |  51 |  28 |
| Fcgamma receptor (FCGR) dependent phagocytosis | -1.726647 |  71 |  19 |
| GLI3 is processed to GLI3R by the proteasome | -1.722373 |  51 |  27 |
| Ovarian tumor domain proteases | -1.708164 |  25 |  10 |
| Selenocysteine synthesis | -1.702643 |  74 |  47 |
| Response of EIF2AK4 (GCN2) to amino acid deficiency | -1.702423 |  82 |  52 |
| Downstream TCR signaling | -1.702062 |  66 |  35 |
| Regulation of NF-kappa B signaling | -1.701554 |  13 |   9 |
| tRNA Aminoacylation | -1.698055 |  37 |  12 |
| SARS-CoV-2-host interactions | -1.694829 | 133 |  49 |
| Cytosolic sulfonation of small molecules | -1.693484 |  12 |   7 |
| Toll-like Receptor Cascades | -1.685631 |  96 |  31 |
| Translation of Structural Proteins | -1.683454 |  21 |  12 |
| Viral mRNA Translation | -1.681921 |  70 |  45 |
| Interleukin-1 family signaling | -1.677088 |  95 |  45 |
| Regulation of actin dynamics for phagocytic cup formation | -1.676820 |  48 |  17 |
| Signaling by Interleukins | -1.667910 | 242 |  76 |
| Membrane binding and targetting of GAG proteins | -1.667162 |  11 |   7 |
| Synthesis And Processing Of GAG, GAGPOL Polyproteins | -1.667162 |  11 |   7 |
| IKK complex recruitment mediated by RIP1 | -1.663769 |  13 |   5 |
| Host Interactions of HIV factors | -1.646932 | 102 |  43 |
| Formation of a pool of free 40S subunits | -1.645196 |  83 |  51 |
| Downstream signaling events of B Cell Receptor (BCR) | -1.642879 |  67 |  36 |
| Regulation of RAS by GAPs | -1.633428 |  56 |  30 |
| TICAM1, RIP1-mediated IKK complex recruitment | -1.625165 |  12 |   5 |
| Chondroitin sulfate biosynthesis | -1.619294 |  10 |   5 |
| SRP-dependent cotranslational protein targeting to membrane | -1.618303 |  89 |  56 |
| Regulation of ornithine decarboxylase (ODC) | -1.617211 |  46 |  26 |
| ABC transporter disorders | -1.605174 |  59 |  33 |
| Antigen Presentation: Folding, assembly and peptide loading of class I MHC | -1.603848 |  19 |   9 |
| Eukaryotic Translation Termination | -1.602546 |  73 |  46 |
| Ubiquitin-dependent degradation of Cyclin D | -1.602206 |  46 |  26 |
| Eukaryotic Translation Elongation | -1.594347 |  73 |  47 |
| Ribosomal scanning and start codon recognition | -1.589249 |  51 |  29 |
| Nucleotide salvage | -1.588633 |  20 |   7 |
| C-type lectin receptors (CLRs) | -1.586923 |  88 |  39 |
| Late SARS-CoV-2 Infection Events | -1.583797 |  48 |  24 |
| NOD1/2 Signaling Pathway | -1.579994 |  25 |   9 |
| Translation initiation complex formation | -1.579983 |  50 |  29 |
| RAC2 GTPase cycle | -1.577927 |  75 |  17 |
| Neutrophil degranulation | -1.574600 | 281 |  81 |
| Eukaryotic Translation Initiation | -1.567530 |  98 |  60 |
| Cap-dependent Translation Initiation | -1.567530 |  98 |  60 |
| Nuclear events stimulated by ALK signaling in cancer | -1.563483 |  26 |   6 |
| Peptide chain elongation | -1.557095 |  71 |  45 |
| IRE1alpha activates chaperones | -1.554711 |  39 |  22 |
| Dectin-1 mediated noncanonical NF-kB signaling | -1.552256 |  52 |  28 |
| L13a-mediated translational silencing of Ceruloplasmin expression | -1.549401 |  91 |  56 |
| RSV-host interactions | -1.547365 |  55 |  17 |
| GTP hydrolysis and joining of the 60S ribosomal subunit | -1.544596 |  92 |  56 |
| Influenza Viral RNA Transcription and Replication | -1.538636 | 109 |  67 |
| Plasma lipoprotein assembly, remodeling, and clearance | -1.531178 |  39 |  11 |
| Regulation of expression of SLITs and ROBOs | -1.520219 | 134 |  83 |
| Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC) | -1.518961 |  75 |  47 |
| Selenoamino acid metabolism | -1.513412 |  95 |  44 |
| Activation of the mRNA upon binding of the cap-binding complex and eIFs, and subsequent binding to 43S | -1.512898 |  51 |  29 |
| Metabolism of polyamines | -1.510906 |  53 |  28 |
| Signaling by SCF-KIT | -1.510103 |  34 |  11 |
| Influenza Infection | -1.501986 | 123 |  72 |
| N-glycan trimming in the ER and Calnexin/Calreticulin cycle | -1.490813 |  30 |  12 |
| XBP1(S) activates chaperone genes | -1.489123 |  38 |  21 |
| RAC1 GTPase cycle | -1.488345 | 143 |  33 |
| RHO GTPases Activate WASPs and WAVEs | -1.488342 |  30 |   8 |
| Signaling by ERBB4 | -1.477576 |  31 |   9 |
| Toll Like Receptor 4 (TLR4) Cascade | -1.476812 |  82 |  24 |
| CLEC7A (Dectin-1) signaling | -1.470010 |  76 |  36 |
| p75 NTR receptor-mediated signalling | -1.452464 |  59 |  16 |
| SARS-CoV-2 Infection | -1.448607 | 197 |  69 |
| Somitogenesis | -1.447189 |  44 |  24 |
| G2/M Checkpoints | -1.446871 | 109 |  43 |
| rRNA modification in the nucleus and cytosol | -1.437797 |  55 |  30 |
| Oxygen-dependent proline hydroxylation of Hypoxia-inducible Factor Alpha | -1.437746 |  59 |  29 |
| Respiratory Syncytial Virus Infection Pathway | -1.418738 |  73 |  22 |
| rRNA processing | -1.417312 | 166 |  96 |
| rRNA processing in the nucleus and cytosol | -1.399510 | 158 |  91 |
| Ion channel transport | -1.390974 |  87 |  21 |
| Major pathway of rRNA processing in the nucleolus and cytosol | -1.389697 | 148 |  87 |
| Degradation of beta-catenin by the destruction complex | -1.389099 |  68 |  32 |
| Adaptive Immune System | -1.352343 | 455 |  97 |
| Disorders of transmembrane transporters | -1.348539 | 109 |  46 |
| Asparagine N-linked glycosylation | -1.258391 | 190 |  79 |
| Response to elevated platelet cytosolic Ca2+ |  1.397186 |  79 |  14 |
| Protein localization |  1.428193 | 122 |  56 |
| Platelet degranulation |  1.436491 |  75 |  14 |
| Transmission across Chemical Synapses |  1.486063 | 103 |  34 |
| Extra-nuclear estrogen signaling |  1.519005 |  44 |  22 |
| Diseases associated with O-glycosylation of proteins |  1.548805 |  24 |  11 |
| Regulation of pyruvate dehydrogenase (PDH) complex |  1.554923 |  12 |   7 |
| MTOR signalling |  1.560298 |  36 |  19 |
| Retinoid metabolism and transport |  1.571927 |  17 |   8 |
| Organelle biogenesis and maintenance |  1.579733 | 180 |  65 |
| Signaling by Activin |  1.584837 |  11 |   5 |
| Maturation of TCA enzymes and regulation of TCA cycle |  1.589006 |  17 |  10 |
| Glyoxylate metabolism and glycine degradation |  1.592081 |  13 |   7 |
| Metabolism of nitric oxide: NOS3 activation and regulation |  1.603027 |  10 |  10 |
| Formation of Fibrin Clot (Clotting Cascade) |  1.620177 |  12 |   4 |
| FOXO-mediated transcription of cell cycle genes |  1.629367 |  12 |   7 |
| Striated Muscle Contraction |  1.631289 |  22 |   6 |
| Ion homeostasis |  1.634688 |  28 |  15 |
| Regulation of TP53 Activity through Association with Co-factors |  1.636267 |  10 |   7 |
| Signaling by Retinoic Acid |  1.649813 |  23 |  13 |
| Sensory Perception |  1.652458 |  68 |  26 |
| Cardiogenesis |  1.667700 |  20 |   9 |
| Regulation of Complement cascade |  1.676912 |  12 |   5 |
| Defective B3GALTL causes PpS |  1.685302 |  15 |   9 |
| Protein-protein interactions at synapses |  1.685851 |  39 |  15 |
| Signaling by BMP |  1.701077 |  19 |   7 |
| Nuclear Receptor transcription pathway |  1.705083 |  30 |  14 |
| Heme signaling |  1.711499 |  37 |  18 |
| Neuronal System |  1.716657 | 147 |  52 |
| Mitochondrial protein degradation |  1.735200 |  78 |  38 |
| Lysine catabolism |  1.741260 |  12 |   6 |
| Aberrant regulation of mitotic cell cycle due to RB1 defects |  1.742286 |  28 |  13 |
| Diseases of mitotic cell cycle |  1.744056 |  30 |  14 |
| O-glycosylation of TSR domain-containing proteins |  1.744456 |  16 |  10 |
| Energy dependent regulation of mTOR by LKB1-AMPK |  1.744984 |  24 |  15 |
| Citric acid cycle (TCA cycle) |  1.746457 |  30 |  18 |
| Complex I biogenesis |  1.785015 |  44 |  32 |
| FOXO-mediated transcription of oxidative stress, metabolic and neuronal genes |  1.791713 |  15 |   9 |
| Cardiac conduction |  1.797730 |  54 |  22 |
| Cristae formation |  1.835000 |  24 |  17 |
| Transcriptional activation of mitochondrial biogenesis |  1.852857 |  42 |  28 |
| Formation of ATP by chemiosmotic coupling |  1.871819 |  12 |  11 |
| Muscle contraction |  1.889581 |  99 |  42 |
| Respiratory electron transport |  1.950513 |  83 |  58 |
| FOXO-mediated transcription |  1.970064 |  44 |  22 |
| Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins. |  2.018368 |  98 |  70 |
| Aerobic respiration and respiratory electron transport |  2.033819 | 160 | 106 |
| Mitochondrial biogenesis |  2.128263 |  73 |  46 |



























### Interpretation of GSEA Results



From ChatGPT when prompted to explain a negative *Normalised Enrichment Score* (NES):

*In summary, a negative NES in GSEA indicates that the gene set is  significantly downregulated in the condition of interest compared to the control, reflecting a negative correlation with the phenotype being  studied. This information can provide valuable insights into the  biological processes that are less active or suppressed in the  experimental condition.*


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
