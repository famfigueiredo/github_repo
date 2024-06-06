There are no explicitly immune related pathways highly ranked in GSEA at 1WPC. However, they show up prominently at 4WPC.

Can this be related to the fact that the heart is not an 'immune organ'? What is the case in the liver, head-kidney, and spleen?

However, in EOMES, more general GO terms like *negative regulation of response to external stimulus*, *positive regulation of response to external stimulus*, and *regulation of response to external stimulus* are present, with relatively large gene sets, generally upregulated. Is this evidence of an early antiviral response?

Responses to external stimuli show up in the DNA vaccine, EOMES, and GATA3, but not in IV-HD? Evidence of viral infection recognition?



# Heart

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



#### Summary

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


#### Overall Interpretation

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



## 10WPC

### DNA vaccine

ReactomePA results

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| Plasma lipoprotein clearance | -2.134091 |  14 |  3 |
| Plasma lipoprotein assembly, remodeling, and clearance | -2.014837 |  18 |  3 |
| E3 ubiquitin ligases ubiquitinate target proteins | -2.003184 |  11 |  5 |
| Protein ubiquitination | -1.963409 |  14 |  6 |
| Complex I biogenesis | -1.898460 |  35 | 15 |
| Prefoldin mediated transfer of substrate  to CCT/TriC | -1.826903 |  15 |  7 |
| Cooperation of Prefoldin and TriC/CCT  in actin and tubulin folding | -1.826903 |  15 |  7 |
| Aerobic respiration and respiratory electron transport | -1.781939 | 112 | 41 |
| Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins. | -1.750500 |  76 | 37 |
| Macroautophagy | -1.684341 |  44 |  8 |
| Respiratory electron transport | -1.652091 |  64 | 30 |
| Signaling by Receptor Tyrosine Kinases |  1.579934 | 118 | 30 |
| Metabolism of water-soluble vitamins and cofactors |  1.735677 |  17 | 10 |
| Non-integrin membrane-ECM interactions |  1.769223 |  20 |  6 |
| RHOQ GTPase cycle |  1.769246 |  15 |  7 |
| SLC transporter disorders |  1.786018 |  10 |  5 |
| MET promotes cell motility |  1.789446 |  11 |  5 |
| ECM proteoglycans |  1.796718 |  25 |  9 |
| RHOG GTPase cycle |  1.849696 |  22 | 11 |
| Signaling by PDGF |  1.851149 |  12 |  7 |
| Assembly of collagen fibrils and other multimeric structures |  1.862744 |  16 | 10 |
| MET activates PTK2 signaling |  1.864245 |  10 |  5 |
| RHOJ GTPase cycle |  1.981134 |  13 |  7 |
| Cell surface interactions at the vascular wall |  1.988497 |  31 | 14 |
| RAC3 GTPase cycle |  2.011709 |  18 | 11 |
| RAC2 GTPase cycle |  2.035045 |  21 | 12 |
| Syndecan interactions |  2.045470 |  13 |  7 |

#### Overall Interpretation

- **Downregulation**: The downregulated pathways are primarily involved in lipid metabolism, protein ubiquitination, mitochondrial function, and autophagy. This suggests a suppression of energy production, lipid clearance, protein turnover, and cellular recycling processes.
- **Upregulation**: The upregulated pathways involve cell signaling, particularly through receptor tyrosine kinases, GTPase cycles, extracellular matrix interactions, and vitamin metabolism. This indicates enhanced cellular communication, structural dynamics, and metabolic activities.

The combination of these changes suggests a cellular state where there's a shift away from energy production and basic metabolic processes towards more complex signaling and structural functions. This could be indicative of a response to external stimuli, tissue repair, or a shift towards a more specialized cellular function.



### EOMES

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| Post-translational protein phosphorylation |  1.858368 |  72 |  23 |
| Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs) |  1.851716 |  82 |  27 |
| Formation of Fibrin Clot (Clotting Cascade) |  1.833651 |  22 |  11 |
| Common Pathway of Fibrin Clot Formation |  1.810430 |  11 |   8 |
| Regulation of Complement cascade |  1.739036 |  23 |  11 |
| Platelet degranulation |  1.703373 |  96 |  23 |
| Intrinsic Pathway of Fibrin Clot Formation |  1.684758 |  13 |   5 |
| Response to elevated platelet cytosolic Ca2+ |  1.672309 | 100 |  23 |
| Plasma lipoprotein remodeling |  1.652534 |  22 |  10 |
| Gluconeogenesis |  1.650880 |  20 |   9 |
| Nuclear Receptor transcription pathway |  1.639433 |  41 |  18 |
| Surfactant metabolism |  1.633069 |  14 |   8 |
| Complement cascade |  1.620782 |  29 |  11 |
| Platelet activation, signaling and aggregation |  1.435253 | 200 |  62 |
| G alpha (i) signalling events |  1.394216 | 151 |  51 |
| GPCR ligand binding |  1.376630 | 197 |  66 |
| Hemostasis |  1.341971 | 429 | 168 |
| Signaling by GPCR |  1.293262 | 385 | 118 |
| Chromosome Maintenance | -1.362768 |  83 |  28 |
| Dual incision in TC-NER | -1.480243 |  54 |  19 |
| Gap-filling DNA repair synthesis and ligation in TC-NER | -1.494346 |  54 |  20 |
| Homology Directed Repair | -1.523449 |  92 |  24 |
| DNA strand elongation | -1.716627 |  30 |  12 |
| Activation of ATR in response to replication stress | -1.716692 |  35 |  14 |
| Cytoprotection by HMOX1 | -1.881835 |  47 |   2 |
| Triglyceride metabolism | -2.211085 |  25 |   4 |
| Triglyceride catabolism | -2.249671 |  15 |   4 |
| Metabolism of porphyrins | -2.452706 |  21 |   2 |

#### Overall Interpretation

- **Upregulated Pathways**: The upregulated pathways in Atlantic salmon heart tissue are predominantly related to blood clotting (hemostasis), platelet activation, immune response (complement cascade), and various signaling pathways (GPCR signaling, insulin-like growth factor transport). This suggests an enhanced readiness for hemostasis, immune defense, and cellular communication, which are critical for maintaining cardiovascular function and responding to physiological stress or injury.
- **Downregulated Pathways**: The downregulated pathways are mainly involved in DNA repair and maintenance (chromosome maintenance, nucleotide excision repair, homologous recombination), cytoprotection (HMOX1), and lipid metabolism (triglyceride metabolism and catabolism). This indicates a reduced capacity for DNA repair, potentially leading to genomic instability, and decreased lipid metabolism, which may affect energy balance and cellular health.

This combination of upregulated and downregulated pathways provides a snapshot of the metabolic and regulatory state of the heart tissue, highlighting areas of increased activity in response to physiological demands and potential vulnerabilities in DNA repair and lipid metabolism.



### GATA3

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| ECM proteoglycans |  2.104700 |  25 |  11 |
| Syndecan interactions |  2.099630 |  13 |   7 |
| Signaling by NTRK1 (TRKA) |  2.083449 |  36 |  12 |
| Signaling by NTRKs |  2.026909 |  38 |  12 |
| Nuclear Events (kinase and transcription factor activation) |  2.017561 |  21 |   8 |
| Non-integrin membrane-ECM interactions |  1.996021 |  20 |   7 |
| RHOJ GTPase cycle |  1.900899 |  13 |   6 |
| RHOQ GTPase cycle |  1.847776 |  15 |   6 |
| Signaling by PDGF |  1.813468 |  12 |   8 |
| Elastic fibre formation |  1.782823 |  10 |   8 |
| NGF-stimulated transcription |  1.776159 |  11 |   6 |
| RAC3 GTPase cycle |  1.744556 |  18 |   8 |
| MAPK targets/ Nuclear events mediated by MAP kinases |  1.735685 |  12 |   3 |
| Striated Muscle Contraction |  1.731935 |  13 |   6 |
| MET activates PTK2 signaling |  1.729049 |  10 |   4 |
| Assembly of collagen fibrils and other multimeric structures |  1.724211 |  16 |  10 |
| NCAM signaling for neurite out-growth |  1.723873 |  12 |   9 |
| MET promotes cell motility |  1.618208 |  11 |   4 |
| Signaling by GPCR |  1.613790 |  57 |  27 |
| Cell surface interactions at the vascular wall |  1.613695 |  31 |  11 |
| Signaling by Receptor Tyrosine Kinases |  1.608698 | 118 |  32 |
| Integrin cell surface interactions |  1.588736 |  24 |  12 |
| GPCR downstream signalling |  1.511815 |  53 |  25 |
| Signal Transduction |  1.212872 | 467 | 115 |
| Metabolism of RNA | -1.309920 | 246 |  68 |
| Immune System | -1.401149 | 383 | 125 |
| Infectious disease | -1.410132 | 308 | 113 |
| Aerobic respiration and respiratory electron transport | -1.442175 | 112 |  48 |
| Signaling by Hedgehog | -1.451002 |  53 |  23 |
| Viral Infection Pathways | -1.453647 | 261 |  96 |
| TCR signaling | -1.470209 |  46 |  20 |
| RHO GTPase Effectors | -1.470920 |  71 |  27 |
| C-type lectin receptors (CLRs) | -1.483550 |  44 |  19 |
| Assembly of the pre-replicative complex | -1.486404 |  39 |  20 |
| Macroautophagy | -1.489676 |  44 |  17 |
| S Phase | -1.503583 |  50 |  25 |
| FCERI mediated NF-kB activation | -1.510886 |  39 |  19 |
| Hedgehog ligand biogenesis | -1.524564 |  40 |  20 |
| Innate Immune System | -1.524781 | 232 |  79 |
| Cell Cycle | -1.526362 | 140 |  53 |
| MITF-M-dependent gene expression | -1.529256 |  22 |   8 |
| MITF-M-regulated melanocyte development | -1.541983 |  29 |  11 |
| Cytosolic tRNA aminoacylation | -1.553471 |  11 |   5 |
| Protein folding | -1.556639 |  22 |  11 |
| ABC transporter disorders | -1.556750 |  41 |  19 |
| Translation | -1.559095 | 163 |  73 |
| Switching of origins to a post-replicative state | -1.562362 |  38 |  20 |
| Resolution of Sister Chromatid Cohesion | -1.567546 |  27 |  13 |
| Cell Cycle, Mitotic | -1.569749 | 117 |  46 |
| mRNA Splicing - Minor Pathway | -1.570548 |  21 |  11 |
| Aggrephagy | -1.579937 |  12 |   7 |
| Chaperonin-mediated protein folding | -1.585963 |  21 |  11 |
| DNA Replication | -1.601078 |  45 |  25 |
| Neutrophil degranulation | -1.603845 | 123 |  51 |
| Protein ubiquitination | -1.610171 |  14 |   8 |
| The role of GTSE1 in G2/M progression after G2 checkpoint | -1.613267 |  42 |  20 |
| Cilium Assembly | -1.617508 |  26 |  12 |
| Synthesis of DNA | -1.636671 |  42 |  24 |
| Cell Cycle Checkpoints | -1.641058 |  68 |  31 |
| ABC-family proteins mediated transport | -1.650507 |  46 |  24 |
| Mitotic G2-G2/M phases | -1.655254 |  62 |  28 |
| G2/M Transition | -1.655254 |  62 |  28 |
| M Phase | -1.657101 |  94 |  39 |
| RNA Polymerase II Pre-transcription Events | -1.659920 |  18 |   8 |
| Interferon alpha/beta signaling | -1.663718 |  19 |   8 |
| tRNA Aminoacylation | -1.681849 |  13 |   6 |
| Mitotic Metaphase and Anaphase | -1.689240 |  73 |  32 |
| Mitotic Anaphase | -1.689240 |  73 |  32 |
| Formation of tubulin folding intermediates by CCT/TriC | -1.711465 |  10 |   7 |
| Separation of Sister Chromatids | -1.730907 |  61 |  30 |
| Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins. | -1.757196 |  76 |  31 |
| Host Interactions of HIV factors | -1.793768 |  50 |  27 |
| Amyloid fiber formation | -1.794465 |  15 |   5 |
| Respiratory electron transport | -1.805950 |  64 |  29 |
| HIV Infection | -1.887990 |  75 |  37 |
| Prefoldin mediated transfer of substrate  to CCT/TriC | -1.911113 |  15 |  10 |
| Cooperation of Prefoldin and TriC/CCT  in actin and tubulin folding | -1.911113 |  15 |  10 |
| Mitochondrial translation initiation | -1.946643 |  39 |  23 |
| Mitochondrial translation termination | -1.946643 |  39 |  23 |
| Plasma lipoprotein clearance | -1.955145 |  14 |   5 |
| Mitochondrial translation | -1.978389 |  40 |  24 |
| Mitochondrial translation elongation | -1.978389 |  40 |  24 |
| Complex I biogenesis | -2.008212 |  35 |  19 |



#### Conclusion

The pathways in Atlantic salmon heart tissue exhibit significant upregulation in signaling pathways related to cell adhesion, growth, and cytoskeletal dynamics, which are essential for maintaining cardiac function and responding to stress or damage. Simultaneously, there is downregulation in pathways associated with the immune response, DNA repair, and mitochondrial function, which may impact the heart’s ability to cope with oxidative stress and maintain cellular homeostasis. This combination of pathway activities suggests a complex regulatory environment in the heart tissue, balancing between structural integrity, signaling, and energy metabolism.



### IV-HD

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| FOXO-mediated transcription |  2.027766 |  40 |  21 |
| Nuclear Receptor transcription pathway |  1.902506 |  23 |  14 |
| Response of EIF2AK1 (HRI) to heme deficiency |  1.816329 |  12 |   7 |
| PI3K/AKT Signaling in Cancer |  1.816199 |  44 |  19 |
| Gluconeogenesis |  1.786022 |  13 |   5 |
| IRS-related events triggered by IGF1R |  1.764965 |  19 |  10 |
| SUMOylation of intracellular receptors |  1.715626 |  11 |   7 |
| Myogenesis |  1.709299 |  17 |  10 |
| Signaling by ALK |  1.705864 |  20 |  14 |
| IGF1R signaling cascade |  1.704783 |  20 |  10 |
| TNFR1-induced proapoptotic signaling |  1.702236 |  13 |   6 |
| DAG and IP3 signaling |  1.678151 |  20 |  15 |
| Signaling by NTRK1 (TRKA) |  1.676073 |  67 |  26 |
| Signaling by NTRKs |  1.663244 |  75 |  27 |
| Transcriptional activation of mitochondrial biogenesis |  1.653007 |  36 |  13 |
| Striated Muscle Contraction |  1.630470 |  20 |   6 |
| NGF-stimulated transcription |  1.622448 |  20 |  10 |
| Insulin receptor signalling cascade |  1.609549 |  20 |  10 |
| Sensory processing of sound by outer hair cells of the cochlea |  1.602291 |  17 |   6 |
| Constitutive Signaling by AKT1 E17K in Cancer |  1.601176 |  18 |   9 |
| Transport of bile salts and organic acids, metal ions and amine compounds |  1.595090 |  20 |   9 |
| IRS-mediated signalling |  1.594991 |  16 |   8 |
| Sema4D in semaphorin signaling |  1.587064 |  19 |  11 |
| Transcriptional regulation of brown and beige adipocyte differentiation |  1.586154 |  18 |   7 |
| Transcriptional regulation of brown and beige adipocyte differentiation by EBF2 |  1.586154 |  18 |   7 |
| NCAM1 interactions |  1.585243 |  11 |   4 |
| Circadian Clock |  1.584840 |  49 |  22 |
| PI3K Cascade |  1.584450 |  14 |   7 |
| Constitutive Signaling by Aberrant PI3K in Cancer |  1.583423 |  25 |  10 |
| HSF1-dependent transactivation |  1.577341 |  18 |  13 |
| PI5P, PP2A and IER3 Regulate PI3K/AKT Signaling |  1.576598 |  39 |  14 |
| Negative regulation of the PI3K/AKT network |  1.558053 |  43 |  15 |
| Mitochondrial biogenesis |  1.532060 |  67 |  20 |
| Chromatin modifying enzymes |  1.464375 | 145 |  53 |
| Chromatin organization |  1.464375 | 145 |  53 |
| Adaptive Immune System | -1.279300 | 380 |  90 |
| Class I MHC mediated antigen processing & presentation | -1.310807 | 230 |  56 |
| Late Phase of HIV Life Cycle | -1.455370 | 102 |  22 |
| Toll-like Receptor Cascades | -1.466204 |  79 |  21 |
| Parasite infection | -1.471162 |  45 |  17 |
| Leishmania phagocytosis | -1.471162 |  45 |  17 |
| FCGR3A-mediated phagocytosis | -1.471162 |  45 |  17 |
| Nucleotide Excision Repair | -1.472963 |  73 |  19 |
| G2/M Checkpoints | -1.474537 |  91 |  28 |
| HIV Infection | -1.475500 | 172 |  50 |
| HCMV Early Events | -1.495223 |  45 |  16 |
| Orc1 removal from chromatin | -1.498078 |  59 |  22 |
| APC/C-mediated degradation of cell cycle proteins | -1.503137 |  64 |  17 |
| Regulation of mitotic cell cycle | -1.503137 |  64 |  17 |
| APC/C:Cdc20 mediated degradation of mitotic proteins | -1.513624 |  56 |  16 |
| Disorders of transmembrane transporters | -1.519665 |  96 |  36 |
| Diseases of metabolism | -1.526189 |  99 |  28 |
| Toll Like Receptor 4 (TLR4) Cascade | -1.537714 |  67 |  12 |
| Regulation of TP53 Activity through Phosphorylation | -1.543226 |  56 |  14 |
| tRNA processing in the nucleus | -1.546369 |  36 |  18 |
| Rev-mediated nuclear export of HIV RNA | -1.547519 |  24 |  14 |
| Global Genome Nucleotide Excision Repair (GG-NER) | -1.552979 |  52 |  15 |
| Diseases of DNA repair | -1.556339 |  24 |   8 |
| RHO GTPases Activate Formins | -1.558956 |  69 |  24 |
| RHO GTPase Effectors | -1.559479 | 153 |  33 |
| Interactions of Rev with host cellular proteins | -1.561168 |  26 |  14 |
| Neutrophil degranulation | -1.565604 | 253 |  69 |
| Metabolism of nucleotides | -1.572199 |  64 |  24 |
| Innate Immune System | -1.574096 | 500 | 126 |
| RHO GTPases Activate WASPs and WAVEs | -1.584989 |  28 |  14 |
| M Phase | -1.589299 | 209 |  45 |
| The role of GTSE1 in G2/M progression after G2 checkpoint | -1.593841 |  55 |  16 |
| Antiviral mechanism by IFN-stimulated genes | -1.596708 |  85 |  27 |
| Nuclear import of Rev protein | -1.597922 |  25 |  12 |
| Inactivation of CSF3 (G-CSF) signaling | -1.598513 |  17 |   5 |
| Postmitotic nuclear pore complex (NPC) reformation | -1.606989 |  21 |   7 |
| Gene Silencing by RNA | -1.616104 |  48 |  20 |
| Regulation of APC/C activators between G1/S and early anaphase | -1.617718 |  60 |  17 |
| Interferon Signaling | -1.621163 | 119 |  28 |
| Antigen processing-Cross presentation | -1.622320 |  69 |  23 |
| Export of Viral Ribonucleoproteins from Nucleus | -1.635787 |  23 |  14 |
| NEP/NS2 Interacts with the Cellular Export Machinery | -1.635787 |  23 |  14 |
| ZBP1(DAI) mediated induction of type I IFNs | -1.637511 |  10 |   7 |
| TICAM1, RIP1-mediated IKK complex recruitment | -1.638459 |  11 |   8 |
| IKK complex recruitment mediated by RIP1 | -1.638459 |  11 |   8 |
| RAS processing | -1.642218 |  11 |   6 |
| Late endosomal microautophagy | -1.659282 |  22 |   7 |
| Translesion synthesis by REV1 | -1.666679 |  13 |   7 |
| Translesion synthesis by POLK | -1.666679 |  13 |   7 |
| Translesion synthesis by POLI | -1.666679 |  13 |   7 |
| DNA Replication Pre-Initiation | -1.671594 |  74 |  25 |
| snRNP Assembly | -1.675309 |  36 |  20 |
| Metabolism of non-coding RNA | -1.675309 |  36 |  20 |
| HIV Life Cycle | -1.675813 | 108 |  26 |
| APC-Cdc20 mediated degradation of Nek2A | -1.682318 |  13 |   6 |
| Evasion by RSV of host interferon responses | -1.685048 |  15 |   6 |
| Activation of APC/C and APC/C:Cdc20 mediated degradation of mitotic proteins | -1.685854 |  57 |  17 |
| Transcription-Coupled Nucleotide Excision Repair (TC-NER) | -1.686556 |  55 |  18 |
| Mitotic Metaphase and Anaphase | -1.691937 | 145 |  41 |
| Mitotic Anaphase | -1.691937 | 145 |  41 |
| DNA Repair | -1.694000 | 154 |  50 |
| Vpr-mediated nuclear import of PICs | -1.699164 |  22 |  12 |
| Removal of the Flap Intermediate from the C-strand | -1.699755 |  10 |   5 |
| Transport of Mature mRNA Derived from an Intronless Transcript | -1.699973 |  29 |  15 |
| Intraflagellar transport | -1.700724 |  15 |   7 |
| Mitotic G1 phase and G1/S transition | -1.703992 | 104 |  37 |
| Cell Cycle Checkpoints | -1.706984 | 156 |  47 |
| SUMOylation of SUMOylation proteins | -1.711252 |  24 |  12 |
| Mitotic Spindle Checkpoint | -1.711715 |  58 |  20 |
| S Phase | -1.712102 | 114 |  40 |
| Cell Cycle | -1.717875 | 357 |  83 |
| Mitotic Prometaphase | -1.728689 |  96 |  31 |
| Transport of Mature mRNAs Derived from Intronless Transcripts | -1.739550 |  30 |  16 |
| NS1 Mediated Effects on Host Pathways | -1.739940 |  27 |  13 |
| G0 and Early G1 | -1.740154 |  19 |   8 |
| Regulation of Glucokinase by Glucokinase Regulatory Protein | -1.744939 |  21 |  12 |
| Defective TPR may confer susceptibility towards thyroid papillary carcinoma (TPC) | -1.744939 |  21 |  12 |
| Separation of Sister Chromatids | -1.747329 | 112 |  34 |
| DNA Double-Strand Break Repair | -1.747986 |  70 |  21 |
| Activation of the pre-replicative complex | -1.752334 |  20 |  12 |
| Transport of Ribonucleoproteins into the Host Nucleus | -1.760923 |  22 |  12 |
| ISG15 antiviral mechanism | -1.763340 |  49 |  22 |
| Assembly Of The HIV Virion | -1.764618 |  12 |   4 |
| Mitotic Prophase | -1.767757 |  60 |  25 |
| Interactions of Vpr with host cellular proteins | -1.768833 |  23 |  13 |
| Transport of the SLBP independent Mature mRNA | -1.770231 |  24 |  13 |
| Nuclear Pore Complex (NPC) Disassembly | -1.774378 |  27 |  14 |
| Detoxification of Reactive Oxygen Species | -1.777591 |  20 |  11 |
| Interferon alpha/beta signaling | -1.787480 |  30 |   8 |
| Mismatch Repair | -1.790752 |  10 |   6 |
| Mismatch repair (MMR) directed by MSH2:MSH6 (MutSalpha) | -1.790752 |  10 |   6 |
| Processive synthesis on the C-strand of the telomere | -1.793722 |  11 |   6 |
| Amplification of signal from the kinetochores | -1.798714 |  51 |  19 |
| Amplification  of signal from unattached  kinetochores via a MAD2  inhibitory signal | -1.798714 |  51 |  19 |
| SARS-CoV-2 activates/modulates innate and adaptive immune responses | -1.806068 |  60 |  19 |
| SUMOylation of ubiquitinylation proteins | -1.809136 |  25 |  13 |
| Membrane binding and targetting of GAG proteins | -1.810276 |  10 |   4 |
| Synthesis And Processing Of GAG, GAGPOL Polyproteins | -1.810276 |  10 |   4 |
| Cell Cycle, Mitotic | -1.817094 | 296 |  75 |
| Viral Messenger RNA Synthesis | -1.839346 |  34 |  16 |
| Transport of the SLBP Dependant Mature mRNA | -1.845530 |  25 |  14 |
| Removal of the Flap Intermediate | -1.846852 |  10 |   8 |
| SUMOylation of DNA replication proteins | -1.849425 |  32 |  13 |
| Condensation of Prophase Chromosomes | -1.855634 |  12 |   6 |
| G1/S Transition | -1.856244 |  96 |  35 |
| Resolution of Sister Chromatid Cohesion | -1.867741 |  65 |  25 |
| HDR through Homologous Recombination (HRR) or Single Strand Annealing (SSA) | -1.873296 |  49 |  19 |
| Homology Directed Repair | -1.878242 |  54 |  20 |
| Activation of ATR in response to replication stress | -1.890795 |  22 |  11 |
| Nuclear Envelope Breakdown | -1.891515 |  40 |  19 |
| Aggrephagy | -1.898108 |  19 |   6 |
| Processing of DNA double-strand break ends | -1.900354 |  34 |  16 |
| Complement cascade | -1.901982 |  11 |   5 |
| Fanconi Anemia Pathway | -1.906671 |  14 |   7 |
| EML4 and NUDC in mitotic spindle formation | -1.908350 |  59 |  23 |
| Processive synthesis on the lagging strand | -1.915336 |  11 |   9 |
| Translesion synthesis by Y family DNA polymerases bypasses lesions on DNA template | -1.926414 |  24 |  11 |
| Transcriptional regulation by small RNAs | -1.926959 |  37 |  18 |
| Regulation of Complement cascade | -1.938878 |  10 |   5 |
| Synthesis of DNA | -1.944058 |  85 |  34 |
| Dual incision in TC-NER | -1.944684 |  44 |  20 |
| DNA Replication | -1.948402 |  91 |  35 |
| Gap-filling DNA repair synthesis and ligation in TC-NER | -1.963493 |  46 |  21 |
| G1/S-Specific Transcription | -1.966459 |  17 |   9 |
| DNA Damage Bypass | -1.967414 |  31 |  17 |
| Translesion Synthesis by POLH | -1.979068 |  14 |   9 |
| Polymerase switching on the C-strand of the telomere | -1.988556 |  12 |  10 |
| HDR through Homologous Recombination (HRR) | -1.994096 |  30 |  12 |
| Dual Incision in GG-NER | -2.018387 |  24 |  11 |
| Base Excision Repair | -2.024403 |  27 |  12 |
| Resolution of Abasic Sites (AP sites) | -2.040838 |  24 |  11 |
| Polymerase switching | -2.050733 |  10 |   9 |
| Leading Strand Synthesis | -2.050733 |  10 |   9 |
| Interconversion of nucleotide di- and triphosphates | -2.063412 |  18 |   7 |
| Cytosolic sensors of pathogen-associated DNA | -2.095129 |  33 |   8 |
| DNA strand elongation | -2.108992 |  22 |  15 |
| Telomere C-strand (Lagging Strand) Synthesis | -2.159225 |  18 |  13 |
| Resolution of AP sites via the multiple-nucleotide patch replacement pathway | -2.160075 |  19 |  14 |
| Recognition of DNA damage by PCNA-containing replication complex | -2.172609 |  21 |  14 |
| Termination of translesion DNA synthesis | -2.187390 |  20 |  11 |
| PCNA-Dependent Long Patch Base Excision Repair | -2.196541 |  15 |  11 |
| Lagging Strand Synthesis | -2.208056 |  15 |  12 |
| Gap-filling DNA repair synthesis and ligation in GG-NER | -2.257043 |  17 |  10 |
| Telomere Maintenance | -2.307736 |  45 |  24 |
| Extension of Telomeres | -2.308929 |  29 |  19 |
| Chromosome Maintenance | -2.314235 |  53 |  27 |

#### Conclusion

The combination of upregulated and downregulated pathways in Atlantic salmon heart tissue indicates a complex regulatory environment focused on maintaining cardiac function. Upregulation of pathways involved in stress response, metabolism, and muscle maintenance suggests an adaptive response to ensure energy supply and muscle integrity. In contrast, downregulation of immune responses and cell cycle activities indicates a shift away from proliferation and inflammation towards stability and homeostasis. Understanding these dynamics can provide deeper insights into how salmon hearts adapt to various environmental stresses and maintain function.



### IV-LD

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| Complex I biogenesis | -2.377325 |  44 |  31 |
| Respiratory electron transport | -2.261479 |  83 |  53 |
| Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins. | -2.195094 |  98 |  60 |
| Aerobic respiration and respiratory electron transport | -1.986631 | 162 |  95 |
| Recognition of DNA damage by PCNA-containing replication complex | -1.964693 |  25 |  13 |
| Vesicle-mediated transport |  1.306716 | 455 | 124 |
| Signaling by Rho GTPases, Miro GTPases and RHOBTB3 |  1.477766 | 497 | 214 |
| Signaling by GPCR |  1.483679 | 228 |  75 |
| Signaling by Rho GTPases |  1.484225 | 485 | 212 |
| Hemostasis |  1.497283 | 338 | 135 |
| Platelet activation, signaling and aggregation |  1.506632 | 164 |  69 |
| GPCR downstream signalling |  1.535800 | 208 |  74 |
| Transcriptional regulation of granulopoiesis |  1.580512 |  25 |   8 |
| Chromatin modifying enzymes |  1.596529 | 170 |  83 |
| Chromatin organization |  1.596529 | 170 |  83 |
| Death Receptor Signaling |  1.633401 | 100 |  40 |
| PPARA activates gene expression |  1.656990 |  92 |  32 |
| RHOQ GTPase cycle |  1.666565 |  50 |  34 |
| Regulation of lipid metabolism by PPARalpha |  1.671811 |  94 |  33 |
| Kinesins |  1.680260 |  28 |  11 |
| Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs) |  1.697851 |  54 |  15 |
| RHO GTPase cycle |  1.703439 | 345 | 175 |
| RAC1 GTPase cycle |  1.723827 | 144 |  80 |
| RHOB GTPase cycle |  1.750238 |  57 |  35 |
| Heme signaling |  1.759474 |  37 |  14 |
| RHOC GTPase cycle |  1.761050 |  61 |  37 |
| NRAGE signals death through JNK |  1.769988 |  38 |  22 |
| Post-translational protein phosphorylation |  1.779441 |  49 |  12 |
| G alpha (12/13) signalling events |  1.841983 |  47 |  26 |
| RHOA GTPase cycle |  1.895649 | 117 |  71 |
| Binding and Uptake of Ligands by Scavenger Receptors |  1.925324 |  23 |   6 |
| Factors involved in megakaryocyte development and platelet production |  1.941952 |  94 |  39 |
| Regulation of Complement cascade |  1.959766 |  12 |   4 |
| CDC42 GTPase cycle |  1.971144 | 119 |  69 |

#### Summary

The combination of these pathways suggests that in Atlantic salmon heart tissue:

- **Positive Enrichment**: There is a significant enrichment in pathways related to signal transduction (GPCR and Rho GTPases), vesicle transport, chromatin organization, lipid metabolism, and platelet production. These processes are essential for cell communication, structural integrity, energy metabolism, and cardiovascular function.
- **Negative Enrichment**: There is a downregulation of pathways involved in mitochondrial respiration and DNA damage recognition. This may indicate a reduced emphasis on oxidative phosphorylation and a potential adaptation in energy metabolism specific to the heart tissue in salmon, possibly reflecting unique metabolic demands or environmental adaptations.



## 10WPI

### DNA vaccine

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| NR1H3 & NR1H2 regulate gene expression linked to cholesterol transport and efflux | -1.875194 |  17 |   8 |
| Peroxisomal lipid metabolism | -1.861756 |  10 |   8 |
| Metabolism of porphyrins | -1.849627 |  15 |   8 |
| Retrograde transport at the Trans-Golgi-Network | -1.680493 |  25 |   7 |
| Biosynthesis of the N-glycan precursor (dolichol lipid-linked oligosaccharide, LLO) and transfer to a nascent protein | -1.660422 |  14 |   5 |
| Synthesis of substrates in N-glycan biosythesis | -1.648169 |  12 |   5 |
| Transcriptional Regulation by MECP2 | -1.644640 |  26 |   9 |
| Iron uptake and transport | -1.481796 |  28 |  10 |
| Intra-Golgi traffic | -1.463751 |  19 |   5 |
| Amino acids regulate mTORC1 | -1.460133 |  26 |   9 |
| NR1H2 and NR1H3-mediated signaling | -1.398787 |  20 |   8 |
| Rab regulation of trafficking | -1.398765 |  53 |  12 |
| Cytoprotection by HMOX1 | -1.380078 |  29 |   9 |
| Fatty acid metabolism | -1.377250 |  68 |  19 |
| Innate Immune System |  1.201911 | 406 | 102 |
| Adaptive Immune System |  1.275866 | 298 | 117 |
| Developmental Biology |  1.288259 | 489 | 178 |
| Translation |  1.326355 | 216 |  73 |
| Organelle biogenesis and maintenance |  1.327841 | 117 |  40 |
| Signaling by WNT |  1.341702 | 135 |  59 |
| Generic Transcription Pathway |  1.341753 | 416 | 117 |
| RNA Polymerase II Transcription |  1.354023 | 489 | 136 |
| Axon guidance |  1.363787 | 272 |  98 |
| Nervous system development |  1.366108 | 279 | 100 |
| Interleukin-1 family signaling |  1.369676 |  69 |  33 |
| PTEN Regulation |  1.375740 |  91 |  46 |
| HIV Infection |  1.379078 | 141 |  58 |
| C-type lectin receptors (CLRs) |  1.380200 |  70 |  32 |
| Oxygen-dependent proline hydroxylation of Hypoxia-inducible Factor Alpha |  1.402031 |  54 |  27 |
| PPARA activates gene expression |  1.405988 |  58 |  16 |
| Regulation of lipid metabolism by PPARalpha |  1.405988 |  58 |  16 |
| Muscle contraction |  1.406642 |  67 |  20 |
| CLEC7A (Dectin-1) signaling |  1.408456 |  64 |  37 |
| Dectin-1 mediated noncanonical NF-kB signaling |  1.409906 |  51 |  25 |
| RAF/MAP kinase cascade |  1.413687 | 122 |  36 |
| Asparagine N-linked glycosylation |  1.414279 | 117 |  38 |
| Cross-presentation of soluble exogenous antigens (endosomes) |  1.416515 |  40 |  22 |
| Transcriptional regulation by RUNX1 |  1.417350 | 103 |  48 |
| Maternal to zygotic transition (MZT) |  1.418396 |  29 |  11 |
| Diseases of signal transduction by growth factor receptors and second messengers |  1.419982 | 252 |  67 |
| Signaling by Receptor Tyrosine Kinases |  1.420946 | 224 |  64 |
| Cell-Cell communication |  1.421635 |  50 |  18 |
| NIK-->noncanonical NF-kB signaling |  1.426048 |  50 |  25 |
| FCERI mediated NF-kB activation |  1.429243 |  53 |  32 |
| Somitogenesis |  1.433963 |  41 |  21 |
| HCMV Early Events |  1.438138 |  34 |  22 |
| Ovarian tumor domain proteases |  1.441190 |  17 |   7 |
| Activation of NF-kappaB in B cells |  1.441425 |  52 |  25 |
| Deubiquitination |  1.441902 | 131 |  66 |
| Signaling by BRAF and RAF1 fusions |  1.441959 |  30 |   8 |
| Infectious disease |  1.444793 | 495 | 189 |
| Signaling by moderate kinase activity BRAF mutants |  1.446979 |  24 |   7 |
| Signaling by RAS mutants |  1.446979 |  24 |   7 |
| Paradoxical activation of RAF signaling by kinase inactive BRAF |  1.446979 |  24 |   7 |
| Signaling downstream of RAS mutants |  1.446979 |  24 |   7 |
| Regulation of RUNX3 expression and activity |  1.452362 |  45 |  25 |
| Gene and protein expression by JAK-STAT signaling after Interleukin-12 stimulation |  1.456245 |  23 |  11 |
| Cell junction organization |  1.458881 |  36 |  14 |
| Striated Muscle Contraction |  1.459084 |  15 |   5 |
| Viral Infection Pathways |  1.461613 | 413 | 167 |
| Transcriptional and post-translational regulation of MITF-M expression and activity |  1.461784 |  17 |   8 |
| Degradation of beta-catenin by the destruction complex |  1.462974 |  63 |  34 |
| Late SARS-CoV-2 Infection Events |  1.464714 |  32 |  17 |
| Signaling by FGFR1 |  1.464931 |  21 |   8 |
| XBP1(S) activates chaperone genes |  1.466037 |  29 |  11 |
| Cytosolic sensors of pathogen-associated DNA |  1.466218 |  27 |   9 |
| tRNA processing |  1.466390 |  30 |  21 |
| Downstream TCR signaling |  1.467353 |  57 |  36 |
| TNFR2 non-canonical NF-kB pathway |  1.468366 |  52 |  27 |
| Metabolism of RNA |  1.469617 | 440 | 200 |
| Transcriptional regulation of brown and beige adipocyte differentiation |  1.473896 |  11 |   7 |
| Transcriptional regulation of brown and beige adipocyte differentiation by EBF2 |  1.473896 |  11 |   7 |
| Synthesis of active ubiquitin: roles of E1 and E2 enzymes |  1.475095 |  19 |   6 |
| Downregulation of ERBB2 signaling |  1.475369 |  12 |   4 |
| Signaling by Interleukins |  1.479267 | 177 |  79 |
| Downregulation of TGF-beta receptor signaling |  1.479837 |  16 |   7 |
| TCR signaling |  1.482881 |  66 |  39 |
| Respiratory Syncytial Virus Infection Pathway |  1.483024 |  51 |  17 |
| RHOV GTPase cycle |  1.486042 |  17 |  14 |
| Regulation of signaling by CBL |  1.487300 |  13 |   8 |
| MAPK1/MAPK3 signaling |  1.488863 | 124 |  37 |
| tRNA processing in the nucleus |  1.489286 |  20 |  14 |
| DNA Damage Recognition in GG-NER |  1.489538 |  24 |   6 |
| Extra-nuclear estrogen signaling |  1.489627 |  33 |  10 |
| Defective CFTR causes cystic fibrosis |  1.491180 |  50 |  27 |
| p53-Dependent G1 DNA Damage Response |  1.493831 |  50 |  29 |
| p53-Dependent G1/S DNA damage checkpoint |  1.493831 |  50 |  29 |
| G1/S DNA Damage Checkpoints |  1.493831 |  50 |  29 |
| Transcriptional regulation by RUNX2 |  1.494496 |  76 |  36 |
| FBXL7 down-regulates AURKA during mitotic entry and in early mitosis |  1.494784 |  46 |  25 |
| Recruitment and ATM-mediated phosphorylation of repair and signaling proteins at DNA double strand breaks |  1.496392 |  23 |  12 |
| DNA Double Strand Break Response |  1.496392 |  23 |  12 |
| RSV-host interactions |  1.498827 |  35 |  14 |
| Degradation of GLI1 by the proteasome |  1.499899 |  48 |  26 |
| Disorders of transmembrane transporters |  1.500319 |  80 |  41 |
| Transcriptional Regulation by TP53 |  1.501665 | 181 |  55 |
| Oncogenic MAPK signaling |  1.501734 |  40 |  10 |
| Antigen processing-Cross presentation |  1.504948 |  61 |  29 |
| Diseases of glycosylation |  1.505632 |  36 |  11 |
| Regulation of RAS by GAPs |  1.507103 |  50 |  26 |
| RHO GTPases Activate ROCKs |  1.510670 |  13 |   6 |
| Negative regulation of MAPK pathway |  1.512410 |  22 |   8 |
| Activation of BH3-only proteins |  1.512440 |  15 |   9 |
| Degradation of DVL |  1.515973 |  46 |  25 |
| Regulation of TP53 Activity |  1.516938 |  82 |  32 |
| MAPK family signaling cascades |  1.519287 | 146 |  41 |
| SCF(Skp2)-mediated degradation of p27/p21 |  1.520072 |  47 |  27 |
| MHC class II antigen presentation |  1.520683 |  56 |  21 |
| EPH-Ephrin signaling |  1.522294 |  40 |  15 |
| Golgi-to-ER retrograde transport |  1.522548 |  65 |  31 |
| IRE1alpha activates chaperones |  1.522732 |  30 |  12 |
| Regulation of RUNX2 expression and activity |  1.523176 |  51 |  27 |
| Autodegradation of Cdh1 by Cdh1:APC/C |  1.523870 |  46 |  25 |
| Cyclin E associated events during G1/S transition |  1.527174 |  56 |  30 |
| Ub-specific processing proteases |  1.529103 | 101 |  57 |
| Nuclear events stimulated by ALK signaling in cancer |  1.529185 |  18 |   4 |
| Signaling by Rho GTPases, Miro GTPases and RHOBTB3 |  1.532654 | 320 |  96 |
| RAF-independent MAPK1/3 activation |  1.533953 |  10 |   4 |
| Regulation of TP53 Activity through Phosphorylation |  1.535663 |  42 |  21 |
| Regulation of TP53 Degradation |  1.536500 |  23 |   7 |
| Regulation of TP53 Expression and Degradation |  1.536500 |  23 |   7 |
| Unfolded Protein Response (UPR) |  1.538254 |  52 |  18 |
| Downstream signaling events of B Cell Receptor (BCR) |  1.539050 |  58 |  33 |
| Signaling by FGFR |  1.539195 |  41 |  19 |
| Vif-mediated degradation of APOBEC3G |  1.539460 |  45 |  25 |
| Vpu mediated degradation of CD4 |  1.540976 |  44 |  25 |
| Hedgehog 'on' state |  1.543479 |  51 |  28 |
| HDR through Single Strand Annealing (SSA) |  1.544441 |  10 |   7 |
| Signaling by Rho GTPases |  1.544952 | 312 |  95 |
| ABC transporter disorders |  1.547030 |  56 |  29 |
| APC/C:Cdc20 mediated degradation of Securin |  1.548144 |  48 |  26 |
| Fc epsilon receptor (FCERI) signaling |  1.549192 |  77 |  45 |
| Regulation of PTEN stability and activity |  1.550553 |  53 |  31 |
| Cell-extracellular matrix interactions |  1.551846 |  10 |   5 |
| Signaling by the B Cell Receptor (BCR) |  1.555098 |  70 |  43 |
| Cyclin A:Cdk2-associated events at S phase entry |  1.555133 |  57 |  31 |
| Regulation of Apoptosis |  1.555829 |  45 |  28 |
| mRNA Splicing |  1.556819 | 156 |  68 |
| Semaphorin interactions |  1.556961 |  35 |  13 |
| Ubiquitin Mediated Degradation of Phosphorylated Cdc25A |  1.557155 |  43 |  25 |
| p53-Independent DNA Damage Response |  1.557155 |  43 |  25 |
| p53-Independent G1/S DNA damage checkpoint |  1.557155 |  43 |  25 |
| Apoptotic cleavage of cellular proteins |  1.557588 |  19 |   7 |
| CDK-mediated phosphorylation and removal of Cdc6 |  1.562933 |  49 |  27 |
| Stabilization of p53 |  1.563927 |  47 |  27 |
| Signaling by MET |  1.564324 |  39 |  11 |
| ABC-family proteins mediated transport |  1.564600 |  63 |  32 |
| Defective pyroptosis |  1.567603 |  13 |   6 |
| GSK3B and BTRC:CUL1-mediated-degradation of NFE2L2 |  1.567627 |  46 |  25 |
| Meiotic recombination |  1.570290 |  12 |   7 |
| mRNA Splicing - Major Pathway |  1.570960 | 155 |  68 |
| Degradation of AXIN |  1.571070 |  46 |  26 |
| Negative regulation of NOTCH4 signaling |  1.573798 |  46 |  25 |
| UCH proteinases |  1.579481 |  60 |  32 |
| Integrin cell surface interactions |  1.580383 |  37 |  13 |
| Sema4D in semaphorin signaling |  1.580631 |  19 |   8 |
| Scavenging by Class A Receptors |  1.580917 |  11 |   4 |
| Hh mutants are degraded by ERAD |  1.586187 |  48 |  27 |
| Hh mutants abrogate ligand secretion |  1.586187 |  48 |  27 |
| ER-Phagosome pathway |  1.586237 |  56 |  29 |
| Regulation of activated PAK-2p34 by proteasome mediated degradation |  1.586941 |  44 |  26 |
| SARS-CoV-1 Infection |  1.587453 |  86 |  27 |
| Nuclear Pore Complex (NPC) Disassembly |  1.589285 |  17 |   9 |
| AUF1 (hnRNP D0) binds and destabilizes mRNA |  1.589295 |  46 |  26 |
| Regulation of mRNA stability by proteins that bind AU-rich elements |  1.589529 |  63 |  33 |
| Autodegradation of the E3 ubiquitin ligase COP1 |  1.592145 |  46 |  27 |
| Degradation of GLI2 by the proteasome |  1.594191 |  49 |  28 |
| GLI3 is processed to GLI3R by the proteasome |  1.594191 |  49 |  28 |
| Host Interactions of HIV factors |  1.596458 |  81 |  44 |
| Potential therapeutics for SARS |  1.597710 |  54 |  29 |
| Formation of Incision Complex in GG-NER |  1.599244 |  22 |   9 |
| Amyloid fiber formation |  1.602731 |  26 |   5 |
| MAPK6/MAPK4 signaling |  1.607457 |  66 |  34 |
| Mitochondrial calcium ion transport |  1.610278 |  19 |   8 |
| Homologous DNA Pairing and Strand Exchange |  1.614107 |  10 |   9 |
| Presynaptic phase of homologous DNA pairing and strand exchange |  1.614107 |  10 |   9 |
| Diseases of DNA Double-Strand Break Repair |  1.614107 |  10 |   9 |
| Defective homologous recombination repair (HRR) due to BRCA2 loss of function |  1.614107 |  10 |   9 |
| SCF-beta-TrCP mediated degradation of Emi1 |  1.615640 |  47 |  27 |
| Sema4D induced cell migration and growth-cone collapse |  1.619142 |  15 |   7 |
| Hedgehog 'off' state |  1.624340 |  61 |  34 |
| Nuclear Envelope (NE) Reassembly |  1.625253 |  35 |  12 |
| RHO GTPases activate PKNs |  1.626478 |  24 |  10 |
| Plasma lipoprotein clearance |  1.626860 |  19 |   7 |
| SARS-CoV Infections |  1.627337 | 203 |  82 |
| APC/C:Cdh1 mediated degradation of Cdc20 and other APC/C:Cdh1 targeted proteins in late mitosis/early G1 |  1.627489 |  50 |  28 |
| SUMOylation |  1.627917 |  71 |  34 |
| Depolymerization of the Nuclear Lamina |  1.628143 |  10 |   4 |
| Hedgehog ligand biogenesis |  1.629750 |  50 |  28 |
| Interleukin-12 signaling |  1.631015 |  27 |  14 |
| Cilium Assembly |  1.633841 |  57 |  26 |
| Nuclear Envelope Breakdown |  1.637358 |  27 |  12 |
| Interconversion of nucleotide di- and triphosphates |  1.639861 |  15 |   5 |
| SUMO E3 ligases SUMOylate target proteins |  1.643844 |  70 |  34 |
| Beta-catenin independent WNT signaling |  1.643980 |  78 |  40 |
| SARS-CoV-2-host interactions |  1.646402 |  89 |  26 |
| Cytokine Signaling in Immune system |  1.646422 | 272 | 105 |
| Asymmetric localization of PCP proteins |  1.647368 |  48 |  29 |
| HSF1 activation |  1.647702 |  10 |   8 |
| Programmed Cell Death |  1.649376 | 117 |  61 |
| Collagen degradation |  1.651167 |  24 |  11 |
| Translation of Structural Proteins |  1.651474 |  29 |  17 |
| Transport of Mature mRNA derived from an Intron-Containing Transcript |  1.657034 |  45 |  28 |
| Base Excision Repair |  1.658025 |  22 |  12 |
| Resolution of Abasic Sites (AP sites) |  1.659551 |  20 |  11 |
| SARS-CoV-2 Infection |  1.661052 | 136 |  44 |
| Collagen chain trimerization |  1.662029 |  12 |   9 |
| Ubiquitin-dependent degradation of Cyclin D |  1.662795 |  46 |  28 |
| SARS-CoV-1 activates/modulates innate immune responses |  1.664464 |  19 |  11 |
| Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs) |  1.666022 |  44 |  11 |
| Signaling by Hedgehog |  1.670115 |  72 |  38 |
| HDR through Homologous Recombination (HRR) or Single Strand Annealing (SSA) |  1.670966 |  29 |  13 |
| Regulation of HSF1-mediated heat shock response |  1.671310 |  36 |  22 |
| Interleukin-12 family signaling |  1.671664 |  29 |  15 |
| Diseases of programmed cell death |  1.672776 |  29 |  17 |
| RHO GTPases activate IQGAPs |  1.675864 |  10 |   7 |
| MET promotes cell motility |  1.678117 |  22 |   8 |
| Post-translational protein phosphorylation |  1.681680 |  41 |  10 |
| RMTs methylate histone arginines |  1.682722 |  20 |   9 |
| Processing of SMDT1 |  1.685182 |  13 |   6 |
| Hemostasis |  1.686948 | 229 |  72 |
| Evasion by RSV of host interferon responses |  1.687224 |  12 |   7 |
| Interferon alpha/beta signaling |  1.687884 |  23 |  10 |
| SARS-CoV-2 activates/modulates innate and adaptive immune responses |  1.688558 |  37 |  19 |
| VEGFR2 mediated vascular permeability |  1.690740 |  15 |   6 |
| SARS-CoV-1 targets host intracellular signalling and regulatory pathways |  1.691042 |  12 |   8 |
| Apoptotic execution phase |  1.693896 |  23 |  10 |
| ISG15 antiviral mechanism |  1.694347 |  38 |  25 |
| HDR through Homologous Recombination (HRR) |  1.694722 |  18 |   9 |
| Bacterial Infection Pathways |  1.695498 |  36 |  10 |
| RHO GTPases activate PAKs |  1.695609 |  13 |   7 |
| Negative regulators of DDX58/IFIH1 signaling |  1.696613 |  18 |  11 |
| Transport of Mature Transcript to Cytoplasm |  1.697061 |  51 |  30 |
| Processing of Capped Intron-Containing Pre-mRNA |  1.697898 | 193 |  91 |
| COPI-dependent Golgi-to-ER retrograde traffic |  1.699115 |  42 |  21 |
| Degradation of the extracellular matrix |  1.699132 |  45 |  23 |
| TP53 Regulates Transcription of Genes Involved in G2 Cell Cycle Arrest |  1.705549 |  11 |   5 |
| Signaling by PDGF |  1.706033 |  29 |  12 |
| Cdc20:Phospho-APC/C mediated degradation of Cyclin A |  1.710573 |  51 |  29 |
| APC:Cdc20 mediated degradation of cell cycle proteins prior to satisfation of the cell cycle checkpoint |  1.710573 |  51 |  29 |
| Cytosolic tRNA aminoacylation |  1.711575 |  17 |  10 |
| Apoptosis |  1.716386 |  96 |  53 |
| Anchoring of the basal body to the plasma membrane |  1.718040 |  33 |  18 |
| PCNA-Dependent Long Patch Base Excision Repair |  1.721329 |  13 |   8 |
| PCP/CE pathway |  1.722499 |  62 |  36 |
| Cellular response to heat stress |  1.727544 |  52 |  29 |
| mRNA 3'-end processing |  1.732851 |  43 |  29 |
| RNA Polymerase II Transcription Termination |  1.733893 |  47 |  32 |
| Transcription-Coupled Nucleotide Excision Repair (TC-NER) |  1.736483 |  45 |  14 |
| Dual incision in TC-NER |  1.739729 |  34 |  14 |
| Diseases associated with glycosaminoglycan metabolism |  1.742922 |  13 |   6 |
| Meiotic synapsis |  1.744094 |  13 |  10 |
| Regulation of APC/C activators between G1/S and early anaphase |  1.747809 |  55 |  31 |
| Aggrephagy |  1.748793 |  18 |  14 |
| APC/C:Cdc20 mediated degradation of mitotic proteins |  1.750323 |  52 |  30 |
| Activation of APC/C and APC/C:Cdc20 mediated degradation of mitotic proteins |  1.750323 |  52 |  30 |
| tRNA Aminoacylation |  1.751536 |  24 |  14 |
| Recruitment of mitotic centrosome proteins and complexes |  1.752162 |  27 |  16 |
| Centrosome maturation |  1.752162 |  27 |  16 |
| Binding and Uptake of Ligands by Scavenger Receptors |  1.756294 |  16 |   8 |
| Assembly of the pre-replicative complex |  1.763883 |  58 |  32 |
| Loss of Nlp from mitotic centrosomes |  1.765790 |  26 |  16 |
| Loss of proteins required for interphase microtubule organization from the centrosome |  1.765790 |  26 |  16 |
| AURKA Activation by TPX2 |  1.765790 |  26 |  16 |
| Switching of origins to a post-replicative state |  1.766310 |  60 |  33 |
| Processing of DNA double-strand break ends |  1.766419 |  20 |  12 |
| APC/C-mediated degradation of cell cycle proteins |  1.775820 |  57 |  32 |
| Regulation of mitotic cell cycle |  1.775820 |  57 |  32 |
| Homology Directed Repair |  1.781699 |  34 |  16 |
| Recognition of DNA damage by PCNA-containing replication complex |  1.786301 |  17 |   8 |
| Resolution of AP sites via the multiple-nucleotide patch replacement pathway |  1.786342 |  16 |  10 |
| Mitotic Prophase |  1.787350 |  44 |  17 |
| Interferon Signaling |  1.790467 |  93 |  43 |
| SARS-CoV-1-host interactions |  1.791140 |  62 |  21 |
| Kinesins |  1.791811 |  19 |  10 |
| Orc1 removal from chromatin |  1.792007 |  55 |  32 |
| Infection with Mycobacterium tuberculosis |  1.796039 |  13 |   3 |
| Recruitment of NuMA to mitotic centrosomes |  1.802060 |  27 |  17 |
| Dual Incision in GG-NER |  1.803412 |  19 |   9 |
| MET activates PTK2 signaling |  1.813860 |  14 |   7 |
| SUMOylation of DNA replication proteins |  1.814669 |  22 |  13 |
| The role of GTSE1 in G2/M progression after G2 checkpoint |  1.815821 |  52 |  33 |
| Translesion synthesis by REV1 |  1.817679 |  10 |   8 |
| Translesion synthesis by POLK |  1.817679 |  10 |   8 |
| Translesion synthesis by POLI |  1.817679 |  10 |   8 |
| Translesion Synthesis by POLH |  1.822036 |  12 |   9 |
| Initiation of Nuclear Envelope (NE) Reformation |  1.822563 |  10 |   5 |
| Syndecan interactions |  1.823063 |  15 |   9 |
| APC/C:Cdc20 mediated degradation of Cyclin B |  1.829894 |  10 |   4 |
| G0 and Early G1 |  1.830428 |  12 |   5 |
| DDX58/IFIH1-mediated induction of interferon-alpha/beta |  1.836886 |  31 |  15 |
| Regulation of PLK1 Activity at G2/M Transition |  1.841462 |  36 |  17 |
| Platelet activation, signaling and aggregation |  1.852929 | 109 |  37 |
| Gap-filling DNA repair synthesis and ligation in TC-NER |  1.862584 |  36 |  14 |
| M Phase |  1.866915 | 169 |  91 |
| Removal of the Flap Intermediate |  1.867534 |  10 |   7 |
| DNA Damage Bypass |  1.869623 |  25 |  11 |
| Condensation of Prophase Chromosomes |  1.870740 |  10 |   6 |
| Prefoldin mediated transfer of substrate  to CCT/TriC |  1.874904 |  16 |  12 |
| Cooperation of Prefoldin and TriC/CCT  in actin and tubulin folding |  1.874904 |  16 |  12 |
| Smooth Muscle Contraction |  1.875071 |  19 |  10 |
| Nucleotide Excision Repair |  1.878546 |  60 |  20 |
| Mitotic G2-G2/M phases |  1.878659 |  91 |  53 |
| G2/M Transition |  1.878659 |  91 |  53 |
| Translesion synthesis by Y family DNA polymerases bypasses lesions on DNA template |  1.878682 |  20 |  11 |
| DNA Double-Strand Break Repair |  1.881783 |  47 |  22 |
| Termination of translesion DNA synthesis |  1.885894 |  17 |  10 |
| Fanconi Anemia Pathway |  1.887161 |  10 |   8 |
| Mitotic Metaphase and Anaphase |  1.890600 | 120 |  45 |
| Mitotic Anaphase |  1.890600 | 120 |  45 |
| Formation of tubulin folding intermediates by CCT/TriC |  1.893730 |  10 |  10 |
| Signaling by ALK in cancer |  1.910403 |  58 |  22 |
| Signaling by ALK fusions and activated point mutants |  1.910403 |  58 |  22 |
| PKR-mediated signaling |  1.911955 |  28 |  17 |
| Gap-filling DNA repair synthesis and ligation in GG-NER |  1.914112 |  15 |  10 |
| Global Genome Nucleotide Excision Repair (GG-NER) |  1.918114 |  43 |  16 |
| Non-integrin membrane-ECM interactions |  1.923018 |  26 |  13 |
| Meiosis |  1.924771 |  24 |  16 |
| Telomere C-strand (Lagging Strand) Synthesis |  1.926484 |  14 |   9 |
| DNA Replication Pre-Initiation |  1.931193 |  66 |  29 |
| Antiviral mechanism by IFN-stimulated genes |  1.948182 |  66 |  36 |
| Processive synthesis on the lagging strand |  1.956764 |  11 |   8 |
| Reproduction |  1.958722 |  26 |  17 |
| Mitotic Spindle Checkpoint |  1.967847 |  43 |  15 |
| Lagging Strand Synthesis |  1.968093 |  13 |   9 |
| DNA Repair |  1.969906 | 106 |  42 |
| Separation of Sister Chromatids |  1.972554 |  94 |  58 |
| Assembly of collagen fibrils and other multimeric structures |  1.976660 |  25 |  14 |
| Amplification of signal from the kinetochores |  1.979575 |  39 |  15 |
| Amplification  of signal from unattached  kinetochores via a MAD2  inhibitory signal |  1.979575 |  39 |  15 |
| E2F mediated regulation of DNA replication |  1.983529 |  11 |   6 |
| G2/M DNA damage checkpoint |  1.989389 |  22 |  16 |
| RHO GTPases Activate Formins |  1.999289 |  53 |  24 |
| Collagen biosynthesis and modifying enzymes |  2.010105 |  21 |  15 |
| DNA Replication |  2.019335 |  79 |  35 |
| Collagen formation |  2.021605 |  34 |  19 |
| EML4 and NUDC in mitotic spindle formation |  2.023676 |  46 |  21 |
| S Phase |  2.023840 |  93 |  41 |
| Synthesis of DNA |  2.028324 |  74 |  34 |
| RHO GTPase Effectors |  2.030571 | 125 |  51 |
| Platelet degranulation |  2.042536 |  57 |  28 |
| Telomere Maintenance |  2.050957 |  37 |  19 |
| Extension of Telomeres |  2.054685 |  22 |  16 |
| Mitotic Prometaphase |  2.062503 |  75 |  50 |
| Cell Cycle, Mitotic |  2.065679 | 231 |  89 |
| Response to elevated platelet cytosolic Ca2+ |  2.069555 |  61 |  29 |
| Cell Cycle |  2.079722 | 276 | 106 |
| Activation of ATR in response to replication stress |  2.083787 |  14 |   9 |
| G2/M Checkpoints |  2.098970 |  74 |  43 |
| Extracellular matrix organization |  2.111018 | 111 |  59 |
| Mitotic G1 phase and G1/S transition |  2.113026 |  87 |  27 |
| G1/S Transition |  2.141832 |  81 |  39 |
| Resolution of Sister Chromatid Cohesion |  2.143426 |  50 |  24 |
| Chromosome Maintenance |  2.156192 |  42 |  22 |
| Cell Cycle Checkpoints |  2.172180 | 122 |  56 |
| ECM proteoglycans |  2.270187 |  31 |  20 |
| DNA strand elongation |  2.289468 |  18 |  14 |
| Activation of the pre-replicative complex |  2.292947 |  16 |  11 |



#### Pathways with Significant NES

**Negative NES (Under-represented Pathways)**:

1. **NR1H3 & NR1H2 regulate gene expression linked to cholesterol transport and efflux**:
   - NES: -1.875194
   - Set Size: 17
   - Count: 8
2. **Peroxisomal lipid metabolism**:
   - NES: -1.861756
   - Set Size: 10
   - Count: 8
3. **Metabolism of porphyrins**:
   - NES: -1.849627
   - Set Size: 15
   - Count: 8

These pathways have significant negative NES, indicating they are under-represented in the dataset compared to other pathways.

**Positive NES (Over-represented Pathways)**:

1. **Innate Immune System**:
   - NES: 1.201911
   - Set Size: 406
   - Count: 102
2. **Adaptive Immune System**:
   - NES: 1.275866
   - Set Size: 298
   - Count: 117
3. **Developmental Biology**:
   - NES: 1.288259
   - Set Size: 489
   - Count: 178

These pathways have significant positive NES, indicating they are over-represented in the dataset.

#### Key Observations

1. **Immune System Pathways**:
   - Both innate and adaptive immune system pathways are highly represented, reflecting their crucial role in response mechanisms.
2. **Developmental Biology**:
   - This pathway shows high representation, highlighting its importance in various biological processes and potentially indicating a focus on developmental studies.
3. **Under-represented Pathways**:
   - Pathways related to lipid metabolism, cholesterol transport, and efflux are under-represented, suggesting less focus or a negative impact in these areas.

#### Interpretation and Usage

- **NES Interpretation**:
  - Positive NES indicates over-representation, suggesting these pathways are more active or prevalent in the given condition or dataset.
  - Negative NES indicates under-representation, suggesting these pathways are less active or prevalent.
- **Set Size and Count**:
  - Set Size refers to the total number of genes or components in the pathway.
  - Count refers to the number of genes or components from the dataset that fall into the pathway.



#### Negative NES (Under-represented Pathways):

1. **NR1H3 & NR1H2 Regulation of Gene Expression Linked to Cholesterol Transport and Efflux:**
   - **NES: -1.875194**
   - **Set Size: 17**
   - **Count: 8**
   - **Reason for Focus:** This pathway is significantly under-represented, suggesting that genes involved in cholesterol transport and efflux are less active in the dataset. Given the critical role of cholesterol metabolism in cellular processes and disease, this under-representation could point to important biological insights or potential dysregulation in the condition studied.
2. **Peroxisomal Lipid Metabolism:**
   - **NES: -1.861756**
   - **Set Size: 10**
   - **Count: 8**
   - **Reason for Focus:** Peroxisomes play a vital role in lipid metabolism, including the breakdown of very long-chain fatty acids. An under-representation of this pathway might indicate metabolic disruptions or alterations in lipid processing, which can have significant implications for metabolic health and diseases.
3. **Metabolism of Porphyrins:**
   - **NES: -1.849627**
   - **Set Size: 15**
   - **Count: 8**
   - **Reason for Focus:** Porphyrins are essential components of hemoglobin, and disruptions in their metabolism can lead to disorders such as porphyrias. The under-representation of this pathway could suggest potential issues in hemoglobin synthesis or related metabolic processes, which are critical for understanding various hematological conditions.

#### Positive NES (Over-represented Pathways):

1. **Innate Immune System:**
   - **NES: 1.201911**
   - **Set Size: 406**
   - **Count: 102**
   - **Reason for Focus:** The innate immune system is the body's first line of defense against pathogens. Its over-representation indicates an elevated immune response, which could be associated with inflammation, infection, or other immune-related conditions.
2. **Adaptive Immune System:**
   - **NES: 1.275866**
   - **Set Size: 298**
   - **Count: 117**
   - **Reason for Focus:** The adaptive immune system is responsible for targeted and long-lasting immune responses. Over-representation here suggests an active immune response possibly due to chronic inflammation, autoimmune conditions, or response to infection.
3. **Developmental Biology:**
   - **NES: 1.288259**
   - **Set Size: 489**
   - **Count: 178**
   - **Reason for Focus:** Developmental biology pathways are crucial for the formation and differentiation of tissues and organs. Over-representation in these pathways could indicate active developmental processes, regeneration, or even pathological conditions like cancer, where developmental pathways can be aberrantly activated.

#### Summary

These pathways were selected based on their extreme NES values, which highlight significant deviations from the norm in the dataset. By focusing on both under-represented and over-represented pathways, we can gain insights into potential metabolic disruptions, immune responses, and developmental changes that may be pertinent to the biological condition under investigation. These insights are valuable for forming hypotheses about underlying mechanisms and potential therapeutic targets.



### EOMES

These are the differentially regulated pathways from ReactomePA.

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| Potassium Channels | 1.603122 |  69 |  14 |
| Neuronal System | 1.456785 | 313 | 131 |
| Muscle contraction | 1.437219 | 150 |  62 |
| Transmission across Chemical Synapses | 1.419551 | 208 |  89 |



The plot below shows only upregulated pathways, since none of the enriched pathways was downregulated. These we obtained by running gseGO against the *org.Hs.eg.db*.

![eomes_10wpi](/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_10wpi/gsea_definitive_plots/eomes_10wpi.png)

#### Summary of Findings

The pathways with positive NES values in the Atlantic salmon heart tissue, obtained through ReactomePA and human orthologs conversion, emphasize several key biological processes:

- **Electrophysiology and Ion Transport (Potassium Channels):** Essential for maintaining cardiac rhythm and function.
- **Neuronal Regulation (Neuronal System):** Suggests significant neural control and integration in cardiac function.
- **Cardiac Muscle Activity (Muscle Contraction):** Reflects the fundamental role of muscle contraction in heart mechanics.
- **Neuro-Cardiac Interaction (Transmission Across Chemical Synapses):** Indicates significant synaptic communication affecting cardiac physiology.

These pathways collectively highlight the intricate regulation and active processes within the salmon heart, underscoring the importance of both muscular and neural components in maintaining cardiac health and performance. This information can be crucial for understanding cardiac physiology in salmon, with potential implications for aquaculture and conservation efforts.



























# Spleen

## 10WPI

### DNA vaccine

ReactomePA differentially regulated pathways.

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| Regulation of Complement cascade | -2.241270 |  15 |  11 |
| Amyloid fiber formation | -2.117998 |  27 |   5 |
| Complement cascade | -1.984553 |  19 |  11 |
| Striated Muscle Contraction | -1.884199 |  14 |   5 |
| Neurotransmitter release cycle | -1.783909 |  14 |   6 |
| Peroxisomal protein import | -1.779507 |  27 |  11 |
| GPCR ligand binding | -1.777843 |  68 |  14 |
| Post-translational modification: synthesis of GPI-anchored proteins | -1.766703 |  17 |  12 |
| Metabolic disorders of biological oxidation enzymes | -1.735140 |  10 |   6 |
| Insulin processing | -1.633392 |  15 |   3 |
| Diseases of metabolism | -1.561942 |  94 |  35 |
| Biosynthesis of the N-glycan precursor (dolichol lipid-linked oligosaccharide, LLO) and transfer to a nascent protein | -1.507741 |  36 |  18 |
| Peptide hormone metabolism | -1.503343 |  31 |   6 |
| Platelet degranulation | -1.494170 |  72 |  14 |
| Class A/1 (Rhodopsin-like receptors) | -1.486493 |  51 |  10 |
| Biological oxidations | -1.477513 |  62 |  23 |
| Response to elevated platelet cytosolic Ca2+ | -1.466929 |  76 |  17 |
| Axon guidance |  1.359716 | 303 | 124 |
| Nervous system development |  1.370982 | 309 |  58 |
| Mitotic Metaphase and Anaphase |  1.400823 | 160 |  73 |
| Mitotic Anaphase |  1.400823 | 160 |  73 |
| Cellular Senescence |  1.413021 |  86 |  34 |
| Mitotic G2-G2/M phases |  1.426075 | 121 |  59 |
| G2/M Transition |  1.426075 | 121 |  59 |
| Cell Cycle Checkpoints |  1.428235 | 176 |  77 |
| Cytokine Signaling in Immune system |  1.438643 | 359 | 113 |
| G1/S Transition |  1.448713 | 102 |  52 |
| Death Receptor Signaling |  1.457112 |  84 |  26 |
| Mitotic Prometaphase |  1.479377 | 120 |  48 |
| Signaling by NTRK1 (TRKA) |  1.481876 |  65 |  31 |
| TCF dependent signaling in response to WNT |  1.485156 | 115 |  33 |
| Ub-specific processing proteases |  1.487941 | 115 |  59 |
| Signaling by VEGF |  1.488078 |  77 |  20 |
| Adipogenesis |  1.504584 |  66 |  24 |
| SARS-CoV-2 activates/modulates innate and adaptive immune responses |  1.506024 |  65 |  22 |
| Signaling by the B Cell Receptor (BCR) |  1.509223 |  87 |  47 |
| SARS-CoV-2-host interactions |  1.511130 | 121 |  32 |
| Apoptosis |  1.513401 | 112 |  53 |
| Activation of anterior HOX genes in hindbrain development during early embryogenesis |  1.513965 |  43 |  16 |
| Activation of HOX genes during differentiation |  1.513965 |  43 |  16 |
| Interferon Signaling |  1.517474 | 126 |  40 |
| Mitotic Prophase |  1.517506 |  61 |  21 |
| Hedgehog 'off' state |  1.519717 |  69 |  38 |
| Fc epsilon receptor (FCERI) signaling |  1.524962 |  99 |  50 |
| RHO GTPase Effectors |  1.526557 | 166 |  73 |
| Signaling by ALK in cancer |  1.528765 |  64 |  21 |
| Signaling by ALK fusions and activated point mutants |  1.528765 |  64 |  21 |
| CLEC7A (Dectin-1) signaling |  1.529449 |  75 |  40 |
| MAPK1/MAPK3 signaling |  1.529760 | 143 |  39 |
| S Phase |  1.532588 | 126 |  66 |
| Costimulation by the CD28 family |  1.533874 |  39 |  15 |
| G alpha (12/13) signalling events |  1.533885 |  38 |  17 |
| Respiratory Syncytial Virus Infection Pathway |  1.537826 |  67 |  25 |
| RAC1 GTPase cycle |  1.539850 | 126 |  52 |
| RAF/MAP kinase cascade |  1.540032 | 140 |  39 |
| Estrogen-dependent gene expression |  1.541615 |  67 |  28 |
| M Phase |  1.542473 | 234 | 104 |
| Signaling by Nuclear Receptors |  1.543015 | 126 |  48 |
| Anchoring of the basal body to the plasma membrane |  1.546421 |  48 |  19 |
| Signaling by NOTCH1 |  1.548265 |  43 |  16 |
| MAPK family signaling cascades |  1.549645 | 167 |  79 |
| Regulation of MITF-M-dependent genes involved in apoptosis |  1.552282 |  11 |   7 |
| Adherens junctions interactions |  1.553615 |  25 |   9 |
| Cell Cycle, Mitotic |  1.555357 | 330 | 111 |
| Deubiquitination |  1.558352 | 161 |  50 |
| RHOJ GTPase cycle |  1.564518 |  40 |  19 |
| Cell Cycle |  1.566349 | 398 | 135 |
| C-type lectin receptors (CLRs) |  1.566550 |  86 |  45 |
| RAC3 GTPase cycle |  1.568830 |  61 |  25 |
| Constitutive Signaling by EGFRvIII |  1.569113 |  11 |   4 |
| Signaling by EGFRvIII in Cancer |  1.569113 |  11 |   4 |
| Transcriptional regulation by RUNX3 |  1.569290 |  67 |  35 |
| ESR-mediated signaling |  1.573844 |  99 |  38 |
| Signaling by WNT |  1.574110 | 166 |  82 |
| VEGFA-VEGFR2 Pathway |  1.575299 |  70 |  19 |
| Constitutive Signaling by Ligand-Responsive EGFR Cancer Variants |  1.576634 |  13 |   4 |
| Signaling by EGFR in Cancer |  1.576634 |  13 |   4 |
| Signaling by Ligand-Responsive EGFR Variants in Cancer |  1.576634 |  13 |   4 |
| Interferon gamma signaling |  1.577620 |  35 |  17 |
| Oxidative Stress Induced Senescence |  1.579621 |  44 |  21 |
| RIPK1-mediated regulated necrosis |  1.582296 |  21 |   6 |
| Regulation of necroptotic cell death |  1.582296 |  21 |   6 |
| SARS-CoV-2 Infection |  1.583791 | 178 |  50 |
| Degradation of beta-catenin by the destruction complex |  1.598609 |  65 |  38 |
| Interferon alpha/beta signaling |  1.606583 |  30 |  13 |
| Beta-catenin independent WNT signaling |  1.607205 |  90 |  50 |
| SARS-CoV-1 targets host intracellular signalling and regulatory pathways |  1.610306 |  11 |   7 |
| RHO GTPase cycle |  1.612038 | 291 | 115 |
| Transcriptional regulation by RUNX1 |  1.615748 | 133 |  40 |
| Platelet calcium homeostasis |  1.621076 |  12 |   5 |
| Diseases of signal transduction by growth factor receptors and second messengers |  1.625372 | 284 |  98 |
| Potential therapeutics for SARS |  1.626597 |  68 |  35 |
| SARS-CoV Infections |  1.632015 | 259 |  81 |
| Signaling by Rho GTPases, Miro GTPases and RHOBTB3 |  1.634491 | 422 | 164 |
| Programmed Cell Death |  1.635565 | 132 |  51 |
| RHOC GTPase cycle |  1.636843 |  55 |  26 |
| Opioid Signalling |  1.637278 |  38 |  20 |
| NOTCH1 Intracellular Domain Regulates Transcription |  1.638939 |  29 |  12 |
| RHOF GTPase cycle |  1.641598 |  30 |  16 |
| Signaling by Rho GTPases |  1.645400 | 413 | 161 |
| Condensation of Prophase Chromosomes |  1.649244 |  14 |   7 |
| RND2 GTPase cycle |  1.659249 |  31 |  17 |
| Maternal to zygotic transition (MZT) |  1.660534 |  33 |  18 |
| RAF activation |  1.662111 |  26 |  12 |
| HSF1-dependent transactivation |  1.663105 |  15 |   8 |
| RHOBTB GTPase Cycle |  1.663225 |  25 |   9 |
| RAC2 GTPase cycle |  1.667253 |  66 |  33 |
| Regulation of RUNX1 Expression and Activity |  1.669452 |  13 |   9 |
| HDMs demethylate histones |  1.672443 |  16 |   9 |
| RHOA GTPase cycle |  1.673734 | 101 |  47 |
| Sema4D induced cell migration and growth-cone collapse |  1.681423 |  13 |   6 |
| Mitotic Telophase/Cytokinesis |  1.683839 |  10 |   7 |
| Regulation of MECP2 expression and activity |  1.687468 |  20 |  12 |
| Signalling to ERKs |  1.688122 |  22 |  13 |
| RND1 GTPase cycle |  1.688583 |  28 |  15 |
| Evasion by RSV of host interferon responses |  1.695392 |  13 |   5 |
| Ca2+ pathway |  1.696858 |  26 |  16 |
| Initiation of Nuclear Envelope (NE) Reformation |  1.707072 |  13 |   7 |
| Interleukin-3, Interleukin-5 and GM-CSF signaling |  1.710003 |  31 |   8 |
| RUNX1 regulates genes involved in megakaryocyte differentiation and platelet function |  1.711621 |  32 |  14 |
| RSV-host interactions |  1.712072 |  51 |  23 |
| CDC42 GTPase cycle |  1.712230 | 102 |  44 |
| Trafficking of AMPA receptors |  1.714581 |  15 |   5 |
| Glutamate binding, activation of AMPA receptors and synaptic plasticity |  1.714581 |  15 |   5 |
| Translation of Structural Proteins |  1.717780 |  18 |   8 |
| RHOB GTPase cycle |  1.718858 |  49 |  26 |
| Chromatin modifying enzymes |  1.734738 | 152 |  58 |
| Chromatin organization |  1.734738 | 152 |  58 |
| Signaling by WNT in cancer |  1.737845 |  17 |  11 |
| Apoptotic cleavage of cellular proteins |  1.751430 |  23 |  14 |
| DDX58/IFIH1-mediated induction of interferon-alpha/beta |  1.754826 |  39 |  11 |
| RHO GTPases activate PAKs |  1.755418 |  15 |   6 |
| Maturation of nucleoprotein |  1.755831 |  13 |   5 |
| Effects of PIP2 hydrolysis |  1.757041 |  16 |   5 |
| Apoptotic execution phase |  1.758126 |  30 |  15 |
| Negative regulators of DDX58/IFIH1 signaling |  1.759331 |  20 |   6 |
| Platelet homeostasis |  1.762752 |  32 |  19 |
| RHOBTB2 GTPase cycle |  1.768215 |  16 |   7 |
| E2F mediated regulation of DNA replication |  1.774100 |  16 |   8 |
| Meiotic synapsis |  1.776616 |  18 |   9 |
| RHO GTPases activate CIT |  1.801145 |  13 |   6 |
| SARS-CoV-1 Infection |  1.802731 |  95 |  27 |
| NR1H2 and NR1H3-mediated signaling |  1.813403 |  22 |  14 |
| RND3 GTPase cycle |  1.814530 |  31 |  18 |
| RHOD GTPase cycle |  1.822609 |  39 |  21 |
| Regulation of signaling by CBL |  1.829728 |  19 |   7 |
| Beta-catenin phosphorylation cascade |  1.843554 |  11 |   8 |
| Signaling by AXIN mutants |  1.843554 |  11 |   8 |
| Signaling by CTNNB1 phospho-site mutants |  1.843554 |  11 |   8 |
| Signaling by APC mutants |  1.843554 |  11 |   8 |
| Signaling by AMER1 mutants |  1.843554 |  11 |   8 |
| Signaling by GSK3beta mutants |  1.843554 |  11 |   8 |
| CTNNB1 S33 mutants aren't phosphorylated |  1.843554 |  11 |   8 |
| CTNNB1 S37 mutants aren't phosphorylated |  1.843554 |  11 |   8 |
| CTNNB1 S45 mutants aren't phosphorylated |  1.843554 |  11 |   8 |
| CTNNB1 T41 mutants aren't phosphorylated |  1.843554 |  11 |   8 |
| APC truncation mutants have impaired AXIN binding |  1.843554 |  11 |   8 |
| AXIN missense mutants destabilize the destruction complex |  1.843554 |  11 |   8 |
| Truncations of AMER1 destabilize the destruction complex |  1.843554 |  11 |   8 |
| PKMTs methylate histone lysines |  1.847820 |  30 |  15 |
| NR1H3 & NR1H2 regulate gene expression linked to cholesterol transport and efflux |  1.848806 |  20 |  13 |
| Disassembly of the destruction complex and recruitment of AXIN to the membrane |  1.872158 |  17 |  11 |
| Negative regulation of MAPK pathway |  1.885807 |  26 |  13 |
| TRAF6 mediated IRF7 activation |  1.982506 |  11 |   6 |







### EOMES

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| Signaling by ALK |  1.649581 |  26 |  10 |
| Voltage gated Potassium channels |  1.566912 |  29 |  14 |
| TNFR1-induced NF-kappa-B signaling pathway |  1.552213 |  28 |  13 |
| Sensory processing of sound |  1.545197 |  57 |  29 |
| Signaling by ERBB4 |  1.533551 |  47 |  22 |
| Sensory processing of sound by inner hair cells of the cochlea |  1.530252 |  53 |  26 |
| Neurexins and neuroligins |  1.515727 |  47 |  15 |
| Neuronal System |  1.427153 | 299 | 104 |
| Transmission across Chemical Synapses |  1.360267 | 195 |  74 |
| Signaling by Receptor Tyrosine Kinases |  1.215053 | 432 | 192 |
| Response to elevated platelet cytosolic Ca2+ | -1.818670 | 100 |  14 |
| Platelet degranulation | -1.853988 |  96 |  14 |
| Post-translational protein phosphorylation | -1.883851 |  76 |  13 |
| Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs) | -1.922559 |  86 |  15 |
| Synthesis of glycosylphosphatidylinositol (GPI) | -1.976571 |  16 |   7 |
| Initial triggering of complement | -1.987510 |  13 |   5 |
| Pyrimidine catabolism | -2.036910 |  11 |   5 |
| Diseases of hemostasis | -2.163460 |  12 |   6 |
| Common Pathway of Fibrin Clot Formation | -2.331327 |  11 |   7 |
| Intrinsic Pathway of Fibrin Clot Formation | -2.350957 |  13 |   8 |
| Formation of Fibrin Clot (Clotting Cascade) | -2.661372 |  22 |  10 |
| Regulation of Complement cascade | -2.773996 |  24 |  13 |
| Complement cascade | -2.832269 |  31 |  14 |





### GATA3

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| CD209 (DC-SIGN) signaling |  1.695850 |  17 |   6 |
| Signaling by Hippo |  1.656248 |  16 |  10 |
| ECM proteoglycans |  1.614106 |  47 |  26 |
| MITF-M-dependent gene expression |  1.606173 |  64 |  31 |
| RHOH GTPase cycle |  1.595730 |  30 |  13 |
| MITF-M-regulated melanocyte development |  1.533134 |  88 |  38 |
| RHOB GTPase cycle |  1.532649 |  61 |  35 |
| L1CAM interactions |  1.510738 |  68 |  35 |
| RHOA GTPase cycle |  1.457084 | 125 |  60 |
| RHO GTPase cycle |  1.441012 | 363 | 155 |
| Chromatin modifying enzymes |  1.433822 | 175 |  75 |
| Chromatin organization |  1.433822 | 175 |  75 |
| RAC1 GTPase cycle |  1.418491 | 150 |  68 |
| Nervous system development |  1.393143 | 393 | 132 |
| Axon guidance |  1.375631 | 377 | 123 |
| Signaling by Receptor Tyrosine Kinases |  1.355121 | 374 | 164 |
| Metabolism of nucleotides | -1.619210 |  76 |  19 |
| Diseases of carbohydrate metabolism | -2.005272 |  26 |   5 |
| Post-translational protein phosphorylation | -2.012732 |  60 |  10 |
| Cholesterol biosynthesis | -2.041372 |  24 |   9 |
| Tryptophan catabolism | -2.055674 |  10 |   3 |
| Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs) | -2.057657 |  69 |  13 |
| Response to elevated platelet cytosolic Ca2+ | -2.143473 |  90 |  14 |
| Platelet degranulation | -2.204138 |  86 |  14 |
| Intrinsic Pathway of Fibrin Clot Formation | -2.463928 |  12 |   6 |
| Common Pathway of Fibrin Clot Formation | -2.465193 |  10 |   8 |
| Regulation of Complement cascade | -2.562415 |  21 |   9 |
| Complement cascade | -2.603283 |  26 |   9 |
| Formation of Fibrin Clot (Clotting Cascade) | -2.750159 |  21 |  11 |











































# Interpretation of GSEA Results



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
