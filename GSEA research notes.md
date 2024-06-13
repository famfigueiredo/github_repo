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



### IV-HD

nMDS plots? Can I use them?



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



### IV-LD

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| Striated Muscle Contraction | -2.560859 |  14 |   2 |
| Complement cascade | -2.413080 |  18 |   8 |
| Platelet degranulation | -2.101032 |  70 |   9 |
| Response to elevated platelet cytosolic Ca2+ | -2.082864 |  74 |   9 |
| SARS-CoV Infections |  1.450543 | 255 |  88 |
| Signaling by Rho GTPases, Miro GTPases and RHOBTB3 |  1.545675 | 402 | 183 |
| Signaling by Rho GTPases |  1.549295 | 393 | 181 |
| RHO GTPase cycle |  1.621215 | 278 |  96 |
| CDC42 GTPase cycle |  1.632554 |  97 |  43 |
| RHOBTB2 GTPase cycle |  1.753969 |  16 |   7 |
| RHOB GTPase cycle |  1.757849 |  48 |  30 |
| RHOA GTPase cycle |  1.784755 |  94 |  50 |
| RHOBTB1 GTPase cycle |  1.813043 |  16 |   8 |
| RND3 GTPase cycle |  1.857320 |  30 |  16 |
| RHOBTB GTPase Cycle |  1.877992 |  25 |  12 |





## 4WPC

### DNA vaccine

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| Metabolism of porphyrins |  2.068871 |  16 |   9 |
| L1CAM interactions |  1.995875 |  56 |  21 |
| Heme biosynthesis |  1.961682 |  11 |   6 |
| Mitochondrial iron-sulfur cluster biogenesis |  1.877300 |  11 |   2 |
| The phototransduction cascade |  1.875759 |  11 |   2 |
| Inactivation, recovery and regulation of the phototransduction cascade |  1.875759 |  11 |   2 |
| Regulated Necrosis |  1.839602 |  37 |  16 |
| Iron uptake and transport |  1.771661 |  35 |  10 |
| Programmed Cell Death |  1.726049 | 138 |  37 |
| Apoptosis |  1.607348 | 116 |  27 |
| Axon guidance |  1.457884 | 318 |  80 |
| Metabolism of lipids | -1.412703 | 395 | 100 |
| Diseases of metabolism | -1.645383 | 108 |  42 |
| Diseases of glycosylation | -1.701822 |  59 |  25 |
| Metabolism of steroids | -1.725580 |  69 |  23 |
| Common Pathway of Fibrin Clot Formation | -1.871216 |  10 |   7 |
| Formation of Fibrin Clot (Clotting Cascade) | -1.932278 |  17 |  10 |
| Complement cascade | -1.984173 |  20 |   5 |
| Regulation of Complement cascade | -2.004786 |  16 |   8 |
| Response to elevated platelet cytosolic Ca2+ | -2.028299 |  80 |  17 |
| Insulin processing | -2.043490 |  16 |   6 |
| Platelet degranulation | -2.059360 |  76 |  17 |
| Post-translational protein phosphorylation | -2.217165 |  53 |  13 |
| Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs) | -2.225907 |  60 |  23 |





### EOMES

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| Metabolism of porphyrins |  2.148982 |  16 |   8 |
| Heme biosynthesis |  2.022551 |  11 |   6 |
| Keratinization |  1.996972 |  17 |   4 |
| Formation of the cornified envelope |  1.996972 |  17 |   4 |
| CD209 (DC-SIGN) signaling |  1.912676 |  14 |   5 |
| Generation of second messenger molecules |  1.906672 |  16 |   6 |
| Intrinsic Pathway for Apoptosis |  1.894452 |  31 |  12 |
| Cobalamin (Cbl, vitamin B12) transport and metabolism |  1.838560 |  10 |   2 |
| Transcriptional regulation by the AP-2 (TFAP2) family of transcription factors |  1.802697 |  17 |  12 |
| The phototransduction cascade |  1.802022 |  11 |   4 |
| Inactivation, recovery and regulation of the phototransduction cascade |  1.802022 |  11 |   4 |
| Muscle contraction |  1.779787 |  72 |  31 |
| Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC) |  1.779769 |  75 |  47 |
| L1CAM interactions |  1.758115 |  56 |  17 |
| Formation of a pool of free 40S subunits |  1.724827 |  83 |  50 |
| GTP hydrolysis and joining of the 60S ribosomal subunit |  1.715622 |  92 |  53 |
| FCERI mediated Ca+2 mobilization |  1.709010 |  24 |  11 |
| Eukaryotic Translation Initiation |  1.708851 |  98 |  56 |
| Cap-dependent Translation Initiation |  1.708851 |  98 |  56 |
| Striated Muscle Contraction |  1.706148 |  14 |   7 |
| L13a-mediated translational silencing of Ceruloplasmin expression |  1.695828 |  91 |  52 |
| Eukaryotic Translation Termination |  1.687454 |  73 |  45 |
| Eukaryotic Translation Elongation |  1.685674 |  72 |  53 |
| Peptide chain elongation |  1.683744 |  70 |  44 |
| rRNA processing in the nucleus and cytosol |  1.679949 | 157 | 103 |
| Visual phototransduction |  1.678981 |  29 |   5 |
| Estrogen-dependent gene expression |  1.660740 |  70 |  28 |
| Response of EIF2AK4 (GCN2) to amino acid deficiency |  1.658338 |  81 |  50 |
| FOXO-mediated transcription |  1.652512 |  36 |  19 |
| Major pathway of rRNA processing in the nucleolus and cytosol |  1.651185 | 147 |  95 |
| Selenoamino acid metabolism |  1.649109 |  92 |  46 |
| Death Receptor Signaling |  1.646330 |  95 |  38 |
| NR1H2 and NR1H3-mediated signaling |  1.639777 |  25 |   6 |
| Selenocysteine synthesis |  1.636396 |  73 |  45 |
| Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell |  1.633773 |  25 |  11 |
| Nonsense-Mediated Decay (NMD) |  1.632377 |  89 |  51 |
| Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC) |  1.632377 |  89 |  51 |
| rRNA processing |  1.631788 | 165 | 106 |
| ESR-mediated signaling |  1.589157 | 104 |  41 |
| Signaling by Nuclear Receptors |  1.503294 | 136 |  44 |
| Axon guidance |  1.357830 | 318 |  70 |
| Asparagine N-linked glycosylation | -1.350545 | 190 |  87 |
| DNA Double-Strand Break Repair | -1.415577 |  90 |  30 |
| RHO GTPases Activate Formins | -1.433743 |  88 |  22 |
| Diseases of metabolism | -1.474885 | 108 |  42 |
| Cell Cycle, Mitotic | -1.478818 | 359 |  80 |
| Homology Directed Repair | -1.486808 |  71 |  26 |
| G2/M Checkpoints | -1.497634 | 105 |  35 |
| M Phase | -1.511043 | 255 |  58 |
| DNA strand elongation | -1.522111 |  28 |  11 |
| Cell Cycle | -1.534260 | 434 |  98 |
| Unfolded Protein Response (UPR) | -1.541637 |  70 |  27 |
| Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs) | -1.553914 |  60 |  22 |
| APC/C:Cdc20 mediated degradation of mitotic proteins | -1.554022 |  63 |  27 |
| COPI-dependent Golgi-to-ER retrograde traffic | -1.594287 |  59 |  21 |
| APC/C-mediated degradation of cell cycle proteins | -1.594544 |  73 |  32 |
| Regulation of mitotic cell cycle | -1.594544 |  73 |  32 |
| Separation of Sister Chromatids | -1.600186 | 137 |  32 |
| Diseases of glycosylation | -1.609731 |  59 |  20 |
| Activation of APC/C and APC/C:Cdc20 mediated degradation of mitotic proteins | -1.631899 |  64 |  28 |
| Telomere C-strand (Lagging Strand) Synthesis | -1.634018 |  24 |  12 |
| HDR through Homologous Recombination (HRR) or Single Strand Annealing (SSA) | -1.639731 |  66 |  25 |
| Mitotic Anaphase | -1.644432 | 171 |  41 |
| Synthesis of very long-chain fatty acyl-CoAs | -1.670293 |  14 |   7 |
| Mitotic Metaphase and Anaphase | -1.679851 | 172 |  42 |
| Nuclear Pore Complex (NPC) Disassembly | -1.699523 |  29 |  11 |
| Resolution of D-Loop Structures | -1.710262 |  19 |   8 |
| Diseases of DNA repair | -1.761204 |  32 |  16 |
| Mitotic Prometaphase | -1.763232 | 132 |  39 |
| Golgi Cisternae Pericentriolar Stack Reorganization | -1.766836 |  10 |   9 |
| Resolution of Sister Chromatid Cohesion | -1.777846 |  85 |  26 |
| Impaired BRCA2 binding to RAD51 | -1.781533 |  23 |  10 |
| Cyclin A/B1/B2 associated events during G2/M transition | -1.783024 |  16 |   8 |
| Phosphorylation of the APC/C | -1.787385 |  15 |   6 |
| Resolution of D-loop Structures through Holliday Junction Intermediates | -1.790687 |  18 |   8 |
| HDR through Homologous Recombination (HRR) | -1.793318 |  42 |  17 |
| XBP1(S) activates chaperone genes | -1.798648 |  37 |  12 |
| Glycosaminoglycan metabolism | -1.799723 |  58 |  17 |
| TP53 regulates transcription of several additional cell death genes whose specific roles in p53-dependent apoptosis remain uncertain | -1.801198 |  11 |   4 |
| Presynaptic phase of homologous DNA pairing and strand exchange | -1.803251 |  25 |  11 |
| G2/M DNA damage checkpoint | -1.818246 |  42 |  15 |
| Resolution of D-loop Structures through Synthesis-Dependent Strand Annealing (SDSA) | -1.822299 |  15 |   8 |
| Diseases of DNA Double-Strand Break Repair | -1.828098 |  26 |  13 |
| Defective homologous recombination repair (HRR) due to BRCA2 loss of function | -1.828098 |  26 |  13 |
| Regulation of Complement cascade | -1.840737 |  16 |   9 |
| EML4 and NUDC in mitotic spindle formation | -1.851267 |  77 |  24 |
| Complement cascade | -1.858660 |  20 |  11 |
| IRE1alpha activates chaperones | -1.859516 |  38 |  13 |
| Homologous DNA Pairing and Strand Exchange | -1.870689 |  27 |  14 |
| HDR through Single Strand Annealing (SSA) | -1.894345 |  26 |  12 |
| Diseases associated with glycosaminoglycan metabolism | -1.925171 |  16 |   4 |
| Defective homologous recombination repair (HRR) due to BRCA1 loss of function | -1.931827 |  14 |   8 |
| Defective homologous recombination repair (HRR) due to PALB2 loss of function | -1.931827 |  14 |   8 |
| Defective HDR through Homologous Recombination Repair (HRR) due to PALB2 loss of BRCA1 binding function | -1.931827 |  14 |   8 |
| Defective HDR through Homologous Recombination Repair (HRR) due to PALB2 loss of BRCA2/RAD51/RAD51C binding function | -1.931827 |  14 |   8 |
| Impaired BRCA2 binding to PALB2 | -1.931827 |  14 |   8 |
| Cell Cycle Checkpoints | -1.995867 | 192 |  49 |
| Mitotic Spindle Checkpoint | -2.074298 |  80 |  25 |
| Amplification of signal from the kinetochores | -2.078635 |  68 |  22 |
| Amplification  of signal from unattached  kinetochores via a MAD2  inhibitory signal | -2.078635 |  68 |  22 |



### GATA3

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| Metabolism of porphyrins |  2.202361 | 16 |  9 |
| Heme biosynthesis |  2.074192 | 11 |  7 |
| Unfolded Protein Response (UPR) | -1.716388 | 70 | 26 |
| Glycosaminoglycan metabolism | -1.740617 | 57 | 23 |
| Post-translational protein phosphorylation | -1.781209 | 50 | 19 |
| Sensory processing of sound | -1.877618 | 25 |  5 |
| Plasma lipoprotein remodeling | -1.948396 | 16 |  5 |
| Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs) | -1.950295 | 57 | 18 |
| XBP1(S) activates chaperone genes | -2.023765 | 37 | 11 |
| IRE1alpha activates chaperones | -2.034111 | 38 | 13 |



### IV-HD

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| Metabolism of porphyrins |  2.025030 |  17 |  9 |
| Heme biosynthesis |  1.956092 |  11 |  6 |
| CD209 (DC-SIGN) signaling |  1.905132 |  16 |  4 |
| Keratinization |  1.866491 |  20 |  5 |
| Formation of the cornified envelope |  1.866491 |  20 |  5 |
| rRNA modification in the nucleus and cytosol |  1.790658 |  55 | 34 |
| rRNA processing in the nucleus and cytosol |  1.698858 | 160 | 84 |
| Major pathway of rRNA processing in the nucleolus and cytosol |  1.655784 | 150 | 77 |
| rRNA processing |  1.612426 | 169 | 86 |
| Platelet activation, signaling and aggregation | -1.792341 | 169 | 34 |
| XBP1(S) activates chaperone genes | -1.853295 |  39 | 16 |
| Sensory processing of sound by inner hair cells of the cochlea | -1.864196 |  25 |  8 |
| IRE1alpha activates chaperones | -1.908621 |  40 | 17 |
| Cholesterol biosynthesis | -1.928222 |  19 | 13 |
| Post-translational protein phosphorylation | -2.089087 |  56 | 11 |
| Intrinsic Pathway of Fibrin Clot Formation | -2.117373 |  11 |  5 |
| Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs) | -2.172118 |  63 | 13 |
| Response to elevated platelet cytosolic Ca2+ | -2.179590 |  85 | 25 |
| Common Pathway of Fibrin Clot Formation | -2.191657 |  10 |  8 |
| Platelet degranulation | -2.214696 |  81 | 25 |
| Regulation of Complement cascade | -2.402575 |  19 |  8 |
| Formation of Fibrin Clot (Clotting Cascade) | -2.411237 |  19 | 10 |
| Complement cascade | -2.440427 |  23 |  9 |



### IV-LD

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| Striated Muscle Contraction |  2.125201 |  18 |   7 |
| Muscle contraction |  2.011388 |  87 |  24 |
| Heme biosynthesis |  1.951132 |  11 |   7 |
| Metabolism of porphyrins |  1.946705 |  17 |   9 |
| Transcription of E2F targets under negative control by DREAM complex |  1.837145 |  16 |   9 |
| G0 and Early G1 |  1.825592 |  23 |  13 |
| TP53 Regulates Transcription of Cell Death Genes |  1.822671 |  25 |   9 |
| Aberrant regulation of mitotic cell cycle due to RB1 defects |  1.819579 |  28 |  13 |
| Formation of tubulin folding intermediates by CCT/TriC |  1.800396 |  13 |   8 |
| snRNP Assembly |  1.797102 |  42 |  19 |
| Metabolism of non-coding RNA |  1.797102 |  42 |  19 |
| Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell |  1.791607 |  27 |  12 |
| HCMV Early Events |  1.775093 |  50 |  19 |
| Prefoldin mediated transfer of substrate  to CCT/TriC |  1.771623 |  19 |  12 |
| Cooperation of Prefoldin and TriC/CCT  in actin and tubulin folding |  1.771623 |  19 |  12 |
| Transport of the SLBP Dependant Mature mRNA |  1.763561 |  27 |  14 |
| Transport of the SLBP independent Mature mRNA |  1.760120 |  26 |  14 |
| Diseases of mitotic cell cycle |  1.756363 |  30 |  10 |
| tRNA processing in the nucleus |  1.745351 |  45 |  19 |
| Polymerase switching |  1.737348 |  12 |   6 |
| Leading Strand Synthesis |  1.737348 |  12 |   6 |
| Transcription of E2F targets under negative control by p107 (RBL1) and p130 (RBL2) in complex with HDAC1 |  1.727841 |  14 |   9 |
| SUMOylation of SUMOylation proteins |  1.706052 |  27 |  13 |
| Extension of Telomeres |  1.699124 |  41 |  18 |
| S Phase |  1.697397 | 142 |  49 |
| tRNA processing |  1.694374 |  85 |  38 |
| Synthesis of DNA |  1.684453 | 105 |  54 |
| Switching of origins to a post-replicative state |  1.661922 |  80 |  42 |
| CDK-mediated phosphorylation and removal of Cdc6 |  1.641596 |  63 |  31 |
| APC/C:Cdh1 mediated degradation of Cdc20 and other APC/C:Cdh1 targeted proteins in late mitosis/early G1 |  1.628545 |  64 |  32 |
| Mitochondrial protein degradation |  1.616809 |  80 |  33 |
| Cdc20:Phospho-APC/C mediated degradation of Cyclin A |  1.599286 |  63 |  35 |
| Respiratory electron transport |  1.566409 |  82 |  31 |
| Factors involved in megakaryocyte development and platelet production |  1.549917 |  95 |  28 |
| G1/S Transition |  1.540709 | 112 |  57 |
| Separation of Sister Chromatids |  1.508940 | 144 |  44 |
| Mitotic G1 phase and G1/S transition |  1.505042 | 123 |  42 |
| Aerobic respiration and respiratory electron transport |  1.459795 | 157 |  57 |
| Metabolism of lipids | -1.545778 | 490 | 129 |
| Metabolism of vitamins and cofactors | -1.624853 | 131 |  36 |
| Diseases of hemostasis | -1.672322 |  11 |   5 |
| Response to elevated platelet cytosolic Ca2+ | -1.708887 |  88 |  21 |
| ABC transporters in lipid homeostasis | -1.723269 |  10 |   4 |
| Post-translational protein phosphorylation | -1.725481 |  56 |  12 |
| Retinoid metabolism and transport | -1.731038 |  22 |   8 |
| Platelet degranulation | -1.741427 |  84 |  21 |
| Diseases of Immune System | -1.772153 |  14 |   4 |
| Diseases associated with the TLR signaling cascade | -1.772153 |  14 |   4 |
| IRS-related events triggered by IGF1R | -1.773082 |  25 |   9 |
| Sensory processing of sound by outer hair cells of the cochlea | -1.775328 |  22 |  10 |
| Signaling by Type 1 Insulin-like Growth Factor 1 Receptor (IGF1R) | -1.786118 |  26 |   9 |
| IGF1R signaling cascade | -1.786118 |  26 |   9 |
| Metabolism of fat-soluble vitamins | -1.795160 |  26 |  10 |
| Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs) | -1.801806 |  63 |  14 |
| Intrinsic Pathway of Fibrin Clot Formation | -1.845032 |  12 |   6 |
| Common Pathway of Fibrin Clot Formation | -1.897044 |  10 |   6 |
| Formation of Fibrin Clot (Clotting Cascade) | -2.077895 |  20 |  10 |
| Regulation of Complement cascade | -2.146441 |  21 |  10 |
| Complement cascade | -2.181194 |  25 |  11 |



Need to do spleen 6WPC, liver (10WPI and 4WPC), and head-kidney (10WPI and 4WPC).

Scripts in use are *vaccine_contrasts.R*, *deseq_modelling_lymphoid.R*, and *gsea_new_and_definitive_4wpc.R*.



# Liver

10/06/2024

Created sampleTables and DGE objects for liver and head-kidney.

Wrote vaccine constrasts and extracted results tables for 10wpi and 4wpc.



*There were 53 pathways for which P-values were not calculated properly due to unbalanced (positive and negative) gene-level statistic values. For such pathways pval, padj, NES, log2err are set to NA.*

This kept happening, so I added an *nPermSimple = 100000* flag to *gseGO* and *nPermSimple = 100000* to *gsePathway*.

## 10WPI

### DNA vaccine

Signaling by interleukins shows upregulated.

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| rRNA modification in the nucleus and cytosol | -2.772010 |  54 |  39 |
| rRNA processing | -2.453593 | 160 |  65 |
| rRNA processing in the nucleus and cytosol | -2.407239 | 154 |  61 |
| Major pathway of rRNA processing in the nucleolus and cytosol | -2.239763 | 145 |  53 |
| tRNA processing | -2.218107 |  64 |  29 |
| tRNA modification in the nucleus and cytosol | -2.091645 |  21 |  15 |
| PERK regulates gene expression | -2.048805 |  24 |  11 |
| ATF4 activates genes in response to endoplasmic reticulum  stress | -1.987790 |  19 |  12 |
| Formation of tubulin folding intermediates by CCT/TriC | -1.887469 |  10 |   4 |
| Transcriptional regulation by small RNAs | -1.805579 |  40 |  20 |
| RNA Polymerase I Promoter Escape | -1.793334 |  23 |  15 |
| RNA Pol II CTD phosphorylation and interaction with CE during HIV infection | -1.781392 |  23 |  13 |
| RNA Pol II CTD phosphorylation and interaction with CE | -1.781392 |  23 |  13 |
| mRNA Capping | -1.750957 |  25 |  14 |
| Unfolded Protein Response (UPR) | -1.743634 |  70 |  28 |
| TAK1-dependent IKK and NF-kappa-B activation | -1.741067 |  18 |   7 |
| Gene Silencing by RNA | -1.733720 |  52 |  24 |
| snRNP Assembly | -1.727473 |  36 |  18 |
| Metabolism of non-coding RNA | -1.727473 |  36 |  18 |
| KSRP (KHSRP) binds and destabilizes mRNA | -1.727124 |  11 |   8 |
| Formation of TC-NER Pre-Incision Complex | -1.720563 |  41 |  22 |
| HIV Life Cycle | -1.711677 | 107 |  47 |
| Late Phase of HIV Life Cycle | -1.710021 | 101 |  45 |
| mRNA decay by 3' to 5' exoribonuclease | -1.707305 |  13 |   9 |
| Nucleotide-binding domain, leucine rich repeat containing receptor (NLR) signaling pathways | -1.669651 |  30 |   9 |
| tRNA processing in the nucleus | -1.667210 |  38 |  16 |
| RNA Polymerase III Chain Elongation | -1.663220 |  12 |   7 |
| Viral Messenger RNA Synthesis | -1.653304 |  34 |  17 |
| Inflammasomes | -1.653136 |  12 |   4 |
| SUMOylation of RNA binding proteins | -1.653020 |  29 |  15 |
| Mitochondrial translation elongation | -1.648479 |  70 |  33 |
| Nuclear import of Rev protein | -1.645935 |  24 |  13 |
| RNA polymerase II transcribes snRNA genes | -1.624966 |  58 |  27 |
| Mitochondrial translation | -1.612842 |  75 |  34 |
| Interleukin-4 and Interleukin-13 signaling | -1.607890 |  39 |  19 |
| Formation of the Early Elongation Complex | -1.596029 |  28 |  14 |
| Formation of the HIV-1 Early Elongation Complex | -1.596029 |  28 |  14 |
| tRNA Aminoacylation | -1.587923 |  33 |  18 |
| Mitochondrial translation initiation | -1.573332 |  70 |  32 |
| Mitochondrial translation termination | -1.563183 |  69 |  31 |
| HCMV Late Events | -1.551581 |  40 |  19 |
| HIV Infection | -1.361084 | 170 |  60 |
| Processing of Capped Intron-Containing Pre-mRNA | -1.346116 | 222 |  69 |
| Cell Cycle |  1.314339 | 363 | 127 |
| Diseases of signal transduction by growth factor receptors and second messengers |  1.316558 | 266 |  77 |
| Axon guidance |  1.344675 | 289 | 118 |
| SARS-CoV Infections |  1.349494 | 247 |  74 |
| Signaling by Receptor Tyrosine Kinases |  1.349672 | 257 |  71 |
| Nervous system development |  1.390382 | 298 | 122 |
| Cell Cycle, Mitotic |  1.405877 | 306 | 110 |
| Signaling by Nuclear Receptors |  1.426774 | 132 |  33 |
| Hemostasis |  1.432707 | 270 | 101 |
| MAPK family signaling cascades |  1.442332 | 159 |  41 |
| Signaling by WNT |  1.476225 | 150 |  54 |
| Transmission across Chemical Synapses |  1.478456 |  70 |  24 |
| Mitotic Metaphase and Anaphase |  1.501343 | 146 |  58 |
| Mitotic Anaphase |  1.501343 | 146 |  58 |
| Cell Cycle Checkpoints |  1.507114 | 156 |  61 |
| GTP hydrolysis and joining of the 60S ribosomal subunit |  1.514329 |  91 |  47 |
| L13a-mediated translational silencing of Ceruloplasmin expression |  1.517851 |  90 |  44 |
| Signaling by GPCR |  1.521175 | 154 |  54 |
| GPCR downstream signalling |  1.529661 | 147 |  52 |
| Beta-catenin independent WNT signaling |  1.541891 |  80 |  33 |
| Muscle contraction |  1.559529 |  51 |  15 |
| Response of EIF2AK4 (GCN2) to amino acid deficiency |  1.562418 |  80 |  42 |
| Formation of a pool of free 40S subunits |  1.583754 |  82 |  41 |
| Integration of energy metabolism |  1.584293 |  38 |  14 |
| RHOG GTPase cycle |  1.586506 |  51 |  16 |
| Platelet activation, signaling and aggregation |  1.590911 | 125 |  39 |
| RAC2 GTPase cycle |  1.594233 |  59 |  19 |
| RND1 GTPase cycle |  1.594742 |  27 |   9 |
| Parasite infection |  1.607800 |  42 |  22 |
| Leishmania phagocytosis |  1.607800 |  42 |  22 |
| FCGR3A-mediated phagocytosis |  1.607800 |  42 |  22 |
| ECM proteoglycans |  1.610333 |  30 |  10 |
| Neurotransmitter receptors and postsynaptic signal transmission |  1.624353 |  53 |  20 |
| Nonsense-Mediated Decay (NMD) |  1.626122 |  87 |  48 |
| Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC) |  1.626122 |  87 |  48 |
| Semaphorin interactions |  1.627524 |  35 |  15 |
| GPVI-mediated activation cascade |  1.628625 |  20 |   9 |
| PLC beta mediated events |  1.629346 |  16 |   7 |
| RAF activation |  1.631989 |  21 |  10 |
| Myogenesis |  1.640567 |  13 |   7 |
| Eukaryotic Translation Termination |  1.642670 |  72 |  40 |
| Separation of Sister Chromatids |  1.648325 | 115 |  48 |
| MAP2K and MAPK activation |  1.650966 |  23 |  11 |
| Meiotic synapsis |  1.654659 |  16 |   7 |
| Oncogenic MAPK signaling |  1.662560 |  44 |  18 |
| Potential therapeutics for SARS |  1.665109 |  63 |  27 |
| Gene and protein expression by JAK-STAT signaling after Interleukin-12 stimulation |  1.671597 |  26 |  12 |
| CDC42 GTPase cycle |  1.674033 |  86 |  22 |
| Sensory Perception |  1.690622 |  60 |  28 |
| Collagen biosynthesis and modifying enzymes |  1.693640 |  23 |   7 |
| Degradation of cysteine and homocysteine |  1.694713 |  10 |   4 |
| RHO GTPase cycle |  1.708945 | 264 |  84 |
| TRAF6 mediated IRF7 activation |  1.709489 |  11 |   6 |
| NR1H2 and NR1H3-mediated signaling |  1.711277 |  26 |   7 |
| Viral mRNA Translation |  1.717562 |  69 |  38 |
| Mitotic Prometaphase |  1.720485 | 107 |  42 |
| Extracellular matrix organization |  1.730960 | 117 |  37 |
| SRP-dependent cotranslational protein targeting to membrane |  1.731785 |  88 |  48 |
| Signaling by high-kinase activity BRAF mutants |  1.732551 |  20 |  11 |
| Notch-HLH transcription pathway |  1.733203 |  18 |   9 |
| Interleukin-12 signaling |  1.735706 |  32 |  14 |
| Cholesterol biosynthesis |  1.735860 |  23 |  14 |
| Mitotic Spindle Checkpoint |  1.748624 |  62 |  28 |
| RHO GTPases Activate Formins |  1.748758 |  68 |  32 |
| EPHB-mediated forward signaling |  1.754624 |  23 |  13 |
| G alpha (12/13) signalling events |  1.760218 |  36 |  20 |
| Interleukin-12 family signaling |  1.761629 |  35 |  15 |
| G-protein beta:gamma signalling |  1.770287 |  14 |   8 |
| EML4 and NUDC in mitotic spindle formation |  1.788346 |  61 |  29 |
| L1CAM interactions |  1.793282 |  53 |  18 |
| Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC) |  1.799826 |  74 |  43 |
| RHOD GTPase cycle |  1.802760 |  39 |  19 |
| Signaling by Rho GTPases, Miro GTPases and RHOBTB3 |  1.818813 | 382 | 120 |
| RHOA GTPase cycle |  1.821890 |  87 |  26 |
| Peptide chain elongation |  1.824288 |  69 |  38 |
| Signaling by Rho GTPases |  1.832108 | 372 | 117 |
| Selenocysteine synthesis |  1.833169 |  72 |  39 |
| RHOC GTPase cycle |  1.847258 |  50 |  17 |
| RND3 GTPase cycle |  1.859686 |  30 |  15 |
| Signaling by RAF1 mutants |  1.860056 |  23 |  13 |
| Eukaryotic Translation Elongation |  1.861776 |  71 |  40 |
| Amplification of signal from the kinetochores |  1.872992 |  54 |  26 |
| Amplification  of signal from unattached  kinetochores via a MAD2  inhibitory signal |  1.872992 |  54 |  26 |
| Cardiac conduction |  1.921692 |  26 |   8 |
| Signaling by moderate kinase activity BRAF mutants |  1.930290 |  25 |  14 |
| Signaling by RAS mutants |  1.930290 |  25 |  14 |
| Paradoxical activation of RAF signaling by kinase inactive BRAF |  1.930290 |  25 |  14 |
| Signaling downstream of RAS mutants |  1.930290 |  25 |  14 |
| RHO GTPase Effectors |  1.934879 | 146 |  72 |
| RHOQ GTPase cycle |  1.973296 |  38 |  14 |
| Resolution of Sister Chromatid Cohesion |  1.978613 |  69 |  36 |
| RHOB GTPase cycle |  2.086523 |  42 |  16 |

### EOMES

Signaling by interleukins shows upregulated.

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| rRNA modification in the nucleus and cytosol | -2.629134 |  57 |  38 |
| Mitochondrial translation | -2.510041 |  82 |  52 |
| Mitochondrial translation elongation | -2.498719 |  76 |  49 |
| Mitochondrial translation termination | -2.465083 |  76 |  48 |
| Mitochondrial translation initiation | -2.452399 |  76 |  47 |
| Cobalamin (Cbl, vitamin B12) transport and metabolism | -2.190745 |  15 |   4 |
| Activation of Matrix Metalloproteinases | -2.169753 |  15 |   2 |
| Mitochondrial tRNA aminoacylation | -2.086415 |  20 |  14 |
| tRNA Aminoacylation | -2.004145 |  40 |  28 |
| rRNA processing | -1.981503 | 174 |  60 |
| rRNA processing in the mitochondrion | -1.974662 |  10 |   8 |
| Removal of the Flap Intermediate | -1.917535 |  13 |   8 |
| Processive synthesis on the lagging strand | -1.902146 |  14 |   8 |
| Kidney development | -1.896871 |  19 |   7 |
| Extension of Telomeres | -1.887146 |  45 |  20 |
| rRNA processing in the nucleus and cytosol | -1.882236 | 164 |  54 |
| Translation | -1.855510 | 250 |  87 |
| Chromosome Maintenance | -1.846192 |  80 |  32 |
| Formation of ATP by chemiosmotic coupling | -1.827672 |  12 |  10 |
| tRNA processing | -1.795627 |  89 |  35 |
| RNA Polymerase III Chain Elongation | -1.783728 |  16 |  10 |
| Telomere Maintenance | -1.765739 |  62 |  30 |
| Major pathway of rRNA processing in the nucleolus and cytosol | -1.724958 | 154 |  46 |
| E2F mediated regulation of DNA replication | -1.701744 |  16 |   8 |
| DNA strand elongation | -1.700657 |  30 |  17 |
| Mitochondrial protein import | -1.698303 |  54 |  24 |
| Cytosolic tRNA aminoacylation | -1.680270 |  23 |  13 |
| Cristae formation | -1.672193 |  24 |  14 |
| Hemostasis |  1.246532 | 399 | 152 |
| Signaling by Interleukins |  1.278670 | 277 | 102 |
| Diseases of signal transduction by growth factor receptors and second messengers |  1.299682 | 349 | 123 |
| MAPK family signaling cascades |  1.305292 | 225 |  78 |
| Nervous system development |  1.321692 | 403 | 130 |
| Axon guidance |  1.324522 | 388 |  91 |
| Intracellular signaling by second messengers |  1.348471 | 230 |  72 |
| MAPK1/MAPK3 signaling |  1.354531 | 193 |  68 |
| RAF/MAP kinase cascade |  1.362480 | 190 |  67 |
| Signaling by Nuclear Receptors |  1.368205 | 181 |  61 |
| ESR-mediated signaling |  1.369432 | 129 |  43 |
| Signaling by Receptor Tyrosine Kinases |  1.370449 | 383 | 184 |
| Chromatin modifying enzymes |  1.398007 | 174 |  77 |
| Chromatin organization |  1.398007 | 174 |  77 |
| MITF-M-regulated melanocyte development |  1.419168 |  85 |  34 |
| Sensory Perception |  1.426483 | 112 |  56 |
| Epigenetic regulation of gene expression |  1.428471 | 116 |  37 |
| Adipogenesis |  1.431171 |  88 |  26 |
| RHO GTPase cycle |  1.437374 | 374 | 144 |
| Signaling by NTRKs |  1.447580 |  99 |  42 |
| Platelet activation, signaling and aggregation |  1.450350 | 186 |  81 |
| Cell-Cell communication |  1.457617 |  96 |  50 |
| Activation of NMDA receptors and postsynaptic events |  1.462013 |  55 |  22 |
| Leishmania infection |  1.466123 | 117 |  48 |
| Parasitic Infection Pathways |  1.466123 | 117 |  48 |
| Sphingolipid metabolism |  1.469900 |  77 |  23 |
| RAC1 GTPase cycle |  1.470595 | 157 |  80 |
| Signaling by TGF-beta Receptor Complex |  1.474899 |  77 |  25 |
| Sensory processing of sound by inner hair cells of the cochlea |  1.475068 |  45 |  27 |
| Sensory processing of sound |  1.475908 |  48 |  28 |
| Death Receptor Signaling |  1.479910 | 118 |  42 |
| Inositol phosphate metabolism |  1.484833 |  40 |  20 |
| Regulation of insulin secretion |  1.497956 |  44 |  17 |
| RHOC GTPase cycle |  1.503305 |  64 |  27 |
| Signaling by TGFB family members |  1.509704 |  97 |  32 |
| NRAGE signals death through JNK |  1.510095 |  45 |  21 |
| Signaling by NTRK1 (TRKA) |  1.515012 |  85 |  37 |
| CDC42 GTPase cycle |  1.520473 | 130 |  73 |
| Ca2+ pathway |  1.523515 |  42 |  25 |
| Platelet Aggregation (Plug Formation) |  1.524694 |  32 |  18 |
| Constitutive Signaling by Aberrant PI3K in Cancer |  1.526696 |  42 |  22 |
| Pre-NOTCH Transcription and Translation |  1.527756 |  30 |  12 |
| G alpha (s) signalling events |  1.529277 |  69 |  36 |
| TGF-beta receptor signaling activates SMADs |  1.532312 |  38 |  15 |
| Potential therapeutics for SARS |  1.537806 |  77 |  42 |
| RUNX1 regulates genes involved in megakaryocyte differentiation and platelet function |  1.538548 |  33 |  16 |
| Nephrin family interactions |  1.539485 |  19 |   9 |
| NOTCH1 Intracellular Domain Regulates Transcription |  1.542501 |  38 |  13 |
| Sema4D induced cell migration and growth-cone collapse |  1.543528 |  17 |   6 |
| Downregulation of SMAD2/3:SMAD4 transcriptional activity |  1.544392 |  25 |   8 |
| Transcriptional regulation by RUNX3 |  1.544888 |  79 |  24 |
| Collagen degradation |  1.545320 |  35 |  18 |
| Post NMDA receptor activation events |  1.552244 |  46 |  19 |
| Neurotransmitter release cycle |  1.552616 |  31 |  15 |
| MAP2K and MAPK activation |  1.556544 |  28 |  16 |
| Signaling by GPCR |  1.556555 | 302 | 133 |
| MET promotes cell motility |  1.557673 |  27 |  17 |
| PI3K events in ERBB2 signaling |  1.560929 |  11 |   6 |
| HDMs demethylate histones |  1.562580 |  17 |   8 |
| NR1H2 and NR1H3-mediated signaling |  1.567284 |  36 |  17 |
| Regulation of CDH11 Expression and Function |  1.570134 |  21 |  13 |
| Adrenaline,noradrenaline inhibits insulin secretion |  1.571358 |  14 |   6 |
| G alpha (q) signalling events |  1.575934 | 100 |  41 |
| Adenylate cyclase inhibitory pathway |  1.576767 |  11 |   6 |
| G alpha (12/13) signalling events |  1.579723 |  55 |  25 |
| Oncogenic MAPK signaling |  1.583526 |  60 |  28 |
| Formation of WDR5-containing histone-modifying complexes |  1.584590 |  35 |  17 |
| MITF-M-dependent gene expression |  1.585433 |  63 |  27 |
| Synthesis of IP3 and IP4 in the cytosol |  1.587776 |  21 |  12 |
| L1CAM interactions |  1.590207 |  72 |  39 |
| Uptake and actions of bacterial toxins |  1.590709 |  18 |   9 |
| Interleukin-4 and Interleukin-13 signaling |  1.593440 |  58 |  27 |
| SUMOylation of intracellular receptors |  1.595570 |  25 |  14 |
| Elastic fibre formation |  1.597056 |  35 |  13 |
| Regulation of beta-cell development |  1.600646 |  24 |  10 |
| G alpha (i) signalling events |  1.602068 | 115 |  51 |
| GPVI-mediated activation cascade |  1.608847 |  27 |  13 |
| Polo-like kinase mediated events |  1.609342 |  14 |   6 |
| RHOA GTPase cycle |  1.615091 | 127 |  51 |
| Glucagon-like Peptide-1 (GLP1) regulates insulin secretion |  1.617192 |  21 |  12 |
| Non-integrin membrane-ECM interactions |  1.622412 |  42 |  17 |
| Interferon gamma signaling |  1.626992 |  40 |  24 |
| FOXO-mediated transcription of cell cycle genes |  1.638358 |  13 |   8 |
| Sphingolipid catabolism |  1.638525 |  10 |   5 |
| GPCR downstream signalling |  1.639542 | 275 | 125 |
| Signaling by NODAL |  1.641624 |  11 |   8 |
| Neurotransmitter receptors and postsynaptic signal transmission |  1.642256 | 107 |  52 |
| Signaling by high-kinase activity BRAF mutants |  1.643647 |  25 |  16 |
| Integrin signaling |  1.645587 |  24 |  17 |
| ECM proteoglycans |  1.645736 |  46 |  22 |
| Bacterial Infection Pathways |  1.646065 |  50 |  18 |
| HSF1-dependent transactivation |  1.648852 |  17 |   8 |
| Anti-inflammatory response favouring Leishmania parasite infection |  1.655132 |  49 |  21 |
| Leishmania parasite growth and survival |  1.655132 |  49 |  21 |
| Interaction between L1 and Ankyrins |  1.657753 |  13 |   7 |
| Integration of energy metabolism |  1.660360 |  71 |  29 |
| RHOQ GTPase cycle |  1.660770 |  50 |  27 |
| NR1H3 & NR1H2 regulate gene expression linked to cholesterol transport and efflux |  1.665302 |  28 |  14 |
| Aquaporin-mediated transport |  1.665910 |  29 |  14 |
| RHO GTPases activate CIT |  1.666394 |  15 |   8 |
| G alpha (z) signalling events |  1.667774 |  29 |  16 |
| FCGR3A-mediated IL10 synthesis |  1.671625 |  29 |  13 |
| Long-term potentiation |  1.675547 |  10 |   7 |
| Molecules associated with elastic fibres |  1.678294 |  28 |  13 |
| Signaling by BMP |  1.682198 |  22 |  10 |
| ADORA2B mediated anti-inflammatory cytokines production |  1.682606 |  24 |  13 |
| RUNX3 regulates NOTCH signaling |  1.688974 |  11 |   7 |
| RHOJ GTPase cycle |  1.690497 |  47 |  26 |
| Neuronal System |  1.697738 | 206 |  98 |
| Regulation of MITF-M-dependent genes involved in cell cycle and proliferation |  1.704705 |  13 |   7 |
| Voltage gated Potassium channels |  1.706622 |  13 |   9 |
| NCAM signaling for neurite out-growth |  1.708120 |  40 |  17 |
| Notch-HLH transcription pathway |  1.712344 |  22 |   9 |
| PKA activation in glucagon signalling |  1.717508 |  13 |   9 |
| Signaling by RAF1 mutants |  1.722014 |  29 |  17 |
| Transmission across Chemical Synapses |  1.727897 | 144 |  68 |
| Glucagon signaling in metabolic regulation |  1.729771 |  20 |  12 |
| TP53 Regulates Transcription of Genes Involved in G1 Cell Cycle Arrest |  1.731729 |  12 |   5 |
| Regulation of gene expression in late stage (branching morphogenesis) pancreatic bud precursor cells |  1.741200 |  11 |   7 |
| PKA-mediated phosphorylation of CREB |  1.744077 |  16 |   8 |
| Myogenesis |  1.750783 |  20 |  13 |
| GPER1 signaling |  1.764746 |  27 |  15 |
| Signaling by moderate kinase activity BRAF mutants |  1.768421 |  32 |  19 |
| Signaling by RAS mutants |  1.768421 |  32 |  19 |
| Paradoxical activation of RAF signaling by kinase inactive BRAF |  1.768421 |  32 |  19 |
| Signaling downstream of RAS mutants |  1.768421 |  32 |  19 |
| Vasopressin regulates renal water homeostasis via Aquaporins |  1.769831 |  24 |  13 |
| PKA activation |  1.771608 |  14 |   8 |
| Phase 0 - rapid depolarisation |  1.777203 |  12 |   7 |
| G-protein mediated events |  1.881901 |  39 |  24 |
| Opioid Signalling |  1.900240 |  55 |  30 |
| DAG and IP3 signaling |  1.907237 |  30 |  18 |
| PLC beta mediated events |  1.910717 |  36 |  23 |
| Calmodulin induced events |  1.970152 |  24 |  16 |
| CaM pathway |  1.970152 |  24 |  16 |
| Ca-dependent events |  2.019479 |  26 |  18 |

### GATA3

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| rRNA processing in the mitochondrion | -2.379128 |  10 |   9 |
| Metabolism of Angiotensinogen to Angiotensins | -2.372706 |  10 |   3 |
| rRNA modification in the nucleus and cytosol | -2.293831 |  57 |  34 |
| Mitochondrial tRNA aminoacylation | -2.253875 |  20 |   9 |
| Mitochondrial translation | -1.936104 |  82 |  41 |
| Mitochondrial translation initiation | -1.887474 |  76 |  38 |
| Mitochondrial translation elongation | -1.875701 |  76 |  39 |
| Mitochondrial translation termination | -1.836147 |  76 |  38 |
| rRNA processing | -1.637892 | 174 |  62 |
| rRNA processing in the nucleus and cytosol | -1.436040 | 164 |  55 |
| Signaling by Receptor Tyrosine Kinases |  1.237291 | 397 | 172 |
| RHO GTPase cycle |  1.256677 | 382 | 199 |
| GPCR downstream signalling |  1.330100 | 309 | 135 |
| Extracellular matrix organization |  1.330252 | 196 |  93 |
| Transmission across Chemical Synapses |  1.330932 | 159 |  75 |
| Signaling by GPCR |  1.346442 | 339 | 135 |
| RAC1 GTPase cycle |  1.362416 | 159 | 104 |
| Neurotransmitter receptors and postsynaptic signal transmission |  1.406526 | 116 |  62 |
| PLC beta mediated events |  1.549178 |  36 |  29 |
| ECM proteoglycans |  1.578233 |  50 |  18 |
| Non-integrin membrane-ECM interactions |  1.581483 |  44 |  17 |
| Myogenesis |  1.661284 |  23 |  13 |

### IV-HD

### IV-LD






## 4WPC

### DNA vaccine

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| SRP-dependent cotranslational protein targeting to membrane |  2.746160 |  88 |  49 |
| Viral mRNA Translation |  2.562936 |  69 |  43 |
| Eukaryotic Translation Elongation |  2.560727 |  71 |  44 |
| Peptide chain elongation |  2.526258 |  69 |  42 |
| Selenocysteine synthesis |  2.465764 |  73 |  44 |
| Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC) |  2.459142 |  74 |  43 |
| Formation of a pool of free 40S subunits |  2.438645 |  82 |  43 |
| L13a-mediated translational silencing of Ceruloplasmin expression |  2.412538 |  90 |  52 |
| GTP hydrolysis and joining of the 60S ribosomal subunit |  2.389440 |  91 |  52 |
| Eukaryotic Translation Initiation |  2.260711 |  97 |  47 |
| Cap-dependent Translation Initiation |  2.260711 |  97 |  47 |
| Eukaryotic Translation Termination |  2.230120 |  73 |  39 |
| Response of EIF2AK4 (GCN2) to amino acid deficiency |  2.215357 |  80 |  42 |
| Formation of the ternary complex, and subsequently, the 43S complex |  2.175515 |  44 |  23 |
| SARS-CoV-1 modulates host translation machinery |  2.158552 |  29 |  15 |
| Nonsense-Mediated Decay (NMD) |  2.128311 |  89 |  43 |
| Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC) |  2.128311 |  89 |  43 |
| Formation of ATP by chemiosmotic coupling |  2.069048 |  12 |   7 |
| Ribosomal scanning and start codon recognition |  2.001157 |  50 |  23 |
| Translation initiation complex formation |  1.997164 |  49 |  23 |
| Activation of the mRNA upon binding of the cap-binding complex and eIFs, and subsequent binding to 43S |  1.985995 |  50 |  23 |
| SARS-CoV-2 modulates host translation machinery |  1.942443 |  38 |  16 |
| Selenoamino acid metabolism |  1.923872 |  93 |  45 |
| Mitochondrial Fatty Acid Beta-Oxidation |  1.897481 |  27 |  10 |
| Cellular response to starvation |  1.895605 | 117 |  51 |
| Defects in vitamin and cofactor metabolism |  1.882877 |  17 |  11 |
| Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins. |  1.826944 |  96 |  45 |
| Regulation of expression of SLITs and ROBOs |  1.720593 | 132 |  49 |
| Fatty acid metabolism |  1.692700 | 108 |  36 |
| Respiratory electron transport |  1.650109 |  80 |  36 |
| Signaling by ROBO receptors |  1.603262 | 161 |  54 |
| Translation |  1.598214 | 237 |  84 |
| Aerobic respiration and respiratory electron transport |  1.559828 | 151 |  61 |
| Metabolism of amino acids and derivatives |  1.373951 | 264 |  70 |
| Cytokine Signaling in Immune system | -1.333362 | 352 | 103 |
| RHO GTPase cycle | -1.341345 | 305 |  81 |
| Signaling by Interleukins | -1.450596 | 224 |  63 |
| Hemostasis | -1.510786 | 305 | 103 |
| RAC1 GTPase cycle | -1.542343 | 128 |  36 |
| Platelet activation, signaling and aggregation | -1.544752 | 145 |  46 |
| SLC-mediated transmembrane transport | -1.547072 | 104 |  36 |
| Chromosome Maintenance | -1.593103 |  66 |  29 |
| VEGFA-VEGFR2 Pathway | -1.608699 |  70 |  26 |
| RAC3 GTPase cycle | -1.653358 |  64 |  20 |
| Mitochondrial tRNA aminoacylation | -1.655229 |  16 |  11 |
| Transport of inorganic cations/anions and amino acids/oligopeptides | -1.657739 |  45 |  15 |
| RHO GTPases Activate NADPH Oxidases | -1.741948 |  17 |   8 |
| Cyclin D associated events in G1 | -1.765330 |  28 |  14 |
| G1 Phase | -1.765330 |  28 |  14 |
| Cell surface interactions at the vascular wall | -1.768807 |  66 |  22 |
| RAF-independent MAPK1/3 activation | -1.791937 |  14 |   6 |
| Amino acid transport across the plasma membrane | -1.807066 |  18 |   8 |
| Interleukin-4 and Interleukin-13 signaling | -1.848273 |  43 |  24 |



### EOMES

There is a marked difference between the number of enriched GO terms and enriched pathways from ReactomePA.

While 257 GO terms were enriched with *gseGO*, 400 pathways were tagged in *ReactomePA*.

In this case, I arranged them in ascending NES, as the downregulated pathways were more interesting.



| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| Interleukin-4 and Interleukin-13 signaling | -2.262873 |  38 |  22 |
| GPER1 signaling | -2.096039 |  15 |   7 |
| Fcgamma receptor (FCGR) dependent phagocytosis | -2.089789 |  52 |  19 |
| FCGR3A-mediated IL10 synthesis | -2.083127 |  16 |   9 |
| Integrin signaling | -2.057657 |  19 |   9 |
| G alpha (s) signalling events | -2.040387 |  24 |  10 |
| Platelet Aggregation (Plug Formation) | -2.038182 |  22 |  11 |
| RHO GTPase Effectors | -2.018924 | 135 |  49 |
| Cell surface interactions at the vascular wall | -2.016042 |  53 |  20 |
| GPVI-mediated activation cascade | -1.987619 |  19 |   8 |
| Hemostasis | -1.979233 | 261 |  87 |
| Interferon gamma signaling | -1.945199 |  29 |  18 |
| Integrin cell surface interactions | -1.937112 |  36 |  15 |
| Leishmania infection | -1.929806 |  79 |  32 |
| Parasitic Infection Pathways | -1.929806 |  79 |  32 |
| Signaling by CSF1 (M-CSF) in myeloid cells | -1.928144 |  18 |  10 |
| Signal transduction by L1 | -1.918460 |  15 |   8 |
| Anti-inflammatory response favouring Leishmania parasite infection | -1.912667 |  29 |  17 |
| Leishmania parasite growth and survival | -1.912667 |  29 |  17 |
| Parasite infection | -1.911090 |  39 |  13 |
| Leishmania phagocytosis | -1.911090 |  39 |  13 |
| FCGR3A-mediated phagocytosis | -1.911090 |  39 |  13 |
| RHO GTPases Activate Formins | -1.883883 |  61 |  24 |
| RAC1 GTPase cycle | -1.871119 | 102 |  38 |
| GPCR downstream signalling | -1.864692 | 128 |  54 |
| Signaling by GPCR | -1.854200 | 134 |  57 |
| Semaphorin interactions | -1.840478 |  35 |  17 |
| RAC2 GTPase cycle | -1.835605 |  54 |  20 |
| EML4 and NUDC in mitotic spindle formation | -1.834083 |  56 |  23 |
| FLT3 Signaling | -1.832020 |  22 |  11 |
| Processive synthesis on the C-strand of the telomere | -1.823964 |  12 |   7 |
| Regulation of actin dynamics for phagocytic cup formation | -1.822127 |  39 |  12 |
| MAPK targets/ Nuclear events mediated by MAP kinases | -1.804437 |  19 |  12 |
| Metabolism of nitric oxide: NOS3 activation and regulation | -1.803939 |  11 |   4 |
| Platelet activation, signaling and aggregation | -1.802865 | 121 |  37 |
| Extra-nuclear estrogen signaling | -1.801008 |  29 |   8 |
| VEGFA-VEGFR2 Pathway | -1.793506 |  57 |  28 |
| DNA strand elongation | -1.792021 |  24 |  15 |
| Collagen degradation | -1.783196 |  28 |   9 |
| Amplification of signal from the kinetochores | -1.782581 |  49 |  19 |
| Amplification  of signal from unattached  kinetochores via a MAD2  inhibitory signal | -1.782581 |  49 |  19 |
| Interleukin-3, Interleukin-5 and GM-CSF signaling | -1.781255 |  26 |  13 |
| Mitotic Prometaphase | -1.780742 |  96 |  43 |
| RAC3 GTPase cycle | -1.776373 |  53 |  18 |
| Signaling by VEGF | -1.771072 |  61 |  30 |
| Toll Like Receptor 7/8 (TLR7/8) Cascade | -1.770673 |  46 |  25 |
| MyD88 dependent cascade initiated on endosome | -1.770673 |  46 |  25 |
| G1/S-Specific Transcription | -1.764554 |  14 |   9 |
| Toll Like Receptor 9 (TLR9) Cascade | -1.764530 |  49 |  26 |
| Muscle contraction | -1.761417 |  48 |  20 |
| Extracellular matrix organization | -1.760416 | 108 |  32 |
| Signaling by Interleukins | -1.759556 | 197 |  71 |
| Cytokine Signaling in Immune system | -1.757916 | 308 | 105 |
| Antigen activates B Cell Receptor (BCR) leading to generation of second messengers | -1.756051 |  16 |   9 |
| EPHA-mediated growth cone collapse | -1.753683 |  10 |   8 |
| MAP kinase activation | -1.752945 |  36 |  19 |
| Processive synthesis on the lagging strand | -1.749877 |  12 |   8 |
| Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell | -1.744003 |  19 |  11 |
| TRAF6 mediated induction of NFkB and MAP kinases upon TLR7/8 or 9 activation | -1.736590 |  45 |  24 |
| Mitotic Spindle Checkpoint | -1.735430 |  57 |  21 |
| SMAD2/SMAD3:SMAD4 heterotrimer regulates transcription | -1.728680 |  19 |   9 |
| Cell-extracellular matrix interactions | -1.727461 |  10 |   5 |
| MyD88-independent TLR4 cascade | -1.725437 |  54 |  26 |
| TRIF (TICAM1)-mediated TLR4 signaling | -1.725437 |  54 |  26 |
| Lagging Strand Synthesis | -1.721041 |  16 |  11 |
| Telomere C-strand (Lagging Strand) Synthesis | -1.719281 |  19 |  11 |
| G-protein beta:gamma signalling | -1.709217 |  12 |   5 |
| Fanconi Anemia Pathway | -1.699367 |  17 |  11 |
| FLT3 signaling in disease | -1.698379 |  17 |   8 |
| G alpha (12/13) signalling events | -1.691303 |  32 |  20 |
| Toll Like Receptor 10 (TLR10) Cascade | -1.687917 |  45 |  23 |
| Toll Like Receptor 5 (TLR5) Cascade | -1.687917 |  45 |  23 |
| MyD88 cascade initiated on plasma membrane | -1.687917 |  45 |  23 |
| RHO GTPases Activate NADPH Oxidases | -1.686677 |  12 |   7 |
| Removal of the Flap Intermediate from the C-strand | -1.684042 |  11 |   6 |
| Regulation of CDH11 Expression and Function | -1.682777 |  14 |   8 |
| Regulation of Homotypic Cell-Cell Adhesion | -1.682777 |  14 |   8 |
| Regulation of Expression and Function of Type II Classical Cadherins | -1.682777 |  14 |   8 |
| Translesion synthesis by Y family DNA polymerases bypasses lesions on DNA template | -1.681754 |  23 |  12 |
| Toll Like Receptor 4 (TLR4) Cascade | -1.681057 |  64 |  29 |
| Smooth Muscle Contraction | -1.679834 |  20 |  10 |
| ERK/MAPK targets | -1.679352 |  13 |   8 |
| Toll Like Receptor 3 (TLR3) Cascade | -1.676690 |  52 |  25 |
| Resolution of Sister Chromatid Cohesion | -1.674593 |  63 |  24 |
| Signaling by TGFB family members | -1.672626 |  63 |  23 |
| Sema4D induced cell migration and growth-cone collapse | -1.670284 |  11 |   7 |
| Signaling by RAF1 mutants | -1.668716 |  22 |   8 |
| Signaling by high-kinase activity BRAF mutants | -1.668466 |  19 |   8 |
| Nuclear Envelope Breakdown | -1.667025 |  39 |  11 |
| Regulation of signaling by CBL | -1.666557 |  15 |   9 |
| DAP12 interactions | -1.664445 |  14 |   5 |
| DAP12 signaling | -1.664445 |  14 |   5 |
| Removal of the Flap Intermediate | -1.664370 |  11 |   8 |
| Signaling by TGF-beta Receptor Complex | -1.663848 |  53 |  21 |
| Uptake and actions of bacterial toxins | -1.659376 |  12 |   2 |
| Signaling by ERBB2 | -1.658931 |  24 |  15 |
| Interleukin-17 signaling | -1.646827 |  37 |  19 |
| ECM proteoglycans | -1.645292 |  29 |   9 |
| Signaling by moderate kinase activity BRAF mutants | -1.635207 |  24 |   8 |
| Signaling by RAS mutants | -1.635207 |  24 |   8 |
| Paradoxical activation of RAF signaling by kinase inactive BRAF | -1.635207 |  24 |   8 |
| Signaling downstream of RAS mutants | -1.635207 |  24 |   8 |
| Depolymerization of the Nuclear Lamina | -1.629480 |  12 |   5 |
| Activation of ATR in response to replication stress | -1.626161 |  23 |  17 |
| Signaling by Rho GTPases | -1.625659 | 346 | 101 |
| Syndecan interactions | -1.622241 |  17 |   7 |
| Elastic fibre formation | -1.617773 |  17 |   9 |
| Role of phospholipids in phagocytosis | -1.613648 |  10 |   5 |
| FCERI mediated MAPK activation | -1.610961 |  18 |  12 |
| Cell Cycle, Mitotic | -1.604524 | 282 | 106 |
| Activation of Matrix Metalloproteinases | -1.604274 |  12 |   4 |
| FCERI mediated Ca+2 mobilization | -1.601557 |  13 |   9 |
| Constitutive Signaling by Ligand-Responsive EGFR Cancer Variants | -1.597643 |  13 |   6 |
| Signaling by EGFR in Cancer | -1.597643 |  13 |   6 |
| Signaling by Ligand-Responsive EGFR Variants in Cancer | -1.597643 |  13 |   6 |
| Basigin interactions | -1.588490 |  13 |   5 |
| Potential therapeutics for SARS | -1.587549 |  60 |  20 |
| Termination of translesion DNA synthesis | -1.586378 |  18 |  10 |
| Signaling by Rho GTPases, Miro GTPases and RHOBTB3 | -1.582732 | 354 | 101 |
| Signaling by FGFR1 | -1.581086 |  17 |   8 |
| HDR through Homologous Recombination (HRR) | -1.579300 |  27 |  13 |
| ADORA2B mediated anti-inflammatory cytokines production | -1.578305 |  12 |   6 |
| Notch-HLH transcription pathway | -1.577985 |  17 |   7 |
| Miscellaneous transport and binding events | -1.573612 |  12 |   6 |
| Interleukin-2 family signaling | -1.572867 |  16 |   7 |
| Signaling by Receptor Tyrosine Kinases | -1.571384 | 243 |  89 |
| Cell junction organization | -1.568458 |  43 |  17 |
| Gap-filling DNA repair synthesis and ligation in GG-NER | -1.566696 |  16 |   9 |
| MAP2K and MAPK activation | -1.565525 |  22 |   8 |
| G0 and Early G1 | -1.564954 |  15 |  10 |
| Mitotic Prophase | -1.562224 |  56 |  13 |
| DDX58/IFIH1-mediated induction of interferon-alpha/beta | -1.560684 |  36 |  13 |
| Signaling by CSF3 (G-CSF) | -1.556781 |  20 |   8 |
| Translesion Synthesis by POLH | -1.555846 |  14 |   8 |
| L1CAM interactions | -1.555558 |  50 |  17 |
| Signaling by FGFR3 | -1.551163 |  13 |   7 |
| Signaling by FGFR4 | -1.551163 |  13 |   7 |
| Death Receptor Signaling | -1.548836 |  77 |  30 |
| MyD88:MAL(TIRAP) cascade initiated on plasma membrane | -1.545719 |  48 |  24 |
| Toll Like Receptor TLR1:TLR2 Cascade | -1.545719 |  48 |  24 |
| Toll Like Receptor TLR6:TLR2 Cascade | -1.545719 |  48 |  24 |
| Toll Like Receptor 2 (TLR2) Cascade | -1.545719 |  48 |  24 |
| DARPP-32 events | -1.545263 |  13 |   5 |
| TNFR1-induced NF-kappa-B signaling pathway | -1.544483 |  16 |   8 |
| Constitutive Signaling by EGFRvIII | -1.543304 |  11 |   5 |
| Signaling by EGFRvIII in Cancer | -1.543304 |  11 |   5 |
| G-protein mediated events | -1.542183 |  17 |  11 |
| Cell-Cell communication | -1.540384 |  57 |  20 |
| Signaling by NTRK2 (TRKB) | -1.537126 |  11 |   7 |
| Downregulation of ERBB2 signaling | -1.536510 |  12 |   7 |
| Amino acid transport across the plasma membrane | -1.536036 |  15 |   6 |
| Signaling by FGFR | -1.533343 |  39 |  16 |
| RHOB GTPase cycle | -1.529351 |  40 |  22 |
| NRAGE signals death through JNK | -1.527677 |  27 |  17 |
| Signaling by FLT3 fusion proteins | -1.527151 |  12 |   8 |
| Chromosome Maintenance | -1.526339 |  55 |  19 |
| Extension of Telomeres | -1.525929 |  31 |  19 |
| VEGFR2 mediated vascular permeability | -1.525305 |  18 |  12 |
| Resolution of Abasic Sites (AP sites) | -1.520702 |  23 |  11 |
| Signaling by ALK in cancer | -1.514907 |  61 |  26 |
| Signaling by ALK fusions and activated point mutants | -1.514907 |  61 |  26 |
| Homology Directed Repair | -1.507547 |  47 |  19 |
| Nuclear events stimulated by ALK signaling in cancer | -1.506139 |  20 |  11 |
| Inactivation of CSF3 (G-CSF) signaling | -1.506066 |  16 |   7 |
| Degradation of the extracellular matrix | -1.505435 |  54 |  13 |
| Platelet degranulation | -1.503440 |  65 |  14 |
| Cell Cycle | -1.495341 | 334 | 121 |
| p75 NTR receptor-mediated signalling | -1.492332 |  47 |  24 |
| Toll-like Receptor Cascades | -1.490329 |  76 |  31 |
| Assembly of collagen fibrils and other multimeric structures | -1.488964 |  22 |   7 |
| Interferon alpha/beta signaling | -1.487745 |  26 |  14 |
| Response to elevated platelet cytosolic Ca2+ | -1.486235 |  69 |  14 |
| Deactivation of the beta-catenin transactivating complex | -1.485505 |  19 |   9 |
| DNA Damage Bypass | -1.485166 |  29 |  13 |
| Factors involved in megakaryocyte development and platelet production | -1.484846 |  71 |  22 |
| Opioid Signalling | -1.479913 |  28 |  15 |
| Signaling by SCF-KIT | -1.476386 |  25 |  11 |
| DNA Repair | -1.465392 | 150 |  40 |
| Transcriptional activity of SMAD2/SMAD3:SMAD4 heterotrimer | -1.464608 |  29 |  10 |
| Negative regulation of the PI3K/AKT network | -1.464134 |  39 |  23 |
| Signaling by BRAF and RAF1 fusions | -1.461940 |  32 |   9 |
| AURKA Activation by TPX2 | -1.459428 |  35 |  10 |
| Cell death signalling via NRAGE, NRIF and NADE | -1.457175 |  38 |  20 |
| TGF-beta receptor signaling activates SMADs | -1.456626 |  25 |  12 |
| RHO GTPase cycle | -1.450703 | 245 |  64 |
| Pre-NOTCH Transcription and Translation | -1.448507 |  22 |  15 |
| Base Excision Repair | -1.445906 |  27 |  11 |
| Interferon Signaling | -1.445318 | 110 |  39 |
| Negative regulation of MAPK pathway | -1.444211 |  22 |  13 |
| NOTCH1 Intracellular Domain Regulates Transcription | -1.439979 |  27 |  10 |
| SUMOylation of DNA damage response and repair proteins | -1.439858 |  43 |   9 |
| Non-integrin membrane-ECM interactions | -1.438765 |  27 |  11 |
| Rab regulation of trafficking | -1.436890 |  67 |  17 |
| Signaling by FGFR2 | -1.436754 |  35 |  14 |
| G alpha (q) signalling events | -1.435915 |  48 |  21 |
| PI5P, PP2A and IER3 Regulate PI3K/AKT Signaling | -1.424922 |  37 |  22 |
| Signaling by MET | -1.414098 |  39 |  20 |
| SARS-CoV-2 activates/modulates innate and adaptive immune responses | -1.406048 |  60 |  12 |
| DNA Double-Strand Break Repair | -1.399285 |  63 |  22 |
| Generic Transcription Pathway | -1.396730 | 473 | 164 |
| CDC42 GTPase cycle | -1.393908 |  80 |  33 |
| Diseases of signal transduction by growth factor receptors and second messengers | -1.391798 | 256 |  71 |
| M Phase | -1.377932 | 201 |  66 |
| SUMO E3 ligases SUMOylate target proteins | -1.372799 |  92 |  23 |
| Cilium Assembly | -1.357340 |  71 |  15 |
| ESR-mediated signaling | -1.340824 |  91 |  38 |
| SUMOylation | -1.336488 |  95 |  23 |
| Adaptive Immune System | -1.335264 | 341 | 112 |
| Mitotic G1 phase and G1/S transition | -1.324119 |  96 |  40 |
| Transcriptional Regulation by TP53 | -1.301197 | 208 |  78 |
| Innate Immune System | -1.296696 | 480 | 125 |
| Diseases of metabolism |  1.316776 | 115 |  44 |
| Respiratory electron transport |  1.353684 |  74 |  34 |
| Aerobic respiration and respiratory electron transport |  1.380650 | 140 |  52 |
| Mitochondrial protein degradation |  1.384518 |  73 |  28 |
| Autodegradation of the E3 ubiquitin ligase COP1 |  1.412268 |  43 |  24 |
| Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs) |  1.419517 |  62 |  31 |
| Stabilization of p53 |  1.443104 |  44 |  25 |
| Nonsense-Mediated Decay (NMD) |  1.459617 |  86 |  45 |
| Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC) |  1.459617 |  86 |  45 |
| Metabolism of lipids |  1.463174 | 381 |  95 |
| Signaling by ROBO receptors |  1.470582 | 150 |  76 |
| Golgi-to-ER retrograde transport |  1.476496 |  79 |  39 |
| SARS-CoV-1 modulates host translation machinery |  1.485222 |  29 |  15 |
| ABC transporter disorders |  1.489243 |  59 |  30 |
| COPI-dependent Golgi-to-ER retrograde traffic |  1.490854 |  52 |  17 |
| Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins. |  1.495957 |  90 |  41 |
| Protein localization |  1.498337 | 106 |  33 |
| Disorders of transmembrane transporters |  1.502778 | 104 |  48 |
| Cellular response to starvation |  1.505566 | 112 |  54 |
| Hh mutants are degraded by ERAD |  1.510001 |  46 |  25 |
| Biosynthesis of the N-glycan precursor (dolichol lipid-linked oligosaccharide, LLO) and transfer to a nascent protein |  1.511547 |  39 |  18 |
| Complement cascade |  1.512167 |  26 |   5 |
| Plasma lipoprotein assembly, remodeling, and clearance |  1.524190 |  43 |  18 |
| rRNA processing in the nucleus and cytosol |  1.534583 | 151 |  75 |
| rRNA processing |  1.535307 | 157 |  61 |
| Defective CFTR causes cystic fibrosis |  1.541313 |  50 |  27 |
| Fatty acid metabolism |  1.545336 |  92 |  18 |
| ABC-family proteins mediated transport |  1.549515 |  66 |  32 |
| Post-translational protein phosphorylation |  1.554096 |  54 |  27 |
| COPII-mediated vesicle transport |  1.562566 |  38 |  16 |
| Lysine catabolism |  1.565716 |  12 |   7 |
| Regulation of Glucokinase by Glucokinase Regulatory Protein |  1.567169 |  22 |   2 |
| Defective TPR may confer susceptibility towards thyroid papillary carcinoma (TPC) |  1.567169 |  22 |   2 |
| Cytochrome P450 - arranged by substrate type |  1.575572 |  21 |  12 |
| Metabolism of vitamins and cofactors |  1.578185 | 103 |  37 |
| Hh mutants abrogate ligand secretion |  1.578914 |  47 |  26 |
| Major pathway of rRNA processing in the nucleolus and cytosol |  1.581543 | 142 |  72 |
| Intra-Golgi and retrograde Golgi-to-ER traffic |  1.592137 | 124 |  57 |
| Activation of the mRNA upon binding of the cap-binding complex and eIFs, and subsequent binding to 43S |  1.603380 |  49 |  27 |
| Translation initiation complex formation |  1.618414 |  48 |  27 |
| Maturation of spike protein |  1.626331 |  19 |   7 |
| Glutathione conjugation |  1.627009 |  13 |   6 |
| Retrograde transport at the Trans-Golgi-Network |  1.627279 |  35 |  14 |
| Hedgehog ligand biogenesis |  1.632777 |  49 |  28 |
| Sensory Perception |  1.633976 |  57 |  17 |
| Ribosomal scanning and start codon recognition |  1.645786 |  49 |  28 |
| Visual phototransduction |  1.679929 |  34 |  12 |
| Regulation of expression of SLITs and ROBOs |  1.681093 | 128 |  71 |
| Fatty acyl-CoA biosynthesis |  1.689937 |  23 |   7 |
| Response of EIF2AK4 (GCN2) to amino acid deficiency |  1.703857 |  78 |  46 |
| Phase II - Conjugation of compounds |  1.721793 |  38 |  14 |
| Metabolism of water-soluble vitamins and cofactors |  1.722364 |  63 |  28 |
| Formation of ATP by chemiosmotic coupling |  1.740552 |  12 |  11 |
| Mitochondrial translation |  1.757240 |  71 |  33 |
| Metabolism of amino acids and derivatives |  1.792006 | 246 | 120 |
| Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC) |  1.796156 |  73 |  42 |
| Transport to the Golgi and subsequent modification |  1.823875 | 101 |  44 |
| Mitochondrial translation elongation |  1.828171 |  67 |  33 |
| Mitochondrial translation initiation |  1.830089 |  66 |  19 |
| Mitochondrial translation termination |  1.839537 |  66 |  19 |
| Nucleotide biosynthesis |  1.841813 |  13 |   7 |
| Cytosolic tRNA aminoacylation |  1.842301 |  21 |  11 |
| Phase I - Functionalization of compounds |  1.842301 |  34 |  19 |
| Asparagine N-linked glycosylation |  1.864055 | 172 |  74 |
| Biological oxidations |  1.886675 |  77 |  32 |
| Eukaryotic Translation Initiation |  1.913272 |  95 |  57 |
| Cap-dependent Translation Initiation |  1.913272 |  95 |  57 |
| Formation of the ternary complex, and subsequently, the 43S complex |  1.917351 |  43 |  21 |
| Selenoamino acid metabolism |  1.919174 |  90 |  49 |
| Defects in vitamin and cofactor metabolism |  1.919440 |  14 |  11 |
| Eukaryotic Translation Elongation |  1.920686 |  70 |  42 |
| Cristae formation |  1.928859 |  24 |  13 |
| COPI-mediated anterograde transport |  1.929901 |  60 |  31 |
| Peptide chain elongation |  1.931277 |  68 |  40 |
| Eukaryotic Translation Termination |  1.957626 |  71 |  47 |
| The canonical retinoid cycle in rods (twilight vision) |  1.978377 |  10 |   8 |
| L13a-mediated translational silencing of Ceruloplasmin expression |  2.003959 |  88 |  53 |
| ER to Golgi Anterograde Transport |  2.008324 |  88 |  42 |
| GTP hydrolysis and joining of the 60S ribosomal subunit |  2.015188 |  89 |  54 |
| Selenocysteine synthesis |  2.036777 |  70 |  41 |
| Purine ribonucleoside monophosphate biosynthesis |  2.048663 |  10 |   7 |
| Viral mRNA Translation |  2.074703 |  68 |  41 |
| Formation of a pool of free 40S subunits |  2.149929 |  80 |  49 |
| Translation |  2.201753 | 223 | 121 |
| Regulation of cholesterol biosynthesis by SREBP (SREBF) |  2.294680 |  45 |  17 |
| Metabolism of steroids |  2.298662 |  86 |  32 |
| SRP-dependent cotranslational protein targeting to membrane |  2.401320 |  87 |  58 |
| Activation of gene expression by SREBF (SREBP) |  2.409997 |  35 |  17 |
| Cholesterol biosynthesis |  2.573984 |  23 |  18 |



### GATA3

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| Transport of inorganic cations/anions and amino acids/oligopeptides | -1.929944 |  44 |  16 |
| RHO GTPases Activate NADPH Oxidases | -1.885335 |  15 |   8 |
| Cell surface interactions at the vascular wall | -1.877235 |  65 |  30 |
| Amino acid transport across the plasma membrane | -1.855301 |  18 |  10 |
| Platelet Aggregation (Plug Formation) | -1.799747 |  25 |  12 |
| Interleukin-4 and Interleukin-13 signaling | -1.785542 |  42 |  25 |
| Signaling by FGFR1 | -1.783707 |  24 |  11 |
| Binding and Uptake of Ligands by Scavenger Receptors | -1.728694 |  23 |   5 |
| Negative regulation of FGFR1 signaling | -1.724222 |  13 |   5 |
| Basigin interactions | -1.718560 |  15 |   7 |
| Signaling by GPCR | -1.712427 | 167 |  56 |
| GPCR downstream signalling | -1.708012 | 157 |  48 |
| Integrin signaling | -1.704380 |  21 |  10 |
| Signaling by CSF1 (M-CSF) in myeloid cells | -1.702810 |  21 |  11 |
| RAC1 GTPase cycle | -1.702467 | 120 |  46 |
| Cyclin D associated events in G1 | -1.698183 |  27 |  12 |
| G1 Phase | -1.698183 |  27 |  12 |
| Signaling by FLT3 fusion proteins | -1.687607 |  13 |   9 |
| Signaling by ALK in cancer | -1.674926 |  65 |  37 |
| Signaling by ALK fusions and activated point mutants | -1.674926 |  65 |  37 |
| PRC2 methylates histones and DNA | -1.667802 |  16 |   7 |
| SLC-mediated transmembrane transport | -1.667667 | 100 |  26 |
| Class A/1 (Rhodopsin-like receptors) | -1.664491 |  47 |  13 |
| Role of LAT2/NTAL/LAB on calcium mobilization | -1.661828 |  10 |   8 |
| RAC3 GTPase cycle | -1.660552 |  60 |  23 |
| HSF1 activation | -1.657144 |  11 |   2 |
| FLT3 signaling in disease | -1.651469 |  18 |  10 |
| Signaling by FGFR3 | -1.638074 |  18 |   7 |
| Signaling by FGFR4 | -1.638074 |  18 |   7 |
| Transcriptional activity of SMAD2/SMAD3:SMAD4 heterotrimer | -1.629110 |  37 |  13 |
| RHOB GTPase cycle | -1.623881 |  45 |  22 |
| Downstream signaling of activated FGFR2 | -1.620808 |  11 |   7 |
| Signaling by FGFR1 in disease | -1.619984 |  18 |  13 |
| Hemostasis | -1.618392 | 289 |  96 |
| RAF-independent MAPK1/3 activation | -1.615505 |  13 |  10 |
| Signaling by ERBB2 ECD mutants | -1.615229 |  13 |   8 |
| Syndecan interactions | -1.613136 |  17 |   6 |
| Signaling by ERBB2 in Cancer | -1.607487 |  14 |   8 |
| Signaling by ERBB2 KD Mutants | -1.607487 |  14 |   8 |
| G alpha (q) signalling events | -1.603283 |  58 |  23 |
| Platelet activation, signaling and aggregation | -1.603022 | 134 |  47 |
| Signaling by ALK | -1.599328 |  17 |  11 |
| GPVI-mediated activation cascade | -1.595708 |  22 |  11 |
| ROS and RNS production in phagocytes | -1.590936 |  18 |   9 |
| Signaling by FGFR | -1.581584 |  46 |  17 |
| Integrin cell surface interactions | -1.576134 |  39 |  14 |
| Downstream signaling of activated FGFR3 | -1.575990 |  10 |   6 |
| Downstream signaling of activated FGFR4 | -1.575990 |  10 |   6 |
| Signaling by TGF-beta Receptor Complex | -1.571323 |  63 |  25 |
| Epigenetic regulation of gene expression | -1.565754 | 104 |  38 |
| Signaling by Type 1 Insulin-like Growth Factor 1 Receptor (IGF1R) | -1.565482 |  20 |  10 |
| IGF1R signaling cascade | -1.565482 |  20 |  10 |
| Interferon gamma signaling | -1.549480 |  30 |  15 |
| SMAD2/SMAD3:SMAD4 heterotrimer regulates transcription | -1.540638 |  25 |  12 |
| Interleukin-37 signaling | -1.539840 |  10 |   8 |
| G alpha (12/13) signalling events | -1.535679 |  39 |  17 |
| Signaling by CSF3 (G-CSF) | -1.529843 |  22 |  12 |
| GPCR ligand binding | -1.529233 |  61 |  14 |
| G alpha (i) signalling events | -1.520189 |  61 |  13 |
| rRNA modification in the nucleus and cytosol | -1.512065 |  54 |  26 |
| Fcgamma receptor (FCGR) dependent phagocytosis | -1.502643 |  59 |  18 |
| Signaling by Interleukins | -1.501554 | 216 |  86 |
| Diseases of signal transduction by growth factor receptors and second messengers | -1.492363 | 279 |  82 |
| Extracellular matrix organization | -1.475897 | 124 |  30 |
| Signaling by Receptor Tyrosine Kinases | -1.437949 | 277 |  76 |
| RHO GTPase cycle | -1.426424 | 289 |  78 |
| CDC42 GTPase cycle | -1.417491 |  97 |  34 |
| Cytokine Signaling in Immune system | -1.333779 | 338 | 112 |
| Metabolism of lipids |  1.362471 | 436 | 114 |
| SLC transporter disorders |  1.500487 |  53 |  13 |
| Diseases of metabolism |  1.515116 | 131 |  36 |
| Complex I biogenesis |  1.525126 |  42 |  17 |
| Metabolism of amino acids and derivatives |  1.530165 | 260 |  99 |
| Ubiquitin Mediated Degradation of Phosphorylated Cdc25A |  1.530984 |  42 |  25 |
| p53-Independent DNA Damage Response |  1.530984 |  42 |  25 |
| p53-Independent G1/S DNA damage checkpoint |  1.530984 |  42 |  25 |
| AUF1 (hnRNP D0) binds and destabilizes mRNA |  1.554179 |  44 |  27 |
| Regulation of activated PAK-2p34 by proteasome mediated degradation |  1.559095 |  42 |  25 |
| Biological oxidations |  1.559514 |  85 |  43 |
| Somitogenesis |  1.577357 |  40 |  24 |
| Autodegradation of the E3 ubiquitin ligase COP1 |  1.579467 |  44 |  27 |
| Disorders of transmembrane transporters |  1.580149 | 114 |  22 |
| Regulation of Apoptosis |  1.587880 |  43 |  26 |
| Mitochondrial protein degradation |  1.594150 |  75 |  30 |
| Sulfur amino acid metabolism |  1.619725 |  22 |  10 |
| ABC transporter disorders |  1.625230 |  61 |  31 |
| Regulation of beta-cell development |  1.625786 |  19 |   4 |
| Respiratory electron transport |  1.626783 |  78 |  33 |
| Metabolism of nucleotides |  1.651396 |  69 |  18 |
| Asymmetric localization of PCP proteins |  1.651439 |  47 |  23 |
| Aquaporin-mediated transport |  1.658361 |  17 |   5 |
| Pyruvate metabolism |  1.667465 |  31 |  14 |
| RHO GTPases activate CIT |  1.667608 |  10 |   3 |
| Degradation of DVL |  1.667783 |  45 |  27 |
| Diseases of carbohydrate metabolism |  1.681234 |  22 |   8 |
| Drug ADME |  1.682511 |  36 |  10 |
| Degradation of AXIN |  1.686945 |  45 |  27 |
| Regulation of Glucokinase by Glucokinase Regulatory Protein |  1.693556 |  24 |   6 |
| Defective TPR may confer susceptibility towards thyroid papillary carcinoma (TPC) |  1.693556 |  24 |   6 |
| Nucleotide catabolism |  1.696699 |  22 |   8 |
| RHO GTPases activate PAKs |  1.697307 |  13 |   4 |
| Triglyceride catabolism |  1.699835 |  10 |   4 |
| Lysine catabolism |  1.713094 |  12 |   8 |
| Hedgehog ligand biogenesis |  1.729542 |  49 |  30 |
| Regulation of pyruvate dehydrogenase (PDH) complex |  1.729713 |  10 |   5 |
| Cristae formation |  1.742941 |  24 |  13 |
| Visual phototransduction |  1.748230 |  35 |  14 |
| Hh mutants abrogate ligand secretion |  1.749674 |  47 |  29 |
| Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins. |  1.751881 |  94 |  41 |
| Aerobic respiration and respiratory electron transport |  1.753255 | 147 |  58 |
| Hh mutants are degraded by ERAD |  1.758956 |  46 |  29 |
| Peroxisomal protein import |  1.759622 |  44 |  28 |
| Formation of ATP by chemiosmotic coupling |  1.788182 |  12 |   7 |
| Branched-chain amino acid catabolism |  1.846272 |  17 |   9 |
| Defective CFTR causes cystic fibrosis |  1.891005 |  51 |  32 |
| Degradation of cysteine and homocysteine |  1.894221 |  11 |   6 |
| Triglyceride metabolism |  1.898246 |  18 |   7 |
| Initial triggering of complement |  1.912913 |  12 |   3 |
| Complement cascade |  1.960510 |  28 |   5 |
| The canonical retinoid cycle in rods (twilight vision) |  2.335860 |  10 |   8 |



### IV-HD

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| tRNA modification in the nucleus and cytosol | -2.177561 |  21 |  14 |
| tRNA processing | -2.114216 |  64 |  33 |
| Activation of the pre-replicative complex | -2.092571 |  23 |  13 |
| Interleukin-4 and Interleukin-13 signaling | -2.024606 |  39 |  13 |
| DNA strand elongation | -1.991399 |  25 |  13 |
| rRNA modification in the nucleus and cytosol | -1.951206 |  54 |  23 |
| Activation of ATR in response to replication stress | -1.881421 |  26 |  13 |
| Synthesis of DNA | -1.771895 |  89 |  30 |
| Regulation of actin dynamics for phagocytic cup formation | -1.769697 |  41 |  11 |
| Fcgamma receptor (FCGR) dependent phagocytosis | -1.733032 |  57 |  15 |
| Amplification of signal from the kinetochores | -1.721040 |  54 |  22 |
| Amplification  of signal from unattached  kinetochores via a MAD2  inhibitory signal | -1.721040 |  54 |  22 |
| Recruitment of mitotic centrosome proteins and complexes | -1.718997 |  41 |  13 |
| Centrosome maturation | -1.718997 |  41 |  13 |
| Cell Cycle, Mitotic | -1.711076 | 306 | 100 |
| RHO GTPases Activate Formins | -1.688922 |  68 |  29 |
| DNA Replication | -1.674303 |  95 |  30 |
| S Phase | -1.672746 | 116 |  37 |
| Mitotic Prometaphase | -1.671646 | 107 |  38 |
| Response to elevated platelet cytosolic Ca2+ | -1.664125 |  70 |  17 |
| Cell Cycle | -1.606754 | 363 | 110 |
| Platelet degranulation | -1.604579 |  66 |  15 |
| G1/S Transition | -1.598515 |  94 |  30 |
| Mitotic G1 phase and G1/S transition | -1.573187 | 103 |  32 |
| Cell Cycle Checkpoints | -1.571586 | 156 |  52 |
| M Phase | -1.538726 | 216 |  67 |
| Sphingolipid metabolism |  1.648678 |  51 |  20 |
| NCAM signaling for neurite out-growth |  1.803092 |  21 |  10 |
| Regulation of beta-cell development |  1.899808 |  17 |   5 |
| FOXO-mediated transcription of oxidative stress, metabolic and neuronal genes |  1.952205 |  14 |   5 |
| FOXO-mediated transcription |  2.146193 |  33 |  14 |



### IV-LD

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| Antigen activates B Cell Receptor (BCR) leading to generation of second messengers | -1.963323 |  21 |  10 |
| Downregulation of SMAD2/3:SMAD4 transcriptional activity | -1.922826 |  21 |   6 |
| RHO GTPases Activate NADPH Oxidases | -1.921551 |  17 |   8 |
| Transcriptional activity of SMAD2/SMAD3:SMAD4 heterotrimer | -1.917286 |  37 |   9 |
| Fcgamma receptor (FCGR) dependent phagocytosis | -1.880702 |  64 |  28 |
| Amino acid transport across the plasma membrane | -1.879409 |  18 |   7 |
| Transport of inorganic cations/anions and amino acids/oligopeptides | -1.856819 |  45 |  15 |
| Signaling by ERBB2 in Cancer | -1.838310 |  14 |   6 |
| Signaling by ERBB2 KD Mutants | -1.838310 |  14 |   6 |
| FCGR3A-mediated IL10 synthesis | -1.836638 |  21 |   9 |
| Platelet activation, signaling and aggregation | -1.835522 | 145 |  52 |
| Cell surface interactions at the vascular wall | -1.833250 |  66 |  26 |
| Signaling by ERBB2 ECD mutants | -1.804882 |  13 |   6 |
| Regulation of insulin secretion | -1.802926 |  26 |   9 |
| FCERI mediated Ca+2 mobilization | -1.802784 |  19 |  10 |
| ROS and RNS production in phagocytes | -1.792430 |  19 |   7 |
| Signaling by ERBB2 TMD/JMD mutants | -1.765167 |  11 |   5 |
| VEGFA-VEGFR2 Pathway | -1.754364 |  70 |  18 |
| RHO GTPases Activate ROCKs | -1.741422 |  12 |   9 |
| PI3K/AKT Signaling in Cancer | -1.737552 |  46 |  26 |
| G alpha (q) signalling events | -1.729742 |  62 |  19 |
| Netrin-1 signaling | -1.724029 |  17 |   6 |
| Signaling by VEGF | -1.718570 |  74 |  18 |
| SMAD2/SMAD3:SMAD4 heterotrimer regulates transcription | -1.714377 |  25 |   7 |
| Semaphorin interactions | -1.709944 |  40 |  20 |
| G alpha (i) signalling events | -1.706816 |  65 |  19 |
| Platelet calcium homeostasis | -1.702921 |  11 |   5 |
| Effects of PIP2 hydrolysis | -1.700834 |  11 |   6 |
| Cholesterol biosynthesis | -1.695521 |  23 |  11 |
| EPHA-mediated growth cone collapse | -1.688049 |  14 |  11 |
| Interleukin-4 and Interleukin-13 signaling | -1.686119 |  43 |  25 |
| Hemostasis | -1.685229 | 305 | 119 |
| RAF-independent MAPK1/3 activation | -1.677323 |  14 |  10 |
| G-protein mediated events | -1.672555 |  22 |   6 |
| Sema4D induced cell migration and growth-cone collapse | -1.668424 |  13 |   9 |
| Signaling by FGFR1 | -1.667585 |  25 |  14 |
| Role of phospholipids in phagocytosis | -1.666463 |  16 |  11 |
| Nuclear Receptor transcription pathway | -1.663007 |  32 |  11 |
| Sema4D in semaphorin signaling | -1.658031 |  17 |  10 |
| Smooth Muscle Contraction | -1.646729 |  24 |  14 |
| RAC3 GTPase cycle | -1.644313 |  64 |  24 |
| Signaling by ERBB2 | -1.644090 |  29 |  12 |
| Signaling by CSF1 (M-CSF) in myeloid cells | -1.643930 |  22 |   8 |
| PLC beta mediated events | -1.640126 |  19 |   5 |
| GPCR downstream signalling | -1.636954 | 170 |  49 |
| Constitutive Signaling by Aberrant PI3K in Cancer | -1.634688 |  28 |  14 |
| Signaling by GPCR | -1.634416 | 184 |  53 |
| Muscle contraction | -1.633358 |  64 |  24 |
| Negative regulation of the PI3K/AKT network | -1.617022 |  49 |  25 |
| PI5P, PP2A and IER3 Regulate PI3K/AKT Signaling | -1.615472 |  46 |  23 |
| SLC-mediated transmembrane transport | -1.615031 | 104 |  29 |
| Signaling by ALK in cancer | -1.614893 |  66 |  26 |
| Signaling by ALK fusions and activated point mutants | -1.614893 |  66 |  26 |
| Activation of gene expression by SREBF (SREBP) | -1.612002 |  38 |  17 |
| Class A/1 (Rhodopsin-like receptors) | -1.609446 |  51 |  14 |
| Platelet degranulation | -1.603596 |  74 |  22 |
| HSP90 chaperone cycle for steroid hormone receptors (SHR) in the presence of ligand | -1.599137 |  31 |   5 |
| Negative regulation of MAPK pathway | -1.593495 |  27 |  14 |
| Parasite infection | -1.591010 |  46 |  18 |
| Leishmania phagocytosis | -1.591010 |  46 |  18 |
| FCGR3A-mediated phagocytosis | -1.591010 |  46 |  18 |
| Fc epsilon receptor (FCERI) signaling | -1.573421 |  93 |  24 |
| Anti-inflammatory response favouring Leishmania parasite infection | -1.571348 |  36 |  17 |
| Leishmania parasite growth and survival | -1.571348 |  36 |  17 |
| Response to elevated platelet cytosolic Ca2+ | -1.569407 |  78 |  31 |
| Ca2+ pathway | -1.567243 |  27 |  10 |
| Signaling by Interleukins | -1.561478 | 224 |  82 |
| Costimulation by the CD28 family | -1.547630 |  35 |  16 |
| RAC2 GTPase cycle | -1.542744 |  66 |  16 |
| Signaling by Receptor Tyrosine Kinases | -1.523483 | 296 |  95 |
| Transcriptional regulation by RUNX1 | -1.500291 | 122 |  41 |
| RHO GTPase Effectors | -1.491957 | 170 |  58 |
| RAC1 GTPase cycle | -1.474237 | 128 |  49 |
| Diseases of signal transduction by growth factor receptors and second messengers | -1.456128 | 289 | 105 |
| Cytokine Signaling in Immune system | -1.443442 | 352 | 125 |
| Adaptive Immune System | -1.396534 | 409 | 120 |
| RHO GTPase cycle | -1.353014 | 305 | 101 |
| Signaling by Rho GTPases | -1.336077 | 429 | 134 |
| Signaling by Rho GTPases, Miro GTPases and RHOBTB3 | -1.332690 | 439 | 139 |
| rRNA processing in the nucleus and cytosol |  1.276087 | 157 |  75 |
| Major pathway of rRNA processing in the nucleolus and cytosol |  1.422418 | 147 |  77 |
| Metabolism of amino acids and derivatives |  1.426549 | 264 | 122 |
| Metabolism of nucleotides |  1.469857 |  72 |  22 |
| Translation |  1.475773 | 237 | 117 |
| Cellular response to starvation |  1.561049 | 117 |  45 |
| Aerobic respiration and respiratory electron transport |  1.664542 | 151 |  77 |
| Regulation of expression of SLITs and ROBOs |  1.668917 | 132 |  50 |
| Selenoamino acid metabolism |  1.730473 |  93 |  56 |
| Retrograde neurotrophin signalling |  1.807523 |  11 |   6 |
| Response of EIF2AK4 (GCN2) to amino acid deficiency |  1.812629 |  80 |  54 |
| Respiratory electron transport |  1.818833 |  80 |  44 |
| The canonical retinoid cycle in rods (twilight vision) |  1.819757 |  10 |   8 |
| Formation of ATP by chemiosmotic coupling |  1.848615 |  12 |   9 |
| SARS-CoV-2 modulates host translation machinery |  1.861811 |  38 |  25 |
| SARS-CoV-1 modulates host translation machinery |  1.885818 |  29 |  21 |
| Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins. |  1.982287 |  96 |  54 |
| Activation of the mRNA upon binding of the cap-binding complex and eIFs, and subsequent binding to 43S |  2.004224 |  50 |  31 |
| Ribosomal scanning and start codon recognition |  2.013422 |  50 |  31 |
| Nonsense-Mediated Decay (NMD) |  2.013986 |  89 |  57 |
| Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC) |  2.013986 |  89 |  57 |
| Cristae formation |  2.016761 |  24 |  16 |
| Translation initiation complex formation |  2.024032 |  49 |  31 |
| Eukaryotic Translation Termination |  2.154457 |  73 |  51 |
| Formation of the ternary complex, and subsequently, the 43S complex |  2.157697 |  44 |  30 |
| SRP-dependent cotranslational protein targeting to membrane |  2.237501 |  88 |  58 |
| Eukaryotic Translation Initiation |  2.271027 |  97 |  63 |
| Cap-dependent Translation Initiation |  2.271027 |  97 |  63 |
| Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC) |  2.283325 |  74 |  52 |
| Selenocysteine synthesis |  2.335805 |  73 |  52 |
| L13a-mediated translational silencing of Ceruloplasmin expression |  2.347708 |  90 |  61 |
| GTP hydrolysis and joining of the 60S ribosomal subunit |  2.348593 |  91 |  61 |
| Peptide chain elongation |  2.375799 |  69 |  51 |
| Eukaryotic Translation Elongation |  2.390348 |  71 |  53 |
| Formation of a pool of free 40S subunits |  2.400952 |  82 |  58 |
| Viral mRNA Translation |  2.415355 |  69 |  51 |





# Head-kidney

## 4 WPC



### DNA vaccine

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| Peptide ligand-binding receptors | -1.898369 |  37 |  16 |
| Class A/1 (Rhodopsin-like receptors) | -1.868515 |  79 |  31 |
| GPCR ligand binding | -1.851556 | 101 |  37 |
| Cellular hexose transport | -1.717279 |  11 |   5 |
| Sensory processing of sound by inner hair cells of the cochlea | -1.714307 |  33 |  15 |
| Signaling by GPCR | -1.686628 | 258 |  85 |
| Sensory processing of sound | -1.674182 |  36 |  17 |
| RORA activates gene expression | -1.668492 |  16 |   5 |
| Constitutive Signaling by Aberrant PI3K in Cancer | -1.667587 |  35 |  14 |
| Nitric oxide stimulates guanylate cyclase | -1.665922 |  10 |   5 |
| Glutamate Neurotransmitter Release Cycle | -1.641438 |  15 |  10 |
| HS-GAG biosynthesis | -1.634413 |  18 |   8 |
| Long-term potentiation | -1.619079 |  10 |   5 |
| IRS-related events triggered by IGF1R | -1.601692 |  28 |  12 |
| GPCR downstream signalling | -1.596801 | 241 |  78 |
| Receptor-type tyrosine-protein phosphatases | -1.588085 |  12 |   5 |
| Sensory processing of sound by outer hair cells of the cochlea | -1.587560 |  25 |  14 |
| NCAM1 interactions | -1.586362 |  16 |   6 |
| G alpha (12/13) signalling events | -1.580873 |  53 |  34 |
| Metal ion SLC transporters | -1.576531 |  14 |   7 |
| NRAGE signals death through JNK | -1.572897 |  43 |  16 |
| Class B/2 (Secretin family receptors) | -1.551199 |  21 |   6 |
| Signaling by Type 1 Insulin-like Growth Factor 1 Receptor (IGF1R) | -1.550559 |  29 |  14 |
| IGF1R signaling cascade | -1.550559 |  29 |  14 |
| G alpha (i) signalling events | -1.548316 | 100 |  35 |
| Negative regulation of the PI3K/AKT network | -1.546372 |  62 |  22 |
| Potassium Channels | -1.545506 |  23 |   9 |
| Dopamine Neurotransmitter Release Cycle | -1.544575 |  14 |   6 |
| Metabolism of fat-soluble vitamins | -1.543244 |  25 |  12 |
| Activation of SMO | -1.540024 |  11 |   6 |
| Signaling by WNT in cancer | -1.515090 |  21 |  10 |
| Laminin interactions | -1.507777 |  18 |   5 |
| RHOB GTPase cycle | -1.502878 |  64 |  27 |
| NCAM signaling for neurite out-growth | -1.499037 |  33 |  10 |
| NGF-stimulated transcription | -1.496005 |  24 |  11 |
| Disassembly of the destruction complex and recruitment of AXIN to the membrane | -1.492778 |  20 |  10 |
| Transport of bile salts and organic acids, metal ions and amine compounds | -1.485597 |  37 |  17 |
| Defective B3GALT6 causes EDSP2 and SEMDJL1 | -1.481914 |  10 |   4 |
| CDC42 GTPase cycle | -1.480730 | 125 |  60 |
| SLC-mediated transmembrane transport | -1.470681 | 129 |  53 |
| Sensory Perception | -1.469184 |  80 |  30 |
| Ca2+ pathway | -1.459038 |  34 |  11 |
| Protein-protein interactions at synapses | -1.457523 |  43 |  14 |
| Muscle contraction | -1.443707 |  90 |  24 |
| Signaling by NTRKs | -1.421794 |  98 |  38 |
| Sphingolipid metabolism | -1.417192 |  70 |  11 |
| G alpha (q) signalling events | -1.414677 |  88 |  44 |
| RHOA GTPase cycle | -1.409192 | 127 |  48 |
| PI3K/AKT Signaling in Cancer | -1.406538 |  57 |  20 |
| Signaling by NTRK1 (TRKA) | -1.380517 |  84 |  34 |
| MyD88:MAL(TIRAP) cascade initiated on plasma membrane | -1.379232 |  78 |  34 |
| Toll Like Receptor TLR1:TLR2 Cascade | -1.379232 |  78 |  34 |
| Toll Like Receptor TLR6:TLR2 Cascade | -1.379232 |  78 |  34 |
| Toll Like Receptor 2 (TLR2) Cascade | -1.379232 |  78 |  34 |
| Neuronal System | -1.361198 | 164 |  57 |
| Transmission across Chemical Synapses | -1.353636 | 117 |  40 |
| Axon guidance |  1.138802 | 361 |  94 |
| Adaptive Immune System |  1.201348 | 488 | 108 |
| Neutrophil degranulation |  1.208304 | 296 | 116 |
| SARS-CoV Infections |  1.311292 | 305 |  94 |
| RHO GTPases Activate Formins |  1.336346 |  97 |  45 |
| SARS-CoV-2 Infection |  1.346521 | 214 |  68 |
| UCH proteinases |  1.378611 |  75 |  24 |
| SARS-CoV-1 Infection |  1.391961 | 111 |  44 |
| CLEC7A (Dectin-1) signaling |  1.407009 |  81 |  20 |
| Downstream TCR signaling |  1.412984 |  70 |  22 |
| KEAP1-NFE2L2 pathway |  1.416750 | 102 |  23 |
| RNA Polymerase II Pre-transcription Events |  1.441413 |  71 |  20 |
| Global Genome Nucleotide Excision Repair (GG-NER) |  1.450598 |  70 |  18 |
| ER-Phagosome pathway |  1.450672 |  63 |  23 |
| Transport of Mature mRNAs Derived from Intronless Transcripts |  1.454109 |  35 |   9 |
| Mitotic Prometaphase |  1.456858 | 153 |  61 |
| Protein localization |  1.459801 | 121 |  34 |
| Neddylation |  1.460006 | 193 |  46 |
| SUMOylation of chromatin organization proteins |  1.470953 |  47 |  12 |
| Formation of paraxial mesoderm |  1.473598 |  54 |  17 |
| PTEN Regulation |  1.476454 | 117 |  28 |
| RHO GTPase Effectors |  1.481098 | 198 |  76 |
| Antigen processing-Cross presentation |  1.481510 |  75 |  29 |
| C-type lectin receptors (CLRs) |  1.483194 |  95 |  23 |
| EML4 and NUDC in mitotic spindle formation |  1.487638 |  84 |  42 |
| Transport of Mature mRNA Derived from an Intronless Transcript |  1.488959 |  34 |   9 |
| Nucleotide Excision Repair |  1.490695 |  91 |  23 |
| ISG15 antiviral mechanism |  1.498335 |  55 |  12 |
| Export of Viral Ribonucleoproteins from Nucleus |  1.502679 |  27 |   8 |
| NEP/NS2 Interacts with the Cellular Export Machinery |  1.502679 |  27 |   8 |
| Translesion synthesis by Y family DNA polymerases bypasses lesions on DNA template |  1.503871 |  32 |   8 |
| Transport of the SLBP independent Mature mRNA |  1.506275 |  28 |   8 |
| Glycolysis |  1.512998 |  51 |  15 |
| HCMV Early Events |  1.516398 |  54 |  15 |
| Mitotic Spindle Checkpoint |  1.521303 |  90 |  40 |
| Amplification of signal from the kinetochores |  1.522766 |  76 |  36 |
| Amplification  of signal from unattached  kinetochores via a MAD2  inhibitory signal |  1.522766 |  76 |  36 |
| ATF4 activates genes in response to endoplasmic reticulum  stress |  1.529761 |  22 |   9 |
| Nuclear Envelope (NE) Reassembly |  1.534652 |  51 |  24 |
| TCR signaling |  1.538692 |  87 |  28 |
| SARS-CoV-2-host interactions |  1.540857 | 139 |  58 |
| Cellular response to chemical stress |  1.542184 | 159 |  54 |
| mRNA Splicing |  1.543553 | 180 |  73 |
| BBSome-mediated cargo-targeting to cilium |  1.546516 |  15 |   7 |
| Mitotic G2-G2/M phases |  1.546591 | 150 |  48 |
| Signaling by the B Cell Receptor (BCR) |  1.547450 |  93 |  25 |
| Transcription of the HIV genome |  1.547464 |  61 |  19 |
| Transcriptional regulation by RUNX1 |  1.549167 | 142 |  32 |
| G2/M Transition |  1.552014 | 149 |  48 |
| Vpr-mediated nuclear import of PICs |  1.552132 |  26 |   8 |
| Somitogenesis |  1.554455 |  46 |  17 |
| Deposition of new CENPA-containing nucleosomes at the centromere |  1.554581 |  21 |   6 |
| Nucleosome assembly |  1.554581 |  21 |   6 |
| Regulation of Glucokinase by Glucokinase Regulatory Protein |  1.563185 |  25 |   8 |
| Defective TPR may confer susceptibility towards thyroid papillary carcinoma (TPC) |  1.563185 |  25 |   8 |
| Metabolism of amino acids and derivatives |  1.563329 | 267 | 112 |
| mRNA Splicing - Minor Pathway |  1.563560 |  44 |  12 |
| Gene and protein expression by JAK-STAT signaling after Interleukin-12 stimulation |  1.563972 |  27 |  12 |
| Meiotic synapsis |  1.564023 |  23 |   5 |
| SUMOylation of RNA binding proteins |  1.565725 |  38 |  11 |
| mRNA 3'-end processing |  1.570448 |  51 |  27 |
| HIV Transcription Initiation |  1.570483 |  40 |  11 |
| RNA Polymerase II HIV Promoter Escape |  1.570483 |  40 |  11 |
| RNA Polymerase II Promoter Escape |  1.570483 |  40 |  11 |
| RNA Polymerase II Transcription Pre-Initiation And Promoter Opening |  1.570483 |  40 |  11 |
| RNA Polymerase II Transcription Initiation |  1.570483 |  40 |  11 |
| RNA Polymerase II Transcription Initiation And Promoter Clearance |  1.570483 |  40 |  11 |
| Signaling by FGFR2 IIIa TM |  1.576421 |  16 |   6 |
| TRAF6 mediated IRF7 activation in TLR7/8 or 9 signaling |  1.580540 |  10 |   8 |
| mRNA Splicing - Major Pathway |  1.581028 | 173 |  71 |
| Termination of translesion DNA synthesis |  1.582349 |  25 |   7 |
| Transport of Mature mRNA derived from an Intron-Containing Transcript |  1.587545 |  61 |  29 |
| mRNA decay by 3' to 5' exoribonuclease |  1.596144 |  14 |   7 |
| Telomere Maintenance |  1.596961 |  62 |  16 |
| APC/C:Cdc20 mediated degradation of Cyclin B |  1.602951 |  20 |   6 |
| Resolution of Sister Chromatid Cohesion |  1.603552 |  94 |  45 |
| FCERI mediated NF-kB activation |  1.610498 |  65 |  21 |
| Mitochondrial protein import |  1.611419 |  50 |  24 |
| Dual Incision in GG-NER |  1.612127 |  33 |   8 |
| Removal of the Flap Intermediate from the C-strand |  1.617365 |  14 |  10 |
| SUMOylation of DNA replication proteins |  1.617774 |  39 |  13 |
| Downstream signaling events of B Cell Receptor (BCR) |  1.621439 |  68 |  23 |
| Transport of Mature Transcript to Cytoplasm |  1.626752 |  69 |  33 |
| MHC class II antigen presentation |  1.629213 |  78 |  19 |
| Processing of Capped Intron-Containing Pre-mRNA |  1.631444 | 239 |  99 |
| NS1 Mediated Effects on Host Pathways |  1.632671 |  32 |  10 |
| MAPK6/MAPK4 signaling |  1.639220 |  76 |  22 |
| SUMOylation of SUMOylation proteins |  1.640954 |  29 |   9 |
| Gap-filling DNA repair synthesis and ligation in TC-NER |  1.641227 |  53 |  13 |
| Hedgehog 'off' state |  1.643923 |  81 |  27 |
| HIV elongation arrest and recovery |  1.644714 |  30 |   9 |
| Pausing and recovery of HIV elongation |  1.644714 |  30 |   9 |
| Chromosome Maintenance |  1.645756 |  80 |  28 |
| PERK regulates gene expression |  1.652403 |  27 |  11 |
| Pausing and recovery of Tat-mediated HIV elongation |  1.653135 |  29 |   9 |
| Tat-mediated HIV elongation arrest and recovery |  1.653135 |  29 |   9 |
| Interconversion of nucleotide di- and triphosphates |  1.655394 |  19 |   8 |
| Regulation of PTEN stability and activity |  1.660899 |  59 |  19 |
| tRNA processing in the nucleus |  1.662955 |  48 |  13 |
| DNA Damage Bypass |  1.663316 |  41 |  12 |
| Base Excision Repair |  1.671815 |  37 |  22 |
| Cell Cycle, Mitotic |  1.674286 | 407 | 163 |
| Formation of TC-NER Pre-Incision Complex |  1.675714 |  45 |  11 |
| Dual incision in TC-NER |  1.678389 |  53 |  13 |
| SARS-CoV-1-host interactions |  1.679226 |  74 |  34 |
| Interactions of Vpr with host cellular proteins |  1.679614 |  27 |   9 |
| Formation of tubulin folding intermediates by CCT/TriC |  1.680516 |  12 |  10 |
| Late Phase of HIV Life Cycle |  1.683110 | 115 |  37 |
| Prefoldin mediated transfer of substrate  to CCT/TriC |  1.688304 |  18 |  11 |
| Cooperation of Prefoldin and TriC/CCT  in actin and tubulin folding |  1.688304 |  18 |  11 |
| rRNA processing |  1.688988 | 172 |  58 |
| Condensation of Prometaphase Chromosomes |  1.691717 |  11 |   5 |
| Translesion synthesis by POLI |  1.693945 |  15 |   4 |
| M Phase |  1.695185 | 287 | 119 |
| RNA Polymerase II Transcription Termination |  1.703843 |  56 |  30 |
| Translesion Synthesis by POLH |  1.709257 |  17 |   5 |
| Disorders of transmembrane transporters |  1.709713 | 116 |  32 |
| SUMOylation of ubiquitinylation proteins |  1.720380 |  31 |  10 |
| SUMOylation of DNA damage response and repair proteins |  1.721382 |  61 |  16 |
| Regulation of mRNA stability by proteins that bind AU-rich elements |  1.721465 |  75 |  26 |
| Collagen degradation |  1.722414 |  31 |  10 |
| snRNP Assembly |  1.726941 |  44 |  13 |
| Metabolism of non-coding RNA |  1.726941 |  44 |  13 |
| Gap-filling DNA repair synthesis and ligation in GG-NER |  1.728069 |  20 |   8 |
| Unwinding of DNA |  1.730626 |  11 |   3 |
| Nuclear events mediated by NFE2L2 |  1.732387 |  76 |  31 |
| ABC-family proteins mediated transport |  1.735827 |  76 |  25 |
| Regulation of RUNX2 expression and activity |  1.739000 |  56 |  19 |
| Cell Cycle Checkpoints |  1.741388 | 215 |  98 |
| Translesion synthesis by REV1 |  1.741781 |  14 |   4 |
| Condensation of Prophase Chromosomes |  1.742712 |  16 |  10 |
| Maternal to zygotic transition (MZT) |  1.746421 |  42 |   3 |
| Activation of NF-kappaB in B cells |  1.752844 |  57 |  20 |
| TNFR2 non-canonical NF-kB pathway |  1.755333 |  60 |  21 |
| Lagging Strand Synthesis |  1.757589 |  18 |   8 |
| Metabolic disorders of biological oxidation enzymes |  1.763286 |  17 |   6 |
| HIV Life Cycle |  1.764963 | 123 |  41 |
| Abortive elongation of HIV-1 transcript in the absence of Tat |  1.767302 |  22 |  12 |
| Transcription-Coupled Nucleotide Excision Repair (TC-NER) |  1.769134 |  65 |  19 |
| Resolution of AP sites via the multiple-nucleotide patch replacement pathway |  1.772279 |  21 |  11 |
| Glutathione conjugation |  1.792277 |  17 |  10 |
| G1/S DNA Damage Checkpoints |  1.793583 |  58 |  21 |
| ATF6 (ATF6-alpha) activates chaperone genes |  1.802080 |  10 |   7 |
| Hh mutants abrogate ligand secretion |  1.803321 |  50 |  20 |
| p53-Dependent G1 DNA Damage Response |  1.806012 |  57 |  21 |
| p53-Dependent G1/S DNA damage checkpoint |  1.806012 |  57 |  21 |
| Signaling by NOTCH4 |  1.807572 |  66 |  23 |
| rRNA processing in the nucleus and cytosol |  1.808687 | 162 |  57 |
| Cyclin A:Cdk2-associated events at S phase entry |  1.813637 |  76 |  25 |
| Inhibition of DNA recombination at telomere |  1.818177 |  21 |   7 |
| Mitochondrial translation termination |  1.832632 |  76 |  43 |
| Cyclin E associated events during G1/S transition |  1.835220 |  74 |  25 |
| G2/M Checkpoints |  1.836895 | 115 |  55 |
| Maturation of TCA enzymes and regulation of TCA cycle |  1.837153 |  17 |  12 |
| Transcriptional regulation by small RNAs |  1.842060 |  45 |  14 |
| Mitochondrial translation initiation |  1.842166 |  76 |  44 |
| Mitotic Anaphase |  1.844995 | 184 |  80 |
| Mitotic Metaphase and Anaphase |  1.849734 | 185 |  81 |
| Hedgehog ligand biogenesis |  1.852418 |  52 |  21 |
| PCNA-Dependent Long Patch Base Excision Repair |  1.871621 |  17 |   8 |
| Host Interactions of HIV factors |  1.874978 | 107 |  36 |
| Complex I biogenesis |  1.875018 |  44 |  22 |
| Mitochondrial translation |  1.876080 |  82 |  46 |
| Mitochondrial translation elongation |  1.887312 |  76 |  45 |
| Dectin-1 mediated noncanonical NF-kB signaling |  1.895239 |  52 |  19 |
| HIV Infection |  1.896702 | 191 |  64 |
| Heme biosynthesis |  1.899694 |  11 |   7 |
| RUNX1 regulates transcription of genes involved in differentiation of HSCs |  1.901211 |  65 |  22 |
| Viral Messenger RNA Synthesis |  1.905425 |  38 |  14 |
| ATF6 (ATF6-alpha) activates chaperones |  1.907050 |  12 |   8 |
| Major pathway of rRNA processing in the nucleolus and cytosol |  1.907547 | 152 |  55 |
| AUF1 (hnRNP D0) binds and destabilizes mRNA |  1.912492 |  47 |  18 |
| FBXL7 down-regulates AURKA during mitotic entry and in early mitosis |  1.921740 |  48 |  19 |
| FGFR2 alternative splicing |  1.923406 |  21 |  10 |
| Cellular response to starvation |  1.923815 | 125 |  41 |
| Degradation of GLI1 by the proteasome |  1.925014 |  50 |  20 |
| PCP/CE pathway |  1.925900 |  70 |  37 |
| APC/C-mediated degradation of cell cycle proteins |  1.927213 |  77 |  28 |
| Regulation of mitotic cell cycle |  1.927213 |  77 |  28 |
| Metabolism of polyamines |  1.938700 |  53 |  27 |
| Activation of the pre-replicative complex |  1.941105 |  28 |  18 |
| Cellular response to hypoxia |  1.941721 |  65 |  23 |
| NIK-->noncanonical NF-kB signaling |  1.945799 |  51 |  19 |
| Activation of Matrix Metalloproteinases |  1.954516 |  11 |   4 |
| GLI3 is processed to GLI3R by the proteasome |  1.968317 |  51 |  20 |
| Activation of ATR in response to replication stress |  1.970960 |  33 |  21 |
| Switching of origins to a post-replicative state |  1.978475 |  81 |  41 |
| Citric acid cycle (TCA cycle) |  1.981681 |  30 |  18 |
| SCF-beta-TrCP mediated degradation of Emi1 |  1.984259 |  48 |  19 |
| Asymmetric localization of PCP proteins |  1.984539 |  52 |  28 |
| Mitotic G1 phase and G1/S transition |  1.988280 | 125 |  62 |
| Separation of Sister Chromatids |  1.988696 | 148 |  74 |
| DNA strand elongation |  1.994820 |  29 |  11 |
| Degradation of GLI2 by the proteasome |  1.997713 |  50 |  20 |
| Stabilization of p53 |  2.004577 |  48 |  18 |
| ABC transporter disorders |  2.006541 |  61 |  24 |
| Cristae formation |  2.019671 |  24 |  17 |
| Regulation of APC/C activators between G1/S and early anaphase |  2.029764 |  71 |  27 |
| Nonsense-Mediated Decay (NMD) |  2.031480 |  92 |  48 |
| Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC) |  2.031480 |  92 |  48 |
| Vpu mediated degradation of CD4 |  2.031771 |  45 |  19 |
| RHO GTPases activate PKNs |  2.032475 |  30 |  18 |
| SCF(Skp2)-mediated degradation of p27/p21 |  2.045367 |  54 |  22 |
| Orc1 removal from chromatin |  2.046792 |  64 |  36 |
| Negative regulation of NOTCH4 signaling |  2.048350 |  47 |  20 |
| Recognition of DNA damage by PCNA-containing replication complex |  2.062587 |  25 |  11 |
| GSK3B and BTRC:CUL1-mediated-degradation of NFE2L2 |  2.062993 |  46 |  19 |
| Degradation of DVL |  2.076107 |  50 |  29 |
| G1/S Transition |  2.078545 | 113 |  60 |
| S Phase |  2.082303 | 144 |  46 |
| Respiratory electron transport |  2.088733 |  83 |  41 |
| APC:Cdc20 mediated degradation of cell cycle proteins prior to satisfation of the cell cycle checkpoint |  2.090415 |  65 |  24 |
| Autodegradation of the E3 ubiquitin ligase COP1 |  2.097057 |  46 |  18 |
| Aerobic respiration and respiratory electron transport |  2.113940 | 161 |  74 |
| Formation of ATP by chemiosmotic coupling |  2.114518 |  12 |   7 |
| APC/C:Cdc20 mediated degradation of mitotic proteins |  2.121036 |  66 |  25 |
| Selenoamino acid metabolism |  2.125163 |  95 |  52 |
| Regulation of ornithine decarboxylase (ODC) |  2.130723 |  46 |  25 |
| Translation |  2.132676 | 249 | 124 |
| Oxygen-dependent proline hydroxylation of Hypoxia-inducible Factor Alpha |  2.149238 |  58 |  21 |
| Activation of APC/C and APC/C:Cdc20 mediated degradation of mitotic proteins |  2.152249 |  67 |  26 |
| Degradation of AXIN |  2.168457 |  48 |  19 |
| Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins. |  2.169438 |  98 |  48 |
| Cross-presentation of soluble exogenous antigens (endosomes) |  2.169850 |  42 |  25 |
| SARS-CoV-1 modulates host translation machinery |  2.171885 |  30 |  17 |
| Cdc20:Phospho-APC/C mediated degradation of Cyclin A |  2.179124 |  64 |  33 |
| CDK-mediated phosphorylation and removal of Cdc6 |  2.189034 |  64 |  33 |
| Ubiquitin Mediated Degradation of Phosphorylated Cdc25A |  2.195298 |  45 |  18 |
| p53-Independent DNA Damage Response |  2.195298 |  45 |  18 |
| p53-Independent G1/S DNA damage checkpoint |  2.195298 |  45 |  18 |
| The role of GTSE1 in G2/M progression after G2 checkpoint |  2.200940 |  58 |  25 |
| APC/C:Cdc20 mediated degradation of Securin |  2.201093 |  58 |  22 |
| APC/C:Cdh1 mediated degradation of Cdc20 and other APC/C:Cdh1 targeted proteins in late mitosis/early G1 |  2.205780 |  64 |  24 |
| Assembly of the pre-replicative complex |  2.206283 |  76 |  27 |
| Regulation of activated PAK-2p34 by proteasome mediated degradation |  2.211546 |  44 |  18 |
| Selenocysteine synthesis |  2.220933 |  75 |  52 |
| Ubiquitin-dependent degradation of Cyclin D |  2.222888 |  46 |  19 |
| Metabolism of porphyrins |  2.223918 |  17 |  11 |
| Regulation of Apoptosis |  2.232589 |  47 |  27 |
| Autodegradation of Cdh1 by Cdh1:APC/C |  2.244973 |  56 |  21 |
| Eukaryotic Translation Termination |  2.252142 |  74 |  53 |
| Regulation of RUNX3 expression and activity |  2.253519 |  47 |  18 |
| Synthesis of DNA |  2.253831 | 106 |  40 |
| Influenza Infection |  2.259532 | 126 |  70 |
| Vif-mediated degradation of APOBEC3G |  2.266014 |  47 |  20 |
| Detoxification of Reactive Oxygen Species |  2.266782 |  22 |  18 |
| DNA Replication |  2.269363 | 115 |  62 |
| SARS-CoV-2 modulates host translation machinery |  2.287955 |  39 |  21 |
| Influenza Viral RNA Transcription and Replication |  2.298649 | 111 |  61 |
| Hh mutants are degraded by ERAD |  2.304187 |  48 |  28 |
| Peptide chain elongation |  2.307492 |  71 |  52 |
| DNA Replication Pre-Initiation |  2.321567 |  90 |  49 |
| Activation of the mRNA upon binding of the cap-binding complex and eIFs, and subsequent binding to 43S |  2.323919 |  51 |  24 |
| Eukaryotic Translation Elongation |  2.335797 |  73 |  53 |
| Translation initiation complex formation |  2.343264 |  50 |  24 |
| Response of EIF2AK4 (GCN2) to amino acid deficiency |  2.344863 |  82 |  46 |
| Defective CFTR causes cystic fibrosis |  2.364787 |  52 |  31 |
| Ribosomal scanning and start codon recognition |  2.389005 |  51 |  25 |
| Formation of the ternary complex, and subsequently, the 43S complex |  2.391335 |  45 |  24 |
| Signaling by ROBO receptors |  2.439710 | 166 |  95 |
| Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC) |  2.440381 |  75 |  53 |
| Viral mRNA Translation |  2.456060 |  70 |  52 |
| SRP-dependent cotranslational protein targeting to membrane |  2.459457 |  89 |  52 |
| Eukaryotic Translation Initiation |  2.477005 |  98 |  66 |
| Cap-dependent Translation Initiation |  2.477005 |  98 |  66 |
| L13a-mediated translational silencing of Ceruloplasmin expression |  2.513564 |  91 |  55 |
| Formation of a pool of free 40S subunits |  2.519407 |  83 |  61 |
| GTP hydrolysis and joining of the 60S ribosomal subunit |  2.535719 |  92 |  56 |
| Regulation of expression of SLITs and ROBOs |  2.575939 | 136 |  83 |
| Mitochondrial protein degradation |  2.577030 |  81 |  50 |



### EOMES

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| Sensory processing of sound by inner hair cells of the cochlea | -1.748837 |  31 |  17 |
| Long-term potentiation | -1.688954 |  10 |   5 |
| Dopamine Neurotransmitter Release Cycle | -1.664041 |  13 |   7 |
| Sensory processing of sound | -1.663711 |  34 |  16 |
| Protein-protein interactions at synapses | -1.581498 |  42 |  15 |
| Transmission across Chemical Synapses | -1.479524 | 115 |  34 |
| Metabolism of carbohydrates |  1.578805 | 189 |  49 |
| Response to elevated platelet cytosolic Ca2+ |  1.581291 |  78 |  21 |
| Metabolism of steroids |  1.584888 |  97 |  24 |
| rRNA processing |  1.593360 | 172 |  64 |
| Diseases of metabolism |  1.595154 | 144 |  35 |
| Platelet degranulation |  1.596438 |  74 |  20 |
| Extracellular matrix organization |  1.600723 | 155 |  29 |
| Glucose metabolism |  1.634594 |  54 |  14 |
| Disorders of transmembrane transporters |  1.636374 | 114 |  33 |
| rRNA processing in the nucleus and cytosol |  1.664360 | 162 |  75 |
| Translation |  1.667454 | 249 | 105 |
| Maternal to zygotic transition (MZT) |  1.700397 |  40 |   4 |
| Pyruvate metabolism |  1.705549 |  37 |  14 |
| Cellular response to starvation |  1.706554 | 124 |  57 |
| Degradation of the extracellular matrix |  1.718058 |  66 |  14 |
| Glycogen metabolism |  1.732991 |  20 |   8 |
| Plasma lipoprotein assembly, remodeling, and clearance |  1.733256 |  46 |  15 |
| Collagen biosynthesis and modifying enzymes |  1.735078 |  30 |  15 |
| Metabolism of nucleotides |  1.742821 |  74 |  21 |
| Nucleotide catabolism |  1.743800 |  24 |   7 |
| Formation of ATP by chemiosmotic coupling |  1.747353 |  12 |  11 |
| Amyloid fiber formation |  1.750055 |  34 |  13 |
| Lysine catabolism |  1.756481 |  11 |   7 |
| Peroxisomal lipid metabolism |  1.757965 |  19 |  12 |
| Metabolism of vitamins and cofactors |  1.758470 | 131 |  35 |
| Activation of Matrix Metalloproteinases |  1.758975 |  11 |   5 |
| Arachidonic acid metabolism |  1.759124 |  30 |  14 |
| Phase II - Conjugation of compounds |  1.784450 |  48 |  16 |
| Assembly of active LPL and LIPC lipase complexes |  1.800929 |  11 |   5 |
| Major pathway of rRNA processing in the nucleolus and cytosol |  1.801319 | 152 |  73 |
| Formation of Fibrin Clot (Clotting Cascade) |  1.802891 |  13 |   3 |
| Fatty acid metabolism |  1.807310 | 118 |  43 |
| Integrin cell surface interactions |  1.809870 |  47 |  11 |
| Metabolism of folate and pterines |  1.831727 |  11 |   6 |
| Degradation of cysteine and homocysteine |  1.840771 |  13 |   8 |
| SARS-CoV-2 modulates host translation machinery |  1.845231 |  39 |  26 |
| SARS-CoV-1 modulates host translation machinery |  1.851809 |  30 |  21 |
| Attachment and Entry |  1.864010 |  11 |   7 |
| Protein localization |  1.884595 | 121 |  39 |
| Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell |  1.888422 |  27 |  10 |
| Respiratory electron transport |  1.896364 |  83 |  41 |
| Metabolism of porphyrins |  1.897188 |  17 |  11 |
| Sulfur amino acid metabolism |  1.897366 |  24 |  12 |
| Binding and Uptake of Ligands by Scavenger Receptors |  1.907952 |  23 |   7 |
| Metabolic disorders of biological oxidation enzymes |  1.925527 |  17 |   7 |
| Nonsense-Mediated Decay (NMD) |  1.928087 |  92 |  56 |
| Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC) |  1.928087 |  92 |  56 |
| Heme biosynthesis |  1.951164 |  11 |   8 |
| Plasma lipoprotein remodeling |  1.955769 |  19 |   9 |
| Translation initiation complex formation |  1.964003 |  50 |  26 |
| Post-translational protein phosphorylation |  1.966010 |  50 |  18 |
| Scavenging by Class A Receptors |  1.980234 |  14 |   5 |
| Activation of the mRNA upon binding of the cap-binding complex and eIFs, and subsequent binding to 43S |  1.992614 |  51 |  27 |
| Cytochrome P450 - arranged by substrate type |  1.993397 |  22 |  10 |
| Peroxisomal protein import |  1.994364 |  44 |  18 |
| Detoxification of Reactive Oxygen Species |  2.005465 |  22 |  13 |
| Ribosomal scanning and start codon recognition |  2.007552 |  51 |  27 |
| Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins. |  2.015502 |  98 |  50 |
| Formation of the ternary complex, and subsequently, the 43S complex |  2.025754 |  45 |  26 |
| Metabolism of fat-soluble vitamins |  2.035588 |  25 |  12 |
| Influenza Infection |  2.037178 | 125 |  67 |
| Influenza Viral RNA Transcription and Replication |  2.066656 | 110 |  64 |
| Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs) |  2.092458 |  56 |  21 |
| Aerobic respiration and respiratory electron transport |  2.110440 | 160 |  78 |
| Retinoid metabolism and transport |  2.147202 |  21 |   9 |
| Glutathione conjugation |  2.149258 |  17 |   7 |
| Collagen degradation |  2.163898 |  30 |  11 |
| Signaling by ROBO receptors |  2.178178 | 165 |  91 |
| Eukaryotic Translation Termination |  2.188410 |  74 |  51 |
| Selenoamino acid metabolism |  2.193095 |  95 |  61 |
| Regulation of expression of SLITs and ROBOs |  2.201981 | 135 |  78 |
| Selenocysteine synthesis |  2.230963 |  75 |  50 |
| Visual phototransduction |  2.236597 |  36 |  15 |
| Peptide chain elongation |  2.250277 |  71 |  50 |
| SRP-dependent cotranslational protein targeting to membrane |  2.268162 |  89 |  60 |
| Eukaryotic Translation Elongation |  2.279521 |  73 |  51 |
| Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC) |  2.287740 |  75 |  51 |
| Response of EIF2AK4 (GCN2) to amino acid deficiency |  2.295444 |  82 |  53 |
| Mitochondrial protein degradation |  2.316663 |  81 |  46 |
| Eukaryotic Translation Initiation |  2.333654 |  98 |  59 |
| Cap-dependent Translation Initiation |  2.333654 |  98 |  59 |
| Viral mRNA Translation |  2.342171 |  70 |  50 |
| L13a-mediated translational silencing of Ceruloplasmin expression |  2.359260 |  91 |  56 |
| GTP hydrolysis and joining of the 60S ribosomal subunit |  2.378197 |  92 |  57 |
| Metabolism of amino acids and derivatives |  2.381380 | 265 | 125 |
| Formation of a pool of free 40S subunits |  2.401536 |  83 |  56 |
| Phase I - Functionalization of compounds |  2.407338 |  41 |  18 |
| Biological oxidations |  2.417072 |  94 |  36 |



### GATA3

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| Synthesis of very long-chain fatty acyl-CoAs | -1.684479 |  19 |   9 |
| Sensory processing of sound | -1.668433 |  32 |  14 |
| Sensory processing of sound by inner hair cells of the cochlea | -1.629033 |  29 |  15 |
| HS-GAG biosynthesis | -1.626025 |  17 |   8 |
| Sensory processing of sound by outer hair cells of the cochlea | -1.608837 |  21 |   9 |
| Metabolism of fat-soluble vitamins | -1.606457 |  25 |  16 |
| Nicotinamide salvaging | -1.595758 |  13 |   7 |
| Plasma lipoprotein remodeling | -1.594360 |  19 |  11 |
| Bile acid and bile salt metabolism | -1.591882 |  22 |   9 |
| Synthesis of PIPs at the late endosome membrane | -1.553295 |  11 |   6 |
| RHOB GTPase cycle | -1.541174 |  63 |  21 |
| NRAGE signals death through JNK | -1.535744 |  42 |  19 |
| Peptide ligand-binding receptors | -1.517745 |  34 |  12 |
| Class A/1 (Rhodopsin-like receptors) | -1.513317 |  74 |  31 |
| G alpha (12/13) signalling events | -1.508970 |  52 |  33 |
| Muscle contraction | -1.501913 |  88 |  21 |
| Signaling by GPCR | -1.499479 | 245 | 105 |
| GPCR ligand binding | -1.484694 |  91 |  33 |
| GPCR downstream signalling | -1.476687 | 230 |  98 |
| Metabolism of vitamins and cofactors | -1.474144 | 131 |  48 |
| SLC-mediated transmembrane transport | -1.441127 | 127 |  34 |
| G alpha (q) signalling events | -1.429814 |  86 |  40 |
| Metabolism of lipids | -1.257625 | 500 | 139 |
| Neutrophil degranulation |  1.283297 | 293 |  96 |
| Nervous system development |  1.287317 | 363 | 112 |
| SARS-CoV Infections |  1.300935 | 303 |  82 |
| Axon guidance |  1.317153 | 350 | 110 |
| SARS-CoV-2-host interactions |  1.342103 | 138 |  40 |
| Signaling by NOTCH |  1.349235 | 139 |  43 |
| Adaptive Immune System |  1.351711 | 481 | 135 |
| Mitotic Prometaphase |  1.363991 | 151 |  46 |
| Cellular response to chemical stress |  1.377313 | 159 |  52 |
| G2/M Checkpoints |  1.403760 | 115 |  35 |
| mRNA Splicing |  1.416873 | 179 |  64 |
| CLEC7A (Dectin-1) signaling |  1.428114 |  80 |  26 |
| Cell Cycle Checkpoints |  1.435963 | 215 |  70 |
| Mitotic G2-G2/M phases |  1.443498 | 148 |  54 |
| G2/M Transition |  1.450079 | 147 |  54 |
| Amplification of signal from the kinetochores |  1.461076 |  76 |  27 |
| Amplification  of signal from unattached  kinetochores via a MAD2  inhibitory signal |  1.461076 |  76 |  27 |
| Cell Cycle |  1.462536 | 497 | 149 |
| Metabolism of amino acids and derivatives |  1.476012 | 264 | 106 |
| Switching of origins to a post-replicative state |  1.480159 |  81 |  31 |
| ABC-family proteins mediated transport |  1.482829 |  74 |  27 |
| rRNA processing |  1.489591 | 172 |  63 |
| mRNA Splicing - Major Pathway |  1.505290 | 172 |  64 |
| Regulation of TNFR1 signaling |  1.508306 |  37 |  11 |
| FCERI mediated NF-kB activation |  1.511784 |  65 |  24 |
| UCH proteinases |  1.515172 |  75 |  28 |
| Processing of Capped Intron-Containing Pre-mRNA |  1.522370 | 236 |  84 |
| Hedgehog 'off' state |  1.526536 |  80 |  30 |
| HCMV Late Events |  1.533330 |  45 |  12 |
| SUMOylation of DNA replication proteins |  1.535361 |  38 |  15 |
| Transport of the SLBP independent Mature mRNA |  1.536768 |  27 |   7 |
| Citric acid cycle (TCA cycle) |  1.545413 |  30 |  15 |
| Postmitotic nuclear pore complex (NPC) reformation |  1.545697 |  23 |   9 |
| Chromosome Maintenance |  1.546066 |  79 |  33 |
| SARS-CoV-1-host interactions |  1.547713 |  74 |  32 |
| Regulation of PTEN stability and activity |  1.551231 |  59 |  22 |
| Late Phase of HIV Life Cycle |  1.551941 | 114 |  46 |
| Nuclear Envelope (NE) Reassembly |  1.552482 |  50 |  20 |
| Transport of Ribonucleoproteins into the Host Nucleus |  1.555139 |  25 |   7 |
| rRNA processing in the nucleus and cytosol |  1.557445 | 162 |  58 |
| MAPK6/MAPK4 signaling |  1.560974 |  75 |  27 |
| Base Excision Repair |  1.561452 |  37 |  13 |
| Transport of Mature mRNA derived from an Intron-Containing Transcript |  1.569210 |  60 |  25 |
| Resolution of Sister Chromatid Cohesion |  1.569948 |  93 |  33 |
| APC/C-mediated degradation of cell cycle proteins |  1.577802 |  77 |  32 |
| Regulation of mitotic cell cycle |  1.577802 |  77 |  32 |
| Nuclear signaling by ERBB4 |  1.580528 |  18 |   6 |
| Vpr-mediated nuclear import of PICs |  1.583347 |  25 |   7 |
| RNA Polymerase II Transcription Termination |  1.584898 |  55 |  24 |
| Export of Viral Ribonucleoproteins from Nucleus |  1.586452 |  26 |   7 |
| NEP/NS2 Interacts with the Cellular Export Machinery |  1.586452 |  26 |   7 |
| RHO GTPase Effectors |  1.587712 | 196 |  71 |
| Cyclin A/B1/B2 associated events during G2/M transition |  1.590258 |  18 |  13 |
| Resolution of AP sites via the multiple-nucleotide patch replacement pathway |  1.592961 |  21 |   8 |
| MHC class II antigen presentation |  1.593630 |  76 |  20 |
| Disorders of transmembrane transporters |  1.593917 | 114 |  40 |
| Transport of Mature mRNAs Derived from Intronless Transcripts |  1.595779 |  34 |  11 |
| Abortive elongation of HIV-1 transcript in the absence of Tat |  1.603469 |  22 |  12 |
| Respiratory electron transport |  1.604512 |  83 |  41 |
| TCR signaling |  1.605577 |  86 |  33 |
| SUMOylation of SUMOylation proteins |  1.612490 |  28 |   8 |
| Transport of Mature Transcript to Cytoplasm |  1.613353 |  68 |  28 |
| Signaling by the B Cell Receptor (BCR) |  1.615003 |  92 |  29 |
| Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell |  1.617468 |  27 |   8 |
| Regulation of Glucokinase by Glucokinase Regulatory Protein |  1.618194 |  24 |   7 |
| Defective TPR may confer susceptibility towards thyroid papillary carcinoma (TPC) |  1.618194 |  24 |   7 |
| TP53 regulates transcription of additional cell cycle genes whose exact role in the p53 pathway remain uncertain |  1.620695 |  17 |   9 |
| Cell Cycle, Mitotic |  1.621489 | 403 | 131 |
| Translation |  1.623650 | 249 |  89 |
| Meiotic synapsis |  1.625306 |  23 |   6 |
| Transcriptional regulation by RUNX1 |  1.630654 | 142 |  45 |
| Transport of Mature mRNA Derived from an Intronless Transcript |  1.632642 |  33 |  11 |
| Regulation of RAS by GAPs |  1.633814 |  56 |  21 |
| Metabolic disorders of biological oxidation enzymes |  1.637941 |  16 |   5 |
| Degradation of the extracellular matrix |  1.638710 |  63 |  17 |
| Antigen processing-Cross presentation |  1.640411 |  75 |  30 |
| M Phase |  1.641700 | 284 |  97 |
| ER-Phagosome pathway |  1.641933 |  63 |  25 |
| PCNA-Dependent Long Patch Base Excision Repair |  1.642133 |  17 |   8 |
| SUMOylation of RNA binding proteins |  1.647227 |  37 |  13 |
| HIV Life Cycle |  1.657428 | 121 |  51 |
| Recognition of DNA damage by PCNA-containing replication complex |  1.662703 |  25 |   8 |
| Collagen chain trimerization |  1.663230 |  14 |   6 |
| Major pathway of rRNA processing in the nucleolus and cytosol |  1.667105 | 152 |  56 |
| Activation of NF-kappaB in B cells |  1.668749 |  57 |  22 |
| Inhibition of DNA recombination at telomere |  1.676651 |  21 |  12 |
| SUMOylation of ubiquitinylation proteins |  1.677816 |  30 |   9 |
| Assembly of the pre-replicative complex |  1.678062 |  76 |  31 |
| Cellular response to starvation |  1.680600 | 124 |  53 |
| Regulation of APC/C activators between G1/S and early anaphase |  1.685531 |  71 |  31 |
| RHO GTPases activate IQGAPs |  1.689130 |  12 |   7 |
| Synthesis of DNA |  1.692688 | 106 |  39 |
| RHO GTPases activate PKNs |  1.694726 |  30 |  18 |
| FGFR2 alternative splicing |  1.695097 |  20 |  12 |
| G1/S DNA Damage Checkpoints |  1.701501 |  58 |  23 |
| Lysine catabolism |  1.702366 |  11 |   4 |
| p53-Dependent G1 DNA Damage Response |  1.703777 |  57 |  23 |
| p53-Dependent G1/S DNA damage checkpoint |  1.703777 |  57 |  23 |
| Synthesis of active ubiquitin: roles of E1 and E2 enzymes |  1.707925 |  22 |  13 |
| APC/C:Cdc20 mediated degradation of Securin |  1.714550 |  58 |  24 |
| Condensation of Prophase Chromosomes |  1.715310 |  16 |   8 |
| CDK-mediated phosphorylation and removal of Cdc6 |  1.715389 |  64 |  29 |
| Deposition of new CENPA-containing nucleosomes at the centromere |  1.716501 |  21 |   9 |
| Nucleosome assembly |  1.716501 |  21 |   9 |
| Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins. |  1.717929 |  98 |  49 |
| Downstream signaling events of B Cell Receptor (BCR) |  1.718758 |  67 |  25 |
| Mitotic G1 phase and G1/S transition |  1.726049 | 125 |  44 |
| snRNP Assembly |  1.728076 |  43 |  13 |
| Metabolism of non-coding RNA |  1.728076 |  43 |  13 |
| Hedgehog ligand biogenesis |  1.728994 |  51 |  21 |
| Autodegradation of Cdh1 by Cdh1:APC/C |  1.729243 |  56 |  23 |
| S Phase |  1.730262 | 143 |  49 |
| APC/C:Cdh1 mediated degradation of Cdc20 and other APC/C:Cdh1 targeted proteins in late mitosis/early G1 |  1.732901 |  64 |  27 |
| Orc1 removal from chromatin |  1.734179 |  64 |  26 |
| SUMOylation of DNA damage response and repair proteins |  1.736060 |  60 |  19 |
| APC/C:Cdc20 mediated degradation of mitotic proteins |  1.738294 |  66 |  29 |
| Nuclear events mediated by NFE2L2 |  1.740403 |  76 |  30 |
| Hh mutants abrogate ligand secretion |  1.743268 |  49 |  21 |
| GLI3 is processed to GLI3R by the proteasome |  1.743908 |  51 |  21 |
| Regulation of RUNX2 expression and activity |  1.750065 |  56 |  21 |
| Interactions of Vpr with host cellular proteins |  1.752998 |  26 |   8 |
| APC:Cdc20 mediated degradation of cell cycle proteins prior to satisfation of the cell cycle checkpoint |  1.755729 |  65 |  29 |
| Aerobic respiration and respiratory electron transport |  1.756970 | 160 |  74 |
| Degradation of GLI1 by the proteasome |  1.758707 |  50 |  21 |
| Activation of APC/C and APC/C:Cdc20 mediated degradation of mitotic proteins |  1.764399 |  67 |  30 |
| Translesion Synthesis by POLH |  1.767317 |  16 |   5 |
| FBXL7 down-regulates AURKA during mitotic entry and in early mitosis |  1.773473 |  48 |  21 |
| tRNA processing in the nucleus |  1.775782 |  47 |  12 |
| Dectin-1 mediated noncanonical NF-kB signaling |  1.780723 |  52 |  21 |
| Degradation of GLI2 by the proteasome |  1.790817 |  50 |  21 |
| Detoxification of Reactive Oxygen Species |  1.792689 |  22 |  16 |
| Host Interactions of HIV factors |  1.800631 | 105 |  40 |
| Metabolism of polyamines |  1.802143 |  53 |  22 |
| Transcriptional regulation by small RNAs |  1.806275 |  44 |  20 |
| G1/S Transition |  1.806739 | 113 |  46 |
| Mitotic Anaphase |  1.808689 | 182 |  70 |
| Mitotic Metaphase and Anaphase |  1.809050 | 183 |  70 |
| DNA Replication |  1.812655 | 115 |  44 |
| Cdc20:Phospho-APC/C mediated degradation of Cyclin A |  1.817971 |  64 |  29 |
| Cellular response to hypoxia |  1.823028 |  65 |  26 |
| Cyclin A:Cdk2-associated events at S phase entry |  1.828468 |  76 |  29 |
| Separation of Sister Chromatids |  1.829491 | 147 |  57 |
| NS1 Mediated Effects on Host Pathways |  1.830066 |  31 |  10 |
| Viral Messenger RNA Synthesis |  1.832078 |  37 |  12 |
| HIV Infection |  1.838232 | 188 |  71 |
| Signaling by NOTCH4 |  1.838687 |  66 |  25 |
| DNA Replication Pre-Initiation |  1.842178 |  90 |  36 |
| Cyclin E associated events during G1/S transition |  1.852647 |  74 |  29 |
| Stabilization of p53 |  1.856818 |  48 |  20 |
| NIK-->noncanonical NF-kB signaling |  1.864619 |  51 |  21 |
| SCF-beta-TrCP mediated degradation of Emi1 |  1.866171 |  48 |  21 |
| Regulation of mRNA stability by proteins that bind AU-rich elements |  1.870693 |  74 |  32 |
| ABC transporter disorders |  1.870972 |  60 |  26 |
| PCP/CE pathway |  1.872663 |  69 |  27 |
| The role of GTSE1 in G2/M progression after G2 checkpoint |  1.877632 |  58 |  26 |
| The phototransduction cascade |  1.880133 |  11 |   4 |
| Inactivation, recovery and regulation of the phototransduction cascade |  1.880133 |  11 |   4 |
| SARS-CoV-1 modulates host translation machinery |  1.891711 |  30 |  20 |
| Autodegradation of the E3 ubiquitin ligase COP1 |  1.897373 |  46 |  20 |
| RUNX1 regulates transcription of genes involved in differentiation of HSCs |  1.905666 |  65 |  26 |
| Formation of ATP by chemiosmotic coupling |  1.906061 |  12 |   7 |
| Activation of Matrix Metalloproteinases |  1.908897 |  11 |   3 |
| Selenoamino acid metabolism |  1.910420 |  95 |  59 |
| Negative regulation of NOTCH4 signaling |  1.911847 |  47 |  23 |
| TNFR2 non-canonical NF-kB pathway |  1.913089 |  58 |  24 |
| SCF(Skp2)-mediated degradation of p27/p21 |  1.914750 |  54 |  24 |
| GSK3B and BTRC:CUL1-mediated-degradation of NFE2L2 |  1.920061 |  46 |  21 |
| Vpu mediated degradation of CD4 |  1.923179 |  45 |  21 |
| Degradation of DVL |  1.942347 |  50 |  22 |
| Asymmetric localization of PCP proteins |  1.959794 |  51 |  22 |
| Regulation of Apoptosis |  1.969835 |  47 |  23 |
| Cristae formation |  1.981851 |  24 |  11 |
| Heme biosynthesis |  1.986583 |  11 |   7 |
| ATF6 (ATF6-alpha) activates chaperone genes |  1.992543 |  10 |   8 |
| SARS-CoV-2 modulates host translation machinery |  1.993772 |  39 |  23 |
| Translation initiation complex formation |  1.993920 |  50 |  32 |
| Activation of the mRNA upon binding of the cap-binding complex and eIFs, and subsequent binding to 43S |  2.011257 |  51 |  33 |
| Oxygen-dependent proline hydroxylation of Hypoxia-inducible Factor Alpha |  2.014222 |  58 |  24 |
| Nonsense-Mediated Decay (NMD) |  2.020740 |  92 |  59 |
| Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC) |  2.020740 |  92 |  59 |
| Cross-presentation of soluble exogenous antigens (endosomes) |  2.021573 |  42 |  20 |
| Ribosomal scanning and start codon recognition |  2.032561 |  51 |  33 |
| Hh mutants are degraded by ERAD |  2.037604 |  48 |  21 |
| Formation of the ternary complex, and subsequently, the 43S complex |  2.037911 |  45 |  21 |
| Regulation of activated PAK-2p34 by proteasome mediated degradation |  2.050765 |  44 |  20 |
| Ubiquitin-dependent degradation of Cyclin D |  2.065511 |  46 |  21 |
| Degradation of AXIN |  2.067437 |  48 |  22 |
| Regulation of ornithine decarboxylase (ODC) |  2.069842 |  46 |  21 |
| Ubiquitin Mediated Degradation of Phosphorylated Cdc25A |  2.075990 |  45 |  20 |
| p53-Independent DNA Damage Response |  2.075990 |  45 |  20 |
| p53-Independent G1/S DNA damage checkpoint |  2.075990 |  45 |  20 |
| Selenocysteine synthesis |  2.096711 |  75 |  51 |
| Mitochondrial protein degradation |  2.108916 |  81 |  39 |
| AUF1 (hnRNP D0) binds and destabilizes mRNA |  2.115587 |  46 |  21 |
| Eukaryotic Translation Termination |  2.123653 |  74 |  52 |
| Collagen degradation |  2.124877 |  30 |  11 |
| Vif-mediated degradation of APOBEC3G |  2.128310 |  47 |  22 |
| ATF6 (ATF6-alpha) activates chaperones |  2.144804 |  12 |   9 |
| Response of EIF2AK4 (GCN2) to amino acid deficiency |  2.146668 |  82 |  56 |
| Maternal to zygotic transition (MZT) |  2.147872 |  40 |   5 |
| Regulation of RUNX3 expression and activity |  2.153411 |  47 |  20 |
| Defective CFTR causes cystic fibrosis |  2.153644 |  52 |  24 |
| Influenza Viral RNA Transcription and Replication |  2.159789 | 110 |  58 |
| Influenza Infection |  2.174912 | 125 |  58 |
| Peptide chain elongation |  2.176142 |  71 |  51 |
| Eukaryotic Translation Elongation |  2.197829 |  73 |  52 |
| Metabolism of porphyrins |  2.211398 |  17 |  10 |
| SRP-dependent cotranslational protein targeting to membrane |  2.220490 |  89 |  61 |
| Signaling by ROBO receptors |  2.257215 | 165 |  76 |
| Viral mRNA Translation |  2.257290 |  70 |  51 |
| Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC) |  2.298739 |  75 |  53 |
| Eukaryotic Translation Initiation |  2.312568 |  98 |  65 |
| Cap-dependent Translation Initiation |  2.312568 |  98 |  65 |
| L13a-mediated translational silencing of Ceruloplasmin expression |  2.371951 |  91 |  62 |
| Regulation of expression of SLITs and ROBOs |  2.388212 | 135 |  67 |
| Formation of a pool of free 40S subunits |  2.389773 |  83 |  59 |
| GTP hydrolysis and joining of the 60S ribosomal subunit |  2.392154 |  92 |  63 |



### IV-HD

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| Signaling by CSF1 (M-CSF) in myeloid cells | -1.869259 |  25 |  13 |
| PPARA activates gene expression | -1.810673 |  75 |  27 |
| Regulation of lipid metabolism by PPARalpha | -1.774851 |  77 |  27 |
| Negative regulation of MET activity | -1.772842 |  13 |  10 |
| Interleukin-6 family signaling | -1.764369 |  11 |   6 |
| Downregulation of SMAD2/3:SMAD4 transcriptional activity | -1.760074 |  19 |   8 |
| Plasma lipoprotein remodeling | -1.753294 |  17 |   8 |
| Signal amplification | -1.749093 |  13 |   6 |
| XBP1(S) activates chaperone genes | -1.720464 |  35 |  19 |
| Signaling by EGFR in Cancer | -1.712812 |  14 |   8 |
| FLT3 signaling in disease | -1.705644 |  19 |  12 |
| SUMOylation of DNA methylation proteins | -1.696502 |  12 |   5 |
| IRE1alpha activates chaperones | -1.688171 |  36 |  19 |
| Signaling by ERBB2 | -1.678333 |  25 |  10 |
| Sensory processing of sound by outer hair cells of the cochlea | -1.669455 |  17 |   9 |
| Signaling by FLT3 fusion proteins | -1.634968 |  14 |   9 |
| Platelet sensitization by LDL | -1.625555 |  13 |   6 |
| Signaling by FGFR3 | -1.609945 |  17 |   8 |
| Signaling by FGFR4 | -1.609945 |  17 |   8 |
| Negative regulation of the PI3K/AKT network | -1.606407 |  45 |  16 |
| Plasma lipoprotein assembly, remodeling, and clearance | -1.602755 |  37 |  12 |
| RHOBTB2 GTPase cycle | -1.595576 |  17 |   3 |
| Listeria monocytogenes entry into host cells | -1.593565 |  11 |   6 |
| PI5P, PP2A and IER3 Regulate PI3K/AKT Signaling | -1.591690 |  42 |  15 |
| Constitutive Signaling by Ligand-Responsive EGFR Cancer Variants | -1.572422 |  13 |   7 |
| Signaling by Ligand-Responsive EGFR Variants in Cancer | -1.572422 |  13 |   7 |
| N-glycan trimming in the ER and Calnexin/Calreticulin cycle | -1.565189 |  28 |  10 |
| Sensory processing of sound by inner hair cells of the cochlea | -1.559412 |  20 |  12 |
| Spry regulation of FGF signaling | -1.556017 |  10 |   4 |
| Negative regulation of FGFR1 signaling | -1.556017 |  10 |   4 |
| Negative regulation of FGFR3 signaling | -1.556017 |  10 |   4 |
| Negative regulation of FGFR4 signaling | -1.556017 |  10 |   4 |
| Sensory processing of sound | -1.551170 |  23 |  13 |
| Chondroitin sulfate/dermatan sulfate metabolism | -1.533951 |  17 |   6 |
| MyD88:MAL(TIRAP) cascade initiated on plasma membrane | -1.521704 |  56 |  22 |
| Toll Like Receptor TLR1:TLR2 Cascade | -1.521704 |  56 |  22 |
| Toll Like Receptor TLR6:TLR2 Cascade | -1.521704 |  56 |  22 |
| Toll Like Receptor 2 (TLR2) Cascade | -1.521704 |  56 |  22 |
| Constitutive Signaling by Aberrant PI3K in Cancer | -1.518995 |  24 |  10 |
| Interferon gamma signaling | -1.510473 |  30 |  14 |
| G alpha (12/13) signalling events | -1.502975 |  38 |  21 |
| Transcriptional regulation of white adipocyte differentiation | -1.487524 |  55 |  16 |
| Platelet homeostasis | -1.483787 |  30 |  14 |
| GPCR downstream signalling | -1.376518 | 160 |  63 |
| Signaling by GPCR | -1.346686 | 168 |  66 |
| Metabolism of lipids | -1.289286 | 384 |  98 |
| Adaptive Immune System |  1.209229 | 414 | 107 |
| Nervous system development |  1.250303 | 301 | 107 |
| Transcriptional Regulation by TP53 |  1.259967 | 242 |  88 |
| mRNA Splicing |  1.283339 | 172 |  71 |
| Processing of Capped Intron-Containing Pre-mRNA |  1.291811 | 226 |  89 |
| Axon guidance |  1.320393 | 292 | 106 |
| mRNA Splicing - Major Pathway |  1.332344 | 168 |  70 |
| Signaling by Nuclear Receptors |  1.339471 | 130 |  27 |
| Late Phase of HIV Life Cycle |  1.357150 | 108 |  38 |
| Resolution of Sister Chromatid Cohesion |  1.360317 |  88 |  35 |
| Signaling by Hedgehog |  1.368366 |  83 |  25 |
| Host Interactions of HIV factors |  1.385170 | 101 |  40 |
| Biological oxidations |  1.398811 |  67 |  21 |
| Mitotic Prophase |  1.403700 |  64 |  23 |
| Transcription-Coupled Nucleotide Excision Repair (TC-NER) |  1.407687 |  61 |  20 |
| Downstream signaling events of B Cell Receptor (BCR) |  1.427761 |  63 |  18 |
| Cellular response to chemical stress |  1.429464 | 145 |  53 |
| Regulation of RUNX2 expression and activity |  1.433935 |  53 |  20 |
| Cell Cycle Checkpoints |  1.437531 | 196 |  63 |
| Nuclear Envelope (NE) Reassembly |  1.442580 |  49 |  16 |
| HIV Life Cycle |  1.455856 | 114 |  42 |
| Degradation of DVL |  1.464981 |  47 |  22 |
| Mitotic G2-G2/M phases |  1.470910 | 134 |  48 |
| G2/M Transition |  1.470910 | 134 |  48 |
| Regulation of RAS by GAPs |  1.473905 |  53 |  21 |
| Complex I biogenesis |  1.478627 |  42 |  19 |
| Respiratory electron transport |  1.485687 |  78 |  30 |
| Diseases of programmed cell death |  1.494375 |  35 |  11 |
| Pre-NOTCH Transcription and Translation |  1.503348 |  24 |   8 |
| ECM proteoglycans |  1.505102 |  29 |  10 |
| Hedgehog 'off' state |  1.510686 |  68 |  29 |
| Deadenylation-dependent mRNA decay |  1.513691 |  40 |  24 |
| Extension of Telomeres |  1.515992 |  39 |  16 |
| Maternal to zygotic transition (MZT) |  1.518402 |  34 |  18 |
| Branched-chain amino acid catabolism |  1.542656 |  15 |   8 |
| Signaling by NOTCH4 |  1.546589 |  60 |  23 |
| Reproduction |  1.546975 |  42 |  14 |
| Polymerase switching on the C-strand of the telomere |  1.547057 |  17 |   8 |
| Degradation of GLI2 by the proteasome |  1.551086 |  48 |  20 |
| GLI3 is processed to GLI3R by the proteasome |  1.551086 |  48 |  20 |
| MAPK6/MAPK4 signaling |  1.557136 |  66 |  27 |
| Chromatin modifications during the maternal to zygotic transition (MZT) |  1.560145 |  10 |   5 |
| G2/M Checkpoints |  1.562135 | 107 |  36 |
| ATF6 (ATF6-alpha) activates chaperones |  1.567495 |  12 |   6 |
| Dectin-1 mediated noncanonical NF-kB signaling |  1.568827 |  50 |  20 |
| HIV Infection |  1.574117 | 179 |  70 |
| Stabilization of p53 |  1.574576 |  46 |  19 |
| Condensation of Prometaphase Chromosomes |  1.580346 |  10 |   6 |
| Diseases of mitotic cell cycle |  1.580636 |  25 |   7 |
| Degradation of GLI1 by the proteasome |  1.581502 |  47 |  20 |
| Meiosis |  1.584745 |  36 |  12 |
| M Phase |  1.584761 | 260 |  67 |
| Regulation of PTEN stability and activity |  1.590671 |  54 |  23 |
| Cell Cycle |  1.592066 | 445 | 110 |
| Regulation of mRNA stability by proteins that bind AU-rich elements |  1.592218 |  68 |  27 |
| TNFR2 non-canonical NF-kB pathway |  1.594814 |  55 |  23 |
| Transcriptional regulation by RUNX1 |  1.595005 | 134 |  32 |
| APC/C-mediated degradation of cell cycle proteins |  1.598256 |  72 |  29 |
| Regulation of mitotic cell cycle |  1.598256 |  72 |  29 |
| Autodegradation of the E3 ubiquitin ligase COP1 |  1.598353 |  45 |  19 |
| Activation of ATR in response to replication stress |  1.598869 |  33 |  13 |
| Hedgehog 'on' state |  1.601294 |  54 |  23 |
| Removal of the Flap Intermediate from the C-strand |  1.602114 |  14 |   7 |
| Neddylation |  1.602366 | 164 |  46 |
| G0 and Early G1 |  1.604597 |  22 |   8 |
| PCP/CE pathway |  1.617148 |  62 |  24 |
| Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins. |  1.619491 |  92 |  37 |
| FBXL7 down-regulates AURKA during mitotic entry and in early mitosis |  1.621595 |  46 |  20 |
| ESR-mediated signaling |  1.622440 | 100 |  25 |
| Recognition of DNA damage by PCNA-containing replication complex |  1.631508 |  24 |   8 |
| SCF-beta-TrCP mediated degradation of Emi1 |  1.631579 |  47 |  20 |
| Cellular response to hypoxia |  1.631995 |  62 |  25 |
| NFE2L2 regulating anti-oxidant/detoxification enzymes |  1.632816 |  13 |   8 |
| Lagging Strand Synthesis |  1.632984 |  18 |   9 |
| Mitochondrial iron-sulfur cluster biogenesis |  1.635088 |  11 |   7 |
| AUF1 (hnRNP D0) binds and destabilizes mRNA |  1.645605 |  45 |  19 |
| RUNX1 regulates transcription of genes involved in differentiation of HSCs |  1.647523 |  63 |  26 |
| Cyclin A/B1/B2 associated events during G2/M transition |  1.649633 |  17 |   8 |
| NIK-->noncanonical NF-kB signaling |  1.650658 |  49 |  20 |
| Orc1 removal from chromatin |  1.653949 |  62 |  24 |
| Telomere C-strand (Lagging Strand) Synthesis |  1.654943 |  25 |  12 |
| GSK3B and BTRC:CUL1-mediated-degradation of NFE2L2 |  1.655468 |  45 |  20 |
| DNA Damage Bypass |  1.655633 |  37 |  11 |
| Activation of the pre-replicative complex |  1.657165 |  28 |  12 |
| Formation of ATP by chemiosmotic coupling |  1.658186 |  12 |  10 |
| Switching of origins to a post-replicative state |  1.660390 |  76 |  25 |
| MET activates PTK2 signaling |  1.662472 |  12 |   4 |
| Mitotic Anaphase |  1.663883 | 172 |  53 |
| Beta-catenin independent WNT signaling |  1.666858 |  84 |  33 |
| Formation of tubulin folding intermediates by CCT/TriC |  1.668526 |  12 |   7 |
| Cellular response to starvation |  1.672384 | 118 |  61 |
| p53-Dependent G1 DNA Damage Response |  1.678593 |  51 |  21 |
| p53-Dependent G1/S DNA damage checkpoint |  1.678593 |  51 |  21 |
| RHO GTPases activate IQGAPs |  1.681201 |  11 |   3 |
| Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell |  1.682871 |  23 |   5 |
| Mitotic Metaphase and Anaphase |  1.683770 | 173 |  54 |
| Resolution of AP sites via the multiple-nucleotide patch replacement pathway |  1.684261 |  20 |   8 |
| Vpu mediated degradation of CD4 |  1.686267 |  44 |  20 |
| Aberrant regulation of mitotic cell cycle due to RB1 defects |  1.688310 |  23 |   7 |
| Resolution of Abasic Sites (AP sites) |  1.690718 |  26 |  11 |
| Telomere Maintenance |  1.691783 |  56 |  23 |
| PCNA-Dependent Long Patch Base Excision Repair |  1.693297 |  17 |   8 |
| Estrogen-dependent gene expression |  1.695618 |  67 |  19 |
| Depolymerization of the Nuclear Lamina |  1.697760 |  11 |   6 |
| Chromosome Maintenance |  1.698144 |  72 |  28 |
| M-decay: degradation of maternal mRNAs by maternally stored factors |  1.703932 |  16 |  11 |
| APC:Cdc20 mediated degradation of cell cycle proteins prior to satisfation of the cell cycle checkpoint |  1.705562 |  60 |  25 |
| Transcriptional regulation by small RNAs |  1.708220 |  42 |  16 |
| G1/S DNA Damage Checkpoints |  1.708472 |  52 |  22 |
| Citric acid cycle (TCA cycle) |  1.709215 |  27 |  14 |
| Hedgehog ligand biogenesis |  1.710239 |  49 |  20 |
| Maturation of TCA enzymes and regulation of TCA cycle |  1.710956 |  15 |  12 |
| Nuclear events mediated by NFE2L2 |  1.713082 |  71 |  32 |
| Hh mutants are degraded by ERAD |  1.717536 |  47 |  20 |
| Hh mutants abrogate ligand secretion |  1.717536 |  47 |  20 |
| Negative regulation of NOTCH4 signaling |  1.721080 |  45 |  21 |
| Separation of Sister Chromatids |  1.726248 | 138 |  44 |
| Aberrant regulation of mitotic G1/S transition in cancer due to RB1 defects |  1.729533 |  10 |   4 |
| Defective binding of RB1 mutants to E2F1,(E2F2, E2F3) |  1.729533 |  10 |   4 |
| Meiotic synapsis |  1.730601 |  20 |   4 |
| Cell Cycle, Mitotic |  1.734849 | 368 |  99 |
| Glutathione conjugation |  1.734893 |  13 |   8 |
| Mitotic Telophase/Cytokinesis |  1.735760 |  11 |   2 |
| APC/C:Cdc20 mediated degradation of mitotic proteins |  1.741056 |  61 |  25 |
| Regulation of APC/C activators between G1/S and early anaphase |  1.748150 |  66 |  28 |
| G1/S-Specific Transcription |  1.754016 |  23 |   9 |
| Gastrulation |  1.757620 |  59 |  25 |
| Translation |  1.762651 | 230 | 130 |
| rRNA processing |  1.768383 | 164 |  85 |
| Activation of APC/C and APC/C:Cdc20 mediated degradation of mitotic proteins |  1.771163 |  62 |  26 |
| APC/C:Cdh1 mediated degradation of Cdc20 and other APC/C:Cdh1 targeted proteins in late mitosis/early G1 |  1.775303 |  60 |  25 |
| APC/C:Cdc20 mediated degradation of Securin |  1.778978 |  54 |  23 |
| The phototransduction cascade |  1.783309 |  11 |   4 |
| Inactivation, recovery and regulation of the phototransduction cascade |  1.783309 |  11 |   4 |
| Aerobic respiration and respiratory electron transport |  1.784765 | 140 |  57 |
| Transcriptional regulation by RUNX3 |  1.790810 |  66 |  25 |
| Metabolism of polyamines |  1.794662 |  51 |  25 |
| Regulation of activated PAK-2p34 by proteasome mediated degradation |  1.794734 |  43 |  19 |
| Cdc20:Phospho-APC/C mediated degradation of Cyclin A |  1.796957 |  59 |  25 |
| Regulation of ornithine decarboxylase (ODC) |  1.804391 |  46 |  23 |
| The role of GTSE1 in G2/M progression after G2 checkpoint |  1.809942 |  56 |  25 |
| SCF(Skp2)-mediated degradation of p27/p21 |  1.810370 |  51 |  22 |
| Autodegradation of Cdh1 by Cdh1:APC/C |  1.817242 |  52 |  23 |
| CDK-mediated phosphorylation and removal of Cdc6 |  1.822288 |  59 |  25 |
| Regulation of Apoptosis |  1.830192 |  44 |  20 |
| TP53 regulates transcription of additional cell cycle genes whose exact role in the p53 pathway remain uncertain |  1.835191 |  15 |  12 |
| Assembly of the pre-replicative complex |  1.841773 |  72 |  26 |
| Cristae formation |  1.843428 |  24 |  14 |
| rRNA processing in the nucleus and cytosol |  1.862370 | 156 |  83 |
| Mitotic G1 phase and G1/S transition |  1.865521 | 116 |  32 |
| Ubiquitin-dependent degradation of Cyclin D |  1.867879 |  45 |  20 |
| SARS-CoV-1 modulates host translation machinery |  1.870299 |  29 |  21 |
| DNA strand elongation |  1.875665 |  29 |  13 |
| Degradation of AXIN |  1.881169 |  47 |  22 |
| Cyclin E associated events during G1/S transition |  1.883151 |  68 |  29 |
| Condensation of Prophase Chromosomes |  1.884843 |  15 |   9 |
| Asymmetric localization of PCP proteins |  1.898340 |  46 |  21 |
| Oxygen-dependent proline hydroxylation of Hypoxia-inducible Factor Alpha |  1.902168 |  55 |  25 |
| Cyclin A:Cdk2-associated events at S phase entry |  1.903061 |  70 |  30 |
| Ubiquitin Mediated Degradation of Phosphorylated Cdc25A |  1.905882 |  43 |  20 |
| p53-Independent DNA Damage Response |  1.905882 |  43 |  20 |
| p53-Independent G1/S DNA damage checkpoint |  1.905882 |  43 |  20 |
| Base Excision Repair |  1.917680 |  32 |  14 |
| Major pathway of rRNA processing in the nucleolus and cytosol |  1.926930 | 146 |  78 |
| Influenza Infection |  1.933389 | 122 |  70 |
| Regulation of RUNX3 expression and activity |  1.948693 |  46 |  20 |
| Influenza Viral RNA Transcription and Replication |  1.967828 | 108 |  66 |
| Synthesis of DNA |  1.972000 | 101 |  37 |
| SARS-CoV-2 modulates host translation machinery |  1.982142 |  37 |  22 |
| Cross-presentation of soluble exogenous antigens (endosomes) |  2.001799 |  41 |  20 |
| Vif-mediated degradation of APOBEC3G |  2.007795 |  45 |  22 |
| Nonsense-Mediated Decay (NMD) |  2.013469 |  89 |  50 |
| Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC) |  2.013469 |  89 |  50 |
| DNA Replication Pre-Initiation |  2.015885 |  86 |  33 |
| SRP-dependent cotranslational protein targeting to membrane |  2.017591 |  88 |  53 |
| G1/S Transition |  2.026620 | 107 |  42 |
| DNA Replication |  2.039751 | 110 |  43 |
| Somitogenesis |  2.052897 |  40 |  20 |
| Mitochondrial protein degradation |  2.057051 |  73 |  37 |
| Detoxification of Reactive Oxygen Species |  2.057222 |  20 |  12 |
| Translation initiation complex formation |  2.087305 |  49 |  31 |
| Activation of the mRNA upon binding of the cap-binding complex and eIFs, and subsequent binding to 43S |  2.088340 |  50 |  32 |
| Formation of paraxial mesoderm |  2.096674 |  46 |  22 |
| Signaling by ROBO receptors |  2.098334 | 156 |  85 |
| Metabolism of amino acids and derivatives |  2.111932 | 224 | 117 |
| Formation of the ternary complex, and subsequently, the 43S complex |  2.126980 |  44 |  30 |
| Heme biosynthesis |  2.128247 |  11 |   7 |
| Ribosomal scanning and start codon recognition |  2.140236 |  50 |  32 |
| S Phase |  2.163457 | 133 |  46 |
| Selenoamino acid metabolism |  2.168885 |  91 |  57 |
| Selenocysteine synthesis |  2.171936 |  72 |  46 |
| Viral mRNA Translation |  2.178346 |  69 |  46 |
| Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC) |  2.206450 |  74 |  52 |
| Eukaryotic Translation Termination |  2.209620 |  72 |  47 |
| Peptide chain elongation |  2.240268 |  69 |  48 |
| Eukaryotic Translation Elongation |  2.250433 |  71 |  47 |
| Metabolism of porphyrins |  2.263622 |  16 |  11 |
| Response of EIF2AK4 (GCN2) to amino acid deficiency |  2.265917 |  80 |  51 |
| Eukaryotic Translation Initiation |  2.325314 |  97 |  64 |
| Cap-dependent Translation Initiation |  2.325314 |  97 |  64 |
| L13a-mediated translational silencing of Ceruloplasmin expression |  2.365928 |  90 |  61 |
| Regulation of expression of SLITs and ROBOs |  2.373838 | 131 |  76 |
| GTP hydrolysis and joining of the 60S ribosomal subunit |  2.394554 |  91 |  57 |
| Formation of a pool of free 40S subunits |  2.394803 |  82 |  55 |



## 10 WPI

### DNA vaccine

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| Striated Muscle Contraction | -2.078391 |  27 |  12 |
| Cell Cycle, Mitotic |  1.325285 | 426 | 288 |
| RHO GTPase Effectors |  1.394363 | 214 | 124 |



### EOMES

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| Neuronal System | 1.394553 | 283 | 99 |

### GATA3

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| rRNA modification in the nucleus and cytosol | -1.761944 |  57 |  55 |
| Striated Muscle Contraction | -1.626816 |  27 |   9 |
| Signaling by GPCR |  1.191060 | 389 | 198 |
| GPCR downstream signalling |  1.193191 | 348 | 168 |
| SLC-mediated transmembrane transport |  1.225929 | 177 |  72 |
| Transmission across Chemical Synapses |  1.241173 | 189 |  89 |
| GPCR ligand binding |  1.250982 | 200 |  99 |
| G alpha (i) signalling events |  1.259988 | 156 |  81 |
| Neuronal System |  1.287592 | 283 | 142 |
| Class A/1 (Rhodopsin-like receptors) |  1.299060 | 139 |  71 |
| Protein-protein interactions at synapses |  1.312714 |  65 |  36 |
| Degradation of the extracellular matrix |  1.321021 |  86 |  32 |
| Potassium Channels |  1.328675 |  65 |  34 |
| Post-translational protein phosphorylation |  1.334892 |  66 |  20 |
| Sensory processing of sound by inner hair cells of the cochlea |  1.353306 |  54 |  26 |
| Sensory processing of sound |  1.355838 |  57 |  27 |
| Visual phototransduction |  1.372641 |  59 |  30 |
| Cell-cell junction organization |  1.374837 |  64 |  30 |
| Sensory Perception |  1.394828 | 134 |  65 |
| Adherens junctions interactions |  1.405572 |  45 |  20 |

### IV-HD

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| Striated Muscle Contraction | -2.643247 | 27 | 15 |



### IV-LD

| Description | NES | setSize | Count |
| --- | --- | --- | --- |
| DNA strand elongation | -2.954534 |  29 |  20 |
| Activation of the pre-replicative complex | -2.823981 |  28 |  20 |
| Activation of ATR in response to replication stress | -2.664441 |  33 |  19 |
| Unwinding of DNA | -2.513751 |  11 |  10 |
| Processive synthesis on the lagging strand | -2.438397 |  14 |  10 |
| Lagging Strand Synthesis | -2.432991 |  18 |  10 |
| Removal of the Flap Intermediate | -2.412073 |  13 |  10 |
| Telomere C-strand (Lagging Strand) Synthesis | -2.368757 |  28 |  13 |
| Polymerase switching on the C-strand of the telomere | -2.212591 |  20 |  10 |
| Polymerase switching | -2.190069 |  12 |   7 |
| Leading Strand Synthesis | -2.190069 |  12 |   7 |
| tRNA processing | -2.178221 |  88 |  31 |
| tRNA modification in the nucleus and cytosol | -2.166711 |  33 |  18 |
| mRNA decay by 3' to 5' exoribonuclease | -2.158842 |  14 |   7 |
| Mismatch repair (MMR) directed by MSH2:MSH6 (MutSalpha) | -2.108978 |  13 |   7 |
| Nuclear import of Rev protein | -2.061305 |  29 |  12 |
| Chromosome Maintenance | -2.025494 |  80 |  36 |
| Mitochondrial translation | -2.003138 |  82 |  40 |
| Synthesis of glycosylphosphatidylinositol (GPI) | -1.990159 |  14 |   7 |
| Interactions of Vpr with host cellular proteins | -1.976322 |  27 |  11 |
| Extension of Telomeres | -1.967876 |  45 |  23 |
| Mitochondrial translation elongation | -1.967231 |  76 |  29 |
| Transport of the SLBP Dependant Mature mRNA | -1.957045 |  29 |  12 |
| Mitochondrial translation initiation | -1.952516 |  76 |  28 |
| Transcriptional regulation by small RNAs | -1.933773 |  45 |  29 |
| Mitochondrial translation termination | -1.923750 |  76 |  37 |
| Postmitotic nuclear pore complex (NPC) reformation | -1.913919 |  24 |  10 |
| tRNA processing in the nucleus | -1.905964 |  48 |  17 |
| HDR through Homologous Recombination (HRR) | -1.900717 |  56 |  20 |
| Mismatch Repair | -1.894024 |  14 |   7 |
| Diseases of DNA repair | -1.889344 |  43 |  13 |
| Regulation of Glucokinase by Glucokinase Regulatory Protein | -1.883906 |  25 |  10 |
| Defective TPR may confer susceptibility towards thyroid papillary carcinoma (TPC) | -1.883906 |  25 |  10 |
| Telomere Maintenance | -1.878765 |  62 |  28 |
| Transport of Ribonucleoproteins into the Host Nucleus | -1.860285 |  26 |  10 |
| RNA Polymerase I Promoter Escape | -1.843944 |  30 |  11 |
| Transport of the SLBP independent Mature mRNA | -1.841741 |  28 |  11 |
| Vpr-mediated nuclear import of PICs | -1.825700 |  26 |  10 |
| Transport of Mature mRNAs Derived from Intronless Transcripts | -1.821860 |  35 |  13 |
| snRNP Assembly | -1.819523 |  44 |  19 |
| Metabolism of non-coding RNA | -1.819523 |  44 |  19 |
| SUMOylation of SUMOylation proteins | -1.814853 |  29 |  11 |
| DNA Replication Pre-Initiation | -1.805919 |  90 |  43 |
| Interactions of Rev with host cellular proteins | -1.801442 |  31 |  12 |
| TRAF6 mediated IRF7 activation in TLR7/8 or 9 signaling | -1.781758 |  10 |   4 |
| Rev-mediated nuclear export of HIV RNA | -1.722768 |  29 |  11 |
| Deposition of new CENPA-containing nucleosomes at the centromere | -1.716512 |  21 |   9 |
| Nucleosome assembly | -1.716512 |  21 |   9 |
| Transport of Mature mRNA Derived from an Intronless Transcript | -1.710378 |  34 |  12 |
| E2F mediated regulation of DNA replication | -1.703348 |  16 |   7 |
| Export of Viral Ribonucleoproteins from Nucleus | -1.701043 |  27 |  10 |
| NEP/NS2 Interacts with the Cellular Export Machinery | -1.701043 |  27 |  10 |
| Mismatch repair (MMR) directed by MSH2:MSH3 (MutSbeta) | -1.691194 |  13 |   6 |
| Resolution of D-loop Structures through Holliday Junction Intermediates | -1.685023 |  30 |  11 |
| Presynaptic phase of homologous DNA pairing and strand exchange | -1.675385 |  32 |  12 |
| Resolution of D-Loop Structures | -1.663293 |  31 |  11 |
| Meiotic recombination | -1.660692 |  22 |   8 |
| SUMOylation of ubiquitinylation proteins | -1.643981 |  31 |  11 |
| Diseases of DNA Double-Strand Break Repair | -1.640567 |  33 |  12 |
| Defective homologous recombination repair (HRR) due to BRCA2 loss of function | -1.640567 |  33 |  12 |
| Homologous DNA Pairing and Strand Exchange | -1.628077 |  35 |  12 |
| Nuclear Pore Complex (NPC) Disassembly | -1.610885 |  31 |  11 |
| Defective homologous recombination repair (HRR) due to BRCA1 loss of function | -1.605874 |  20 |   8 |
| Defective homologous recombination repair (HRR) due to PALB2 loss of function | -1.605874 |  20 |   8 |
| Defective HDR through Homologous Recombination Repair (HRR) due to PALB2 loss of BRCA1 binding function | -1.605874 |  20 |   8 |
| Defective HDR through Homologous Recombination Repair (HRR) due to PALB2 loss of BRCA2/RAD51/RAD51C binding function | -1.605874 |  20 |   8 |
| Impaired BRCA2 binding to PALB2 | -1.605874 |  20 |   8 |
| Gap-filling DNA repair synthesis and ligation in GG-NER | -1.583184 |  20 |   7 |
| DNA Replication | -1.578501 | 115 |  45 |
| Resolution of D-loop Structures through Synthesis-Dependent Strand Annealing (SDSA) | -1.570014 |  22 |   8 |
| RNA Polymerase I Transcription Termination | -1.552286 |  26 |   8 |
| Gap-filling DNA repair synthesis and ligation in TC-NER | -1.539402 |  53 |  29 |
| Dual incision in TC-NER | -1.537704 |  53 |  28 |
| Viral Messenger RNA Synthesis | -1.529689 |  38 |  11 |
| Negative epigenetic regulation of rRNA expression | -1.509520 |  46 |  13 |
| Synthesis of DNA | -1.504299 | 106 |  44 |
| NS1 Mediated Effects on Host Pathways | -1.489298 |  32 |  10 |
| Fanconi Anemia Pathway | -1.487228 |  33 |  12 |
| NoRC negatively regulates rRNA expression | -1.486856 |  43 |  12 |
| SUMOylation of DNA replication proteins | -1.478088 |  39 |  14 |
| RNA Polymerase I Promoter Clearance | -1.458619 |  48 |  17 |
| RNA Polymerase I Transcription | -1.443085 |  49 |  17 |
| Formation of the ternary complex, and subsequently, the 43S complex | -1.421806 |  45 |  44 |
| Transport of Mature mRNA derived from an Intron-Containing Transcript | -1.408246 |  61 |  18 |
| Transport of Mature Transcript to Cytoplasm | -1.403856 |  69 |  20 |
| Assembly of the pre-replicative complex | -1.391016 |  76 |  33 |
| Homology Directed Repair | -1.349732 |  89 |  28 |
| HDR through Homologous Recombination (HRR) or Single Strand Annealing (SSA) | -1.325622 |  83 |  23 |
| DNA Double-Strand Break Repair | -1.314092 | 114 |  34 |
| Axon guidance |  1.239477 | 361 | 103 |
| Hemostasis |  1.249844 | 362 | 101 |
| MAPK1/MAPK3 signaling |  1.251507 | 175 |  64 |
| Nervous system development |  1.264132 | 376 | 110 |
| RAF/MAP kinase cascade |  1.265024 | 172 |  63 |
| Extracellular matrix organization |  1.266893 | 157 |  59 |
| Signaling by GPCR |  1.267244 | 258 | 111 |
| GPCR downstream signalling |  1.267495 | 241 | 104 |
| RAC1 GTPase cycle |  1.286911 | 148 |  79 |
| RHO GTPase cycle |  1.310042 | 363 | 180 |
| Adipogenesis |  1.324017 |  81 |  34 |
| Ion channel transport |  1.333062 |  89 |  23 |
| Cargo recognition for clathrin-mediated endocytosis |  1.351754 |  65 |  18 |
| Factors involved in megakaryocyte development and platelet production |  1.353268 |  98 |  40 |
| Cell-Cell communication |  1.354079 |  87 |  44 |
| G alpha (q) signalling events |  1.356045 |  88 |  32 |
| Transport of small molecules |  1.361181 | 417 | 116 |
| CDC42 GTPase cycle |  1.376493 | 125 |  75 |
| SLC-mediated transmembrane transport |  1.376789 | 129 |  35 |
| Transcriptional regulation of white adipocyte differentiation |  1.387855 |  62 |  25 |
| Cell surface interactions at the vascular wall |  1.396451 |  77 |  39 |
| Fatty acid metabolism |  1.396874 | 119 |  40 |
| Degradation of the extracellular matrix |  1.404099 |  67 |  27 |
| Peptide ligand-binding receptors |  1.413956 |  37 |  15 |
| RHOD GTPase cycle |  1.422377 |  45 |  28 |
| NRAGE signals death through JNK |  1.423426 |  43 |  20 |
| RHOJ GTPase cycle |  1.425257 |  48 |  31 |
| L1CAM interactions |  1.431323 |  67 |  32 |
| Circadian Clock |  1.436316 |  56 |  23 |
| RND3 GTPase cycle |  1.436619 |  36 |  14 |
| TBC/RABGAPs |  1.439978 |  34 |  17 |
| Signaling by NOTCH3 |  1.449675 |  33 |   9 |
| RHOA GTPase cycle |  1.451754 | 127 |  53 |
| Cardiac conduction |  1.459053 |  49 |  15 |
| Arachidonic acid metabolism |  1.462290 |  31 |   5 |
| Metabolism of water-soluble vitamins and cofactors |  1.466094 |  88 |  27 |
| Sensory processing of sound by outer hair cells of the cochlea |  1.478480 |  25 |   9 |
| Signaling by ERBB2 in Cancer |  1.481454 |  17 |   8 |
| Signaling by ERBB2 KD Mutants |  1.481454 |  17 |   8 |
| Uptake and actions of bacterial toxins |  1.484117 |  16 |   6 |
| Nucleotide catabolism |  1.486791 |  24 |   4 |
| EGFR downregulation |  1.490518 |  21 |   7 |
| Constitutive Signaling by Aberrant PI3K in Cancer |  1.491620 |  35 |  17 |
| Regulation of lipid metabolism by PPARalpha |  1.492182 |  93 |  40 |
| Ion transport by P-type ATPases |  1.494421 |  32 |  13 |
| Sensory Perception |  1.494804 |  80 |  33 |
| Glycogen synthesis |  1.494930 |  10 |   2 |
| PPARA activates gene expression |  1.497057 |  91 |  40 |
| Bile acid and bile salt metabolism |  1.497313 |  22 |  10 |
| Signal transduction by L1 |  1.497539 |  20 |   3 |
| Endogenous sterols |  1.500532 |  13 |   7 |
| Assembly of active LPL and LIPC lipase complexes |  1.504366 |  11 |   5 |
| Glyoxylate metabolism and glycine degradation |  1.504830 |  14 |   8 |
| Non-integrin membrane-ECM interactions |  1.508077 |  35 |  26 |
| Syndecan interactions |  1.509090 |  21 |  13 |
| Biological oxidations |  1.511144 |  94 |  27 |
| Interaction between L1 and Ankyrins |  1.514887 |  11 |   6 |
| Attachment and Entry |  1.524825 |  11 |   8 |
| Basigin interactions |  1.525764 |  20 |   7 |
| SUMOylation of intracellular receptors |  1.530752 |  19 |  10 |
| Sulfur amino acid metabolism |  1.533032 |  24 |  10 |
| Glycosphingolipid biosynthesis |  1.533504 |  10 |   4 |
| NR1H2 and NR1H3-mediated signaling |  1.546182 |  32 |  13 |
| Nitric oxide stimulates guanylate cyclase |  1.564377 |  10 |   6 |
| Degradation of cysteine and homocysteine |  1.564391 |  13 |   6 |
| Signaling by ERBB2 TMD/JMD mutants |  1.568626 |  14 |   8 |
| Plasma lipoprotein clearance |  1.569191 |  26 |  11 |
| Ion homeostasis |  1.569326 |  27 |  10 |
| Metabolism of vitamins and cofactors |  1.576338 | 131 |  38 |
| CD209 (DC-SIGN) signaling |  1.588426 |  17 |   9 |
| Erythropoietin activates RAS |  1.590764 |  12 |   8 |
| Metabolism of folate and pterines |  1.597534 |  11 |   5 |
| Visual phototransduction |  1.606751 |  37 |  20 |
| Phase I - Functionalization of compounds |  1.608231 |  41 |  13 |
| Keratinization |  1.608706 |  21 |  11 |
| Formation of the cornified envelope |  1.608706 |  21 |  11 |
| Serotonin Neurotransmitter Release Cycle |  1.610092 |  10 |   2 |
| SHC1 events in ERBB2 signaling |  1.617309 |  13 |  10 |
| NOTCH3 Activation and Transmission of Signal to the Nucleus |  1.625611 |  18 |   5 |
| Activation of IRF3, IRF7 mediated by TBK1, IKKε (IKBKE) |  1.626696 |  11 |   1 |
| Nuclear Receptor transcription pathway |  1.628092 |  32 |  17 |
| Post-translational protein phosphorylation |  1.634458 |  51 |  14 |
| BMAL1:CLOCK,NPAS2 activates circadian gene expression |  1.640093 |  22 |  11 |
| Gluconeogenesis |  1.641454 |  14 |   7 |
| Activation of Matrix Metalloproteinases |  1.648243 |  11 |   4 |
| Triglyceride metabolism |  1.649731 |  20 |  10 |
| Cellular hexose transport |  1.662707 |  11 |   6 |
| Metabolism of fat-soluble vitamins |  1.665754 |  25 |  13 |
| Triglyceride catabolism |  1.666912 |  13 |   6 |
| Cytochrome P450 - arranged by substrate type |  1.671828 |  22 |   7 |
| Regulation of Complement cascade |  1.673959 |  11 |   3 |
| Plasma lipoprotein assembly, remodeling, and clearance |  1.681026 |  46 |  22 |
| Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs) |  1.690045 |  57 |  17 |
| Scavenging by Class A Receptors |  1.691384 |  15 |   7 |
| Binding and Uptake of Ligands by Scavenger Receptors |  1.704992 |  24 |   8 |
| Retinoid metabolism and transport |  1.709875 |  21 |   8 |
| Metabolism of steroids |  1.721594 |  97 |  36 |
| Plasma lipoprotein remodeling |  1.736262 |  19 |   9 |
| Regulation of cholesterol biosynthesis by SREBP (SREBF) |  1.743377 |  48 |  21 |
| Activation of gene expression by SREBF (SREBP) |  1.821012 |  37 |  20 |
| Cholesterol biosynthesis |  1.868195 |  22 |  14 |
| Complement cascade |  1.892080 |  17 |   6 |
| Initial triggering of complement |  1.915240 |  11 |   6 |



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
