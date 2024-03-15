# Investigating Causal Effects of Autoimmune Diseases and Monoclonal Antibody Drug Targets on Pregnancy Outcomes: a Mendelian Randomization study

## BACKGROUND
### AUTOIMMUNE DISEASES AND PREGNANCY OUTCOMES
Chronic autoimmune diseases affect around 13% of women in the UK (1), encompassing conditions such as rheumatic diseases, inflammatory bowel disease, and multiple sclerosis. More common amongst female patients, these diseases are often commonly diagnosed before menopause; in studies of women in high income countries, incidence of IBD is highest amongst those aged 20-29 (2), while women aged 40-49 have the highest incidence of MS (3) and systemic lupus erythematous (4).

Observational evidence suggests an underlying chronic autoimmune disease is associated with an increased odds of adverse pregnancy outcomes. Meta-analyses of observational studies have found that active inflammatory bowel disease is associated with an increased risk of pre-term birth, small-for-gestational-age birthweight, stillbirth, and congenital anomalies (5). Rheumatic diseases, encompassing both rheumatic arthritis and systemic lupus erythematosus, have been associated with increased risk of preterm birth, intrauterine growth restriction (6), and pre-eclampsia (7). Adverse outcomes have similarly been associated with multiple sclerosis (8) and axial spondyloarthritis (9). Though there is comparatively less evidence on psoriasis, a recent population-wide study of Sweden and Denmark found an association between maternal psoriasis and increased risks of gestational diabetes and hypertensive disorders of pregnancy (10). Given that such observational studies may be subject to some residual confounding, the triangulation of these associations with other sources of evidence could provide further useful insights.

### MONOCLONAL ANTIBODY TREATMENTS FOR AUTOIMMUNE DISEASES AND PREGNANCY OUTCOMES
Given that autoimmune conditions are common in women of reproductive age, it is crucial to understand how to manage them during pregnancy. Ensuring that the disease symptoms are adequately treated must be weighed against ensuring that treatments do not increase risk of adverse pregnancy and perinatal outcomes beyond the disease themselves. Current guidelines for best clinical practice with pregnant patients recommend treating active disease and where possible ensuring clinical remission, for IBD (11,12) and rheumatic diseases (13).

Treatment for autoimmune diseases has been revolutionised by the development of novel therapeutic monoclonal antibodies (mAbs). First produced in Nobel-winning experiments in 1975 (14), therapeutic applications to cancer and autoimmune disease swiftly followed, with the first drug receiving approval in the US in 1986 (15,16). Approximately 175 mAbs have now been approved by the US and EU and more are developed each year (17); they now comprise one of the bestselling drug sectors globally (16). Falling within the biological products class of medicines (biologics), their ability to bind highly specific antigen targets enables inhibition of pro-inflammatory signalling molecules such as tumour necrosis factor and interleukins within blood plasma.

Despite these advances, the safety of mAbs during pregnancy is largely unknown (18). Evidence and guidelines can vary between particular drugs and diseases, particularly around whether therapy should be discontinued in the third trimester (13,19). Research suggests that most biologics can cross the placenta (20), with greater rates of transfer occurring later in gestation (18). Anecdotal evidence points to potential risks. Sadly, one infant exposed to infliximab in utero died of disseminated BCG infection after receiving a vaccination at 3 months, suggesting potential for tragic harm (21). Further reports have included infant thrombocytopenia following natalizumab exposure during the third trimester (22). On the other hand, eculizumab has been used anecdotally to treat severe preeclampsia and associated acute kidney injury post-delivery (23,24), prompting speculation about its therapeutic potential here (25). 

The widespread exclusion of pregnant participants from randomised clinical trials (RCTs) (26,27), despite compelling arguments for greater inclusion (28,29), has resulted in a paucity of gold-standard evidence on the safety and efficacy of mAbs during pregnancy. Despite widespread prescription drug use during pregnancy and the needs of patients to continue to manage chronic conditions, very few medicines are explicitly licensed for use during pregnancy in Europe (30). 

In the absence of RCT evidence, various causal inference methodologies have been developed to address causal questions using available observational data. Mendelian randomization (MR) is one causal inference method which leverages genetic data. MR is commonly implemented in instrumental variable analyses which treat genetic variants associated with an exposure of interest as randomly inherited instrumental variables, ‘genetic instruments’, to index an exposures’ causal effect on an outcome (31). Three core instrumental variable assumptions are required for an MR analysis: (1) genetic variants are associated with the exposure (relevance), (2) there are no confounders of the genetic variant-outcome associations (independence), and (3) genetic variants only act on the outcome through their effect on the exposure of interest (the exclusion restriction). Drug target MR, a specialised application of MR, uses genetic variants to estimate the effects of pharmacological modulation of drug targets (32). While this approach has been used successfully within cardiovascular disease contexts (33), such as the use of variants in the CETP gene to mimic CETP inhibitors and investigate their effects on lipoproteins (34), it has only rarely been applied to characterise drug safety and efficacy during pregnancy (35). Monoclonal antibodies are a tractable drug class for drug target MR to proxy via protein quantitative trait loci (pQTLs), since mAbs typically target circulating proteins (32).

This study aims to characterise the safety profile and possible benefits of modulating immune signalling during pregnancy, using drug target MR. First, the causal effects of autoimmune diseases on adverse pregnancy and perinatal outcomes will be ascertained. Second, the causal effects of modulating mAb target proteins will be assessed for the same outcomes. Together, these analyses aim to characterise the safety profile of mAbs and whether these may have beneficial effects mediated through controlling autoimmune disease.

Efforts to redress this evidence gap around mAb safety in pregnancy have recently gained momentum, for instance through the UK obstetric surveillance system (UKOSS) programme to observationally monitor biologic use in hospitals (36). The insights gained here from genetic approaches aim to contribute towards this growing evidence base for clinicians and patients. Additionally, because the investigated drug targets are important mediators of inflammatory responses, these results may also provide broader biological insights into the role of inflammatory signalling in adverse pregnancy outcomes.

## RESEARCH AIMS
Overall, this project aims to characterise the safety profile and potential benefits of using mAbs to manage autoimmune diseases during pregnancy.

First, causal relationships between genetic predisposition to autoimmune diseases and adverse pregnancy outcomes will be characterised in MR.

Second, the causal effects of mAbs on pregnancy outcomes will be evaluated. Evidence on drug safety and efficacy will first be triangulated from across any available existing RCTs including pregnant participants. Drug target MR will then be employed to proxy on-target effects of biologic drugs. Using the MR-PREG collaboration of relevant cohort studies, this will assess the effects of target modulation on a range of adverse pregnancy outcomes.

## METHODS
### OUTCOMES
Pregnancy and perinatal outcomes were selected from those available in the MRPREG collaboration. Seventeen outcomes were selected based on clinical relevance and previously reported observational associations, both with autoimmune diseases (5–10,37) and mAbs (21,23). These were selected by DAL, CMB, CB, KB and LA and finalised on 7/2/24 (see Table 1).


Table 1. Pregnancy and perinatal outcomes.
![](https://github.com/eaiton/autoimmune-drug-pregnancy/blob/main/analysis_plan/table1.png)

All outcomes will be investigated for all autoimmune diseases and drug target exposures, since there is frequent comorbidity and shared molecular pathophysiology among these autoimmune diseases (38), and mAbs are commonly indicated for multiple autoimmune conditions (see Table 3).

### AUTOIMMUNE DISEASES AND PREGNANCY OUTCOMES
#### 1. AUTOIMMUNE DISEASE EXPOSURES
The following conditions of interest were chosen, with selection criteria: (1) autoimmune diseases, (2) which can present in women of reproductive age, and (3) are treated by mAbs, according to NICE clinical guidelines:

* inflammatory bowel disease (IBD; entailing Crohn’s disease and ulcerative colitis)
* rheumatoid arthritis (RA)
* systemic lupus erythematosus (SLE)
* multiple sclerosis (MS)
* psoriasis vulgaris (PV, i.e. plaque psoriasis)
*	psoriatic arthritis (PA)
* axial spondyloarthritis (AS, previously called ankylosing spondylitis)

#### 2. TWO-SAMPLE MENDELIAN RANDOMIZATION 
Two-sample MR will be used to assess causal relationships between the predisposition to each autoimmune diseases and the selected adverse pregnancy outcomes (described above). Though observational studies have established evidence for associations between chronic autoimmune diseases and an increased risk of adverse pregnancy outcomes, they may be subject to residual confounding by factors such as socioeconomic position (1,39) and age. Since MR utilises genetic variation as instruments, these variants are determined at conception and are less likely to relate to socioeconomic and behavioural confounding factors (40), providing a useful basis for comparison against these observational studies.

![](https://github.com/eaiton/autoimmune-drug-pregnancy/blob/main/analysis_plan/fig1.png)
Figure 1. Causal diagram describing causal model underlying the Mendelian Randomization of autoimmune disease exposures on adverse pregnancy outcomes. Produced using Daggity web app (41).

Two-sample MR will be conducting using the TwoSampleMR R package (42), with outcome data obtained from the MR-PREG collaboration, a meta-analysis of pregnancy outcomes across several observational databases (43,44). Since MR-PREG uses predominantly European ancestry cohorts, the largest GWAS in Europeans for each disease will be used to identify genetic variants to proxy the exposure. GWAS were identified in the IEU Open GWAS catalogue or EBI GWAS catalogue (Table 2). Notably, though these GWAS necessarily considered autoimmune diseases as a binary exposure, the genetic instruments are conceptualised as indexing an underlying predisposition to the autoimmune disease (45,46).

Table 2. Autoimmune disease GWAS for Mendelian Randomisation. Largest European GWAS with available summary statistics were selected.
![](https://github.com/eaiton/autoimmune-drug-pregnancy/blob/main/analysis_plan/table2.png)

Inverse variance weighted (IVW) estimates will be presented as the primary analysis. Where instruments are not available in outcome data, LD proxies with r2>0.8 (using a 1000G reference panel) will be used. Where at least 5 SNPs have been selected as genetic instruments, horizontal pleiotropy will be assessed using MR-Egger estimates which allow for an average pleiotropic effect (conceptualised as an intercept term), and weighted median which assumes that only 50% of instruments are valid. Since these pleiotropy robust methods will not be applicable where fewer instruments are available, leave-one-out analyses removing one SNP at a time will be conducted to test for outlying SNPs which may be pleiotropic.

### MONOCLONAL ANTIBODY TREATMENTS FOR AUTOIMMUNE DISEASES AND PREGNANCY OUTCOMES
#### 1.	SELECTION OF DRUGS
All monoclonal antibody drugs indicated for the management of these selected autoimmune diseases listed under Anatomical Therapeutic Chemical (ATC) classification ‘L04 Antineoplastic and Immunomodulating Agents: Immunosuppressants’ were identified as possible exposures for this study. Drug generic names, ATC codes, approval status, targets, and indications were retrieved from DrugBank (54) (https://go.drugbank.com/atc) on 28-29th November 2023. Investigational drugs were included when trials had been conducted or are ongoing for the autoimmune disease indications of interest.

This search identified 25 monoclonal antibody drugs (Table 3), including three currently in trials and/or awaiting approvals (Netakimab, Sirukumab, Opinercept), one withdrawn (Briankinumab) and one currently approved for use by the general population only in Russia (Olokizumab). These have drug target proteins encoded by 12 distinct genes, including several cytokines such as interleukin-6, interleukin-23 and tumour necrosis factor alpha.

Table 3. Selected monoclonal antibody drugs of interest with autoimmune disease indications.
![](https://github.com/eaiton/autoimmune-drug-pregnancy/blob/main/analysis_plan/table3a.png)
![](https://github.com/eaiton/autoimmune-drug-pregnancy/blob/main/analysis_plan/table3b.png)
 
#### 2.	SYSTEMATIC SEARCH FOR RCTS AND PLACENTAL TRANSFER
To identify any RCTs which were conducted to assess the safety or efficacy of the selected mAbs during pregnancy, both the WHO International Clinical Trials Registry Platform (ICTRP) (55) and ClinicalTrails.gov (56) will be searched. Searches will be made for the drug name and the terms “pregnant” or “pregnancy” (e.g. olokizumab AND (pregnant OR pregnancy)). RCTs will be recorded which consider autoimmune disease progression outcomes during pregnancy or pregnancy safety outcomes.

For all identified drugs, current BNF advice on their use in pregnancy will be recorded (57). A systematic literature search will be conducted to summarise any research on placental transfer, searching PubMed (58) for ‘([drug name]) AND (placental transfer)’, recording any relevant results. For all searches, date and search terms will be recorded.

#### 3.	DRUG TARGET MENDELIAN RANDOMIZATION
Drug target MR will be conducted to evaluate whether biologics used to treat autoimmune diseases may cause adverse pregnancy outcomes.

##### 3.A GENETIC INSTRUMENT SELECTION
Genetic variants will be selected for each drug target (see Table 1) using pQTLs recently identified in 54,219 UK Biobank (UKB) participants (59). Instruments selected will be those in cis (defined as those ±1Mb from the gene encoding the protein), to minimise potential for horizontal pleiotropy and since cis variants generally explain more variance in protein levels than trans variants and are more likely to be replicated (59). To ensure independence of instruments, only those in weak LD at the gene loci will be utilised, defined as R<0.02.

If genetic instruments are unavailable in UKB, pQTLs identified by DECODE will be utilised (60). Where instruments for the same drug target are present in both cohorts, results will be replicated for comparison. This may identify epitope effects where protein-coding genetic variants alter the epitope recognised by the assay, thereby biasing the true association with protein levels (61).

Trans-acting instruments will be identified for secondary analyses. Where both cis and trans instruments are available, results for the combined pQTL instruments will be presented. To avoid including trans-acting variants with pleiotropic activity, trans variants will only be included when they have an upstream roles in the signalling pathway of interest. This approach should avoid selecting trans-pQTL that have a very distal relation to the target (for instance, via autoimmune disorders or adiposity), which could be highly pleiotropic and bias our analyses.

Finally, for drug targets comprised of multiple protein subunits, such as Risankizumab which targets both IL-12B and IL-23A subunits comprising IL-23, separate analyses will be conducted for each subunit.


##### 3.B GENETIC INSTRUMENT VALIDATION
The validity of the selected genetic variants will be assessed using proteomic data in the Born in Bradford (BiB) cohort (62). This will provide a general indication of the relevance of these pQTLs as instruments during pregnancy.

The validity of these pQTLs during pregnancy will be assessed using proteomic data measured for BiB participants. Associations between pQTLs and protein concentrations will be tested using conventional multivariable regression in the following model:

	Protein concentration ~ genetic instrument effect alleles + age
	
##### 3.C OBSERVATIONAL ASSOCIATIONS
Conventional multivariable regression will be conducted for comparison against MR estimates. Leveraging within-pregnancy proteomics measured in BiB cohort participants, associations between each drug target protein and pregnancy outcomes will be evaluated using the following model:

	Pregnancy outcome ~ protein concentration during pregnancy + age + parity + offspring sex + best available proxy for socioeconomic position
	
This will be completed for all protein targets which have been measured in BiB.

##### 3.D MENDELIAN RANDOMIZATION
MATERNAL DRUG TARGET EXPOSURE
While high LD between selected instruments would underestimate standard errors and inflate the false positive rate (32), overly stringent pruning thresholds could remove valid instruments and result in an underpowered analysis (63). In order to include several genetic variants at the gene locus while accounting for LD between them, a generalised inverse variance weighted estimator (gIVW) will be used for the primary analysis (64). 
 
 ![](https://github.com/eaiton/autoimmune-drug-pregnancy/blob/main/analysis_plan/fig2.png)
Figure 2. Causal diagram describing causal model underlying the Mendelian Randomization of autoimmune disease exposures on adverse pregnancy outcomes. Produced using Daggity web app (41).

###### FETAL DRUG TARGET EXPOSURE
To explore possible effects of placentally transferred mAbs, a secondary analysis will be conducted to investigate how drug target perturbation in the offspring may causally affect offspring outcomes. For mAbs shown to be transferred across the placenta (see ‘2. Systematic Search for RCTs and placental transfer’), their drug targets will be considered here. This will utilise gIVW estimators with the same genetic instruments identified above, but instrumenting offspring genotypes. Outcomes will be the same as above (see Table 1), excluding maternal outcomes (hypertensive disorders of pregnancy, gestational diabetes, gestational hypertension and preeclampsia) and pregnancy losses (miscarriage and stillbirth).

###### SENSITIVITY ANALYSES
Sensitivity analyses will be conducted for both maternal and fetal drug target exposure analyses.

###### I. CONFOUNDING BY LINKAGE DISEQUILIBRIUM
For each drug target gene locus, colocalization of the genetic associations with drug target protein levels and pregnancy outcomes will be assessed. The colocalization of these associations is necessary for a true causal effect to be present; if instead the variants underlying these associations at the locus are distinct, this indicates that the results may have been confounded by LD.

Colocalization of protein-gene expression data will first be examined to validate the pQTLs. Gene expression data will be obtained for whole blood from eQTLGen (65). Protein-outcome colocalization will then be assessed. All summary statistics at the locus will first be visually inspected LocusZoom (66). Colocalization will then be tested statistically using *coloc*, which assesses the probability of each possible configuration of causal variants using a Bayesian framework, for each cis region ±1Mb from the gene transcription start site. For the prior probabilities, we will assume any SNP is associated with trait 1 (p1 = 1 x 10-4), trait 2 (p2 = 1 x 10-4) or both traits (p12 = 1 x 10-6). A posterior probability for H4 (association with both traits due to a single shared causal variant) of ≥ 70% will be used as a threshold for evidence supporting colocalization, and the same threshold used for H3 (association with both traits due to a distinct causal variants). For tests where the combined posterior probability of association (PPA) for H3 and H4 is lower than 70%, suggesting low power for *coloc*, an additional sensitivity analyses relaxing the prior for both traits (p12 = 1×10−5⁠) will be conducted. This is justified by our candidate gene approach to select specific drug target gene loci, rather than the hypothesis-free genome wide context in which *coloc* was proposed.


###### II. CONFOUNDING BY PLEIOTROPIC EFFECTS
The results of an MR analysis may be biased when genetic instruments affect the outcome via pathways not mediated through the exposure (violation of the exclusion restriction assumption). In this study, this may occur via pleiotropic effects of genetic variants (horizontal pleiotropy) or via direct genetic transmission of instrumented variants from maternal to fetal genotype (fetal genetic effects; see Figure 2).

Horizontal pleiotropy occurs when a variant alters other risk factors associated with the adverse outcome, and may bias effect estimates in either direction. This will be accounted for using MRlink2 (63), an MR method which models pleiotropic variance to separately provide estimates of the causal effect and the horizontal pleiotropic variance. This is implemented through LD clumping variants at varying R2 thresholds, and using an external LD reference panel to model correlation between these selected variants.

To assess possible effects of confounding via direct genetic effects, MR estimates will be replicated using SNP-outcome GWAS adjusted for the fetal/maternal genotype as applicable. For the maternal drug target exposure analysis, this will be maternal SNP-outcome GWAS adjusted for offspring genotype, while for the fetal drug target exposure analysis, this will be fetal SNP-outcome GWAS adjusted for maternal genotype. To derive these conditional GWAS, a weighted linear model (WLM) (67) was implemented in the DONUTS R package (68) using all MR-PREG collaboration studies with both maternal and offspring genotypes available.


## REFERENCES
1.	Conrad N, Misra S, Verbakel JY, Verbeke G, Molenberghs G, Taylor PN, et al. Incidence, prevalence, and co-occurrence of autoimmune disorders over time and by age, sex, and socioeconomic status: a population-based cohort study of 22 million individuals in the UK. The Lancet. 2023 Jun 3;401(10391):1878–90. 
2.	Shivashankar R, Tremaine WJ, Harmsen WS, Loftus EV. Incidence and Prevalence of Crohn’s Disease and Ulcerative Colitis in Olmsted County, Minnesota From 1970 Through 2010. Clinical Gastroenterology and Hepatology. 2017 Jun 1;15(6):857–63. 
3.	Ribbons K, Lea R, Tiedeman C, Mackenzie L, Lechner-Scott J. Ongoing increase in incidence and prevalence of multiple sclerosis in Newcastle, Australia: A 50-year study. Mult Scler. 2017 Jul 1;23(8):1063–71. 
4.	Rees F, Doherty M, Grainge M, Davenport G, Lanyon P, Zhang W. The incidence and prevalence of systemic lupus erythematosus in the UK, 1999–2012. Ann Rheum Dis. 2016 Jan;75(1):136–41. 
5.	O’Toole A, Nwanne O, Tomlinson T. Inflammatory Bowel Disease Increases Risk of Adverse Pregnancy Outcomes: A Meta-Analysis. Dig Dis Sci. 2015 Sep 1;60(9):2750–61. 
6.	Giles I, Yee CS, Gordon C. Stratifying management of rheumatic disease for pregnancy and breastfeeding. Nat Rev Rheumatol. 2019 Jul;15(7):391–402. 
7.	Tian L, Zhang Z, Mao Y, Zong M. Association between pregnant women with rheumatoid arthritis and preeclampsia: A systematic review and meta-analysis. Medicine (Baltimore). 2023 Jun 30;102(26):e34131. 
8.	Fink K, Gorczyca A, Alping P, Englund S, Farmand S, Langer-Gould AM, et al. Multiple sclerosis, disease-modifying drugs and risk for adverse perinatal and pregnancy outcomes: Results from a population-based cohort study. Mult Scler. 2023 May;29(6):731–40. 
9.	Hamroun S, Hamroun A, Bigna JJ, Allado E, Förger F, Molto A. Fertility and pregnancy outcomes in women with spondyloarthritis: a systematic review and meta-analysis. Rheumatology. 2022 Apr 11;61(4):1314–27. 
10.	Bröms G, Haerskjold A, Granath F, Kieler H, Pedersen L, Berglind IA. Effect of Maternal Psoriasis on Pregnancy and Birth Outcomes: A Population-based Cohort Study from Denmark and Sweden. Acta Dermato-Venereologica. 2018 Jun 4;98(8):728–34. 
11.	Bortoli A, Pedersen N, Duricova D, D′Inca R, Gionchetti P, Panelli MR, et al. Pregnancy outcome in inflammatory bowel disease: prospective European case-control ECCO-EpiCom study, 2003–2006. Alimentary Pharmacology & Therapeutics. 2011;34(7):724–34. 
12.	Shand A, Chen J, Selby W, Solomon M, Roberts C. Inflammatory bowel disease in pregnancy: a population-based study of prevalence and pregnancy outcomes. BJOG: An International Journal of Obstetrics & Gynaecology. 2016;123(11):1862–70. 
13.	Skorpen CG, Hoeltzenbein M, Tincani A, Fischer-Betz R, Elefant E, Chambers C, et al. The EULAR points to consider for use of antirheumatic drugs before pregnancy, and during pregnancy and lactation. Annals of the Rheumatic Diseases. 2016 May 1;75(5):795–810. 
14.	Köhler G, Milstein C. Continuous cultures of fused cells secreting antibody of predefined specificity. Nature. 1975 Aug;256(5517):495–7. 
15.	Tsumoto K, Isozaki Y, Yagami H, Tomita M. Future perspectives of therapeutic monoclonal antibodies. Immunotherapy. 2019 Feb;11(2):119–27. 
16.	Lu RM, Hwang YC, Liu IJ, Lee CC, Tsai HZ, Li HJ, et al. Development of therapeutic antibodies for the treatment of diseases. Journal of Biomedical Science. 2020 Jan 2;27(1):1. 
17.	Kaplon H, Crescioli S, Chenoweth A, Visweswaraiah J, Reichert JM. Antibodies to watch in 2023. MAbs. 15(1):2153410. 
18.	Pham-Huy A, Top KA, Constantinescu C, Seow CH, El-Chaâr D. The use and impact of monoclonal antibody biologics during pregnancy. CMAJ. 2021 Jul 26;193(29):E1129–36. 
19.	Nguyen GC, Seow CH, Maxwell C, Huang V, Leung Y, Jones J, et al. The Toronto Consensus Statements for the Management of Inflammatory Bowel Disease in Pregnancy. Gastroenterology. 2016 Mar 1;150(3):734-757.e1. 
20.	Mahadevan U, Wolf DC, Dubinsky M, Cortot A, Lee SD, Siegel CA, et al. Placental Transfer of Anti–Tumor Necrosis Factor Agents in Pregnant Patients With Inflammatory Bowel Disease. Clinical Gastroenterology and Hepatology. 2013 Mar;11(3):286–92. 
21.	Cheent K, Nolan J, Shariq S, Kiho L, Pal A, Arnold J. Case Report: Fatal case of disseminated BCG infection in an infant born to a mother taking infliximab for Crohn’s Disease. Journal of Crohn’s and Colitis. 2010 Nov 1;4(5):603–5. 
22.	Haghikia A, Langer-Gould A, Rellensmann G, Schneider H, Tenenbaum T, Elias-Hamp B, et al. Natalizumab Use During the Third Trimester of Pregnancy. JAMA Neurology. 2014 Jul 1;71(7):891–5. 
23.	Elabd H, Elkholi M, Steinberg L, Acharya A. Eculizumab, a novel potential treatment for acute kidney injury associated with preeclampsia/HELLP syndrome. BMJ Case Rep. 2019 Sep;12(9):e228709. 
24.	Lokki AI, Haapio M, Heikkinen-Eloranta J. Eculizumab Treatment for Postpartum HELLP Syndrome and aHUS—Case Report. Frontiers in Immunology [Internet]. 2020 [cited 2023 Jun 27];11. Available from: https://www.frontiersin.org/articles/10.3389/fimmu.2020.00548
25.	Tong S, Kaitu’u-Lino TJ, Hastie R, Brownfoot F, Cluver C, Hannan N. Pravastatin, proton-pump inhibitors, metformin, micronutrients, and biologics: new horizons for the prevention or treatment of preeclampsia. American Journal of Obstetrics and Gynecology. 2022 Feb 1;226(2, Supplement):S1157–70. 
26.	Heyrana K, Byers HM, Stratton P. Increasing the Participation of Pregnant Women in Clinical Trials. JAMA. 2018 Nov 27;320(20):2077–8. 
27.	Whitelaw S, Sullivan K, Eliya Y, Alruwayeh M, Thabane L, Yancy CW, et al. Trial characteristics associated with under-enrolment of females in randomized controlled trials of heart failure with reduced ejection fraction: a systematic review. European Journal of Heart Failure. 2021;23(1):15–24. 
28.	Assadpour E, Van Spall HGC. Pregnant and lactating women should be included in clinical trials for cardiovascular disease. Nat Med. 2023 Jun 26;1–3. 
29.	Van Spall HGC. Exclusion of pregnant and lactating women from COVID-19 vaccine trials: a missed opportunity. European Heart Journal. 2021 Jul 21;42(28):2724–6. 
30.	Workshop on benefit-risk of medicines used during pregnancy and breastfeeding | European Medicines Agency [Internet]. [cited 2022 Aug 4]. Available from: https://www.ema.europa.eu/en/events/workshop-benefit-risk-medicines-used-during-pregnancy-breastfeeding
31.	Davies NM, Holmes MV, Smith GD. Reading Mendelian randomisation studies: a guide, glossary, and checklist for clinicians. BMJ. 2018 Jul 12;362:k601. 
32.	Gill D, Georgakis MK, Walker VM, Schmidt AF, Gkatzionis A, Freitag DF, et al. Mendelian randomization for studying the effects of perturbing drug targets. Wellcome Open Res. 2021 Feb 10;6:16. 
33.	Interleukin-6 Receptor Mendelian Randomisation Analysis (IL6R MR) Consortium, Swerdlow DI, Holmes MV, Kuchenbaecker KB, Engmann JEL, Shah T, et al. The interleukin-6 receptor as a target for prevention of coronary heart disease: a mendelian randomisation analysis. Lancet. 2012 Mar 31;379(9822):1214–24. 
34.	Ference BA, Kastelein JJP, Ginsberg HN, Chapman MJ, Nicholls SJ, Ray KK, et al. Association of Genetic Variants Related to CETP Inhibitors and Statins With Lipoprotein Levels and Cardiovascular Risk. JAMA. 2017 Sep 12;318(10):947–56. 
35.	Ardissino M, Slob EAW, Rajasundaram S, Reddy RK, Woolf B, Girling J, et al. Safety of beta-blocker and calcium channel blocker antihypertensive drugs in pregnancy: a Mendelian randomization study. BMC Med. 2022 Sep 6;20:288. 
36.	Biological agents in pregnancy | UKOSS | NPEU [Internet]. [cited 2023 Jun 23]. Available from: https://www.npeu.ox.ac.uk/ukoss/current-surveillance/bioagents
37.	Kornfeld D, Cnattingius S, Ekbom A. Pregnancy outcomes in women with inflammatory bowel disease—A population-based cohort study. American Journal of Obstetrics and Gynecology. 1997 Oct 1;177(4):942–6. 
38.	Ellinghaus D, Jostins L, Spain SL, Cortes A, Bethune J, Han B, et al. Analysis of five chronic inflammatory diseases identifies 27 new associations and highlights disease-specific patterns at shared loci. Nat Genet. 2016 May;48(5):510–8. 
39.	Thomson K, Moffat M, Arisa O, Jesurasa A, Richmond C, Odeniyi A, et al. Socioeconomic inequalities and adverse pregnancy outcomes in the UK and Republic of Ireland: a systematic review and meta-analysis. BMJ Open. 2021 Mar 1;11(3):e042753. 
40.	Smith GD, Lawlor DA, Harbord R, Timpson N, Day I, Ebrahim S. Clustered environments and randomized genes: a fundamental distinction between conventional and genetic epidemiology. PLoS Med. 2007 Dec;4(12):e352. 
41.	Textor J, van der Zander B, Gilthorpe MS, Liśkiewicz M, Ellison GT. Robust causal inference using directed acyclic graphs: the R package ‘dagitty’. International Journal of Epidemiology. 2016 Dec 1;45(6):1887–94. 
42.	Hemani G, Zheng J, Elsworth B, Wade KH, Haberland V, Baird D, et al. The MR-Base platform supports systematic causal inference across the human phenome. eLife. 7:e34408. 
43.	Bulik-Sullivan BK, Loh PR, Finucane HK, Ripke S, Yang J, Patterson N, et al. LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. Nat Genet. 2015 Mar;47(3):291–5. 
44.	van Rheenen W, Peyrot WJ, Schork AJ, Lee SH, Wray NR. Genetic correlations of polygenic disease traits: from theory to practice. Nat Rev Genet. 2019 Oct;20(10):567–81. 
45.	Howe LJ, Tudball M, Davey Smith G, Davies NM. Interpreting Mendelian-randomization estimates of the effects of categorical exposures such as disease status and educational attainment. International Journal of Epidemiology. 2022 Jun 1;51(3):948–57. 
46.	Burgess S, Labrecque JA. Mendelian randomization with a binary exposure variable: interpretation and presentation of causal estimates. Eur J Epidemiol. 2018 Oct 1;33(10):947–52. 
47.	de Lange KM, Moutsianas L, Lee JC, Lamb CA, Luo Y, Kennedy NA, et al. Genome-wide association study implicates immune activation of multiple integrin genes in inflammatory bowel disease. Nat Genet. 2017 Feb;49(2):256–61. 
48.	Okada Y, Wu D, Trynka G, Raj T, Terao C, Ikari K, et al. Genetics of rheumatoid arthritis contributes to biology and drug discovery. Nature. 2014 Feb;506(7488):376–81. 
49.	Bentham J, Morris DL, Graham DSC, Pinder CL, Tombleson P, Behrens TW, et al. Genetic association analyses implicate aberrant regulation of innate and adaptive immunity genes in the pathogenesis of systemic lupus erythematosus. Nat Genet. 2015 Dec;47(12):1457–64. 
50.	INTERNATIONAL MULTIPLE SCLEROSIS GENETICS CONSORTIUM. Multiple sclerosis genomic map implicates peripheral immune cells and microglia in susceptibility. Science. 2019 Sep 27;365(6460):eaav7188. 
51.	Sakaue S, Kanai M, Tanigawa Y, Karjalainen J, Kurki M, Koshiba S, et al. A cross-population atlas of genetic associations for 220 human phenotypes. Nat Genet. 2021 Oct;53(10):1415–24. 
52.	Soomro M, Stadler M, Dand N, Bluett J, Jadon D, Jalali-Najafabadi F, et al. Comparative Genetic Analysis of Psoriatic Arthritis and Psoriasis for the Discovery of Genetic Risk Factors and Risk Prediction Modeling. Arthritis Rheumatol. 2022 Sep;74(9):1535–43. 
53.	Cortes A, Hadler J, Pointon JP, Robinson PC, Karaderi T, Leo P, et al. Identification of multiple risk variants for ankylosing spondylitis through high-density genotyping of immune-related loci. Nat Genet. 2013 Jul;45(7):730–8. 
54.	Wishart DS, Feunang YD, Guo AC, Lo EJ, Marcu A, Grant JR, et al. DrugBank 5.0: a major update to the DrugBank database for 2018. Nucleic Acids Res. 2018 Jan 4;46(D1):D1074–82. 
55.	International Clinical Trials Registry Platform (ICTRP) [Internet]. [cited 2023 Jul 5]. Available from: https://www.who.int/clinical-trials-registry-platform
56.	Home | ClinicalTrials.gov [Internet]. [cited 2023 Dec 4]. Available from: https://clinicaltrials.gov/
57.	BNF content published by NICE [Internet]. 2023 [cited 2023 Dec 4]. Available from: https://bnf.nice.org.uk/
58.	PubMed [Internet]. [cited 2023 Dec 4]. PubMed. Available from: https://pubmed.ncbi.nlm.nih.gov/
59.	Sun BB, Chiou J, Traylor M, Benner C, Hsu YH, Richardson TG, et al. Plasma proteomic associations with genetics and health in the UK Biobank. Nature. 2023 Oct;622(7982):329–38. 
60.	Ferkingstad E, Sulem P, Atlason BA, Sveinbjornsson G, Magnusson MI, Styrmisdottir EL, et al. Large-scale integration of the plasma proteome with genetics and disease. Nat Genet. 2021 Dec;53(12):1712–21. 
61.	Suhre K, McCarthy MI, Schwenk JM. Genetics meets proteomics: perspectives for large population-based studies. Nat Rev Genet. 2021 Jan;22(1):19–37. 
62.	Wright J, Small N, Raynor P, Tuffnell D, Bhopal R, Cameron N, et al. Cohort Profile: The Born in Bradford multi-ethnic family cohort study. International Journal of Epidemiology. 2013 Aug 1;42(4):978–91. 
63.	Schmidt AF, Finan C, Gordillo-Marañón M, Asselbergs FW, Freitag DF, Patel RS, et al. Genetic drug target validation using Mendelian randomisation. Nat Commun. 2020 Jun 26;11(1):3255. 
64.	Burgess S, Zuber V, Valdes-Marquez E, Sun BB, Hopewell JC. Mendelian randomization with fine-mapped genetic data: Choosing from large numbers of correlated instrumental variables. Genetic Epidemiology. 2017;41(8):714–25. 
65.	Võsa U, Claringbould A, Westra HJ, Bonder MJ, Deelen P, Zeng B, et al. Large-scale cis- and trans-eQTL analyses identify thousands of genetic loci and polygenic scores that regulate blood gene expression. Nat Genet. 2021 Sep;53(9):1300–10. 
66.	Boughton AP, Welch RP, Flickinger M, VandeHaar P, Taliun D, Abecasis GR, et al. LocusZoom.js: interactive and embeddable visualization of genetic association study results. Bioinformatics. 2021 Sep 29;37(18):3017–8. 
67.	Genome-wide association study of placental weight in 179,025 children and parents reveals distinct and shared genetic influences between placental and fetal growth | medRxiv [Internet]. [cited 2023 Dec 15]. Available from: https://www.medrxiv.org/content/10.1101/2022.11.25.22282723v1
68.	Wu Y, Zhong X, Lin Y, Zhao Z, Chen J, Zheng B, et al. Estimating genetic nurture with summary statistics of multigenerational genome-wide association studies. Proceedings of the National Academy of Sciences. 2021 Jun 22;118(25):e2023184118. 