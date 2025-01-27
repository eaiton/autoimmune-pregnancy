# Investigating Causal Effects of Autoimmune Conditions on Pregnancy Outcomes: a Mendelian Randomization study protocol
	
## BACKGROUND
### AUTOIMMUNE DISEASES AND PREGNANCY OUTCOMES
Chronic autoimmune conditions affect around 13% of women in the UK (1), encompassing conditions such as rheumatic diseases, inflammatory bowel disease, and multiple sclerosis. More common amongst female patients, these conditions are often commonly diagnosed before menopause; in studies of women in high income countries, incidence of IBD is highest amongst those aged 20-29 (2), while women aged 40-49 have the highest incidence of MS (3) and systemic lupus erythematous (4).
	
Observational evidence suggests an underlying autoimmune condition is associated with an increased odds of adverse pregnancy outcomes. Meta-analyses of observational studies have found that active inflammatory bowel disease is associated with an increased risk of pre-term birth, small-for-gestational-age birthweight, stillbirth, and congenital anomalies (5). Rheumatic diseases, encompassing both rheumatic arthritis and systemic lupus erythematosus, have been associated with increased risk of preterm birth, intrauterine growth restriction (6), and pre-eclampsia (7). Adverse outcomes have similarly been associated with multiple sclerosis (8) and axial spondyloarthritis (9). Though there is comparatively less evidence on psoriasis, a recent population-wide study of Sweden and Denmark found an association between maternal psoriasis and increased risks of gestational diabetes and hypertensive disorders of pregnancy (10). There is a need to provide comprehensive evidence since rarer conditions and outcomes are relatively under-studied, and existing observational studies may be subject to residual confounding such as by medication use or socioeconomic position.
	
Mendelian randomization (MR) is a causal inference method which leverages genetic data. MR is commonly implemented in instrumental variable analyses which treat genetic variants associated with an exposure of interest as randomly inherited instrumental variables, ‘genetic instruments’, to index an exposures’ causal effect on an outcome (11). Three core instrumental variable assumptions are required for an MR analysis: (1) genetic variants are associated with the exposure (relevance), (2) there are no confounders of the genetic variant-outcome associations (independence), and (3) genetic variants only act on the outcome through their effect on the exposure of interest (the exclusion restriction).
	
## RESEARCH AIMS
This study aims to characterise the causal effects of autoimmune diseases on adverse pregnancy and perinatal outcomes will be ascertained.

## METHODS
### OUTCOMES
Pregnancy and perinatal outcomes were selected from those available in the MRPREG collaboration. Seventeen outcomes were selected based on clinical relevance and previously reported observational associations with autoimmune diseases (5-10). These were selected by DAL, CMB, CB, KB and LA and finalised on 7/2/24 (see Table 1).
	
Table 1. Pregnancy and perinatal outcomes.
![](https://github.com/eaiton/autoimmune-pregnancy/blob/main/analysis_plan/table1.png)
	
### AUTOIMMUNE DISEASES EXPOSURES
We selected autoimmune conditions which can present in women of reproductive age (15-44 years) (1), for which a publicly available GWAS with at least 5,000 cases had been published:

* axial spondyloarthritis (AS)
* celiac disease (CD)
* Hashimoto's thyroiditis (HT)
* inflammatory bowel disease (IBD)
* multiple sclerosis (MS)
* psoriasis (Ps)
* rheumatoid arthritis (RA)
* systemic lupus erythematosus (SLE)
* systemic sclerosis
* type 1 diabetes (T1D)
	
### TWO-SAMPLE MENDELIAN RANDOMIZATION 
Two-sample MR will be conducting using the TwoSampleMR R package (12), with outcome data obtained from the MR-PREG collaboration, a meta-analysis of pregnancy outcomes across several observational databases. GWAS were identified in the IEU Open GWAS catalogue or EBI GWAS catalogue, and the largest GWAS for each disease will be used to identify genetic instruments.
	
Inverse variance weighted (IVW) estimates will be presented as the primary analysis. Where instruments are not available in outcome data, LD proxies with r2>0.8 (using a 1000G reference panel) will be used. Where at least 5 SNPs have been selected as genetic instruments, horizontal pleiotropy will be assessed using MR-Egger estimates which allow for an average pleiotropic effect (conceptualised as an intercept term), and weighted median which assumes that only 50% of instruments are valid. Since these pleiotropy robust methods will not be applicable where fewer instruments are available, leave-one-out analyses removing one SNP at a time will be conducted to test for outlying SNPs which may be pleiotropic. Instruments within the MHC region will also be excluded in a sensitivity analysis, since this region is highly pleiotropic.
	
Finally, MR estimates will be replicated using a maternal SNP-outcome GWAS adjusted for fetal genotype, to assess robustness of effects to confounding by direct genetic transmission. To derive these conditional outcome GWAS, a weighted linear model (WLM) (13) was implemented in the DONUTS R package (14) using all MR-PREG collaboration studies with both maternal and offspring genotypes available.
	
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
11.	Davies NM, Holmes MV, Smith GD. Reading Mendelian randomisation studies: a guide, glossary, and checklist for clinicians. BMJ. 2018 Jul 12;362:k601.
12.	Hemani G, Zheng J, Elsworth B, Wade KH, Haberland V, Baird D, et al. The MR-Base platform supports systematic causal inference across the human phenome. eLife. 7:e34408. 
13.	Genome-wide association study of placental weight in 179,025 children and parents reveals distinct and shared genetic influences between placental and fetal growth | medRxiv [Internet]. [cited 2023 Dec 15]. Available from: https://www.medrxiv.org/content/10.1101/2022.11.25.22282723v1
14.	Wu Y, Zhong X, Lin Y, Zhao Z, Chen J, Zheng B, et al. Estimating genetic nurture with summary statistics of multigenerational genome-wide association studies. Proceedings of the National Academy of Sciences. 2021 Jun 22;118(25):e2023184118. 
