# MT
A mixture-model-based method to improve the predictive power of polygenic risk score by analyzing multiple genetically correlated traits

Instruction for running the MT pipeline:

Note: All bash files are submission scripts that is tailored towards UCLA Hoffman2 IDRE system. The user may write scripts on their own to suit for other servers/job submission system.


For the whole genome analysis (without stratification of SNPs into bins)

Step 1: Run PLINK for GWAS (See src/PLINK/)

	1.1 Filter individuals for GWAS analysis (e.g. for training and validation purposes)

	Example script see: src/PLINK/runPlink_filterIndivs.sh

	1.2 Filter SNPs for GWAS analysis (e.g. remove SNPs in regions where the effect size might be inflated)

	Example script see: src/PLINK/runPlink_filterSNPs.sh

	1.3 Run PLINK GWAS using linear or logistic regression for multiple traits

	Example script see: src/PLINK/runPlink_gwas_job_array.sh

	Details for running PLINK please check the following website: https://www.cog-genomics.org/plink/2.0/

	1.4 Generate SNP frequency files

	Example script see: src/PLINK/runPlink_generate_frq.sh

Step 2: Reformat GWAS summary statistics to run MT mixture model. 

	2.1 Add MAF and A2 to GWAS output.

	Example script see helpers: add_A2allele_to_new_plinkout.py and addMAF_to_plinkout.py


	2.2 Reformat to suit MT code.

	Reformat script see MT_parse.py. Use -h to check functional arguments. Note that our method does not support case-control (binary) study right now, please only use the script for study with numeric phenotypes. 
	Example reformatted GWAS summary statistics file see Reformatted_GWAS_example.txt

	CHR	SNP	BP	A1	NMISS	BETA	SE	STAT	P	A2	
	1	1:100000	100000	A	50000	-0.002044	0.001671	-1.223	0.2213	C

Step 3: Run LDSC to prepare genetic covariance and heritability information. (See src/LDSC/)
	
	Details for LDSC functions please check: https://github.com/bulik/ldsc. Install of ldsc environment required.
	For detailed tutorial of this step please check: https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation
	All scripts need to be adjusted for running multiple traits and trait sets.

	3.1 Run LD Score Estimation:

	Please check https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial for detailed information. In short, ldsc is used to estimate LD score from PLINK formatted files .bed/.bim/.fam file sets (added MAF). 

	Script see: run_ldsc_w_maf.sh

	Example output LD Score file: LDSC_file_example.txt

	CHR	SNP	BP	MAF	L2
	1	1:897738	897738	0.063	1.121
	
	3.2 Munge GWAS Data:
	
	Script see: src/LDSC/munge_sumstats.py.
	Example bash script see: src/LDSC/munge.sh

	3.3 Estimate heritability and genetic covariance:

	Script see: src/LDSC/heritability.sh. 

	3.4 Parse heritability and genetic covariance results into one file:

	Script see: parse_ldsc.py. Use -h to check functional arguments for input. # Example command

	Example output format from the above steps see (numbers do not necessarily represent the actual relationship between traits):

		Genetic covariance file: Genetic_covariance_file_example.txt

		trait1	trait2	gcov	gcovSE	meanZ1Z2	gcovIntercept	gcovInterceptSE	gcorr	gcorrSE	gcorrZ	gcorrP
		armfat_percent	trunkfat_percent	-0.0367	0.0166	-0.1835	-0.1087	0.0321	-0.1644	0.0721	-2.2807	0.0226

		Heritability file: Heritability_file_example.txt
	
		trait	nSnpMunged	h2	h2SE	lambda_GC	meanChisq	intercept	interceptSE	ratio
		armfat_percent	61215	0.0648	0.0177	1.0895	1.0984	0.978	0.0319	-9

Step 4: Run MT using GWAS summary statistics. (See src/MT/)
	
	Check README.txt in MT for more detailed description of the scripts and how to run.

	Job submission script see: MT_run_bash.sh

	Example output of MT see: MT_estimates_example.MT_estimates (Estimates of effect sizes and corresponding p values of all SNPs using different methods)

	SNP	CHR	BP	A1	NMISS	BETA	SE	STAT	P	A2	BETA_MTAG	BETA_GMMa	BETA_GMMb	BETA_GMMc	BETA_GMMd	PVAL_MTAG	PVAL_GMMa	PVAL_GMMb	PVAL_GMMc	PVAL_GMMd
	1:100001	1	100001	A	196573	-0.002044	0.001671	-1.223	0.2213	C	-0.002213070214195456	-0.0015691694742631956	-0.0015606696870015893	-0.0015719249659124897	-0.0015719901114600806	0.18537066880295794	0.3476998738934918	0.35031757250082696	0.3468539359648244	0.34683395210015167

	SNP_assignment.Omega_0n.gammas (Assignments of each SNP to different components in GMM methods. GMMa corresponds to .Omega_0a.gammas, and so on)

	snp_id	gamma_0	gamma_1
	1:100001	0.001	0.999

	Weights_of_components.Omega_0s.pi (The overall weight of each component in GMM methods.)
	Omega1_example.Omega_1 (Omega 1 estimated in MT script)

	# trait1	trait2	trait3
	1.359e-05	1.362e-05	1.388e-05
	1.362e-05	1.365e-05	1.390e-05
	1.388e-05	1.390e-05	1.416e-05

Step 5: Run LD Pred on GWAS summary statistics as well as MT estimations to compute polygenic risk score. (See src/LDPred/)

	We were using an older version of LD Pred, v. 1.0.10, committed on Oct 21st, 2019.

The package of this version of LDpred is in the folder named ldpred. For details please check the README.md in ldpred folder.

	Example scripts see run_ldpred_coord.sh, run_ldpred_gibbs.sh and run_ldpred_score.sh

	Example of LD Pred results see LDPred_score_example.txt



For stratification of SNPs into 4 bins (see src/Stratify_SNPs/):

Step 1 - Step 3.1 is the same as previous.

Step 3.1 will generate a log file that will contain a table looking like this:

Summary of LD Scores
         MAF        L2
mean  0.1563    3.2521
std   0.1234    3.9089
min   0.0104    1.0004
25%   0.0568    1.9512
50%   0.1096    2.5946
75%   0.2271    3.5550
max   0.5000  109.8475

Step 3.1.1: Split SNPs into 4 bins according to MAF and L2 (LD Score). Example script see stratify_snps_4bins.py and stratify_snps.sh

	After splitting, the results should have the same format as LDSC_file_example.txt
	CHR	SNP	BP	MAF	L2
	1	1:897738	897738	0.063	1.121

Step 3.1.2: Separate (reformatted) GWAS summary statistics into 4 bins according to SNP assignment from 3.1.1. Example script see MT_gwas_separater_4bins.py and MT_gwas_separater.sh

	After this step, the splitted summary statistics should have the same format as Reformatted_GWAS_example.txt
	CHR	SNP	BP	A1	NMISS	BETA	SE	STAT	P	A2
	3	3:87867846	87867846	T	196573	0.00447	0.001671	2.675	0.007473	C

For the rest of the steps in Step 3, the process should be the same although we need to now munge the GWAS data, estimate heritability and genetic covariance, and parse heritability and genetic covariance results into files for each of the sumstats file generated in step 3.1.2 respectively. 

Step 4: Run MT (only MTAG analysis, no GMMs) on splitted summary statistics (see src/MT_stratified) 

	Check README.txt in MT for more detailed description of the scripts and how to run.

	Job submission script see: MT_run_bash.sh

	Example output of MT see: MT_estimates_example_stratified.MT_estimates. Formatted like the following. (MTAG was run on each bin of summary statistics separately and then the results are combined)

	SNP	CHR	BP	A1	NMISS	BETA	SE	STAT	P	A2	BETA_MTAG	Z_MTAG	P_MTAG
	1:1000001	1	1000001	T	199998	0.0009708	0.001671	0.5809	0.5613	C	0.001294804903927894	0.7748682848162143	0.43841749442142186

Step 5 is the same as previous.



Making Incremental R^2 Plots.

Check src/Incremental_R2_plots.

Step 1: Generate R^2 Dataframes. Use code generate_df.sh which calls generate_r2_df.R. The results are rds files that contains dataframes with the following format:

   Trait Method         R2     Adj_R2
1 height    Cov 0.53373774 0.53370245
2 height   GMMa 0.10168125 0.10154842
3 height   GMMb 0.10280894 0.10267663
4 height   GMMc 0.10042015 0.10028674
5 height   GMMd 0.10042016 0.10028674
6 height   MTAG 0.07018063 0.07003328
7 height   GWAS 0.09890004 0.09876593
8 height  4bins 0.01416871 0.01399552

Step 2: Generate Incremental R2 barplots. Use code generate_barplots.sh and barplots.R. 

Example plot should look like Adjusted_R2_plot.png






