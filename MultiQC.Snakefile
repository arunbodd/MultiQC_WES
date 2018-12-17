import os
from os.path import join
dirName = "BATCH5_QC"

if not os.path.exists(dirName):
    os.mkdir(dirName)
    print("Directory " , dirName ,  " Created ")
else:
   print("Directory " , dirName ,  " already exists")

configfile:"/hpcdata/dir/CIDR_DATA_QC_RESULTS/multiqc_run.json"
BAM_DIR = "/hpcdata/dir/CIDR_DATA_RENAMED/BATCH5/BAM/"
RES_FOLDER = "/hpcdata/dir/CIDR_DATA_QC_RESULTS/BATCH5_QC/"
#BAM_SUFFIX = ".bam"
BAM_ID, = glob_wildcards(BAM_DIR + "{ID}" + ".bam")
#CUMULATIVE_REPORT_DIR = "/hpcdata/dir/CIDR_DATA_QC_RESULTS/"


VCF_FILE = "/hpcdata/dir/CIDR_DATA_RENAMED/BATCH5/VCF_rerun2_111618/BATCH5_RERUN2.vcf"



rule all:
    input:fastqc=expand(join(RES_FOLDER,"FastQC/{Sample}"),Sample=BAM_ID),
    	  flagstats=expand(join(RES_FOLDER,"Flagstats/{Sample}.flagstats"),Sample=BAM_ID),
    	  qualimap=expand(join(RES_FOLDER, "{Sample}", "qualimapReport.html"),Sample = BAM_ID),
    	  selectvariants=expand(join(RES_FOLDER,"{Sample}.vcf.gz"),Sample=BAM_ID),
    	  varianteval=expand(join(RES_FOLDER,"VariantEval/{Sample}"),Sample = BAM_ID),
    	  snpeff= expand(join(RES_FOLDER,"SNPeff/{Sample}/{Sample}"),Sample = BAM_ID),
 		  bcftools=expand(join(RES_FOLDER,"BCFStats/{Sample}"),Sample = BAM_ID),
 		  multipqc=join(RES_FOLDER,"Batch5QC_Report.html")

rule fastqc_fastaq:
	input :join(BAM_DIR+ "{Sample}.bam")
	output : join(RES_FOLDER,"FastQC/{Sample}")
        params: adapters=config['references']['fastqc_adapters']
        threads: 10
        shell: "module load fastqc;fastqc -o BATCH5_QC/FastQC -f fastq --threads {threads} -f bam --contaminants {params.adapters} {input}"

rule qualimap:
 	input: join(BAM_DIR+ "{Sample}.bam")
   	output: txt = join(RES_FOLDER,"{Sample}","genome_results.txt"), html = join(RES_FOLDER, "{Sample}", "qualimapReport.html")
	threads:8
   	params:regions=config['references']['REGIONS'], dir = join(RES_FOLDER, "{Sample}")
   	shell: "module load qualimap;unset DISPLAY; qualimap bamqc -bam {input} --java-mem-size=48G -c gd hg19 -ip -outdir {params.dir} -gff {params.regions} -outformat HTML -nt {threads} --skip-duplicated -nw 500 -p NON-STRAND-SPECIFIC"

rule samtools_flagstats:
  	input:bam= join(BAM_DIR+ "{Sample}.bam")
  	output:join(RES_FOLDER,"Flagstats/{Sample}.flagstats")
  	shell: "module load samtools; samtools flagstat {input} > {output}"
  	

rule Gatk_SelectVariants:
        input:selectvariants= VCF_FILE
        output:join(RES_FOLDER,"{Sample}.vcf.gz")          
        params:genome=config['references']['GENOME'], Sname = "{Sample}"
        shell: "module load GATK/3.7-0-Java-1.8.0_92;java -Xmx48g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T SelectVariants -R {params.genome} -o {output} -V {input} --sample_name {params.Sname} --excludeNonVariants"

rule bcftools:
		input:join(RES_FOLDER,"{Sample}.vcf.gz")
  		output:join(RES_FOLDER,"BCFStats/{Sample}")
  		shell: "module load bcftools/1.4.1-goolf-1.7.20; bcftools stats {input} > {output}"
  	
rule varianteval:
  		input:vcf = join(RES_FOLDER,"{Sample}.vcf.gz")
  		output:join(RES_FOLDER,"VariantEval/{Sample}")
   		params:genome=config['references']['GENOME'],vcf=config['references']['DBSNP']	
   		shell:"module load GATK/3.7-0-Java-1.8.0_92;java -Xmx48g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantEval -R {params.genome} -o {output} --dbsnp {params.vcf} --eval {input.vcf}"
   		
rule snpeff:
  		input: join(RES_FOLDER,"{Sample}.vcf.gz")
  		output:vcf=join(RES_FOLDER,"SNPeff/{Sample}/{Sample}_exome.vcf"),
	           csv = join(RES_FOLDER,"SNPeff/{Sample}/{Sample}.stats.csv"),
	           html = join(RES_FOLDER,"SNPeff/{Sample}/{Sample}")
  		params:genome=config['references']['SNPEFF_GENOME'],effconfig=config['references']['SNPEFF_CONFIG']
  		shell: "module load java/1.8.0_92; java -Xmx12g -jar /hpcdata/scratch/lackjb/snpEff/snpEff.jar -v -canon -c {params.effconfig} -csvstats {output.csv} -stats {output.html} {params.genome} {input} > {output.vcf}"
  	
   		
rule multiqc:
		input:expand(join(RES_FOLDER,"FastQC/{Sample}"),Sample=BAM_ID), 
		      expand(join(RES_FOLDER,"Flagstats/{Sample}.flagstats"),Sample=BAM_ID),
		 	  expand(join(RES_FOLDER, "{Sample}", "qualimapReport.html"),Sample = BAM_ID),
			  expand(join(RES_FOLDER,"VariantEval/{Sample}"),Sample = BAM_ID), 
		      expand(join(RES_FOLDER,"SNPeff/{Sample}/{Sample}"),Sample = BAM_ID), 
 		      expand(join(RES_FOLDER,"BCFStats/{Sample}"),Sample = BAM_ID),
		output:out1 = join(RES_FOLDER,"Batch5QC_Report.html")
		shell:"module load multiqc;multiqc -f -n {output.out1} BATCH5_QC/;rm BATCH5_QC/*.vcf.gz*"
