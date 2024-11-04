#!/bin/bash

TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

logfile=$(grep "LOGFILE=" $1 | cut -d "=" -f 2 | sed "s/\r//g")
if [[ -z "${logfile}" ]]
	then
    	printf "\tThe LOGFILE variable does not exist in the parameter file\n"
		exit 1
else
	if [ ! -f "${logfile}" ]
		then
			printf "\tWelcome to FastGBS!\n" | tee -a "${logfile}"
			printf "\tStarting time: ${TIMESTAMP}\n" | tee -a "${logfile}"
	else
		printf "\nWarning! The file ${logfile} already exist! All subsequent information will be append to it\n" | tee -a "${logfile}"
		printf "Starting time: ${TIMESTAMP}\n" | tee -a "${logfile}"
	fi
fi

printf "\nChecking the parameter file\n"  | tee -a "${logfile}"
if [ "$1" != "" ]; then
    printf "\tThe parameter file is:  ${1}\n"  | tee -a "${logfile}"
else
    printf "\tYou have not specified a parameter file\n\n"  | tee -a "${logfile}"
    exit 1
fi

printf "\nChecking the variables in the parameter file\n"  | tee -a "${logfile}"

flowcell=$(grep "FLOWCELL=" $1 | cut -d "=" -f 2 | sed "s/\r//g")
if [[ -z "${flowcell}" ]]
	then
    	printf "\tThe FLOWCELL variable does not exist in the parameter file\n"  | tee -a "${logfile}"
		exit 1
else
	printf  "\tFLOWCELL : ${flowcell}\n"  | tee -a "${logfile}"
fi

techno=$(grep "TECHNOLOGY=" $1 | cut -d "=" -f 2 | sed "s/\r//g")
if [[ -z "${techno}" ]]
	then
    	printf "\tThe TECHNOLOGY variable does not exist in the parameter file\n"  | tee -a "${logfile}"
		exit 1
else
	printf  "\tTECHNOLOGY : ${techno}\n"  | tee -a "${logfile}"
	if [ "${techno}" = "ILLUMINA" ]
		then
			genIndels= 
			minMapQual=20
			minBaseQual=20
	elif [ "${techno}" = "IONTORRENT" ]
		then
			genIndels=--skip-indels
			minMapQual=10
			minBaseQual=10
	fi
fi
printf "\t\tThe following parameters will be used in the bcftools step:\n"  | tee -a "${logfile}"
printf "\t\tcollapse: ${collapse}\n"  | tee -a "${logfile}"
printf "\t\tminMapQual: ${minMapQual}\n"  | tee -a "${logfile}"
printf "\t\tminBaseQual: ${minBaseQual}\n\n"  | tee -a "${logfile}"

seqtype=$(grep "SEQTYPE=" $1 | cut -d "=" -f 2 | sed "s/\r//g")
if [[ -z "${seqtype}" ]]
	then
    	printf "\tThe SEQTYPE variable does not exist in the parameter file\n"  | tee -a "${logfile}"
		exit 1
else
	printf  "\tSEQTYPE : ${seqtype}\n"  | tee -a "${logfile}"
fi

refgen=$(grep "REFGEN=" $1 | cut -d "=" -f 2 | sed "s/\r//g")
if [[ -z "${refgen}" ]]
	then
    	printf "\tThe REFGEN variable does not exist in the parameter file\n"  | tee -a "${logfile}"
		exit 1
else
	printf  "\tREFGEN : ${refgen}\n"
fi

adapfor=$(grep "ADAPFOR=" $1 | cut -d "=" -f 2 | sed "s/\r//g")
if [[ -z "${adapfor}" ]]
	then
    	printf "\tThe ADAPFOR variable does not exist in the parameter file\n"  | tee -a "${logfile}"
		exit 1
else
	printf  "\tADAPFOR : ${adapfor}\n"  | tee -a "${logfile}"
fi

readlen=$(grep "READLEN=" $1 | cut -d "=" -f 2 | sed "s/\r//g")
if [[ -z "${readlen}" ]]
	then
    	printf "\tThe READLEN variable does not exist in the parameter file\n"  | tee -a "${logfile}"
		exit 1
else
	printf  "\tREADLEN : ${readlen}\n"  | tee -a "${logfile}"
fi

nbcor=$(grep "NBCOR=" $1 | cut -d "=" -f 2 | sed "s/\r//g")
if [[ -z "${nbcor}" ]]
	then
    	printf "\tThe NBCOR variable does not exist in the parameter file\n"  | tee -a "${logfile}"
		exit 1
else
	printf  "\tNBCOR : ${nbcor}\n"  | tee -a "${logfile}"
fi

bwapar=$(grep "BWAPAR=" $1 | cut -d "=" -f 2 | sed "s/\r//g")
if [[ -z "${bwapar}" ]]
	then
    	printf "\tThe BWAPAR variable does not exist in the parameter file\n"  | tee -a "${logfile}"
		exit 1
else
	printf  "\tBWAPAR : ${bwapar}\n"  | tee -a "${logfile}"
fi

bwathr=$(grep "BWATHR=" $1 | cut -d "=" -f 2 | sed "s/\r//g")
if [[ -z "${bwathr}" ]]
	then
    	printf "\tThe BWATHR variable does not exist in the parameter file\n"  | tee -a "${logfile}"
		exit 1
else
	printf  "\tBWATHR : ${bwathr}\n"  | tee -a "${logfile}"
fi

delfiles=$(grep "DELFILES=" $1 | cut -d "=" -f 2 | sed "s/\r//g")
if [[ -z "${delfiles}" ]]
	then
    	printf "\tThe DELFILES variable does not exist in the parameter file\n"  | tee -a "${logfile}"
		exit 1
else
	printf  "\tDELFILES : ${delfiles}\n"  | tee -a "${logfile}"
fi

dellist=$(grep "DELLIST=" $1 | cut -d "=" -f 2 | sed "s/\r//g")
if [[ -z "${dellist}" ]]
	then
    	printf "\tThe DELLIST variable does not exist in the parameter file\n"  | tee -a "${logfile}"
		exit 1
else
	printf  "\tDELLIST : ${dellist}\n"  | tee -a "${logfile}"
fi

mvsamples=$(grep "MVSAMPLES=" $1 | cut -d "=" -f 2 | sed "s/\r//g")
if [[ -z "${mvsamples}" ]]
	then
    	printf "\tThe MVSAMPLES variable does not exist in the parameter file\n"  | tee -a "${logfile}"
		exit 1
else
	printf  "\tMVSAMPLES : ${mvsamples}\n"  | tee -a "${logfile}"
fi

bamlist=$(grep "BAMLIST=" $1 | cut -d "=" -f 2 | sed "s/\r//g")
if [[ -z "${bamlist}" ]]
	then
    	printf "\tThe BAMLIST variable does not exist in the parameter file\n"  | tee -a "${logfile}"
		exit 1
else
	printf  "\tBAMLIST : ${bamlist}\n"  | tee -a "${logfile}"
fi

minreads=$(grep "MINREADS=" $1 | cut -d "=" -f 2 | sed "s/\r//g")
if [[ -z "${minreads}" ]]
	then
    	printf "\tThe MINREADS variable does not exist in the parameter file\n"  | tee -a "${logfile}"
		exit 1
else
	printf  "\tMINREADS : ${minreads}\n"  | tee -a "${logfile}"
fi

logplat=$(grep "LOGPLAT=" $1 | cut -d "=" -f 2 | sed "s/\r//g")
if [[ -z "${logplat}" ]]
	then
    	printf "\tThe LOGPLAT variable does not exist in the parameter file\n"  | tee -a "${logfile}"
		exit 1
else
	printf  "\tLOGPLAT : ${logplat}\n"  | tee -a "${logfile}"
fi

source=$(grep "SOURCE=" $1 | cut -d "=" -f 2 | sed "s/\r//g")
if [[ -z "${source}" ]]
	then
    	printf "\tThe SOURCE variable does not exist in the parameter file\n"  | tee -a "${logfile}"
		exit 1
else
	printf  "\tSOURCE : ${source}\n"  | tee -a "${logfile}"
fi

outplat=$(grep "OUTPLAT=" $1 | cut -d "=" -f 2 | sed "s/\r//g")
if [[ -z "${outplat}" ]]
	then
    	printf "\tThe OUTPLAT variable does not exist in the parameter file\n"  | tee -a "${logfile}"
		exit 1
else
	printf  "\tOUTPLAT : ${outplat}\n"  | tee -a "${logfile}"
fi

maxmis=$(grep "MAX-MISSING=" $1 | cut -d "=" -f 2 | sed "s/\r//g")
if [[ -z "${maxmis}" ]]
	then
    	printf "\tThe MAX-MISSING variable does not exist in the parameter file\n"  | tee -a "${logfile}"
		exit 1
else
	printf  "\tMAX-MISSING : ${maxmis}\n"  | tee -a "${logfile}"
fi

maf=$(grep "MAF=" $1 | cut -d "=" -f 2 | sed "s/\r//g")
if [[ -z "${maf}" ]]
	then
    	printf "\tThe MAF variable does not exist in the parameter file\n"  | tee -a "${logfile}"
		exit 1
else
	printf  "\tMAF : ${maf}\n"  | tee -a "${logfile}"
fi

printf "\nChecking for the presence of the checkpoint file\n" | tee -a "${logfile}"
if [ ! -f "checkpoint_${1}" ]
	then
		printf "\tThe file checkpoint_${1} does not exist, we will create it.\n" | tee -a "${logfile}"
    	touch "checkpoint_${1}"
else
	printf "\tThe file checkpoint_${1} already exist\n"  | tee -a "${logfile}"
fi

if [ "${seqtype}" = "SE" ]
	then
	printf "\nSequence type: Single End"

	printf "\nSE: Production of specific barcodes files\n" | tee -a "${logfile}"
	Step=$(grep "BARCODES" checkpoint_${1})
	if [ "${Step}" != "BARCODES" ]
		then
		printf "\tList of flowcell and their lanes" | tee -a "${logfile}"
		for f in ${flowcell}
		do
			printf "\n\t${f}\n" | tee -a "${logfile}"
			lanes=$(grep "${f}_LANES=" ${1} | cut -d "=" -f 2 | sed "s/\r//g")
			for l in ${lanes}
			do
				printf "\t\t${l}\n" | tee -a "${logfile}"
				./makeBarcodeSabre_V2.py ${f} ${l} 'SE'
				if [ $? -ne 0 ]
					then 
						printf "\t!!! There is a problem in step makeBarcode with ${f} ${l}\n" | tee -a "${logfile}"
						printf "\t!!! Check the names given to the barcode files in the barcode directory"  | tee -a "${logfile}"
						exit 1
				fi
			done
		done
	
		printf "BARCODES\n" >> checkpoint_${1}
	else
		printf  "\tThe variable BARCODES is in the checkpoint file. This step will then be passed\n"| tee -a "${logfile}"
	fi

	printf "\nSE: Demultiplex with sabre\n" | tee -a "${logfile}"
	Step=$(grep "SABRE" checkpoint_${1})
	if [ "${Step}" != "SABRE" ]
		then
			cd data
			parallel -j "${nbcor}" sabre se -f {}.fq.gz -b {}_SE -u {}.unknown.barcodes ::: $(ls -1 *.fq.gz | sed 's/.fq.gz//')
			if [ $? -ne 0 ]
				then 
					printf "\t!!! There is a problem in the demultiplex step with sabre for ${f} ${l}\n" | tee -a ../"${logfile}"
					printf "\t!!! Check the names given to the sequences and barcode files\n" | tee -a ../"${logfile}"
					exit 1
			fi
			cd ..
			printf "SABRE\n" >> checkpoint_${1}
	else
		printf  "\tThe variable SABRE is in the checkpoint file. This step will be passed\n" | tee -a "${logfile}"
	fi

	printf "\nSE: Removing adaptor with cutadapt\n" | tee -a "${logfile}"
	Step=$(grep "CUTADAP" checkpoint_${1})
	if [ "${Step}" != "CUTADAP" ]
		then
			cd data
			parallel -j "${nbcor}" cutadapt -a ${adapfor} -m ${readlen} -o {}.fastq {}.fq ::: $(ls -1 *.fq | sed 's/.fq//')
			if [ $? -ne 0 ]
				then 
					printf "\t!!! There is a problem in the cutadapt step\n" | tee -a ../"${logfile}"
					exit 1
			fi

			cd ..
			printf "CUTADAP\n" >> checkpoint_${1}
	else
		printf  "\tThe variable CUTADAP is in the checkpoint file. This step will be passed\n" | tee -a "${logfile}"
	fi

	printf "\nSE: Moving samples files with less than 10 percent of the average number of read in the pool\n" | tee -a "${logfile}"
	Step=$(grep "MVSAMPLES" checkpoint_${1})
	if [ "${Step}" != "MVSAMPLES" ]
		then
		if [ ${mvsamples} = YES ]
			then
				./count_nbseq_V2.sh			
				if [ $? -ne 0 ]
					then 
						printf "\t!!! There is a problem with the moving samples step\n" | tee -a ../"${logfile}"
						exit 1
				fi
		else
			printf "No file move"
			printf "MVSAMPLES\n" >> checkpoint_${1}
		fi
			printf "MVSAMPLES\n" >> checkpoint_${1}
	else
		printf  "\tThe variable MVSAMPLES is in the checkpoint file. This step will be passed\n" | tee -a "${logfile}"
	fi

	printf "\nSE: Alignment of reads with BWA-MEM\n" | tee -a "${logfile}"
	Step=$(grep "ALIGN" checkpoint_${1})
	if [ "${Step}" != "ALIGN" ]
		then printf "\tStep no. 8: Alignment of single-end reads\n\n" | tee -a "${logfile}"
		cd data
		parallel -j "${bwapar}" bwa mem -t "${bwathr}" ../refgenome/"${refgen}" {}.fastq ">" {}.sam ::: $(ls -1 *.fastq | sed 's/.fastq//')
		if [ $? -ne 0 ]
			then 
				printf "\t!!! There is a problem in the alignment step\n" | tee -a ../"${logfile}"
				exit 1
		fi

		cd ..
		printf "ALIGN\n" >> checkpoint_${1}
	else
		printf  "\tThe variable ALIGN is in the checkpoint file. This step will be passed\n" | tee -a "${logfile}"
	fi

elif [ "${seqtype}" = "PE" ]
	then
	printf "\nSequence type: Paired Ends\n"
	
	printf "\nPE: Production of specific barcodes files\n" | tee -a "${logfile}"
	Step=$(grep "BARCODES" checkpoint_${1})
	if [ "${Step}" != "BARCODES" ]
		then
		printf "\tList of flowcell and their lanes" | tee -a "${logfile}"
		for f in ${flowcell}
		do
			printf "\n\t${f}\n" | tee -a "${logfile}"
			lanes=$(grep "${f}_LANES=" ${1} | cut -d "=" -f 2 | sed "s/\r//g")
			for l in ${lanes}
			do
				printf "\t\t${l}\n" | tee -a "${logfile}"
				./makeBarcodeSabre_V2.py ${f} ${l} 'PE'
				if [ $? -ne 0 ]
					then 
						printf "\t!!! There is a problem in step makeBarcode with ${f} ${l}\n" | tee -a "${logfile}"
						printf "\t!!! Check the names given to the barcode files in the barcode directory"  | tee -a "${logfile}"
						exit 1
				fi
			done
		done
	
		printf "BARCODES\n" >> checkpoint_${1}
	else
		printf  "\tThe variable BARCODES is in the checkpoint file. This step will then be passed\n"| tee -a "${logfile}"
	fi

	printf "\nPE: Demultiplex with sabre\n" | tee -a "${logfile}"
	Step=$(grep "SABRE" checkpoint_${1})
	if [ "${Step}" != "SABRE" ]
		then
			cd data
			parallel -j "${nbcor}" sabre pe -c -f {}_R1.fq.gz -r {}_R2.fq.gz -b {}_PE -u {}_R1.unknown.barcodes -w {}_R2.unknown.barcodes ::: $(ls -1 *_R1.fq.gz | sed 's/_R1.fq.gz//')
			if [ $? -ne 0 ]
				then 
					printf "\t!!! There is a problem in the demultiplex step with sabre for ${f} ${l}\n" | tee -a ../"${logfile}"
					printf "\t!!! Check the names given to the sequences and barcode files\n" | tee -a ../"${logfile}"
					exit 1
			fi
			cd ..
			printf "SABRE\n" >> checkpoint_${1}
	else
		printf  "\tThe variable SABRE is in the checkpoint file. This step will be passed\n" | tee -a "${logfile}"
	fi

	adaprev=$(grep "ADAPREV=" $1 | cut -d "=" -f 2 | sed "s/\r//g")
	if [[ -z "${adaprev}" ]]
		then
		printf "\tThe ADAPREV variable does not exist in the parameter file\n"  | tee -a "${logfile}"
			exit 1
	else
		printf  "\tADAPREV : ${adaprev}\n"  | tee -a "${logfile}"
	fi

	printf "\nPE: Removing adaptor with cutadapt\n" | tee -a "${logfile}"
	Step=$(grep "CUTADAP" checkpoint_${1})
	if [ "${Step}" != "CUTADAP" ]
		then
			cd data
			parallel -j "${nbcor}" cutadapt -a ${adapfor} -A ${adaprev} -m ${readlen} -o {}_R1.fastq -p {}_R2.fastq {}_R1.fq {}_R2.fq ::: $(ls -1 *_R1.fq | sed 's/_R1.fq//')
			if [ $? -ne 0 ]
				then 
					printf "\t!!! There is a problem in the cutadapt step\n" | tee -a ../"${logfile}"
					exit 1
			fi

			cd ..
			printf "CUTADAP\n" >> checkpoint_${1}
	else
		printf  "\tThe variable CUTADAP is in the checkpoint file. This step will be passed\n" | tee -a "${logfile}"
	fi


	printf "\nPE: Moving samples files with less than 10 percent of the average number of read in the pool\n" | tee -a "${logfile}"
	Step=$(grep "MVSAMPLES" checkpoint_${1})
	if [ "${Step}" != "MVSAMPLES" ]
		then
		if [ ${mvsamples} = YES ]
			then
				./count_nbseq_V2.sh			
				if [ $? -ne 0 ]
					then 
						printf "\t!!! There is a problem with the moving samples step\n" | tee -a ../"${logfile}"
						exit 1
				fi
		else
			printf "No file move"
			printf "MVSAMPLES\n" >> checkpoint_${1}
		fi
			printf "MVSAMPLES\n" >> checkpoint_${1}
	else
		printf  "\tThe variable MVSAMPLES is in the checkpoint file. This step will be passed\n" | tee -a "${logfile}"
	fi


	printf "\nPE: Alignment of reads with BWA-MEM\n" | tee -a "${logfile}"
	Step=$(grep "ALIGN" checkpoint_${1})
	if [ "${Step}" != "ALIGN" ]
		then printf "\tStep no. 8: Alignment of paired-end reads\n\n" | tee -a "${logfile}"
		cd data
		parallel -j "${bwapar}" bwa mem -t "${bwathr}" ../refgenome/"${refgen}" {}_R1.fastq {}_R2.fastq ">" {}.sam ::: $(ls -1 *_R1.fastq | sed 's/_R1.fastq//')
		if [ $? -ne 0 ]
			then 
				printf "\t!!! There is a problem in the alignment step\n" | tee -a ../"${logfile}"
				exit 1
		fi

		cd ..
		printf "ALIGN\n" >> checkpoint_${1}
	else
		printf  "\tThe variable ALIGN is in the checkpoint file. This step will be passed\n" | tee -a "${logfile}"
	fi

fi

printf "\nConvert sam files into bam files\n" | tee -a "${logfile}"
Step=$(grep "SAM2BAM" checkpoint_${1})
if [ "${Step}" != "SAM2BAM" ]
	then
		cd data
		parallel -j "${nbcor}" samtools view -b -S -h {}.sam ">" {}.temp.bam ::: $(ls -1 *.sam | sed 's/.sam//')
		if [ $? -ne 0 ]
			then 
				printf "\t!!! There is a problem in the samtools-view step\n" | tee -a ../"${logfile}"
				exit 1
		fi

		cd ..
		printf "SAM2BAM\n" >> checkpoint_${1}
else
	printf  "\tThe variable SAM2BAM is in the checkpoint file. This step will be passed\n" | tee -a "${logfile}"
fi

printf "\nSort of the bam files\n" | tee -a "${logfile}"
Step=$(grep "SORTBAM" checkpoint_${1})
if [ "${Step}" != "SORTBAM" ]
	then
		cd data
		parallel -j "${nbcor}" samtools sort {}.temp.bam -o {}.sort.bam ::: $(ls -1 *.temp.bam | sed 's/.temp.bam//')
		if [ $? -ne 0 ]
			then 
				printf "\t!!! There is a problem in the samtools-sort step\n" | tee -a ../"${logfile}"
				exit 1
		fi
		cd ..
		printf "SORTBAM\n" >> checkpoint_${1}
else
	printf  "\tThe variable SORTBAM is in the checkpoint file. This step will be passed\n" | tee -a "${logfile}"
fi

printf "\nIndex of the bam files\n" | tee -a "${logfile}"
Step=$(grep "INDEXBAM" checkpoint_${1})
if [ "${Step}" != "INDEXBAM" ]
	then
		cd data
		parallel -j "${nbcor}" samtools index {} ::: $(ls -1 *.sort.bam)
		if [ $? -ne 0 ]
			then 
				printf "\t!!! There is a problem in the samtools-index step\n" | tee -a ../"${logfile}"
				exit 1
		fi
		cd ..
		printf "INDEXBAM\n" >> checkpoint_${1}
else
	printf  "\tThe variable INDEXBAM is in the checkpoint file. This step will be passed\n" | tee -a "${logfile}"
fi

printf "\nDeletion of intermediary files to save disk space\n" | tee -a "${logfile}"
Step=$(grep "DELFILES" checkpoint_${1})
if [ "${Step}" != "DELFILES" ]
	then
	if [ ${delfiles} = YES ]
		then
			printf "\tDeletion of intermediary files\n"
			cd data
			rm ${dellist}
			cd ..
			if [ $? -ne 0 ]
				then 
					printf "\t!!! There is a problem with the deletion step\n" | tee -a ../"${logfile}"
					exit 1
			fi
	else
		printf "\tNo deletion of intermediary files\n"
	fi
		printf "DELFILES\n" >> checkpoint_${1}
else
	printf  "\tThe variable DELFILES is in the checkpoint file. This step will be passed\n" | tee -a "${logfile}"
fi

printf "\nProduction of the file containing the list of bam files to be process by bcftools\n" | tee -a "${logfile}"
Step=$(grep "BAMLIST" checkpoint_${1})
if [ "${Step}" != "BAMLIST" ]
	then
		cd data
		if [ ! -f "${bamlist}" ]
			then
				printf "\tThe file ${bamlist} does not exist\n" | tee -a ../"${logfile}"
				touch "${bamlist}"
		else
			printf "\tWarning! The file ${bamlist} already exist! It will be deleted.\n" | tee -a ../"${logfile}"
			rm ${bamlist}
		fi

		for i in $(ls -1 *.sort.bam)
			do
				printf "$PWD/${i}\n" >> "${bamlist}"
				if [ $? -ne 0 ]
					then 
					printf "\t!!! There is a problem in the production of the bam file list\n" | tee -a ../"${logfile}"
					exit 1
				fi
			done
			mv "${bamlist}" ../results
			cd ..
			printf "BAMLIST\n" >> checkpoint_${1}
else
	printf  "\tThe variable BAMLIST is in the checkpoint file. This step will be passed\n" | tee -a "${logfile}"
fi

printf "\nSearch the variants with bcftools\n" | tee -a "${logfile}"
Step=$(grep "BCFTOOLS" checkpoint_${1})
if [ "${Step}" != "BCFTOOLS" ]
        then
                cd results

                if [ "${source}" = "None" ]
                        then sopt=""
                else
                        sopt="--source=${source}"
                fi

                #############################################################################################
                # From /home/vtbos/scripts/fast-gbs-v2B/fast-gbs-v2B.sh
                bcftools mpileup --threads="${nbcor}" \
                -Ou -f ../refgenome/"${refgen}" \
                -min-MQ="${minMapQual}" --min-BQ="${minBaseQual}" \
                "${genIndels}" \
                -b ../results/"${bamlist}" ${sopt} --annotate FORMAT/DP,FORMAT/AD | \

                bcftools call --threads="${nbcor}" -mv -Oz -o raw.vcf

                #bcftools view --threads="${nbcor}" -q 0.002:minor raw.vcf -Oz -o "${outplat}".vcf
                bcftools view --threads="${nbcor}" -i 'FORMAT/DP>=2' raw.vcf | \
                bcftools view --threads="${nbcor}" -q 0.002:minor -Oz -o "${outplat}".vcf

                #############################################################################################

                if [ $? -ne 0 ]
                        then 
                                printf "\t!!! There is a problem at the bcftools step\n" | tee -a ../"${logfile}"
                                exit 1
                fi

        # Extract and clean the #CHROM line from the VCF
    	cat "${outplat}".vcf | awk 'BEGIN{OFS="\t"} /^#CHROM/ {
			for(i=1; i<=NF; i++) {
				if(i > 9) {
					sub(/.*\//, "", $i); 
					sub(/\.sort\.bam$/, "", $i); 
				}
			}
			print
		}' > chrom_line_cleaned.txt

    	# Replace the old #CHROM line in the VCF file with the cleaned one
    	cat "${outplat}".vcf | awk -v new_line="$(cat chrom_line_cleaned.txt)" '{
        	if(/^#CHROM/) {
            	print new_line
        	} else {
            	print
        	}
    	}' > temp_"${outplat}".vcf

    	mv temp_"${outplat}".vcf "${outplat}".vcf
    	rm chrom_line_cleaned.txt
		rm raw.vcf

    	if [ $? -ne 0 ]; then 
        	printf "\t!!! There is a problem at the processing step\n" | tee -a ../"${logfile}"
        	exit 1
    	fi

            cd ..
            printf "BCFTOOLS\n" >> checkpoint_${1}
else
        printf  "\tThe variable BCFTools is in the checkpoint file. This step will be passed\n" | tee -a "${logfile}"
fi

printf "\nImputation of the vcf file\n" | tee -a "${logfile}"
Step=$(grep "IMPUTATION" checkpoint_${1})
if [ "${Step}" != "IMPUTATION" ]
	then
		cd results
			vcftools --vcf "${outplat}".vcf --remove-filtered-all --max-missing ${maxmis} --maf ${maf} --remove-indels --mac 1 --min-alleles 2 --max-alleles 2 --recode --out "${outplat}"
		
			java -Xmx25000m -jar /prg/beagle/5.0/beagle.jar gt="${outplat}".recode.vcf out="${outplat}"_recode_imputed

		if [ $? -ne 0 ]
			then 
				printf "\t!!! There is a problem at the imputation step\n" | tee -a ../"${logfile}"
				exit 1
		fi
	    cd ..
	    printf "IMPUTATION\n" >> checkpoint_${1}

else
	printf  "\tThe variable IMPUTATION is in the checkpoint file. This step will be passed\n" | tee -a "${logfile}"
fi

printf "\nSummary from the vcf file\n" | tee -a "${logfile}"
Step=$(grep "SUMMARY" checkpoint_${1})
if [ "${Step}" != "SUMMARY" ]
	then
		cd results
			vcftools --vcf "${outplat}".recode.vcf --extract-FORMAT-info GT --out "${outplat}".recode
			../Summary4VCF.py "${outplat}".recode.GT.FORMAT

		if [ $? -ne 0 ]
			then 
				printf "\t!!! There is a problem at the summary step\n" | tee -a ../"${logfile}"
				exit 1
		fi
	    cd ..
	    printf "SUMMARY\n" >> checkpoint_${1}

else
	printf  "\tThe variable SUMMARY is in the checkpoint file. This step will be passed\n" | tee -a "${logfile}"
fi
