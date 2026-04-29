#######################################################################
### Quality control and trimming of reads (paired end & single end) ###
#######################################################################

## programs needed:
#export PATH=$PATH:/ngs-data/tools/FastQC_v0.11.4/
#export PATH=$PATH:/ngs-data/tools/Trimmomatic-0.33/
#source ~/.bash_profile


iDATE=$(date +%s)


###################################################################################################
###################################################################################################
# PARAMETERS TO SET MANUALLY

CONTROL_SAMPLES="false"   ## set to "true" or "false"
compound="APA"
tissue="Hepatic"       ## set to "Cardiac" or "Hepatic"
compoundfolder="Acetaminophen"
###############################
### Setting new samplenames ###
###############################

timeTHE=("000" "002" "008" "024" "072" "168" "240" "336")  
timeTOX=("002" "008" "024" "072" "168" "240")
# timeTOX=("000" "002" "008" "024" "072" "168" "240")
dose=("The" "Tox")
replicates=("1" "2" "3")
reads=("R1" "R2")
samplename=()

## Paths to directories
maindir="/ngs-data/data/hecatos/${tissue}/${compoundfolder}/TotalRNA"
inputdir="${maindir}/concatenated"
outputdir="/share/data/hecatos/juantxo/${tissue}/${compoundfolder}/TotalRNA/trimmed"
trimmoversion="0.36"
trimmo="/ngs-data/data/Juantxo_tools/Trimmomatic-${trimmoversion}"
trimmomatic="$trimmo/trimmomatic-${trimmoversion}.jar"
## Trimmomatic parameters
trimlog="false"    ## set to "true" if a trimlog is desired for each trimmomatic run (additional output file)
unpaired="false"    ## set to "true" if unpaired reads from the paired end filtering are desired (additional output file)
singleend='false'
## fastqQC parameters
fastqQC="false"    ## set to "true" or "false"
perlsoft="/ngs-data/data/Juantxo_tools/FastQC/fastqc"

###################################################################################################
outputdirFQR="${outputdir}/fastQC_raw"
outputdirFQT="${outputdir}/fastQC_trimmed"
trimlogdir="${outputdir}/trimlogs_trimmomatic"
unpaireDir="${outputdir}/unpaired_SE_reads"
###################################################################################################
###################################################################################################

echo "ATTENTION: Deleting previous folder with the same name"
# rm -r ${outputdir}
mkdir -p ${outputdir}
mkdir -p ${outputdirFQR}
mkdir -p ${outputdirFQT}
mkdir -p ${trimlogdir}
mkdir -p ${unpaireDir}

cd ${inputdir}



# if ${CONTROL_SAMPLES}; then
# 	for i in ${time[@]}; do
# 		for j in ${replicates[@]}; do
# 			samplename+=("${compound}_${i}_${j}");
# 		done;
# 	done;
# else

unset samplename

for h in ${dose[@]}; do
if [ "$h"  ==  "Tox" ]; then
time=${timeTOX[@]}
else
  time=${timeTHE[@]}
fi
for i in ${time[@]}; do
for j in ${replicates[@]}; do
samplename+=("${compound}_${h}_${i}_${j}");
done;
done;
done;
# fi

#          for g in ${reads[@]}; do
#  				done;

samplenumber=${#samplename[@]}
  echo "The number of samples is ${samplenumber} (${#dose[@]} doses x ${#timeTHE[@]}/${#timeTOX[@]} timepoints x ${#replicates[@]} replicates)"
  
  ########################
  ### perl /ngs-data/data/Juantxo/FastQC/fastqcraw reads ###
  ########################
  
  namepos=0
  ##trimm_reads=("R1_trimmed_PE" "R1_unpaired" "R2_trimmed_PE" "R2_unpaired")
  
  ##for i in ${samplename[@]}; do
  for (( i = 0; i < ${samplenumber} ; i+=1 )); do
  thiSample=${samplename[i]}
  DATE1=$(date +%s)
  if $fastqQC; then
  echo "-----------------------------------------------------"
  echo "Executing FastQC on the raw reads of $thiSample"
  echo "-----------------------------------------------------"
  for g in ${reads[@]}; do
  perl ${perlsoft} ${inputdir}/${thiSample}_${g}.fastq.gz -o $outputdirFQR
  echo "FastQC ${thiSample} completed";
  done
  DATE2=$(date +%s)
  echo "FastQC raw reads completed; time needed: $[ $DATE2 - $DATE1 ]s"
  fi
  echo "-----------------------------------------------------"
  echo "Executing Trimmomatic on ${thiSample}"
  echo "-----------------------------------------------------"
  java -Xms2G -Xmx3G -jar $trimmomatic \
  PE -threads 12 -phred33 -trimlog ${trimlogdir}/trimlog_${thiSample}.txt \
  ${inputdir}/${thiSample}_R1.fastq.gz ${inputdir}/${thiSample}_R2.fastq.gz  \
  ${outputdir}/${thiSample}_R1_trimmed_PE.fastq.gz ${unpaireDir}/${thiSample}_R1_unpaired.fastq.gz \
  ${outputdir}/${thiSample}_R2_trimmed_PE.fastq.gz ${unpaireDir}/${thiSample}_R2_unpaired.fastq.gz \
  ILLUMINACLIP:$trimmo/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:12
  
  if ${singleend}; then
  java -Xms2G -Xmx3G -jar $trimmomatic \
  SE -threads 12 -phred33 -trimlog ${trimlogdir}/trimlog_${thiSample}.txt \
  ${inputdir}/${thiSample}_R1.fastq.gz   \
  ${outputdir}/${thiSample}_R1_trimmed_SE.fastq.gz \
  ILLUMINACLIP:$trimmo/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:12
  
  java -Xms2G -Xmx3G -jar $trimmomatic \
  SE -threads 12 -phred33 -trimlog ${trimlogdir}/trimlog_${thiSample}.txt \
  ${inputdir}/${thiSample}_R2.fastq.gz  \
  ${outputdir}/${thiSample}_R2_trimmed_SE.fastq.gz \
  ILLUMINACLIP:$trimmo/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:12
  fi
  
  ## PE: paired end. -phred33 (there is also a second possibility: phred64). trimlog: log of all read trimmings (read name, surviving seq length)
  ## LEADING: # rm low-qual (below x) bases from the start. ## TRAILING: # rm low-qual (below x) bases from the end.
  ## SLIDINGWINDOW: # rm bases once the avg is below a threshold (x (window):y(avg qual)) ## MINLEN: # rm reads shorter than x. ##HEADCROP: removes x bases from the start of the read.
  DATE3=$(date +%s)
  echo "Trimming of ${thiSample} Time needed: $(($DATE3 - $DATE1))s";
  
  # if $fastqQC; then
  #   echo "-----------------------------------------------------"
  #   echo "Executing FastQC on the trimmed reads of $thiSample"
  #   echo "-----------------------------------------------------"
  #   for l in ${reads[@]}; do
  #     perl ${perlsoft} ${outputdir}/${thiSample}_${l}_trimmed_PE.fastq.gz -o $outputdirFQT
  # 		perl ${perlsoft} ${outputdir}/${thiSample}_${l}_trimmed_SE.fastq.gz -o $outputdirFQT
  # 		perl ${perlsoft} ${unpaireDir}/${thiSample}_${l}_unpaired.fastq.gz -o $outputdirFQT
  #     echo "FastQC trimmed reads ${thiSample} completed";
  #   done
  # 	DATE4=$(date +%s)
  # 	echo "FastQC on the trimmed reads of $thiSample completed; time needed: $[ $DATE4 - $DATE3 ]s";
  # fi
  echo "Testing number ${i} total ${#samplename[@]}"
  namepos=$(($namepos+1))
  # DATEB=$(date +%s)
  # sec=$(( $DATEB - $DATE1))
  # h1=$(($sec/3600))
  # m1=$((($sec-$h1*3600)/60))
  # s1=$(($sec-$h1*3600-$m1*60))
  # timeperfolder=$((($DATEB - $iDATE)/($namepos))) 	   ## final time - initial time == duration --> current duration / current numb of folders completed --> time/folder
  # remfol=$((${#samplename[@]} - $namepos))   ## all folders - folders completed = folders remaining
  # remtime=$(($timeperfolder * $remfol))    ## time/folder * folder remaining == time remaining
  # h2=$(($remtime/3600))
  # m2=$((($remtime-$h2*3600)/60))
  # s2=$(($remtime-$h2*3600-$m2*60))
  # echo "${thiSample} took ${h1} hours, ${m1} minutes and ${s1} seconds   |  $(($namepos*100/${#samplename[@]}))% completed     |  Time remaining: ${h2} hours, ${m2} minutes and ${s2} seconds";
  done;
  
  if ! $trimlog; then
  # rm -r $trimlogdir
  fi
  
  if ! $unpaired; then
  # rm -r $unpaireDir
  fi
  DATE5=$(date +%s)
  
  echo "Quality control and trimming of reads done."
  echo "Total time needed: $((($DATE5 - $iDATE)/3600)) hours"
  