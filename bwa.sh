USER=$1
RESOURCE=/projects/micb405/resources/project_2/2018/Metatranscriptomes/
PROJECT2=/home/$USER/project2
NAME=_10m_qtrim.artifact.rRNA.clean.fastq.gz
mkdir sam
SAM=$PROJECT2/sam

 for f in $RESOURCE/Metatrans/*$NAME; do prefix=$( basename $f | sed 's/$NAME//g' );
        echo $prefix;
        bwa mem -t 40 $REF/ref_genome.fasta $RESOURCE/$prefix\$NAME > $SAM/$prefix.sam

done