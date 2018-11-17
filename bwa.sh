USER=$1
RESOURCE=/projects/micb405/resources/project_2/2018/Metatranscriptomes
PROJECT2=/home/$USER/project2
NAME=_10m.qtrim.artifact.rRNA.clean.fastq

SAM=$PROJECT2/sam
REF=$PROJECT2/bwa/*.fasta

if [ ! -d $SAM ]; then
    mkdir sam
fi

for f in $RESOURCE/*$NAME*; do 
    
    prefix=$( basename $f | sed 's/_10m.qtrim.artifact.rRNA.clean.fastq//g' | sed 's/.gz//g');
    echo $prefix;
    bwa mem -t 40 $REF $RESOURCE/$prefix$NAME* > $SAM/$prefix.sam
done

