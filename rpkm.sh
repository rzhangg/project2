USER=$1
RESOURCE=/projects/micb405/resources/project_2/2018
PROJECT2=/home/$USER/project2
cd $PROJECT2
SAM=$PROJECT2/sam
RPKM=$PROJECT2/rpkm_output
BWA=$PROJECT2/bwa
mkdir rpkm_output
for f in $SAM/*.sam; do
    prefix=$( basename $f | sed 's/.sam//g');
    echo $prefix;
    $RESOURCE/rpkm -c $BWA/fat_reference.fasta -a $SAM/$prefix.sam -o $RPKM/$prefix.rpkm.csv
done
