USER=$1
RESOURCE=/projects/micb405/resources/project_2/2018
PROJECT2=/home/$USER/project2

PROKKA_OUTPUT=/home/$USER/project2/prokka_output
KAAS=$PROJECT2/kaas
BWA=$PROJECT2/bwa
NAME=_10m.qtrim.artifact.rRNA.clean.fastq
SAM=$PROJECT2/sam
REF=$PROJECT2/bwa/*.fa
RPKM=$PROJECT2/rpkm_output
mkdir rpkm_output

if [ ! -d $PROKKA_OUTPUT ]; then
    mkdir prokka_output
fi

for f in $RESOURCE/SaanichInlet_10m/MetaBAT2_SaanichInlet_10m/MedQPlus_MAGs/*fa;
    do sid=$( basename $f | sed 's/.fa//g' );
    echo $sid;
    tax=$(grep -w $sid $RESOURCE/SaanichInlet_10m/MetaBAT2_SaanichInlet_10m/gtdbtk_output/gtdbtk.*.classification_pplacer.tsv | awk '{ print $2 }' | awk -F";" '{ print $1 }' | sed 's/d__//g');
    echo $sid,$tax;
    prokka --kingdom $tax --cpus 40 --outdir $PROKKA_OUTPUT/ --prefix $sid --force $f
done

if [ ! -d $KAAS ]; then
    mkdir kaas
fi

for f in $PROKKA_OUTPUT/*; do

    # echo $folder;
    prefix=$( basename $f | sed 's/.faa//g' );
    echo 'copying';
    cp $PROKKA_OUTPUT/$prefix.faa $KAAS/$prefix.faa;
    echo 'done';

done

cat $KAAS/*.faa >$KAAS/kaas.fasta;
echo 'cat done';
rm -r $KAAS/*.faa;


if [ ! -d $BWA ]; then
    mkdir bwa
fi

for f in $PROKKA_OUTPUT/*; do
    prefix=$( basename $f | sed 's/.ffn//g' );
    echo 'copying';
    cp $PROKKA_OUTPUT/$prefix.ffn $BWA/$prefix.ffn;
    echo 'done';
done

cat $BWA/*.ffn >$BWA/reference.fasta;
echo 'cat done';
rm -r $BWA/*.ffn;
echo 'removed ffn files';
bwa index $BWA/reference.fasta;

if [ ! -d $SAM ]; then
    mkdir sam
fi

for f in $RESOURCE/*$NAME*; do 
    
    prefix=$( basename $f | sed 's/_10m.qtrim.artifact.rRNA.clean.fastq//g' | sed 's/.gz//g');
    echo $prefix;
    bwa mem -t 40 $REF $RESOURCE/$prefix$NAME* > $SAM/$prefix.sam
done

mkdir rpkm_output
for f in $SAM/*.sam; do
    prefix=$( basename $f | sed 's/.sam//g');
    echo $prefix;
    $RESOURCE/rpkm -c $BWA/fat_reference.fasta -a $SAM/$prefix.sam -o $RPKM/$prefix.rpkm.csv
done