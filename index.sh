USER=$1
RESOURCE=/projects/micb405/resources/project_2/2018
PROJECT2=/home/$USER/project2
PROKKA_OUTPUT=/home/$USER/project2/prokka_output
BWA=$PROJECT2/bwa
cd $PROJECT2
if [ ! -d $BWA ]; then
    mkdir bwa
fi

for f in $PROKKA_OUTPUT/*; do
    prefix=$( basename $f | sed 's/.ffn//g' );
    echo 'copying';
    cp $PROKKA_OUTPUT/$prefix.ffn $BWA/$prefix.ffn;
    echo 'done';
done

cat $BWA/*.ffn >$BWA/fat_reference.fasta;
echo 'cat done';
rm -r $BWA/*.ffn;
echo 'removed ffn files';
bwa index $BWA/fat_reference.fasta;
