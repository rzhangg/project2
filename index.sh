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
    folder=$f;
    # echo $folder;
    for f in $folder/*.faa; do
        prefix=$( basename $f | sed 's/.faa//g' );
        echo 'copying';
        cp $folder/$prefix.ffn $BWA/${folder: -7}.ffn;
        echo 'done';
    done
done

cat $BWA/*.ffn >$BWA/fat_reference.fasta;
echo 'cat done';

bwa index $BWA/fat_reference.fasta;
