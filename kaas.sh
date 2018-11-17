USER=$1
RESOURCE=/projects/micb405/resources/project_2/2018
PROJECT2=/home/$USER/project2
PROKKA_OUTPUT=/home/$USER/project2/prokka_output
KAAS=$PROJECT2/kaas
cd $PROJECT2
if [ ! -d $KAAS ]; then
    mkdir kaas
fi

for f in $PROKKA_OUTPUT/*; do
    folder=$f;
    # echo $folder;
    for f in $folder/*.faa; do
        prefix=$( basename $f | sed 's/.faa//g' );
        echo 'copying';
        cp $folder/$prefix.faa $KAAS/${folder: -7}.faa;
        echo 'done';
    done
done

cat $KAAS/*.faa >$KAAS/fatfile.fasta;
echo 'cat done';
rm -r $KAAS/*.faa;
