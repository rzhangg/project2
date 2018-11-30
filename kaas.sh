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

    # echo $folder;
    prefix=$( basename $f | sed 's/.faa//g' );
    echo 'copying';
    cp $PROKKA_OUTPUT/$prefix.faa $KAAS/$prefix.faa;
    echo 'done';

done


cat $KAAS/*.faa >$KAAS/fatfile.fasta;
echo 'cat done';
rm -r $KAAS/*.faa;
