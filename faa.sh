USER=$1
RESOURCE=/projects/micb405/resources/project_2/2018
PROJECT2=/home/$USER/project2
cd $PROJECT2
mkdir prokka_output
PROKKA_OUTPUT=/home/$USER/project2/prokka_output

for f in $PROJECT2/*; do
    folder=$( basename $f | sed 's///g' );
    echo $folder;
    for f in $PROJECT2/$folder/*.faa; do
        prefix=$( basename $f | sed 's/.faa//g' );
        echo $prefix;
    done
done