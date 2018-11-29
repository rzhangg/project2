USER=$1
RESOURCE=/projects/micb405/resources/project_2/2018
PROJECT2=/home/$USER/project2
PROKKA_OUTPUT=/home/$USER/project2/prokka_output
MAP=$PROJECT2/map
cd $PROJECT2
mkdir map
for f in $PROKKA_OUTPUT/*; do
    prokka_id=$( head -1 $f | awk -F_ '{ print $1 }' | sed 's/^>//g' )
    mag_id=$( echo $f | sed 's/.faa//g')
    echo $prokka_id,$mag_id
done >Prokka_MAG_map.csv
