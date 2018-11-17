USER=$1
RESOURCE=/projects/micb405/resources/project_2/2018
PROJECT2=/home/$USER/project2
cd $PROJECT2
mkdir prokka_output
PROKKA_OUTPUT=/home/$USER/project2/prokka_output
for f in $RESOURCE/SaanichInlet_10m/MetaBAT2_SaanichInlet_10m/MedQPlus_MAGs/*fa;
    do sid=$( basename $f | sed 's/.fa//g' );
    echo $sid;
    tax=$(grep -w $sid $RESOURCE/SaanichInlet_10m/MetaBAT2_SaanichInlet_10m/gtdbtk_output/gtdbtk.*.classification_pplacer.tsv | awk '{ print $2 }' | awk -F";" '{ print $1 }' | sed 's/d__//g');
    echo $sid,$tax;
    prokka --kingdom $tax --cpus 40 --outdir $PROKKA_OUTPUT/$sid --force $f
done
