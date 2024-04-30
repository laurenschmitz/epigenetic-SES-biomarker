
cd /treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/results/HRS\ Enet/EWAS/logs



for l in 2
do
for k in ALL #M1 M2 M3 F1 F2 F3 ALL 1 2 3 M F
do
for j in rameduc2 rafeduc2 #COLLEDUC EDYEARS HHINC HSEDUC NDS PEDUC POV lnhhincx lnwealth2x
do
for i in 0 100000 200000 300000 400000 500000 600000 700000 800000
do qsub -cwd -l vf=50G -b y -N ${j}_${k}_${i} "Rscript /treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/programs/HRS_socialepi_EWAS_row_col_plate.R ${i} ${j} ${k} ${l}"
done
done
done
done


