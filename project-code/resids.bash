
cd /treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/data/HRS_methy_resid_logs




for i in 0 100000 200000 300000 400000 500000 600000 700000 800000
do qsub -cwd -l vf=50G -b y -N resid4_${i} "Rscript /treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/programs/HRS_methy_resids.R ${i}"
done



