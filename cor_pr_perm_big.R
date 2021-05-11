#=============================================================+
# This is the top-level script.                               |
# This script generates pbs scripts that call cor_pr_perm.R   |
#=============================================================+

# Constants
output_folder = "~/scratch/brain/results/cor_pr_perm/"
n_perm = 10000 # 200 matrices total because bhve and ctrl
n_perm_per_job = 500
pbs_file_name = "~/scratch/brain/brain_scripts/cor_pr_perm.pbs"

# Generate pbs script and submit jobs
system(paste0("mkdir -p ", output_folder))
sapply(1:(n_perm/n_perm_per_job), function(x) {
  fileConn <- file(pbs_file_name)
  writeLines(c("#PBS -A GT-js585-biocluster",
               paste0("#PBS -N cor_pr_perm_", x),
               "#PBS -l mem=128gb",
               "#PBS -l nodes=2:ppn=4",
               "#PBS -l walltime=90:00:00",
               "#PBS -j oe",
               paste0("#PBS -o cor_pr_perm_", x, ".out"),
               "#PBS -m abe",
               "#PBS -M gwg@gatech.edu\n",
               "cd $PBS_O_WORKDIR",
               "module purge",
               "module load anaconda3",
               "module load r",
               "conda activate r4\n",
               paste("python grn.py", x, n_perm_per_job)),
             fileConn)
  close(fileConn)
  system(paste("qsub", pbs_file_name))
})

  
# paste("Rscript cor_pr_perm.R", x, n_perm_per_job)),