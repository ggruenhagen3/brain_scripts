#=============================================================+
# This is the top-level script.                               |
# This script generates pbs scripts that call cor_pr_perm.R   |
#=============================================================+

# Constants
output_folder = "~/scratch/brain/results/cor_pr_perm/"
n_perm = 100 # 200 matrices total because bhve and ctrl
n_perm_per_job = 10
pbs_file_name = "~/scratch/brain/brain_scripts/cor_pr_perm.pbs"

# Generate pbs script and submit jobs
system(paste0("mkdir -p ", output_folder))
sapply(1:1, function(x) {
# sapply(1:(n_perm/n_perm_per_job), function(x) {
  fileConn <- file(pbs_file_name)
  writeLines(c("#PBS -A GT-js585",
               paste0("#PBS -N cor_pr_perm_", x),
               "#PBS -l mem=400gb",
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
               paste("Rscript cor_pr_perm.R", x, n_perm_per_job)),
             fileConn)
  close(fileConn)
  system(paste("qsub", pbs_file_name))
})
