#======================================================+
# This is the top-level script.                        |
# This script generates pbs scripts that call grn.py   |
#======================================================+

# Constants
output_folder = "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/results/py_ns_cluster15/"
n_perm = 1000 # 200 matrices total because bhve and ctrl
n_perm_per_job = 500
pbs_file_name = "~/scratch/brain/brain_scripts/py_ns_cluster.pbs"
cluster = 15

# Generate pbs script and submit jobs
system(paste0("mkdir -p ", output_folder))
sapply(1:(n_perm/n_perm_per_job), function(x) {
  for (cluster in 0:(cluster-1)) {
    fileConn <- file(pbs_file_name)
    writeLines(c("#PBS -A GT-js585-biocluster",
                 paste0("#PBS -N py_ns_cluster", as.character(cluster), "_", as.character(x)),
                 "#PBS -l mem=128gb",
                 "#PBS -l nodes=2:ppn=4",
                 "#PBS -l walltime=90:00:00",
                 "#PBS -j oe",
                 paste0("#PBS -o py_ns_cluster", as.character(cluster), "_", as.character(x), ".out"),
                 "#PBS -m abe",
                 "#PBS -M gwg@gatech.edu\n",
                 "cd $PBS_O_WORKDIR",
                 "module purge",
                 "module load anaconda3",
                 "module load r",
                 "conda activate r4\n",
                 paste("python grn.py", x, n_perm_per_job, "-c", cluster, "-o", output_folder)),
               fileConn)
    close(fileConn)
    system(paste("qsub", pbs_file_name))
  }
})
