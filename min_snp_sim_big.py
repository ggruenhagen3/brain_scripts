import os

for depth in [0.001, 0.01, 0.1]:
    print(depth)
    pbs_script = ["#PBS -A GT-js585",
                 "#PBS -N min_snp_sim_" + str(depth),
                 "#PBS -l mem=64gb",
                 "#PBS -l nodes=2:ppn=4",
                 "#PBS -l walltime=10:00:00",
                 "#PBS -j oe",
                 "#PBS -o min_snp_sim_" + str(depth) + ".out",
                 "#PBS -m abe",
                 "#PBS -M gwg@gatech.edu\n",
                 "cd $PBS_O_WORKDIR",
                 "module purge",
                 "module load anaconda3",
                 "module load r",
                 "conda activate scSplit\n",
                 "python min_snp_sim.py " + str(1) + " " + str(3) + " " + str(depth)]

    f = open("min_snp_sim_depth.pbs", "w")
    f.write('\n'.join(pbs_script))
    f.close()
    os.system("qsub min_snp_sim_depth.pbs")
    print("done")
print("All Done")
