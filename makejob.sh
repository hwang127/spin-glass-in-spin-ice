#!/bin/csh

#foreach h (0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00)
#foreach h (0.05 0.30 0.45)

foreach nsamples (51200000000)
foreach seed ( ` seq 0  9` )#0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99)
#foreach seed (20 21 22 23 24 25 26 27 28 29)
foreach H(0)#0.001 0.002 0.003 0.005 0.01 0.02 0.03 0.04)# 0.07 0.1 0.15 0.2 0.25 0.30)
foreach nwait(10000)
foreach L(40)
foreach T1(0.5)
foreach T2(0.3)
foreach dj(0.6 1.0 0.2 0.15 0.1 0.05 0.02 0.01 0.0   )
# Input file
cat >/home-2/hwang127@jhu.edu/data/3d_Ising/deltaJ/jobs/in_seed_${seed}_T1_${T1}_T2_${T2}_dj_${dj}_L_${L}.txt << EOFm
nsamples=${nsamples}
nwait=0
L=${L}
Jnn=1.8
Dnn=0
dj=${dj}
T1=${T1}
T2=${T2}
H=${H}
seed=${seed}
outfilename=/home-2/hwang127@jhu.edu/data/3d_Ising/deltaJ/runs/seed_${seed}_T1_${T1}_T2_${T2}_dj_${dj}_L_${L}.txt

EOFm

cat >/home-2/hwang127@jhu.edu/data/3d_Ising/deltaJ/jobs/job_seed_${seed}_T1_${T1}_T2_${T2}_dj_${dj}_L_${L}<<EOFm
#!/bin/bash -l
# Batch Queue Script
#SBATCH --time=12:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

#SBATCH --cpus-per-task=1
#SBATCH --account=olegt
/home-2/hwang127@jhu.edu/data/3d_Ising/deltaJ/mcfile /home-2/hwang127@jhu.edu/data/3d_Ising/deltaJ/jobs/in_seed_${seed}_T1_${T1}_T2_${T2}_dj_${dj}_L_${L}.txt 

EOFm

chmod 755 /home-2/hwang127@jhu.edu/data/3d_Ising/deltaJ/jobs/job_seed_${seed}_T1_${T1}_T2_${T2}_dj_${dj}_L_${L}
sbatch /home-2/hwang127@jhu.edu/data/3d_Ising/deltaJ/jobs/job_seed_${seed}_T1_${T1}_T2_${T2}_dj_${dj}_L_${L}

end
end
end
end
end
end
end
end
