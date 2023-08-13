import os
import subprocess

### Define environment
# Path to local neup-ies repo
neup_ies = '/Users/unadavies/Documents/GitHub/neup-ies'
# Path to RAVEN executable
raven_exec = '/Users/unadavies/Documents/Madison_WI/NE2/RAVEN/raven/raven_framework'
# Path to RAVEN-ARMA input file
raven_arma = neup_ies + '/simulations/data/RAVEN_PySAM_coupling/RAVEN_ARMA_coupled.xml'
# Path to RAVEN-PySAM input file
raven_pysam = neup_ies + '/simulations/data/RAVEN_PySAM_coupling/RAVEN_PySAM_coupled.xml'
# Number of perturbations
num_pert = 5

## RAVEN-ARMA
# Generate perturbed datasets
#print('  Running RAVEN-ARMA...')
#subprocess.run([raven_exec, raven_arma])
#print('  Finished RAVEN-ARMA run.')

## Processing perturbed data file names
# convert synData_i.csv outputs to have .0.csv endings for RAVEN
'''
print('  Renaming perturbed data files to allow RAVEN runs...')
for p in range(num_pert):
    mv_comm = ' ARMA/synData_{:d}.csv'.format(p)
    mv_to = ' ARMA/synData_{:d}.0.csv'.format(p)
    os.system('mv' + ' ' + str(mv_comm) + str(mv_to))
print('  Finished renaming data files.')
'''
#
## RAVEN-PySAM
# Run PySAM with synthetic data
print('  Running RAVEN-PySAM with ARMA perturbed files...')
subprocess.run([raven_exec, raven_pysam])
print('  Finished RAVEN-PySAM run.')