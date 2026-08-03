==User guidle==

For current installation and usage, see 
https://wiki.openwfm.org/wiki/Running_WRF-SFIRE_with_real_data_in_the_WRFx_system

==Smoke test==

conda activate wrfx
./simple_forecast.sh
   press enter for all questions - leave at default except as noted
./forecast.sh jobs/experiment.json
```

Customizations currently needed:
Select 1 node, 64 cpu cores or less.
Select queue system alderaan 
After the slurm job is submitted, add to the job any reservation needed: 
scontrol update job nnnnnnnn reservation=<name>
