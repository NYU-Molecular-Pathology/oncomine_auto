## Oncomine/Myeloseq Automation  
### Process Overview  

This workflow performs the following steps in sequence:  

1. **Detect samplesheet drop-off** on the research network drive (performed by the lab).  
2. **Parse the samplesheet** to extract the run ID and sample IDs.  
3. **Download VCF and TSV files** from Ion Reporter via API calls.  
4. **Filter raw VCFs** to generate outputs identical to Ion Reporter (CSV/XLSX).  
5. **Generate QC plots** in PNG and PDF formats.  
6. **Calculate dropouts** and produce a `dropout.html` report.  
7. **Download BAM files** from Ion Reporter.  


## Instructions to Run  
### Watchdog process  
This process is currently on our production VM (@MOLAPPLPDCPVM01). Scripts located under `~/oncomine_automation/`. At the moment, all processes are protected by their `systemctl` services. You can find the details of `systemctl` under `~/.config/systemd/user/`   
```
# log in to production VM
# check if processes are up do
ps axu | grep "myeloseq"
# or
ps axu | grep "oncomine"
# to confirm processes are running for example myeloseqer
ps axu | grep "myeloseq"
# There will be two processes for each assay because one is for IR2 runs, one is for IR3 runs
#yangy15+  384241  0.1  0.1 2135648 161020 ?      Sl   Dec08   5:37 python -u /home/yangy15_adm/oncomine_auto/myeloseq_report/myelo_watchdog.py /home/yangy15_adm/oncomine_auto/ion_config_IR2.conf /mnt/Z_drive/Molecular/IonTorrent/myeloseqer/IR2.dropoff
#yangy15+  384242  0.0  0.0 2051648 79528 ?       Sl   Dec08   1:37 python -u /home/yangy15_adm/oncomine_auto/myeloseq_report/myelo_watchdog.py /home/yangy15_adm/oncomine_auto/ion_config_IR3.conf /mnt/Z_drive/Molecular/IonTorrent/myeloseqer/IR3.dropoff

```
### Manual rerun  
To run a run manually, depending on the run is processed by IR2 or IR3  
__Oncomine__  
```
# log in to the production VM
### Make sure you have conda installed in your home ###
cd /apps/bi_shared/bi_pipelines/scripts/oncomine_automation/ion_report
# load oncomine environment
conda activate /apps/bi_shared/bi_pipelines/envs/ion_env
# if the run is on IR3
python3 ion_runner.py /apps/bi_shared/bi_pipelines/scripts/oncomine_auto/ion_config_IR3.conf /mnt/Z_drive/Molecular/IonTorrent/oncomine/IR3.dropoff/<runid>.xlsm
# if the run is on IR2
python3 ion_runner.py /apps/bi_shared/bi_pipelines/scripts/oncomine_auto/ion_config_IR2.conf /mnt/Z_drive/Molecular/IonTorrent/oncomine/IR2.dropoff/<runid>.xlsm
```
__Myeloseqer__  
```
# log in to the production VM
### Make sure you have conda installed in your home ###
cd /apps/bi_shared/bi_pipelines/scripts/oncomine_automation/myeloseq_report
# load myeloseqer environment
conda activate /apps/bi_shared/bi_pipelines/envs/myeloseq_env
# if the run is on IR3
python3 myelo_runner.py /apps/bi_shared/bi_pipelines/scripts/oncomine_auto/ion_config_IR3.conf /mnt/Z_drive/Molecular/IonTorrent/myeloseqer/IR3.dropoff/<runid>.xlsm
# if the run is on IR2
python3 myelo_runner.py /apps/bi_shared/bi_pipelines/scripts/oncomine_auto/ion_config_IR2.conf /mnt/Z_drive/Molecular/IonTorrent/myeloseqer/IR2.dropoff/<runid>.xlsm
```
