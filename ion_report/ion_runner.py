import sys
from ion_worker import oncomine_solid

if len(sys.argv) != 3:
    print("Usage: ion_runner.py <config_path> <worksheet>")
else:
    config_path = sys.argv[1]
    worksheet = sys.argv[2]

    ion_worker = oncomine_solid(config_path)
    ion_worker.workbook = worksheet
    ion_worker.start()
