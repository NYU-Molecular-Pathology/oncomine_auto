import sys
from ion_worker import oncomine_solid

if len(sys.argv) != 2:
    print("Worksheet was not provided, try again!")
else:
    ion_worker = oncomine_solid("/Users/yangy15/Documents/oncomine_auto_copy/ion_config.conf")
    ion_worker.workbook = sys.argv[1]
    ion_worker.start()
