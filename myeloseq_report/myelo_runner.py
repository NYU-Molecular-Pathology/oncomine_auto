import sys
from myelo_worker import myeloseq

if len(sys.argv) != 2:
    print("Worksheet was not provided, try again!")
else:
    myelo_runner = myeloseq("/Users/yangy15/Documents/oncomine_automation/ion_config.conf")
    myelo_runner.workbook = sys.argv[1]
    myelo_runner.start()
