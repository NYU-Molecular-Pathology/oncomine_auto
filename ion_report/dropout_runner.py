import sys
from dropout_onco_worker import dropout

if len(sys.argv) != 2:
    print("Worksheet was not provided, try again!")
else:
    dropout_runner = dropout("/Users/yangy15/Documents/oncomine_automation/ion_config.conf")
    dropout_runner.workbook = sys.argv[1]
    dropout_runner.start()