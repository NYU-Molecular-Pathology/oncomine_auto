import sys
from dropout_worker import dropout

if len(sys.argv) != 3:
    print("Config or Worksheet was not provided, try again!")
else:
    config_path = sys.argv[1]
    
    dropout_runner = dropout(config_path)
    dropout_runner.workbook = sys.argv[2]
    dropout_runner.start()
