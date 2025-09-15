import sys
from myelo_worker import myeloseq

if len(sys.argv) != 3:
    print("Usage: myelo_runner.py <config_path> <worksheet>")
else:
    config_path = sys.argv[1]
    worksheet = sys.argv[2]

    myelo_runner = myeloseq(config_path)
    myelo_runner.workbook = worksheet
    myelo_runner.start()
