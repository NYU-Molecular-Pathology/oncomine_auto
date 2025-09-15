__author__ = "Kelsey Zhu"
__version__ = "1.0.3"
import time
from watchdog.observers.polling import PollingObserver
from watchdog.events import FileSystemEventHandler
from myelo_worker import myeloseq

class Watcher:
    DIRECTORY_TO_WATCH = "/mnt/Z_drive/Molecular/IonTorrent/myeloseqer_test/worksheet.dropoffs/IR3"
    def __init__(self):
        self.observer = PollingObserver()

    def run(self):
        event_handler = Handler()
        self.observer.schedule(event_handler, self.DIRECTORY_TO_WATCH, recursive=True)
        self.observer.start()
        try:
            while True:
                time.sleep(5)
        except:
            self.observer.stop()

        self.observer.join()


class Handler(FileSystemEventHandler):
    @staticmethod
    def on_any_event(event):
        if event.is_directory:
            return None

        elif event.event_type == 'created':
            # Take any action here when a file is first created.
            print("Received created event - %s." % event.src_path)
            time.sleep(5)
            myelo_worker = myeloseq("/Users/yangy15/Documents/oncomine_auto_copy/ion_config.conf") 
            myelo_worker.workbook = event.src_path
            myelo_worker.start()
            del myelo_worker

if __name__ == '__main__':
    w = Watcher()
    w.run()

