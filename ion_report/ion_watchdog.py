__author__ = "Kelsey Zhu"
__version__ = "1.0.2"
import time
from watchdog.observers.polling import PollingObserver
from watchdog.events import FileSystemEventHandler
from ion_worker import oncomine_solid

class Watcher:
    DIRECTORY_TO_WATCH = "/mnt/Z_drive/Molecular/IonTorrent/oncosolid_autoreport/dropoff"
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
            ion_worker = oncomine_solid("/home/ionadmin/ion_config.conf")
            ion_worker.workbook = event.src_path
            ion_worker.start()
            del ion_worker

if __name__ == '__main__':
    w = Watcher()
    w.run()

