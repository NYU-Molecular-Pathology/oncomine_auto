__author__ = "Kelsey Zhu, Yiying Yang"
__version__ = "1.1"

import sys
import time
from watchdog.observers.polling import PollingObserver
from watchdog.events import FileSystemEventHandler
from ion_worker import oncomine_solid

class Watcher:
    def __init__(self, directory_to_watch, config_path):
        self.DIRECTORY_TO_WATCH = directory_to_watch
        self.config_path = config_path
        self.observer = PollingObserver()

    def run(self):
        event_handler = Handler(self.config_path)
        self.observer.schedule(event_handler, self.DIRECTORY_TO_WATCH, recursive=True)
        self.observer.start()
        try:
            while True:
                time.sleep(5)
        except KeyboardInterrupt:
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
            ion_worker = oncomine_solid(self.config_path)
            ion_worker.workbook = event.src_path
            ion_worker.start()
            del ion_worker

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: ion_watchdog.py <config_path> <directory_to_watch>")
        sys.exit(1)

    config_path = sys.argv[1]
    directory_to_watch = sys.argv[2]

    w = Watcher(directory_to_watch, config_path)
    w.run()

