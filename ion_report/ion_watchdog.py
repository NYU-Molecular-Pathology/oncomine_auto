__author__ = "Kelsey Zhu, Yiying Yang"
__version__ = "1.1"

import os
import sys
import time
import threading
from datetime import datetime
from watchdog.observers.polling import PollingObserver
from watchdog.events import FileSystemEventHandler
from ion_worker import oncomine_solid

# Keep-alive interval (seconds) – hits drive every hour
DRIVE_KEEPALIVE_INTERVAL = 3600


def hit_drive_recursive(directory):
    """
    Recursively hit/awake the mounted drive (os.walk → os.listdir per dir).
    Keeps network mount responsive; returns True if OK, False on error.
    """
    try:
        for _ in os.walk(directory):
            pass
        return True
    except (OSError, PermissionError) as e:
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{ts}] Drive hit FAIL: {e}")
        return False


class Watcher:
    def __init__(self, directory_to_watch, config_path):
        self.DIRECTORY_TO_WATCH = directory_to_watch
        self.config_path = config_path
        self.observer = PollingObserver()

    def _drive_keepalive_loop(self):
        """Hit drive recursively every hour to keep mount awake."""
        while True:
            time.sleep(DRIVE_KEEPALIVE_INTERVAL)
            if hit_drive_recursive(self.DIRECTORY_TO_WATCH):
                ts = datetime.now().strftime("%H:%M")
                print(f"[{ts}] Drive OK")
            # on failure, hit_drive_recursive already logs

    def run(self):
        event_handler = Handler(self.config_path)
        self.observer.schedule(event_handler, self.DIRECTORY_TO_WATCH, recursive=True)
        self.observer.start()

        keepalive = threading.Thread(target=self._drive_keepalive_loop, daemon=True)
        keepalive.start()

        try:
            while True:
                time.sleep(5)
        except KeyboardInterrupt:
            self.observer.stop()
        self.observer.join()


class Handler(FileSystemEventHandler):
    def __init__(self, config_path):
        self.config_path = config_path

    def on_any_event(self, event):
        if event.is_directory:
            return

        if event.event_type == "created":
            print(f"Received created event - {event.src_path}")
            time.sleep(5)
            oncomine_worker = oncomine_solid(self.config_path)
            oncomine_worker.workbook = event.src_path
            oncomine_worker.start()
            del oncomine_worker


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: ion_watchdog.py <config_path> <directory_to_watch>")
        sys.exit(1)

    config_path = sys.argv[1]
    directory_to_watch = sys.argv[2]

    w = Watcher(directory_to_watch, config_path)
    w.run()
