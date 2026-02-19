__author__ = "Kelsey Zhu, Yiying Yang"
__version__ = "1.1"

import os
import sys
import time
import threading
from watchdog.observers.polling import PollingObserver
from watchdog.events import FileSystemEventHandler
from myelo_worker import myeloseq
from datetime import datetime

# Keep-alive interval in seconds (1 hour)
DRIVE_KEEPALIVE_INTERVAL = 3600


def hit_drive_recursive(directory):
    """
    Recursively hit/awake the mounted drive by walking the tree (os.listdir per directory).
    Helps keep the network mount responsive and avoids stale connections.
    Returns True if OK, False on error (and logs the error).
    """
    try:
        for _ in os.walk(directory):
            pass  # os.walk() calls os.listdir() on each dir, which hits the drive
        return True
    except (OSError, PermissionError) as e:
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{timestamp}] Keep-alive: error - {e}")
        return False


class Watcher:
    def __init__(self, directory_to_watch, config_path):
        self.DIRECTORY_TO_WATCH = directory_to_watch
        self.config_path = config_path
        self.observer = PollingObserver()

    def _drive_keepalive_loop(self):
        """Background loop: hit the drive recursively every hour to keep mount awake."""
        while True:
            time.sleep(DRIVE_KEEPALIVE_INTERVAL)
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            if hit_drive_recursive(self.DIRECTORY_TO_WATCH):
                print(f"[{timestamp}] Keep-alive: drive OK")
            # on failure, hit_drive_recursive already logs the error

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
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"[{timestamp}] Received created event - {event.src_path}")
            #print(f"Received created event - {event.src_path}")
            time.sleep(5)
            myelo_worker = myeloseq(self.config_path)
            myelo_worker.workbook = event.src_path
            myelo_worker.start()
            del myelo_worker


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: myelo_watchdog.py <config_path> <directory_to_watch>")
        sys.exit(1)

    config_path = sys.argv[1]
    directory_to_watch = sys.argv[2]

    w = Watcher(directory_to_watch, config_path)
    w.run()
