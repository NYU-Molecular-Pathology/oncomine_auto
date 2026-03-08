__author__ = "Kelsey Zhu, Yiying Yang"
__version__ = "1.1"

import os
import sys
import time
import threading
import smtplib
import traceback
import configparser
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from watchdog.observers.polling import PollingObserver
from watchdog.events import FileSystemEventHandler
from myelo_worker import myeloseq
from datetime import datetime

# Keep-alive interval in seconds (1 hour)
DRIVE_KEEPALIVE_INTERVAL = 3600

def send_email_notification(subject, body, smtp_server, smtp_port, sender_email, receiver_email):
    """Send a plain-text email notification through institutional SMTP relay."""
    if not smtp_server or not sender_email or not receiver_email:
        print(
            f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] "
            f"Email skipped for '{subject}': SMTP settings not configured."
        )
        return

    receiver_list = [email.strip() for email in receiver_email.split(",") if email.strip()]
    if not receiver_list:
        print(
            f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] "
            f"Email skipped for '{subject}': no valid receiver configured."
        )
        return

    message = MIMEMultipart()
    message["From"] = sender_email
    message["To"] = ", ".join(receiver_list)
    message["Subject"] = subject
    message.attach(MIMEText(body, "plain"))

    server = None
    try:
        server = smtplib.SMTP(smtp_server, smtp_port)
        server.sendmail(sender_email, receiver_list, message.as_string())
        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Email sent: {subject}")
    except Exception as e:
        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Error sending email: {e}")
    finally:
        if server is not None:
            server.quit()


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
        config = configparser.ConfigParser()
        config.read(config_path)
        self.smtp_server = config["DEFAULT"].get("SMTP_SERVER")
        self.smtp_port = config["DEFAULT"].getint("SMTP_PORT", 25)
        self.sender_email = config["DEFAULT"].get("SENDER_EMAIL")
        self.receiver_email = config["DEFAULT"].get("RECEIVER_EMAIL")

    def on_any_event(self, event):
        if event.is_directory:
            return

        if event.event_type == "created":
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"[{timestamp}] Received created event - {event.src_path}")
            send_email_notification(
                subject="MyeloSeq watchdog: file detected",
                body=(
                    f"Watchdog detected a new file and will start a run.\n\n"
                    f"Time: {timestamp}\n"
                    f"File: {event.src_path}\n"
                ),
                smtp_server=self.smtp_server,
                smtp_port=self.smtp_port,
                sender_email=self.sender_email,
                receiver_email=self.receiver_email,
            )
            time.sleep(5)
            try:
                myelo_worker = myeloseq(self.config_path)
                myelo_worker.workbook = event.src_path
                myelo_worker.start()
                send_email_notification(
                    subject="MyeloSeq watchdog: run finished successfully",
                    body=(
                        f"MyeloSeq run completed successfully.\n\n"
                        f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
                        f"File: {event.src_path}\n"
                    ),
                    smtp_server=self.smtp_server,
                    smtp_port=self.smtp_port,
                    sender_email=self.sender_email,
                    receiver_email=self.receiver_email,
                )
                del myelo_worker
            except Exception as e:
                send_email_notification(
                    subject="MyeloSeq watchdog: run failed",
                    body=(
                        f"MyeloSeq run failed with an error.\n\n"
                        f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
                        f"File: {event.src_path}\n"
                        f"Error: {e}\n\n"
                        f"Traceback:\n{traceback.format_exc()}"
                    ),
                    smtp_server=self.smtp_server,
                    smtp_port=self.smtp_port,
                    sender_email=self.sender_email,
                    receiver_email=self.receiver_email,
                )
                raise


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: myelo_watchdog.py <config_path> <directory_to_watch>")
        sys.exit(1)

    config_path = sys.argv[1]
    directory_to_watch = sys.argv[2]

    w = Watcher(directory_to_watch, config_path)
    w.run()
