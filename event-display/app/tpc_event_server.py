#!/usr/bin/env python3
import argparse
import json
import os
import re
import subprocess
import threading
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse


class EventDisplayHandler(SimpleHTTPRequestHandler):
    converter_lock = threading.Lock()

    def __init__(self, *args, project_dir=None, web_dir=None, converter=None, **kwargs):
        self.project_dir = Path(project_dir).resolve()
        self.web_dir = Path(web_dir).resolve()
        self.converter = Path(converter).resolve()
        super().__init__(*args, directory=str(self.web_dir), **kwargs)

    def end_headers(self):
        self.send_header("Cache-Control", "no-store")
        super().end_headers()

    def do_GET(self):
        parsed = urlparse(self.path)
        if parsed.path == "/api/render":
            self.render_entry(parsed.query)
            return
        super().do_GET()

    def render_entry(self, query):
        params = parse_qs(query)
        source = self.first_param(params, "input")
        tree = self.first_param(params, "tree")
        entry_text = self.first_param(params, "entry", "0")

        try:
            entry = int(entry_text)
            if entry < 0:
                raise ValueError
        except ValueError:
            self.write_json(400, {"ok": False, "error": "entry must be a non-negative integer"})
            return

        if not source:
            self.write_json(400, {"ok": False, "error": "input ROOT file is required"})
            return

        output = self.web_dir / "display.root"
        cmd = [
            str(self.converter),
            "-i",
            source,
            "-o",
            str(output),
            "-e",
            str(entry),
            "-n",
            "1",
        ]
        if tree:
            cmd.extend(["-t", tree])

        with self.converter_lock:
            result = subprocess.run(
                cmd,
                cwd=str(self.project_dir),
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )

        if result.returncode != 0:
            self.write_json(
                500,
                {
                    "ok": False,
                    "error": result.stderr.strip() or result.stdout.strip() or "converter failed",
                },
            )
            return

        metadata = self.parse_converter_stdout(result.stdout)
        metadata.update(
            {
                "ok": True,
                "entry": entry,
                "file": "display.root",
                "stdout": result.stdout.strip(),
            }
        )

        self.write_json(
            200,
            metadata,
        )

    @staticmethod
    def first_param(params, key, default=""):
        values = params.get(key)
        return values[0] if values else default

    @staticmethod
    def parse_converter_stdout(stdout):
        metadata = {}
        for key in ("tracks", "raw_hits", "clusters", "track_clusters", "vertex_pairs"):
            match = re.search(rf"{key}=([0-9]+)", stdout)
            if match:
                metadata[key] = int(match.group(1))
        track_match = re.search(r"track_info=(\[.*\])", stdout)
        if track_match:
            try:
                metadata["track_info"] = json.loads(track_match.group(1))
            except json.JSONDecodeError:
                metadata["track_info"] = []
        vertex_match = re.search(r"vertex_info=(\[.*?\])", stdout)
        if vertex_match:
            try:
                metadata["vertex_info"] = json.loads(vertex_match.group(1))
            except json.JSONDecodeError:
                metadata["vertex_info"] = []
        return metadata

    def write_json(self, status, payload):
        body = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)


def main():
    parser = argparse.ArgumentParser(description="TPC JSROOT event-display server")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", default=8080, type=int)
    parser.add_argument("--project-dir", default=".")
    parser.add_argument("--web-dir", default="web")
    parser.add_argument("--converter", default="bin/tpc_event_export")
    args = parser.parse_args()

    project_dir = Path(args.project_dir).resolve()
    web_dir = (project_dir / args.web_dir).resolve()
    converter = (project_dir / args.converter).resolve()

    def handler(*handler_args, **handler_kwargs):
        return EventDisplayHandler(
            *handler_args,
            project_dir=project_dir,
            web_dir=web_dir,
            converter=converter,
            **handler_kwargs,
        )

    server = ThreadingHTTPServer((args.host, args.port), handler)
    print(f"Serving TPC event display on http://{args.host}:{args.port}/")
    print("Render API: /api/render?input=/path/to/file.root&tree=tpc&entry=0")
    server.serve_forever()


if __name__ == "__main__":
    main()
