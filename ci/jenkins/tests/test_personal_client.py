"""Exercise personal CI selection through a local Jenkins HTTP fixture."""
from http.server import BaseHTTPRequestHandler, HTTPServer
import json
import os
from pathlib import Path
import subprocess
import tempfile
import threading
import unittest
from urllib.parse import parse_qs


CLIENT = Path(__file__).resolve().parents[1] / 'jenkins-personal.sh'


class PersonalClientTests(unittest.TestCase):
    def test_ci_selector_is_sent_independently_and_unsupported_jobs_are_rejected(self):
        submissions = []
        supported = [True]

        class Handler(BaseHTTPRequestHandler):
            def log_message(self, *args):
                pass

            def do_GET(self):
                self.send_response(200)
                self.end_headers()
                names = ['CI_REF'] if supported[0] else []
                self.wfile.write(json.dumps({'property': [{'parameterDefinitions': [
                    {'name': name} for name in names]}]}).encode())

            def do_POST(self):
                submissions.append(parse_qs(self.rfile.read(int(self.headers['Content-Length'])).decode(),
                                            keep_blank_values=True))
                self.send_response(201)
                self.send_header('Location', '/queue/item/42/')
                self.end_headers()

        with tempfile.TemporaryDirectory() as directory, HTTPServer(('127.0.0.1', 0), Handler) as server:
            worker = threading.Thread(target=server.serve_forever, daemon=True)
            worker.start()
            try:
                auth = Path(directory) / 'personal.auth'
                auth.write_text('test:fixture-token\n')
                auth.chmod(0o600)
                env = dict(os.environ, JENKINS_URL='http://127.0.0.1:{}'.format(server.server_port),
                           JENKINS_CLI_AUTH_FILE=str(auth), JENKINS_JOB='gkeyll-ci-personal')
                command = [str(CLIENT), 'run', '--candidate-ref', 'candidate', '--baseline-ref', 'main',
                           '--ci-ref', 'feature/ci']
                result = subprocess.run(command, env=env, capture_output=True, text=True, timeout=15)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(submissions[0]['CI_REF'], ['feature/ci'])
                self.assertEqual(submissions[0]['CANDIDATE_REF'], ['candidate'])
                self.assertEqual(submissions[0]['BASELINE_REF'], ['main'])
                supported[0] = False
                result = subprocess.run(command, env=env, capture_output=True, text=True, timeout=15)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('does not expose CI_REF', result.stderr)
                self.assertEqual(len(submissions), 1)
            finally:
                server.shutdown()
                worker.join(timeout=5)
