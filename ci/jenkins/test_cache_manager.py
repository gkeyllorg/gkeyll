#!/usr/bin/env python3
"""Focused safety tests for the persistent Jenkins cache manager."""

import io
import json
import os
from pathlib import Path
import subprocess
import tempfile
import types
import unittest
from unittest.mock import patch

import cache_manager


class Response(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


class CacheManagerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="gkeyll-cache-manager-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.states = {}
        self.slurm = []
        self.current = "jenkins-test-2"
        self.platform = "stellar-cpu"
        self.fake_api = patch("cache_manager.urllib.request.urlopen", side_effect=self.api)
        self.fake_queue = patch("cache_manager.subprocess.run", side_effect=self.command)
        self.fake_api.start()
        self.fake_queue.start()
        self.addCleanup(self.fake_api.stop)
        self.addCleanup(self.fake_queue.stop)

    def args(self, tag=None, number="2"):
        return types.SimpleNamespace(root=str(self.root), platform=self.platform,
            tag=tag or self.current, job="test", number=number,
            url=f"http://jenkins.invalid/job/test/{number}/", controller_url="", auth_file="")

    def manager(self, tag=None, number="2"):
        return cache_manager.CacheManager(self.args(tag, number))

    def api(self, request, timeout):
        number = request.full_url.split("/job/test/")[1].split("/")[0]
        if self.states.get(number) == "unavailable":
            raise cache_manager.urllib.error.URLError("controller down")
        if self.states.get(number) == "missing":
            raise cache_manager.urllib.error.HTTPError(request.full_url, 404, "missing", {}, None)
        return Response(json.dumps({"building": self.states.get(number) == "building"}).encode())

    def command(self, argv, **_):
        if argv[0] == "squeue":
            return subprocess.CompletedProcess(argv, 0, "\n".join(self.slurm), "")
        raise AssertionError(f"unexpected command: {argv}")

    def register(self, tag, number, baseline=""):
        runs = self.root / "cache-runs" / self.platform
        runs.mkdir(parents=True, exist_ok=True)
        (runs / f"{tag}.json").write_text(json.dumps({"tag": tag, "job": "test", "number": number,
            "baseline": baseline, "state": "building"}))

    def candidate(self, tag):
        path = self.root / "candidate-cache" / self.platform / tag
        path.mkdir(parents=True)
        (path / "results.txt").write_text("inspect me")
        return path

    def test_begin_preserves_live_candidate_and_slurm_workspace(self):
        active, done = "jenkins-test-1", "jenkins-test-0"
        self.register(active, "1")
        self.register(done, "0")
        active_path = self.candidate(active)
        done_path = self.candidate(done)
        job_root = self.root / "workspaces" / "test"
        for tag in (active, done, self.current, active + "@tmp", done + "@tmp", self.current + "@tmp"):
            (job_root / tag).mkdir(parents=True)
        self.states = {"1": "done", "0": "done"}
        self.slurm = [str(job_root / active)]
        self.manager().begin()
        self.assertTrue(active_path.exists())
        self.assertTrue((job_root / active).exists())
        self.assertTrue((job_root / (active + "@tmp")).exists())
        self.assertFalse(done_path.exists())
        self.assertFalse((job_root / done).exists())
        self.assertFalse((job_root / (done + "@tmp")).exists())
        self.assertTrue((job_root / self.current).exists())

    def test_active_login_build_is_preserved(self):
        tag = "jenkins-test-1"
        self.register(tag, "1")
        path = self.candidate(tag)
        self.states = {"1": "building"}
        self.manager().begin()
        self.assertTrue(path.exists())

    def test_controller_outage_fails_closed(self):
        tag = "jenkins-test-1"
        self.register(tag, "1")
        path = self.candidate(tag)
        self.states = {"1": "unavailable"}
        with self.assertRaisesRegex(RuntimeError, "unavailable"):
            self.manager().begin()
        self.assertTrue(path.exists())

    def test_controller_query_uses_local_url_instead_of_browser_tunnel(self):
        self.register("jenkins-test-1", "1")
        self.candidate("jenkins-test-1")
        args = self.args()
        args.url = "http://127.0.0.1:8083/job/test/2/"
        args.controller_url = "http://127.0.0.1:8080"
        manager = cache_manager.CacheManager(args)
        with patch.object(manager, "active_slurm_directories", return_value=set()), \
             patch("cache_manager.urllib.request.urlopen", wraps=self.api) as query:
            manager.begin()
        self.assertEqual(query.call_args.args[0].full_url,
                         "http://127.0.0.1:8080/job/test/1/api/json?tree=building,result")

    def test_different_runs_of_same_commit_keep_separate_directories(self):
        sha = "a" * 40
        first = self.root / "candidate-cache" / self.platform / "jenkins-test-1" / sha
        second = self.root / "candidate-cache" / self.platform / "jenkins-test-2" / sha
        first.mkdir(parents=True)
        second.mkdir(parents=True)
        self.assertNotEqual(first, second)
        self.assertTrue(first.exists() and second.exists())

    def test_baseline_owner_wait_and_failed_build_cleanup(self):
        sha = "a" * 40
        first = self.manager("jenkins-test-1", "1")
        self.states = {"1": "building"}
        first.begin()
        self.manager().begin()
        self.assertEqual(first.claim_baseline(sha), "build")
        self.assertEqual(self.manager().claim_baseline(sha), "wait")
        self.states["1"] = "done"
        first.abort_baseline(sha)
        self.assertFalse((self.root / "baseline-cache" / self.platform / sha).exists())
        self.assertEqual(self.manager().claim_baseline(sha), "build")

    def test_second_run_waits_then_reuses_published_baseline(self):
        sha = "c" * 40
        first = self.manager("jenkins-test-1", "1")
        second = self.manager("jenkins-test-2", "2")
        self.states = {"1": "building", "2": "building"}
        first.begin()
        second.begin()
        self.assertEqual(first.claim_baseline(sha), "build")
        self.assertEqual(second.claim_baseline(sha), "wait")
        baseline = self.root / "baseline-cache" / self.platform / sha
        with patch.object(cache_manager.CacheManager, "baseline_valid",
                          side_effect=lambda path, commit: path == baseline and
                          (path / "cache-manifest.txt").exists()):
            (baseline / "cache-manifest.txt").write_text("complete")
            first.publish_baseline(sha)
            self.assertEqual(second.claim_baseline(sha), "hit")
            self.assertFalse((baseline / "cache-owner.txt").exists())

    def test_symlink_is_not_followed(self):
        outside = self.root / "outside"
        outside.mkdir()
        (outside / "keep").write_text("yes")
        parent = self.root / "candidate-cache" / self.platform
        parent.mkdir(parents=True)
        (parent / "jenkins-test-1").symlink_to(outside, target_is_directory=True)
        self.manager().begin()
        self.assertEqual((outside / "keep").read_text(), "yes")

    def test_successful_baseline_and_latest_candidate_survive_finish(self):
        sha = "a" * 40
        manager = self.manager()
        manager.begin()
        candidate = self.candidate(self.current)
        self.assertEqual(manager.claim_baseline(sha), "build")
        baseline = self.root / "baseline-cache" / self.platform / sha
        with patch.object(cache_manager.CacheManager, "baseline_valid",
                          side_effect=lambda path, commit: path == baseline and path.exists()):
            manager.publish_baseline(sha)
            manager.finish()
            self.assertTrue(baseline.exists())
            self.assertTrue(candidate.exists())
            next_manager = self.manager("jenkins-test-3", "3")
            next_manager.begin()
            self.assertFalse(candidate.exists())
            self.assertTrue(baseline.exists())
            self.assertEqual(next_manager.claim_baseline(sha), "hit")

    def test_failed_baseline_removed_but_completed_candidate_kept(self):
        sha = "b" * 40
        manager = self.manager()
        manager.begin()
        candidate = self.candidate(self.current)
        self.assertEqual(manager.claim_baseline(sha), "build")
        baseline = self.root / "baseline-cache" / self.platform / sha
        (baseline / "partial-result").write_text("failure evidence copied to Jenkins")
        manager.abort_baseline(sha)
        manager.finish()
        self.assertFalse(baseline.exists())
        self.assertTrue(candidate.exists())

    def test_unknown_slurm_working_directory_stops_cleanup(self):
        previous = "jenkins-test-1"
        self.register(previous, "1")
        candidate = self.candidate(previous)
        self.slurm = ["(null)"]
        with self.assertRaisesRegex(RuntimeError, "unknown Slurm working directory"):
            self.manager().begin()
        self.assertTrue(candidate.exists())

    def test_active_run_keeps_old_baseline_during_new_baseline_build(self):
        old_sha, new_sha = "a" * 40, "b" * 40
        old_tag = "jenkins-test-1"
        self.register(old_tag, "1", old_sha)
        self.states = {"1": "building"}
        old_baseline = self.root / "baseline-cache" / self.platform / old_sha
        old_baseline.mkdir(parents=True)
        with patch.object(cache_manager.CacheManager, "baseline_valid",
                          side_effect=lambda path, sha: path == old_baseline):
            manager = self.manager()
            manager.begin()
            self.assertEqual(manager.claim_baseline(new_sha), "build")
            self.assertTrue(old_baseline.exists())

    def test_newest_started_candidate_wins_when_older_run_finishes_last(self):
        old = self.manager("jenkins-test-1", "1")
        new = self.manager("jenkins-test-2", "2")
        self.states = {"1": "building", "2": "building"}
        old.begin()
        old_candidate = self.candidate(old.tag)
        new.begin()
        new_candidate = self.candidate(new.tag)
        new.finish()
        self.assertTrue(old_candidate.exists() and new_candidate.exists())
        old.finish()
        self.assertFalse(old_candidate.exists())
        self.assertTrue(new_candidate.exists())

    def test_newest_started_candidate_wins_when_newer_run_finishes_last(self):
        old = self.manager("jenkins-test-1", "1")
        new = self.manager("jenkins-test-2", "2")
        self.states = {"1": "building", "2": "building"}
        old.begin()
        old_candidate = self.candidate(old.tag)
        new.begin()
        new_candidate = self.candidate(new.tag)
        old.finish()
        self.assertFalse(old_candidate.exists())
        self.assertTrue(new_candidate.exists())
        new.finish()
        self.assertFalse(old_candidate.exists())
        self.assertTrue(new_candidate.exists())

    def test_new_run_removes_completed_candidate_during_jenkins_tail_steps(self):
        old = self.manager("jenkins-test-1", "1")
        self.states = {"1": "building"}
        old.begin()
        old_candidate = self.candidate(old.tag)
        old.finish()
        self.assertTrue(old_candidate.exists())
        self.manager().begin()
        self.assertFalse(old_candidate.exists())

    def test_newer_run_without_candidate_does_not_retain_older_overlap(self):
        old = self.manager("jenkins-test-1", "1")
        new = self.manager("jenkins-test-2", "2")
        self.states = {"1": "building", "2": "building"}
        old.begin()
        old_candidate = self.candidate(old.tag)
        new.begin()
        self.assertTrue(old_candidate.exists())
        new.finish()
        old.finish()
        self.assertFalse(old_candidate.exists())


if __name__ == "__main__":
    unittest.main()
