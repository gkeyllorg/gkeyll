"""Offline publication, retry, identity and report-format tests."""
import argparse
import contextlib
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from ci.jenkins import github_report as report


class ReportingTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.previous = os.getcwd()
        os.chdir(self.tmp.name)
        self.args = argparse.Namespace(
            platform='personal', context='ci/personal', commit='a' * 40, pr='1157', ref='',
            result='success', description='Passed: unit tests and C regressions.', stage='tests',
            repo='gkeyllorg/gkeyll', report='ci-report.md', queue_id='42')
        self.comments, self.statuses, self.calls = [], [], []
        self.pages(['Report summary'])
        self.environment = patch.dict(os.environ, GITHUB_TOKEN='existing-token', BUILD_NUMBER='11')
        self.environment.start()
        self.http = patch.object(report, 'api', side_effect=self.api)
        self.http.start()
        self.sleep_patch = patch.object(report.time, 'sleep')
        self.sleep = self.sleep_patch.start()

    def tearDown(self):
        self.sleep_patch.stop()
        self.http.stop()
        self.environment.stop()
        os.chdir(self.previous)
        self.tmp.cleanup()

    def pages(self, sections):
        pages = [report.MARKER_FORMAT.format(self.args.context if i == 0 else
                 self.args.context + '/part-{}'.format(i + 1)) + '\n\n' + text
                 for i, text in enumerate(sections)]
        Path('ci-report.md.json').write_text(json.dumps(pages))
        Path('ci-report.md').write_text('\n\n'.join(pages))

    def api(self, method, url, token, payload=None):
        self.assertEqual(token, 'existing-token')
        self.calls.append((method, url, payload))
        base = '{}/repos/{}'.format(report.API, self.args.repo)
        comments_url = base + '/commits/{}/comments'.format(self.args.commit)
        if url.endswith('/user'):
            return {'id': 1}, ''
        if method == 'GET' and '/status?' in url:
            self.assertEqual(url.split('?')[0], base + '/commits/{}/status'.format(self.args.commit))
            return {'statuses': [s for s in self.statuses
                                 if s.get('sha', self.args.commit) == self.args.commit][:1]}, ''
        if method == 'POST' and '/statuses/' in url:
            self.assertEqual(url, base + '/statuses/' + self.args.commit)
            status = dict(payload, id=len(self.statuses) + 100, sha=self.args.commit)
            self.statuses.insert(0, status)
            return status, ''
        if method == 'GET':
            self.assertEqual(url.split('?')[0], comments_url)
            return [c for c in self.comments if c.get('commit_id', self.args.commit) == self.args.commit], ''
        if method == 'POST':
            self.assertEqual(url, comments_url)
            comment = dict(payload, id=len(self.comments) + 1, user={'id': 1},
                           commit_id=self.args.commit,
                           html_url='https://github.com/{}/commit/{}#commitcomment-{}'.format(
                               self.args.repo, self.args.commit, len(self.comments) + 1))
            self.comments.append(comment)
        else:
            self.assertEqual(method, 'PATCH')
            comment = next(c for c in self.comments if c['id'] == int(url.rsplit('/', 1)[-1]))
            self.assertEqual(url, base + '/comments/{}'.format(comment['id']))
            comment.update(payload)
        return dict(comment), ''

    def update(self):
        with contextlib.redirect_stdout(io.StringIO()):
            report.update_report(self.args)
        return json.loads(Path('ci-report-delivery.json').read_text())

    def test_uses_existing_token_for_comments_and_linked_commit_status(self):
        delivery = self.update()
        self.assertTrue(delivery['comment'])
        self.assertTrue(delivery['status'])
        self.assertEqual(self.statuses[0]['target_url'], self.comments[0]['html_url'])
        self.assertEqual(self.statuses[0]['state'], 'success')
        self.assertIn('Jenkins build #11; ID 42', self.statuses[0]['description'])
        self.assertNotIn('queue', self.statuses[0]['description'])
        self.assertTrue(any('/statuses/' + self.args.commit in c[1] for c in self.calls))

    def test_retries_update_existing_comment_and_do_not_repeat_status(self):
        self.update()
        self.update()
        self.assertEqual(len(self.comments), 1)
        self.assertEqual(len(self.statuses), 1)

    def test_commit_runs_publish_commit_comments(self):
        self.args.pr = ''
        self.update()
        self.assertTrue(any('/commits/' + self.args.commit + '/comments' in c[1] for c in self.calls))
        self.assertFalse(any('/issues/' in c[1] for c in self.calls))

    def test_pr_runs_publish_on_the_recorded_commit(self):
        self.args.ref = 'feature-branch'
        self.update()
        self.assertTrue(any(c[:2] == ('POST', '{}/repos/{}/commits/{}/comments'.format(
            report.API, self.args.repo, self.args.commit)) for c in self.calls))
        self.assertFalse(any('/issues/' in c[1] or self.args.ref in c[1] for c in self.calls))
        self.assertEqual(self.comments[0]['commit_id'], self.args.commit)

    def test_older_pr_commit_keeps_its_report_after_a_newer_commit_finishes(self):
        original_commit = self.args.commit
        self.args.commit = 'b' * 40
        self.args.queue_id = '43'
        self.pages(['Newer commit summary', 'Warnings', 'More warnings'])
        self.update()
        newer_comments = [dict(c) for c in self.comments]
        newer_status = dict(self.statuses[0])

        self.args.commit = original_commit
        self.args.queue_id = '42'
        self.pages(['Original commit summary'])
        delivery = self.update()
        self.assertTrue(delivery['comment'])
        self.assertTrue(delivery['status'])
        self.assertEqual(len(self.comments), 4)
        self.assertEqual(self.comments[:3], newer_comments)
        self.assertEqual(self.comments[-1]['commit_id'], original_commit)
        self.assertEqual(report.comment_metadata(self.comments[-1]['body'])['commit'], original_commit)
        self.assertEqual(len(self.statuses), 2)
        self.assertEqual(self.statuses[1], newer_status)
        self.assertEqual(self.statuses[0]['target_url'], self.comments[-1]['html_url'])

    def test_progress_posts_status_without_comment_notifications(self):
        self.args.result = 'pending'
        self.args.description = 'Running: Build candidate.'
        self.update()
        self.assertEqual(self.comments, [])
        self.assertIn('Gkeyll CI running: Build candidate', self.statuses[0]['description'])
        self.assertIn(self.args.commit, self.statuses[0]['target_url'])

    def test_stage_updates_preserve_jenkins_timing_within_status_limit(self):
        self.statuses.append(dict(state='pending', context=self.args.context,
            description='Gkeyll CI running: candidate unit build; elapsed 29 m; ETA ~1 hr 50 m '
                        '(Jenkins build #11; ID 42).'))
        self.args.result = 'pending'
        self.args.description = 'Running: ' + 'candidate regression build ' * 10
        self.update()
        description = self.statuses[0]['description']
        self.assertTrue(description.startswith('Gkeyll CI running: candidate regression build'))
        self.assertIn('; elapsed 29 m; ETA ~1 hr 50 m', description)
        self.assertTrue(description.endswith(' (Jenkins build #11; ID 42).'))
        self.assertLessEqual(len(description), 140)
        self.assertEqual(self.comments, [])

    def test_new_run_does_not_inherit_previous_run_timing(self):
        self.statuses.append(dict(state='pending', context=self.args.context,
            description='Gkeyll CI running: candidate unit build; elapsed 29 m; ETA ~1 hr 50 m '
                        '(Jenkins build #10; ID 41).'))
        self.args.result = 'pending'
        self.args.description = 'Running: Build candidate.'
        self.update()
        self.assertNotIn('elapsed 29 m', self.statuses[0]['description'])

    def test_legacy_statuses_retain_newer_run_and_terminal_protection(self):
        for state, queue_id in [('pending', '43'), ('success', '42')]:
            with self.subTest(state=state, queue_id=queue_id):
                self.statuses[:] = [dict(state=state, context=self.args.context,
                    description='Legacy status (Jenkins queue #{}).'.format(queue_id))]
                self.args.result = 'pending'
                delivery = self.update()
                self.assertTrue(delivery['status_skipped'])
                self.assertEqual(len(self.statuses), 1)

    def test_late_progress_does_not_reopen_terminal_status(self):
        self.update()
        self.args.result = 'pending'
        delivery = self.update()
        self.assertTrue(delivery['status_skipped'])
        self.assertEqual(self.statuses[0]['state'], 'success')
        self.assertEqual(len(self.statuses), 1)

    def test_older_run_cannot_replace_newer_report_or_status(self):
        self.args.queue_id = '43'
        self.update()
        self.args.queue_id = '42'
        delivery = self.update()
        self.assertTrue(delivery['comment_skipped'])
        self.assertTrue(delivery['status_skipped'])
        self.assertEqual(len(self.statuses), 1)
        self.assertEqual(report.comment_metadata(self.comments[0]['body'])['queue_id'], '43')

    def test_build_numbers_protect_reports_without_queue_listener(self):
        self.assertTrue(report.superseded({'job': 'ci/PR-1', 'build_number': '10'},
                                         {'job': 'ci/PR-1', 'build_number': '9'}))
        self.assertFalse(report.superseded({'job': 'ci/PR-2', 'build_number': '10'},
                                          {'job': 'ci/PR-1', 'build_number': '9'}))

    def test_overflow_pages_are_linked_and_old_pages_are_cleared(self):
        self.pages(['Summary', 'Warnings', 'More warnings'])
        self.update()
        self.assertIn('part-3', self.comments[0]['body'])
        self.assertIn('Next diagnostics page', self.comments[1]['body'])
        main = next(c for c in self.comments if 'Summary' in c['body'])
        self.assertIn('[Part 2](', main['body'])
        self.assertIn('[Part 3](', main['body'])
        self.assertEqual(self.statuses[0]['target_url'], main['html_url'])
        self.pages(['Short clean report'])
        self.update()
        self.assertEqual(len(self.comments), 3)
        self.assertTrue(all('no longer needed' in c['body'] for c in self.comments[:2]))

    def test_other_users_and_other_contexts_are_never_edited(self):
        self.comments.append(dict(id=1, user={'id': 2}, body=Path('ci-report.md').read_text()))
        self.comments.append(dict(id=2, user={'id': 1}, body='<!-- gkeyll-ci-report context="ci/another" -->\nOther report'))
        originals = [dict(c) for c in self.comments]
        self.update()
        self.assertEqual(self.comments[:2], originals)
        self.assertEqual(len(self.comments), 3)

    def test_large_reports_keep_navigation_within_the_comment_limit(self):
        self.pages(['Summary\n' + 'x' * 59500] + ['Diagnostics'] * 200)
        self.update()
        self.assertIn('[Continue diagnostics]', self.comments[-1]['body'])
        self.assertTrue(all(len(c['body'].encode('utf-8')) <= 65000 for c in self.comments))
        self.assertTrue(all('Next diagnostics page' in c['body'] for c in self.comments[1:-1]))

    def test_cancellation_timeout_and_errors_use_error_status(self):
        for result in ('cancelled', 'timed_out', 'error'):
            self.args.result = result
            self.args.description = ''
            self.update()
            self.assertEqual(self.statuses[0]['state'], 'error')
        self.assertIn('Errored at stage', self.statuses[0]['description'])

    def test_comment_outage_still_publishes_test_result_and_records_error(self):
        with patch.object(report, 'publish_report', side_effect=report.GitHubError('POST', 403)):
            with self.assertRaisesRegex(RuntimeError, 'Comment: GitHub POST failed: HTTP 403'):
                self.update()
        self.assertEqual(self.statuses[0]['state'], 'success')
        self.assertIn('report unavailable', self.statuses[0]['description'])
        self.assertTrue(json.loads(Path('ci-report-delivery.json').read_text())['status'])
        self.sleep.assert_not_called()

    def test_status_outage_still_publishes_comment_and_records_error(self):
        with patch.object(report, 'publish_status', side_effect=OSError('network unavailable')):
            with self.assertRaisesRegex(RuntimeError, 'Commit status'):
                self.update()
        self.assertEqual(len(self.comments), 1)
        self.assertTrue(json.loads(Path('ci-report-delivery.json').read_text())['comment'])
        self.assertEqual(self.sleep.call_count, 2)

    def test_lost_comment_response_is_retried_without_duplicate_comment(self):
        failed = False
        def flaky(method, url, token, payload=None):
            nonlocal failed
            result = self.api(method, url, token, payload)
            if method == 'POST' and '/comments' in url and not failed:
                failed = True
                raise OSError('lost response after POST')
            return result
        with patch.object(report, 'api', side_effect=flaky):
            self.update()
        self.assertEqual(len(self.comments), 1)
        self.sleep.assert_called_once_with(1)

    def test_missing_token_and_invalid_sha_are_recorded_without_api_calls(self):
        for commit, token in [('a' * 40, ''), ('main', 'existing-token')]:
            self.args.commit = commit
            with patch.dict(os.environ, GITHUB_TOKEN=token), self.assertRaises(RuntimeError):
                self.update()
            self.assertTrue(json.loads(Path('ci-report-delivery.json').read_text())['errors'])
        self.assertEqual(self.calls, [])

    def test_standalone_publication_requires_sha_even_with_pr_or_ref(self):
        for pr in ('', '1157'):
            for commit in ('', 'main', 'a' * 7):
                with self.subTest(pr=pr, commit=commit):
                    self.args.pr = pr
                    self.args.ref = 'main'
                    self.args.commit = commit
                    with self.assertRaisesRegex(ValueError, '--commit must be a full 40-character SHA'):
                        report.publish_report(self.args)
        self.assertEqual(self.calls, [])

    def test_standalone_cli_publishes_pr_report_on_the_commit(self):
        with contextlib.redirect_stdout(io.StringIO()):
            report.main(['publish', '--context', self.args.context, '--commit', self.args.commit,
                         '--pr', self.args.pr, '--ref', 'feature-branch'])
        self.assertEqual(len(self.comments), 1)
        self.assertEqual(self.comments[0]['commit_id'], self.args.commit)
        self.assertFalse(any('/issues/' in c[1] for c in self.calls))

    def test_oversized_and_mismatched_pages_fail_before_network_access(self):
        for page in ['bad marker', report.MARKER_FORMAT.format(self.args.context) + '\n' + '漢' * 30000]:
            Path('ci-report.md.json').write_text(json.dumps([page]))
            with self.assertRaises(ValueError):
                report.publish_report(self.args)
        self.assertEqual(self.calls, [])

    def test_pagination_discovers_existing_report(self):
        self.update()
        saved = dict(self.comments[0])
        def paginated(method, url, token, payload=None):
            if method == 'GET' and '/comments?' in url and 'page=2' not in url:
                return [], '<' + url + '&page=2>; rel="next"'
            return self.api(method, url, token, payload)
        with patch.object(report, 'api', side_effect=paginated):
            self.update()
        self.assertEqual(len(self.comments), 1)
        self.assertEqual(self.comments[0], saved)

    def test_failed_checkout_still_records_the_accepted_candidate(self):
        args = argparse.Namespace(platform='personal', context='ci/test', result='failure', pr='1157', output='ci-report.md')
        with patch.dict(os.environ, CI_REPORT_COMMIT='a' * 40), contextlib.redirect_stdout(io.StringIO()):
            report.build_report(args)
        self.assertIn('PR #1157 @ `aaaaaaa`', Path('ci-report.md').read_text())

    def test_failure_status_and_report_include_execution_queue_and_stage_time(self):
        args = argparse.Namespace(platform='personal', context='ci/test', result='failure', pr='', output='ci-report.md')
        Path('ci-stage-timings.json').write_text(json.dumps({'stages': [dict(
            stage='Build candidate', elapsed_ms=7123, result='failure')]}))
        with patch.dict(os.environ, CI_REPORT_START_MS='1700000000000', CI_REPORT_END_MS='1700000007000',
                        CI_QUEUE_ENQUEUED_MS='1699999990000', CI_BOOTSTRAP_RETRY_WAIT_MS='1500',
                        CI_FAILURE_STAGE='Build candidate'), contextlib.redirect_stdout(io.StringIO()):
            report.build_report(args)
        self.assertIn('Failed after 7 s', Path('ci-status-description.txt').read_text())
        text = Path('ci-report.md').read_text()
        self.assertIn('Queue wait:** 10 s', text)
        self.assertIn('Setup retry waiting:** 1 s', text)
        self.assertIn('7.123 s | failure', text)

    def test_pending_report_does_not_claim_failure(self):
        args = argparse.Namespace(platform='personal', context='ci/test', result='pending', pr='1157', output='ci-report.md')
        with patch.dict(os.environ, CI_FAILURE_STAGE='Build candidate'), contextlib.redirect_stdout(io.StringIO()):
            report.build_report(args)
        text = Path(args.output).read_text()
        self.assertIn('Current stage:** Build candidate', text)
        self.assertNotIn('Failed at stage', text)
        self.assertNotIn('finished ', text)

    def test_utf8_reports_fit_github_limits_without_broken_sections(self):
        sections = report.detail_sections('Warnings', 'warning: λ漢字\n' * 14000)
        pages = report.report_pages('Summary', sections, 'ci/test')
        self.assertGreater(len(pages), 2)
        for page in pages:
            self.assertLessEqual(len(page.encode('utf-8')), report.COMMENT_LIMIT)
            self.assertEqual(page.count('<details>'), page.count('</details>'))
            self.assertEqual(page.count('<div>'), page.count('</div>'))
            self.assertEqual(page.count('<code>'), page.count('</code>'))

    def test_warning_comparison_requires_the_same_completed_step(self):
        Path('candidate-unit-build.log').write_text('core/zero/a.c:50: warning: example\n')
        Path('baseline-unit-build.log').write_text('completed with no warnings\n')
        Path('baseline-unit-build-seconds.txt').write_text('5\n')
        # An identical warning in a different baseline command is not evidence
        # that it already existed in the corresponding candidate command.
        Path('baseline-install.log').write_text('core/zero/a.c:10: warning: example\n')
        Path('baseline-install-seconds.txt').write_text('1\n')
        report.diagnostic_sections('diagnostics.json')
        self.assertEqual(json.loads(Path('diagnostics.json').read_text())['new_warnings'], 1)
        Path('baseline-unit-build.log').write_text('core/zero/a.c:10: warning: example\n')
        report.diagnostic_sections('diagnostics.json')
        self.assertEqual(json.loads(Path('diagnostics.json').read_text())['new_warnings'], 0)

    def test_failed_baseline_does_not_claim_clean_warning_comparison(self):
        Path('candidate-unit-build.log').write_text('core/zero/a.c:50: warning: example\n')
        Path('baseline-unit-build.log').write_text('failed before compiling\n')
        Path('baseline-unit-build.log.exit').write_text('2\n')
        report.diagnostic_sections('diagnostics.json')
        values = json.loads(Path('diagnostics.json').read_text())
        self.assertIsNone(values['new_warnings'])
        self.assertEqual(values['unclassified_warnings'], 1)


if __name__ == '__main__':
    unittest.main()
