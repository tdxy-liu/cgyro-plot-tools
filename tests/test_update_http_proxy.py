"""Offline update proxy tests; never connect to the public default proxy in CI."""
import os
from pathlib import Path
import shutil
import ssl
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import cgyro_update as u


class HttpProxyConfigTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.config = Path(self.tmp.name)/'update_connection.json'
        self.env = patch.dict(os.environ, {u.UPDATE_CONFIG_ENV: str(self.config)})
        self.env.start()
        self.addCleanup(self.env.stop)

    def test_new_machine_default_and_saved_direct_are_distinct(self):
        self.assertEqual(u.load_update_proxy_config().mode, 'http')
        self.assertEqual(u.load_update_proxy_config().http_proxy, 'http://47.102.120.146:18889')
        self.config.write_text('{"mode":"direct"}')
        self.assertEqual(u.load_update_proxy_config().mode, 'direct')
        self.assertEqual(u.load_update_proxy_config().http_proxy, u.DEFAULT_HTTP_PROXY)

    def test_save_load_and_legacy_modes(self):
        for mode in ('http', 'direct', 'socks5', 'ssh-socks'):
            with self.subTest(mode=mode):
                config = u.UpdateProxyConfig(mode=mode, http_proxy='127.0.0.1:8888',
                                            socks_proxy='socks5h://localhost:1080', ssh_host='relay')
                u.save_update_proxy_config(config)
                self.assertEqual(u.load_update_proxy_config(), config)

    def test_proxy_parse_bare_url_ipv6_and_default_port(self):
        for value, expected in (
            ('47.102.120.146:18889', ('47.102.120.146', 18889, u.DEFAULT_HTTP_PROXY)),
            ('http://localhost/', ('localhost', 80, 'http://localhost:80')),
            ('http://[::1]:8888', ('::1', 8888, 'http://[::1]:8888')),
        ):
            self.assertEqual(u._parse_http_proxy(value), expected)

    def test_invalid_proxy_never_stores_credentials_or_falls_back(self):
        for value in ('', 'http://', 'http://host:0', 'http://host:65536',
                      'https://host:8888', 'socks5://host:1080',
                      'http://user:secret@host:80', 'http://host/path',
                      'http://host?x=1', 'http://host#frag', 'host\n:80', 'host:abc'):
            with self.subTest(value=value), self.assertRaises(u.UpdateConnectionError):
                u.save_update_proxy_config(u.UpdateProxyConfig(http_proxy=value))
        self.assertFalse(self.config.exists())
        with self.assertRaises(u.UpdateConnectionError):
            u.normalize_update_proxy_config({'mode': 'typo'})

    def test_cli_modes_and_exclusive_overrides(self):
        parser = u._build_cli_parser()
        for argv, mode in (([], 'http'), (['--direct'], 'direct'),
                           (['--http-proxy'], 'http'),
                           (['--http-proxy', 'localhost:8888'], 'http'),
                           (['--socks5-proxy', 'localhost:1080'], 'socks5'),
                           (['--ssh-relay', 'relay'], 'ssh-socks')):
            self.assertEqual(u._connection_config_from_cli(parser.parse_args(argv)).mode, mode)
        for argv in (['--direct', '--http-proxy'],
                     ['--http-proxy', '--socks5-proxy', 'localhost:1080'],
                     ['--http-proxy', '--ssh-relay', 'relay']):
            with self.assertRaises(u.UpdateConnectionError):
                u._connection_config_from_cli(parser.parse_args(argv))

    def test_session_does_not_mutate_process_environment(self):
        before = dict(os.environ)
        with u._update_proxy_session({'mode': 'http', 'http_proxy': 'localhost:8888'}) as proxy:
            self.assertEqual(proxy, 'http://localhost:8888')
        with u._update_proxy_session({'mode': 'direct'}) as proxy:
            self.assertIsNone(proxy)
        self.assertEqual(dict(os.environ), before)

    def test_git_environment_wins_over_no_proxy_without_global_changes(self):
        with patch.dict(os.environ, {'NO_PROXY': '*', 'no_proxy': 'github.com',
                                     'HTTP_PROXY': 'http://old:1'}):
            before = dict(os.environ)
            env = u._git_environment_for_proxy(u.DEFAULT_HTTP_PROXY)
            for key in ('HTTP_PROXY', 'HTTPS_PROXY', 'ALL_PROXY', 'http_proxy', 'https_proxy'):
                self.assertEqual(env[key], u.DEFAULT_HTTP_PROXY)
            self.assertEqual(env['NO_PROXY'], '')
            self.assertEqual(env['no_proxy'], '')
            self.assertEqual(os.environ['HTTP_PROXY'], before['HTTP_PROXY'])
            self.assertEqual(dict(os.environ), before)
        self.assertIsNone(u._git_environment_for_proxy(None))

    def test_socks_ssh_environment_regression(self):
        env = u._git_environment_for_proxy('socks5h://localhost:1080')
        self.assertIn('--socks-proxy-command', env['GIT_SSH_COMMAND'])
        self.assertIn('ProxyCommand=', env['GIT_SSH_COMMAND'])


class HttpConnectTests(unittest.TestCase):
    def test_connect_target_and_tls_checks_are_preserved(self):
        for target, host, port in (('github.com', 'github.com', 443),
                                   ('github.com:8443', 'github.com', 8443)):
            c = u._HTTPProxyHTTPSConnection(target, '127.0.0.1', 8888)
            self.assertEqual((c.host, c.port), ('127.0.0.1', 8888))
            self.assertEqual((c._tunnel_host, c._tunnel_port), (host, port))
            self.assertTrue(c._context.check_hostname)
            self.assertEqual(c._context.verify_mode, ssl.CERT_REQUIRED)

    def test_fetch_http_handler_is_explicit_even_with_no_proxy(self):
        response = Mock()
        response.__enter__ = Mock(return_value=response)
        response.__exit__ = Mock(return_value=False)
        response.read.return_value = b'0.2.20'
        opener = Mock()
        opener.open.return_value = response
        with patch.dict(os.environ, {'NO_PROXY': '*'}), patch.object(u, 'build_opener', return_value=opener) as build:
            self.assertEqual(u._fetch_text(u.VERSION_URL, proxy_url=u.DEFAULT_HTTP_PROXY), '0.2.20')
            empty_proxy_handler, tunnel_handler = build.call_args.args
            self.assertEqual(empty_proxy_handler.proxies, {})
            self.assertIsInstance(tunnel_handler, u._HTTPProxyHTTPSHandler)

    def test_https_handlers_support_python39_and_312(self):
        for cls in (u._HTTPProxyHTTPSHandler, u._SocksHTTPSHandler):
            handler = cls('localhost', 8888)
            handler.do_open = Mock(return_value='ok')
            request = u.Request(u.VERSION_URL)
            self.assertEqual(handler.https_open(request), 'ok')
            self.assertNotIn('check_hostname', handler.do_open.call_args.kwargs)
            self.assertTrue(handler.do_open.call_args.kwargs['context'].check_hostname)

    def test_connection_failure_does_not_retry_direct(self):
        opener = Mock()
        opener.open.side_effect = u.URLError('proxy unavailable')
        with patch.object(u, 'build_opener', return_value=opener), patch.object(u, 'urlopen') as direct:
            with self.assertRaises(u.UpdateCheckError):
                u._fetch_text(u.VERSION_URL, proxy_url=u.DEFAULT_HTTP_PROXY)
            direct.assert_not_called()

    def test_version_check_uses_same_proxy_for_manifest_fallback(self):
        with patch.object(u, '_fetch_text', side_effect=[None, '0.2.21']) as fetch:
            result = u.check_for_updates('0.2.20', proxy_config={'mode': 'http'})
        self.assertTrue(result.update_available)
        self.assertEqual(result.latest_version, '0.2.21')
        for call in fetch.call_args_list:
            self.assertEqual(call.kwargs['proxy_url'], u.DEFAULT_HTTP_PROXY)


@unittest.skipUnless(shutil.which('git'), 'Git required for local URL rewrite test')
class GitHttpProxyTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.repo = Path(self.tmp.name)
        self.git('init', '-q')
        self.url = 'git@github.com:tdxy-liu/cgyro-plot-tools.git'
        self.git('remote', 'add', 'origin', self.url)

    def git(self, *args):
        return subprocess.run(['git', *args], cwd=self.repo, capture_output=True,
                              text=True, check=True, timeout=10).stdout.strip()

    def test_ssh_to_https_is_effective_and_not_persisted(self):
        before = (self.repo/'.git/config').read_bytes()
        opts = u._http_remote_override(self.repo, 'origin', 10)
        self.assertEqual(self.git(*opts, 'ls-remote', '--get-url', 'origin'),
                         'https://github.com/tdxy-liu/cgyro-plot-tools.git')
        self.assertEqual((self.repo/'.git/config').read_bytes(), before)
        self.assertEqual(self.git('config', '--get', 'remote.origin.url'), self.url)

    def test_https_and_other_remote_types(self):
        for url in ('https://github.com/owner/repo.git', 'https://git.example.org/owner/repo.git'):
            self.git('remote', 'set-url', 'origin', url)
            self.assertEqual(u._http_remote_override(self.repo, 'origin', 10), [])
        for url in ('git@other.example:owner/repo.git', 'http://github.com/owner/repo.git'):
            self.git('remote', 'set-url', 'origin', url)
            with self.assertRaises(u.GitUpdateError):
                u._http_remote_override(self.repo, 'origin', 10)
        self.git('config', '--add', 'remote.origin.url', 'https://github.com/owner/repo.git')
        with self.assertRaisesRegex(u.GitUpdateError, 'multiple URLs'):
            u._http_remote_override(self.repo, 'origin', 10)

    def test_update_keeps_fast_forward_guard_and_proxy_scope(self):
        replies = [str(self.repo), 'refs/heads/main', '', 'oldsha',
                   self.url, 'Already up to date.', 'newsha']
        def fake_git(repo, args, timeout, env=None):
            return SimpleNamespace(returncode=0, stdout=replies.pop(0))
        with patch.object(u, '_run_git', side_effect=fake_git) as git:
            r = u.update_from_git(self.repo, proxy_config={'mode': 'http'})
        self.assertTrue(r.updated)
        pull = next(c for c in git.call_args_list if 'pull' in c.args[1])
        self.assertIn('--ff-only', pull.args[1])
        self.assertIn('http.proxy='+u.DEFAULT_HTTP_PROXY, pull.args[1])
        self.assertTrue(any('.insteadOf=' in a for a in pull.args[1]))
        self.assertEqual(pull.kwargs['env']['HTTPS_PROXY'], u.DEFAULT_HTTP_PROXY)


if __name__ == '__main__':
    unittest.main()
