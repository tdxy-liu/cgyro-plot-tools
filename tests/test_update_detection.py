"""Version discovery must not turn a partial network failure into 'up to date'."""

import base64
import json
import unittest
from unittest.mock import patch

import cgyro_update as updater


class UpdateDetectionTests(unittest.TestCase):
    @staticmethod
    def release(version):
        return json.dumps({"tag_name": "v" + version})

    @staticmethod
    def manifest(version):
        return json.dumps({"encoding": "base64", "content":
                           base64.b64encode((version + "\n").encode()).decode() + "\n"})

    def test_raw_failure_uses_api_with_the_same_proxy_and_no_cache(self):
        with patch.object(updater, "_fetch_text", side_effect=[
                self.release("0.2.22"), updater.UpdateCheckError("raw blocked"),
                self.manifest("0.2.33")]) as fetch:
            result = updater.check_for_updates("0.2.31", proxy_config={"mode": "http"})
        self.assertTrue(result.update_available)
        self.assertEqual(result.latest_version, "0.2.33")
        self.assertEqual(result.source, "GitHub VERSION API")
        self.assertTrue(fetch.call_args_list[-1].args[0].startswith(updater.VERSION_API_URL + "&cb="))
        for call in fetch.call_args_list:
            self.assertIn("cb=", call.args[0])
            self.assertEqual(call.kwargs["proxy_url"], updater.DEFAULT_HTTP_PROXY)

    def test_old_or_equal_release_cannot_confirm_latest_when_manifest_unreachable(self):
        for release in ("0.2.22", "0.2.31"):
            with self.subTest(release=release), patch.object(updater, "_fetch_text", side_effect=[
                    self.release(release), updater.UpdateCheckError("raw blocked"),
                    updater.UpdateCheckError("API blocked")]):
                with self.assertRaisesRegex(updater.UpdateCheckError, "Could not verify.*main/VERSION"):
                    updater.check_for_updates("0.2.31", proxy_config={"mode": "direct"})

    def test_verified_api_can_confirm_no_update(self):
        with patch.object(updater, "_fetch_text", side_effect=[
                self.release("0.2.22"), updater.UpdateCheckError("raw blocked"),
                self.manifest("0.2.31")]):
            result = updater.check_for_updates("0.2.31", proxy_config={"mode": "direct"})
        self.assertFalse(result.update_available)
        self.assertEqual(result.latest_version, "0.2.31")
        self.assertEqual(result.source, "GitHub VERSION API")

    def test_release_failure_does_not_hide_accessible_source_version(self):
        with patch.object(updater, "_fetch_text", side_effect=[
                updater.UpdateCheckError("release unavailable"), "0.2.33"]) as fetch:
            result = updater.check_for_updates("0.2.31", proxy_config={"mode": "direct"})
        self.assertTrue(result.update_available)
        self.assertEqual(fetch.call_count, 2)

    def test_malformed_api_is_not_treated_as_version_text(self):
        for payload in ("not JSON", "[]", '{"encoding":"utf-8","content":"0.2.33"}',
                        '{"encoding":"base64","content":"!invalid!"}',
                        '{"encoding":"base64","content":null}', self.manifest("invalid")):
            with self.subTest(payload=payload), self.assertRaises(updater.UpdateCheckError):
                updater._api_manifest_info(payload, "0.2.31")

    def test_invalid_api_does_not_fall_back_to_old_release_as_latest(self):
        with patch.object(updater, "_fetch_text", side_effect=[
                self.release("0.2.22"), updater.UpdateCheckError("raw blocked"), "{}"]):
            with self.assertRaisesRegex(updater.UpdateCheckError, "Could not verify"):
                updater.check_for_updates("0.2.31", proxy_config={"mode": "direct"})


if __name__ == "__main__":
    unittest.main()
