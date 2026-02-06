"""Tests for dcd_mapping.resource_utils"""

from contextlib import ExitStack
from unittest import mock

import pytest
import requests

from dcd_mapping.resource_utils import request_with_backoff


class _DummyResponse:
    def __init__(self, status_code=200, headers=None):
        self.status_code = status_code
        self.headers = headers or {}

    def raise_for_status(self):
        if not (200 <= self.status_code < 300):
            msg = f"HTTP {self.status_code}"
            raise requests.HTTPError(msg)


def _sequence_side_effect(values):
    """Turn a list of values/exceptions into a side_effect callable."""
    it = iter(values)

    def _next(*args, **kwargs):  # noqa: ANN002
        v = next(it)
        if isinstance(v, BaseException):
            raise v
        return v

    return _next


def test_success_200_returns_response():
    dummy = _DummyResponse(200)
    with mock.patch(
        "dcd_mapping.resource_utils.requests.get", return_value=dummy
    ) as get_mock:
        resp = request_with_backoff("http://example.com/resource")
        assert resp is dummy
        get_mock.assert_called_once()


def test_timeout_then_success_retries_with_backoff():
    first_exc = requests.Timeout("timeout")
    second_resp = _DummyResponse(200)
    with ExitStack() as stack:
        get_mock = stack.enter_context(
            mock.patch(
                "dcd_mapping.resource_utils.requests.get",
                side_effect=_sequence_side_effect([first_exc, second_resp]),
            )
        )
        sleep_mock = stack.enter_context(
            mock.patch("dcd_mapping.resource_utils.time.sleep")
        )
        resp = request_with_backoff("http://example.com/resource", backoff_factor=0.5)
        assert resp is second_resp
        assert get_mock.call_count == 2
        # first backoff attempt uses factor * (2**0) == 0.5
        sleep_mock.assert_called_once_with(0.5)


def test_connection_error_until_max_raises():
    seq = [requests.ConnectionError("conn err")] * 5
    with ExitStack() as stack:
        stack.enter_context(
            mock.patch(
                "dcd_mapping.resource_utils.requests.get",
                side_effect=_sequence_side_effect(seq),
            )
        )
        sleep_mock = stack.enter_context(
            mock.patch("dcd_mapping.resource_utils.time.sleep")
        )
        with pytest.raises(requests.ConnectionError):
            request_with_backoff(
                "http://example.com/resource", max_retries=5, backoff_factor=0.1
            )
        # should have slept 4 times (no sleep on final raise)
        assert sleep_mock.call_count == 4
        # verify exponential values: 0.1, 0.2, 0.4, 0.8
        assert [c.args[0] for c in sleep_mock.mock_calls] == [0.1, 0.2, 0.4, 0.8]


def test_5xx_then_success_retries():
    seq = [_DummyResponse(503), _DummyResponse(200)]
    with ExitStack() as stack:
        get_mock = stack.enter_context(
            mock.patch(
                "dcd_mapping.resource_utils.requests.get",
                side_effect=_sequence_side_effect(seq),
            )
        )
        sleep_mock = stack.enter_context(
            mock.patch("dcd_mapping.resource_utils.time.sleep")
        )
        resp = request_with_backoff("http://example.com/resource", backoff_factor=0.25)
        assert resp.status_code == 200
        assert get_mock.call_count == 2
        sleep_mock.assert_called_once_with(0.25)  # first attempt sleep


def test_5xx_until_max_raises_http_error():
    seq = [_DummyResponse(500), _DummyResponse(500), _DummyResponse(500)]
    with ExitStack() as stack:
        stack.enter_context(
            mock.patch(
                "dcd_mapping.resource_utils.requests.get",
                side_effect=_sequence_side_effect(seq),
            )
        )
        sleep_mock = stack.enter_context(
            mock.patch("dcd_mapping.resource_utils.time.sleep")
        )
        with pytest.raises(requests.HTTPError):
            request_with_backoff(
                "http://example.com/resource", max_retries=3, backoff_factor=0.1
            )
        # slept for first two attempts, then raised on third
        assert [c.args[0] for c in sleep_mock.mock_calls] == [0.1, 0.2]


def test_429_with_retry_after_header_respected():
    resp1 = _DummyResponse(429, headers={"Retry-After": "1.7"})
    resp2 = _DummyResponse(200)
    with ExitStack() as stack:
        get_mock = stack.enter_context(
            mock.patch(
                "dcd_mapping.resource_utils.requests.get",
                side_effect=_sequence_side_effect([resp1, resp2]),
            )
        )
        sleep_mock = stack.enter_context(
            mock.patch("dcd_mapping.resource_utils.time.sleep")
        )
        resp = request_with_backoff("http://example.com/resource", backoff_factor=0.9)
        assert resp.status_code == 200
        assert get_mock.call_count == 2
        sleep_mock.assert_called_once_with(1.7)


def test_429_with_bad_retry_after_falls_back_to_backoff():
    resp1 = _DummyResponse(429, headers={"Retry-After": "not-a-number"})
    resp2 = _DummyResponse(200)
    with ExitStack() as stack:
        stack.enter_context(
            mock.patch(
                "dcd_mapping.resource_utils.requests.get",
                side_effect=_sequence_side_effect([resp1, resp2]),
            )
        )
        sleep_mock = stack.enter_context(
            mock.patch("dcd_mapping.resource_utils.time.sleep")
        )
        resp = request_with_backoff("http://example.com/resource", backoff_factor=0.3)
        assert resp.status_code == 200
        sleep_mock.assert_called_once_with(0.3)  # backoff_factor * (2**0)


def test_non_retryable_4xx_raises_immediately():
    resp = _DummyResponse(404)
    with ExitStack() as stack:
        stack.enter_context(
            mock.patch("dcd_mapping.resource_utils.requests.get", return_value=resp)
        )
        sleep_mock = stack.enter_context(
            mock.patch("dcd_mapping.resource_utils.time.sleep")
        )
        with pytest.raises(requests.HTTPError):
            request_with_backoff("http://example.com/resource")
        # no sleeps for non-retryable errors
        sleep_mock.assert_not_called()


def test_exhausted_retries_without_response_raises_request_exception():
    # The only way to trigger the terminal state in the function is to not even
    # attempt a request (max_retries=0)
    with mock.patch(
        "dcd_mapping.resource_utils.requests.get", return_value=_DummyResponse(500)
    ), mock.patch("dcd_mapping.resource_utils.time.sleep"), pytest.raises(
        Exception  # noqa: PT011
    ) as exc:
        request_with_backoff("http://example.com/resource", max_retries=0)
    assert "Failed to fetch" in str(exc.value)


def test_kwargs_are_passed_through_to_requests_get():
    dummy = _DummyResponse(200)
    with mock.patch(
        "dcd_mapping.resource_utils.requests.get", return_value=dummy
    ) as get_mock:
        request_with_backoff(
            "http://example.com/resource", headers={"X-Test": "1"}, params={"q": "x"}
        )
        called_kwargs = get_mock.call_args.kwargs
        assert called_kwargs["headers"] == {"X-Test": "1"}
        assert called_kwargs["params"] == {"q": "x"}
        assert called_kwargs["timeout"] == 60
