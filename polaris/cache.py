"""On-disk request cache + polite HTTP with retry/backoff.

Every external call is cached by a stable hash of (method,url,params/body) so the
entire pipeline is reproducible and reruns are offline-fast. ASCII-safe logging.
"""
from __future__ import annotations
import hashlib, json, os, time, pathlib, sys
import requests

ROOT = pathlib.Path(__file__).resolve().parents[1]
CACHE_DIR = ROOT / "data" / "cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

_SESSION = requests.Session()
_SESSION.headers.update({"User-Agent": "POLARIS-research/1.0 (non-coding variant pipeline)"})

# polite per-host minimum interval (seconds) applied only on cache miss
_HOST_DELAY = {
    "rest.ensembl.org": 0.30,
    "eutils.ncbi.nlm.nih.gov": 0.34,
    "api.genome.ucsc.edu": 0.20,
    "gtexportal.org": 0.20,
    "gnomad.broadinstitute.org": 0.20,
    "jaspar.elixir.no": 0.15,
    "api.platform.opentargets.org": 0.15,
}
_last_call = {}


def _key(method, url, params, data):
    blob = json.dumps([method, url, params, data], sort_keys=True, default=str)
    return hashlib.sha256(blob.encode()).hexdigest()[:24]


def _host(url):
    return url.split("/")[2] if "://" in url else url


def _throttle(url):
    h = _host(url)
    d = _HOST_DELAY.get(h, 0.0)
    if d:
        dt = time.time() - _last_call.get(h, 0.0)
        if dt < d:
            time.sleep(d - dt)
    _last_call[h] = time.time()


def request(method, url, params=None, json_body=None, headers=None, timeout=60,
            retries=4, cache=True, is_json=True):
    """Cached HTTP request. Returns parsed JSON (is_json) or text."""
    subdir = CACHE_DIR / _host(url).replace(":", "_")
    subdir.mkdir(parents=True, exist_ok=True)
    k = _key(method, url, params, json_body)
    fp = subdir / (k + (".json" if is_json else ".txt"))
    if cache and fp.exists():
        txt = fp.read_text(encoding="utf-8")
        return json.loads(txt) if is_json else txt

    last_err = None
    for attempt in range(retries):
        try:
            _throttle(url)
            r = _SESSION.request(method, url, params=params, json=json_body,
                                 headers=headers, timeout=timeout)
            if r.status_code in (429, 500, 502, 503, 504):
                raise requests.HTTPError(f"HTTP {r.status_code}")  # transient -> retry
            r.raise_for_status()
            out = r.json() if is_json else r.text
            if cache:
                fp.write_text(json.dumps(out) if is_json else out, encoding="utf-8")
            return out
        except requests.HTTPError as e:  # terminal 4xx (e.g. dyneqtl 400) -> do not retry
            resp = getattr(e, "response", None)
            if resp is not None and 400 <= resp.status_code < 500 and resp.status_code != 429:
                raise
            last_err = e
        except Exception as e:  # noqa
            last_err = e
            wait = 1.5 * (2 ** attempt)
            sys.stderr.write(f"[cache] retry {attempt+1}/{retries} {_host(url)} "
                             f"({str(e)[:80]}) wait {wait:.1f}s\n")
            time.sleep(wait)
    raise RuntimeError(f"request failed after {retries} tries: {url} :: {last_err}")


def get_json(url, params=None, **kw):
    return request("GET", url, params=params, is_json=True, **kw)


def get_text(url, params=None, **kw):
    return request("GET", url, params=params, is_json=False, **kw)


def post_json(url, json_body, **kw):
    return request("POST", url, json_body=json_body, is_json=True, **kw)


def graphql(url, query, variables=None, **kw):
    return post_json(url, {"query": query, "variables": variables or {}}, **kw)
