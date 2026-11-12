#!/usr/bin/env python3
"""Add one version of the documentation to the versioned documentation site.

The site has one directory per documented version: "develop" and one per minor
release ("v1.0", "v1.1", ...). This script replaces the directory of one
version with a new Doxygen build and regenerates the files shared by all
versions:

  version-switcher.js  version dropdown and banner, loaded by every page
  index.html           redirects to the newest release, or to develop
  404.html             sends missing pages to the front page of their version,
                       and links from before the site was versioned to the
                       same page in the newest release
"""

import argparse
import html
import json
import os
import pathlib
import re
import shutil
import sys

RELEASE = re.compile(r"v(\d+)\.(\d+)")
SWITCHER = "version-switcher.js"
SITE_DATA = re.compile(r"^(\s*const SITE = ).*(; // @site-data@)$", re.MULTILINE)

REDIRECT_PAGE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<meta http-equiv="refresh" content="0; url={target}">
<link rel="canonical" href="{target}">
</head>
<body>
<p>Redirecting to <a href="{target}">{label}</a>.</p>
</body>
</html>
"""

NOT_FOUND_PAGE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Page not found - {title}</title>
<script>
// Served for every missing path, so all links here are absolute.
(function () {{
  var base = {base};
  var versions = {versions};
  var fallback = {fallback};
  var path = window.location.pathname;
  if (path.indexOf(base) !== 0) return;
  var rest = path.slice(base.length);
  var version = rest.split("/")[0];
  var target = null;
  if (versions.indexOf(version) >= 0) {{
    // A page that does not exist in this version: go to its front page.
    if (rest !== version + "/" && rest !== version + "/index.html") {{
      target = base + version + "/";
    }}
  }} else if (fallback) {{
    // A link from before the site was versioned: try the newest release.
    target = base + fallback + "/" + rest;
  }}
  if (target) window.location.replace(target + window.location.hash);
}})();
</script>
</head>
<body>
<h1>Page not found</h1>
<p>Go to the <a href="{base_html}">documentation front page</a>.</p>
</body>
</html>
"""


def site_versions(site):
    """Return the versions in the site, newest release first, develop last."""
    names = [path.name for path in site.iterdir() if path.is_dir()]
    releases = [name for name in names if RELEASE.fullmatch(name)]
    releases.sort(
        key=lambda name: tuple(map(int, RELEASE.fullmatch(name).groups())),
        reverse=True)
    return releases + (["develop"] if "develop" in names else [])


def load_switcher(version_dir, site):
    """Make every page of one version load the version switcher."""
    for page in version_dir.rglob("*.html"):
        content = page.read_bytes()
        if SWITCHER.encode() in content or b"</head>" not in content:
            continue
        source = os.path.relpath(site / SWITCHER, page.parent)
        tag = f'<script src="{source.replace(os.sep, "/")}" defer></script>\n'
        page.write_bytes(content.replace(b"</head>", tag.encode() + b"</head>", 1))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--site", required=True, type=pathlib.Path,
        help="directory of the site, i.e. a checkout of the gh-pages branch")
    parser.add_argument(
        "--version", required=True,
        help='version to add or replace: "develop" or "vMAJOR.MINOR"')
    parser.add_argument(
        "--html", required=True, type=pathlib.Path,
        help="Doxygen HTML output to publish as that version")
    parser.add_argument(
        "--base-path", default="",
        help='path of the site on its host, e.g. "/neko-top"')
    parser.add_argument(
        "--title", default="Neko-TOP documentation",
        help="title of the shared pages (default: %(default)s)")
    args = parser.parse_args()

    if args.version != "develop" and not RELEASE.fullmatch(args.version):
        sys.exit(f"Invalid version '{args.version}', expected develop or vX.Y.")
    if not (args.html / "index.html").is_file():
        sys.exit(f"No Doxygen output found in '{args.html}'.")

    # Replace the directory of this version.
    args.site.mkdir(parents=True, exist_ok=True)
    target = args.site / args.version
    if target.exists():
        shutil.rmtree(target)
    shutil.copytree(args.html, target)

    versions = site_versions(args.site)
    latest = next((v for v in versions if v != "develop"), None)
    front = latest or "develop"
    for version in versions:
        load_switcher(args.site / version, args.site)

    # The switcher, with the current list of versions.
    template = (pathlib.Path(__file__).parent / SWITCHER).read_text()
    data = json.dumps({"latest": latest, "versions": versions})
    script, count = SITE_DATA.subn(
        lambda match: match.group(1) + data + match.group(2), template)
    if count != 1:
        sys.exit(f"Could not find the site data line in {SWITCHER}.")
    (args.site / SWITCHER).write_text(script)

    # The shared pages.
    title = html.escape(args.title)
    (args.site / "index.html").write_text(REDIRECT_PAGE.format(
        title=title, target=f"{front}/",
        label=html.escape(f"the {front} documentation")))
    base = args.base_path.rstrip("/") + "/"
    (args.site / "404.html").write_text(NOT_FOUND_PAGE.format(
        title=title, base=json.dumps(base), versions=json.dumps(versions),
        fallback=json.dumps(front), base_html=html.escape(base)))

    print(f"Published {args.version}. Versions: {', '.join(versions)}; "
          f"front page: {front}.")


if __name__ == "__main__":
    main()
