/* Version dropdown and banner for the versioned documentation site.
 *
 * The site has one directory per documented version ("develop", "v1.0", ...)
 * and this script at its root. Every page loads it to get a dropdown of the
 * versions next to the project name, and a banner when the page is not part of
 * the newest release. The site data below is filled in when the site is
 * published, see update_site.py.
 */
(function () {
  "use strict";

  const SITE = { latest: null, versions: [] }; // @site-data@

  // The first path component below the site root, where this script lives,
  // is the version of the current page.
  const script = document.currentScript;
  if (!script) return;
  const root = new URL(".", script.src).href;
  const here = window.location.href;
  if (!here.startsWith(root)) return;
  const relative = here.slice(root.length);
  const current = relative.split("/")[0];
  if (!SITE.versions.includes(current)) return;

  // The same page in another version. Pages missing from that version are
  // sent to its front page by the site's 404 page.
  const page = relative.slice(current.length + 1);
  function pageIn(version) {
    return root + version + "/" + page;
  }

  // Colours follow the doxygen-awesome theme, including its dark mode.
  const style = document.createElement("style");
  style.textContent = `
    #version-switcher {
      margin-left: 0.75em;
      padding: 0.1em 0.4em;
      vertical-align: middle;
      font-size: 0.6em;
      color: var(--page-foreground-color, inherit);
      background: var(--page-background-color, #fff);
      border: 1px solid var(--separator-color, #ccc);
      border-radius: var(--border-radius-small, 4px);
    }
    #version-banner {
      /* The theme lays out the header as a wrapping flex row. */
      flex: 0 0 100%;
      box-sizing: border-box;
      padding: 0.5em 1em;
      text-align: center;
      color: var(--warning-color-darker, #5f4204);
      background: var(--warning-color, #faf3d8);
      border-bottom: 1px solid var(--warning-color-dark, #f3a600);
    }
    #version-banner a {
      color: inherit;
      font-weight: bold;
      text-decoration: underline;
    }`;
  document.head.appendChild(style);

  // Dropdown of all versions, next to the project name in the title area.
  const select = document.createElement("select");
  select.id = "version-switcher";
  select.title = "Documentation version";
  select.setAttribute("aria-label", "Documentation version");
  for (const version of SITE.versions) {
    const option = document.createElement("option");
    option.value = version;
    option.textContent =
      version === SITE.latest ? version + " (latest)" : version;
    option.selected = version === current;
    select.appendChild(option);
  }
  select.addEventListener("change", () => {
    window.location.href = pageIn(select.value);
  });
  const title = document.getElementById("projectname");
  if (title) title.appendChild(select);

  // Banner on every version other than the newest release.
  if (current !== SITE.latest) {
    const banner = document.createElement("div");
    banner.id = "version-banner";
    banner.setAttribute("role", "note");
    banner.append(
      current === "develop"
        ? "This is the documentation of the development version, which may " +
            "describe features that are not released yet. "
        : "This is the documentation of " + current + ", an older release. "
    );
    if (SITE.latest) {
      const link = document.createElement("a");
      link.href = pageIn(SITE.latest);
      link.textContent = "Go to the latest release (" + SITE.latest + ")";
      banner.append(link, ".");
    }
    // Inside the header, so Doxygen sizes the page panels around it.
    const header = document.getElementById("top") || document.body;
    header.prepend(banner);
  }
})();
