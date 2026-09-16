(() => {
  const storageKey = "bcmp-r-manual-theme";
  const root = document.documentElement;
  const rCodeMarker = /(<-|::|library\s*\(|install\.packages\s*\(|bcmp(?:_embedding)?\s*\(|\b(?:NULL|TRUE|FALSE)\b)/;
  const rTokenPattern = /(#.*$)|(\"(?:\\.|[^\"])*\"|'(?:\\.|[^'])*')|(\b(?:TRUE|FALSE|NULL|NA|NaN|Inf)\b)|(\b[A-Za-z.][A-Za-z0-9._]*)(?=\s*\()|(<-|::|\$|=)|(\b\d+(?:\.\d+)?L?\b)/gm;

  const escapeHtml = (value) => value.replace(/[&<>\"]/g, (character) => ({
    "&": "&amp;",
    "<": "&lt;",
    ">": "&gt;",
    '"': "&quot;",
  })[character]);

  const highlightRLine = (source) => {
    rTokenPattern.lastIndex = 0;
    let output = "";
    let cursor = 0;
    let match;

    while ((match = rTokenPattern.exec(source)) !== null) {
      output += escapeHtml(source.slice(cursor, match.index));
      const tokenClass = match[1] ? "tok-comment"
        : match[2] ? "tok-string"
        : match[3] ? "tok-literal"
        : match[4] ? "tok-function"
        : match[5] ? "tok-operator"
        : "tok-number";
      output += '<span class="' + tokenClass + '">' + escapeHtml(match[0]) + "</span>";
      cursor = match.index + match[0].length;
    }

    return output + escapeHtml(source.slice(cursor));
  };

  const highlightR = (source) => source.split("\n")
    .map(highlightRLine)
    .join("\n");

  const applyRHighlighting = () => {
    document.querySelectorAll("pre code").forEach((code) => {
      const source = code.textContent;
      if (rCodeMarker.test(source)) {
        code.innerHTML = highlightR(source);
      }
    });
  };

  const applyThemeImages = (theme) => {
    document.querySelectorAll("img[data-theme-src-light]").forEach((image) => {
      image.src = image.dataset[theme === "dark" ? "themeSrcDark" : "themeSrcLight"];
    });
  };

  const applyTheme = (theme, button) => {
    root.dataset.theme = theme;
    applyThemeImages(theme);
    const dark = theme === "dark";
    button.setAttribute(
      "aria-label",
      dark ? "Switch to light mode" : "Switch to dark mode",
    );
    button.setAttribute(
      "title",
      dark ? "Switch to light mode" : "Switch to dark mode",
    );
    button.querySelector(".theme-toggle__icon").textContent = dark ? "☼" : "☾";
    button.querySelector(".theme-toggle__label").textContent = dark ? "Light" : "Dark";
  };

  document.addEventListener("DOMContentLoaded", () => {
    const button = document.querySelector("[data-theme-toggle]");
    if (!button) return;

    const saved = localStorage.getItem(storageKey);
    const initial = saved || (
      window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light"
    );
    applyTheme(initial, button);
    applyRHighlighting();

    button.addEventListener("click", () => {
      const next = root.dataset.theme === "dark" ? "light" : "dark";
      localStorage.setItem(storageKey, next);
      applyTheme(next, button);
    });
  });
})();
