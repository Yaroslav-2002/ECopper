/* CuRemediate – progressive enhancements. Every page works without this file. */
(function () {
  "use strict";

  // ---- Sortable tables: <table data-sortable>, <th class="sortable" data-type="num|text">, cells may carry data-value.
  function cellValue(row, idx, type) {
    var cell = row.cells[idx];
    if (!cell) return type === "num" ? -Infinity : "";
    var raw = cell.hasAttribute("data-value") ? cell.getAttribute("data-value") : cell.textContent.trim();
    if (type === "num") {
      var n = parseFloat(raw);
      return isNaN(n) ? -Infinity : n;
    }
    return raw.toLowerCase();
  }

  document.querySelectorAll("table[data-sortable]").forEach(function (table) {
    var tbody = table.tBodies[0];
    if (!tbody) return;
    var headers = table.querySelectorAll("thead th");
    headers.forEach(function (th, idx) {
      if (!th.classList.contains("sortable")) return;
      th.setAttribute("tabindex", "0");
      th.setAttribute("aria-sort", "none");
      function sort() {
        var type = th.getAttribute("data-type") || "text";
        var asc = !th.classList.contains("sort-asc");
        if (!th.classList.contains("sort-asc") && !th.classList.contains("sort-desc") && type === "num") asc = false;
        headers.forEach(function (h) {
          h.classList.remove("sort-asc", "sort-desc");
          if (h.hasAttribute("aria-sort")) h.setAttribute("aria-sort", "none");
        });
        th.classList.add(asc ? "sort-asc" : "sort-desc");
        th.setAttribute("aria-sort", asc ? "ascending" : "descending");
        var rows = Array.prototype.slice.call(tbody.rows);
        rows.sort(function (a, b) {
          var va = cellValue(a, idx, type), vb = cellValue(b, idx, type);
          if (va < vb) return asc ? -1 : 1;
          if (va > vb) return asc ? 1 : -1;
          return 0;
        });
        rows.forEach(function (r) { tbody.appendChild(r); });
      }
      th.addEventListener("click", sort);
      th.addEventListener("keydown", function (e) {
        if (e.key === "Enter" || e.key === " ") { e.preventDefault(); sort(); }
      });
    });
  });

  // ---- Level filter chips + text search: [data-filter-target="#tableId"]
  document.querySelectorAll("[data-filter-for]").forEach(function (bar) {
    var table = document.querySelector(bar.getAttribute("data-filter-for"));
    if (!table) return;
    var chips = bar.querySelectorAll(".chip");
    var search = bar.querySelector(".search");
    var empty = document.querySelector(bar.getAttribute("data-empty") || "#none");
    var level = "all";
    function apply() {
      var q = search ? search.value.trim().toLowerCase() : "";
      var shown = 0;
      Array.prototype.forEach.call(table.tBodies[0].rows, function (row) {
        var okLevel = level === "all" || row.getAttribute("data-level") === level;
        var okText = !q || row.textContent.toLowerCase().indexOf(q) !== -1;
        var show = okLevel && okText;
        row.classList.toggle("is-hidden", !show);
        if (show) shown++;
      });
      if (empty) empty.style.display = shown ? "none" : "block";
    }
    chips.forEach(function (chip) {
      chip.addEventListener("click", function () {
        chips.forEach(function (c) { c.classList.remove("active"); c.setAttribute("aria-pressed", "false"); });
        chip.classList.add("active");
        chip.setAttribute("aria-pressed", "true");
        level = chip.getAttribute("data-level");
        apply();
      });
    });
    if (search) search.addEventListener("input", apply);
  });

  // ---- Clickable rows: <tr data-href="...">
  document.querySelectorAll("tr[data-href]").forEach(function (row) {
    row.style.cursor = "pointer";
    row.addEventListener("click", function (e) {
      if (e.target.closest("a, button, input")) return;
      window.location.href = row.getAttribute("data-href");
    });
  });

  // ---- Upload form
  var form = document.getElementById("upload-form");
  if (!form) return;

  var dropzone = form.querySelector(".dropzone");
  var fileInput = dropzone && dropzone.querySelector("input[type=file]");
  var fileName = dropzone && dropzone.querySelector(".file-chip .name");
  var clearFile = dropzone && dropzone.querySelector(".clear-file");
  var textarea = form.querySelector("textarea");
  var counter = form.querySelector("[data-seq-counter]");
  var typeHint = form.querySelector("[data-type-hint]");
  var typeRadios = form.querySelectorAll("input[name=seq_type]");
  var submit = form.querySelector("button[type=submit]");
  var overlay = document.querySelector(".loading-overlay");

  function selectedType() {
    var r = form.querySelector("input[name=seq_type]:checked");
    return r ? r.value : "protein";
  }

  function showFile() {
    var has = fileInput.files && fileInput.files.length > 0;
    dropzone.classList.toggle("has-file", has);
    if (has) {
      var f = fileInput.files[0];
      var kb = f.size / 1024;
      fileName.textContent = f.name + " (" + (kb < 1024 ? kb.toFixed(1) + " KB" : (kb / 1024).toFixed(1) + " MB") + ")";
      if (textarea.value.trim()) {
        textarea.value = "";
        updateTextInfo();
      }
      fileInput.files[0].text && fileInput.files[0].slice(0, 200000).text().then(function (t) { guessType(t); });
    }
  }

  if (fileInput) {
    fileInput.addEventListener("change", showFile);
    ["dragenter", "dragover"].forEach(function (ev) {
      dropzone.addEventListener(ev, function () { dropzone.classList.add("dragover"); });
    });
    ["dragleave", "drop"].forEach(function (ev) {
      dropzone.addEventListener(ev, function () { dropzone.classList.remove("dragover"); });
    });
    clearFile.addEventListener("click", function (e) {
      e.preventDefault();
      e.stopPropagation();
      fileInput.value = "";
      showFile();
      if (!textarea.value.trim()) typeHint.textContent = "";
    });
  }

  // Count sequences and residues in the pasted text, and warn when it looks like the other sequence type.
  function parseFasta(text) {
    var seqs = [];
    text.split(/\r?\n/).forEach(function (line) {
      line = line.trim();
      if (!line) return;
      if (line.charAt(0) === ">") { seqs.push(""); return; }
      if (!seqs.length) seqs.push("");
      seqs[seqs.length - 1] += line.replace(/[\s*\-]/g, "").toUpperCase();
    });
    return seqs;
  }

  function guessType(text) {
    var seqs = parseFasta(text).filter(Boolean);
    if (!seqs.length) { typeHint.textContent = ""; return; }
    var joined = seqs.join("").slice(0, 50000);
    var nt = (joined.match(/[ACGTUN]/g) || []).length / joined.length;
    var looks = nt > 0.9 ? "nucleotide" : "protein";
    if (looks !== selectedType()) {
      typeHint.innerHTML = "";
      var span = document.createElement("span");
      span.className = "hint-warn";
      span.textContent = "This looks like a " + looks + " sequence. ";
      var btn = document.createElement("a");
      btn.href = "#";
      btn.textContent = "Switch to " + looks;
      btn.addEventListener("click", function (e) {
        e.preventDefault();
        form.querySelector("input[name=seq_type][value=" + looks + "]").checked = true;
        typeHint.textContent = "";
        updateTextInfo();
      });
      typeHint.appendChild(span);
      typeHint.appendChild(btn);
    } else {
      typeHint.textContent = "";
    }
  }

  function updateTextInfo() {
    var text = textarea.value;
    if (!text.trim()) {
      counter.textContent = "";
      if (!(fileInput && fileInput.files.length)) typeHint.textContent = "";
      return;
    }
    var seqs = parseFasta(text).filter(Boolean);
    var residues = seqs.reduce(function (s, x) { return s + x.length; }, 0);
    var unit = selectedType() === "nucleotide" ? "nt" : "aa";
    counter.textContent = seqs.length + " sequence" + (seqs.length === 1 ? "" : "s") + " · " +
      residues.toLocaleString() + " " + unit + (seqs.length > 1 ? " · batch screen" : "");
    guessType(text);
  }

  textarea.addEventListener("input", function () {
    if (textarea.value.trim() && fileInput && fileInput.files.length) {
      fileInput.value = "";
      showFile();
    }
    updateTextInfo();
  });
  typeRadios.forEach(function (r) { r.addEventListener("change", updateTextInfo); });
  updateTextInfo();

  // Example sequences
  var examplesEl = document.getElementById("examples-data");
  if (examplesEl) {
    var examples = JSON.parse(examplesEl.textContent);
    form.querySelectorAll("[data-example]").forEach(function (btn) {
      btn.addEventListener("click", function () {
        var ex = examples[btn.getAttribute("data-example")];
        if (!ex) return;
        form.querySelector("input[name=seq_type][value=" + ex.seq_type + "]").checked = true;
        if (fileInput) { fileInput.value = ""; showFile(); }
        textarea.value = ex.text.trim();
        updateTextInfo();
        textarea.focus();
        textarea.setSelectionRange(0, 0);
        textarea.scrollTop = 0;
      });
    });
  }

  // Loading state: BLAST plus literature lookups can take a while.
  form.addEventListener("submit", function () {
    submit.disabled = true;
    submit.classList.add("is-loading");
    submit.querySelector(".label").textContent = "Screening…";
    if (overlay) setTimeout(function () { overlay.classList.add("show"); }, 250);
  });
  // Restore the button when the page comes back from the back/forward cache.
  window.addEventListener("pageshow", function (e) {
    if (!e.persisted) return;
    submit.disabled = false;
    submit.classList.remove("is-loading");
    submit.querySelector(".label").textContent = "Analyse";
    if (overlay) overlay.classList.remove("show");
  });
})();
