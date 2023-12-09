sections = {
  "theory": 1,
}


window.MathJax = {
  loader: {load: ['[tex]/tagformat', '[tex]/ams']},
  tex: {
    packages: {'[+]': ['tagformat', 'ams']},
    macros: {
      dd: "{\\, \\mathrm{d}}",
      R: "{\\mathbb{R}}",
    },
    tags: 'ams',
    tagformat: {
      number: (n) => sections[window.location.pathname.split("/").pop().split(".")[0]] + '.' + n,
    },
    ams: {
      multilineWidth: '100%',
      multilineIndent: '50em'
    }
  },
}
