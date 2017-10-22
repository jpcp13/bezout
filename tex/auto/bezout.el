(TeX-add-style-hook
 "bezout"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("hyperref" "linktoc=all")))
   (TeX-run-style-hooks
    "latex2e"
    "abstract"
    "univariate"
    "multivariate"
    "barnett_formula"
    "roots"
    "bibliography"
    "article"
    "art10"
    "standalone"
    "amsmath"
    "amsfonts"
    "amsthm"
    "hyperref"
    "algorithm2e"
    "graphicx"
    "multirow")
   (LaTeX-add-environments
    "prop"
    "defn"
    "exmp"
    "rem")))

