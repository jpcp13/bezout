(TeX-add-style-hook
 "univariate"
 (lambda ()
   (TeX-run-style-hooks
    "latex2e"
    "standalone"
    "standalone10")
   (LaTeX-add-labels
    "compan"
    "compan2roots"
    "g_1=g_2"
    "def_bez"
    "xBg"
    "exmp_1"
    "relations_prop"
    "relations"
    "Bar"
    "Barnett"
    "Bar_gen"
    "GBG")))

