
PROB=$1

glpsol  --lp "$PROB".lp
echo "______________________________________________"
time ../../cpplex/cpplex-read-only/bin/solver "$PROB".cpplp
echo "______________________________________________"
time ../standard_ref/bin/main "$PROB".dlp
