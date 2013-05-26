
PROB=$1
time ../../cpplex/cpplex-read-only/bin/solver "$PROB".cpplp
echo "______________________________________________"
glpsol  --lp "$PROB".lp
echo "______________________________________________"
time ../bin/main "$PROB".dlp
