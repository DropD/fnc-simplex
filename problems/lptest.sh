
PROB=$1

time ../../cpplex/cpplex-read-only/bin/solver "$PROB".cpplp
echo "______________________________________________"
time ../../gurobi/run_gurobi.sh "$PROB".lp
echo "______________________________________________"
time ../../soplex/bin/soplex "$PROB".lp
echo "______________________________________________"
glpsol  --lp "$PROB".lp
echo "______________________________________________"
time ../bin/main "$PROB".dlp
