lp_gen
======

DEPENDENCIES: 

* python
* cheetah (http://www.cheetahtemplate.org/).

Generate a benchmark suite of feasible (easy) linear programs
in both CPLEX ".lp" and a simplified ".dlp" format (and any format you care to write a template for).

The Generator uses cheetah, a python templating engine to make it easy to output to various formats.

To build a complete benchmark with problems of various sizes you may use the included CMake build system.

to generate an lp with specific size:

* compile the templates with cheetah
* then run the gen_lp.py script (run without arguments for usage hints).

to generate an lp with specific size:

* compile the templates with cheetah
* then run the gen_lp.py script (run without arguments for usage hints).
