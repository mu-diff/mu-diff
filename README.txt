This is mu-diff: an open matlab toolbox for solving multiple scattering problems

mu-diff is copyright (C) 2014 X. Antoine and B. Thierry, University of Lorraine, France,
and is distributed under the terms of the GNU General Public License, Version 2
or later. See Doc/LICENSE.txt and Doc/CREDITS.txt for more information.

See the Doc/ and Examples/ directories for documentation and examples. The reference manual is
located in Doc/. See the web site http://mu-diff.math.cnrs.fr for more informations.

The mu-diff toolbox requieres the software Matlab available to purchase on the website
http://www.mathworks.fr/products/matlab/
To install mu-diff, unzip the mu-diff toolbox (available on the website), go to Matlab and add the file directory (and its subdirectories) to the Matlab path file. Then, you can use mu-diff and launch some of the examples.

To verify that mu-diff is well installed, try to launch the Benchmarks:
Examples/Benchmark/BenchmarkDirichlet
Examples/Benchmark/BenchmarkNeumann
Examples/Benchmark/BenchmarkPenetrable

This should result in a GMRES history of convergence, the radar cross section (far field) and the scattered and total fields, computed on a grid.
