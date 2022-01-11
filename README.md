NR-grep: a fast and flexible pattern-matching tool
==================================================

Author: Gonzalo Navarro
Department of Computer Science, University of Chile, Blanco Encalada 2120, Santiago, Chile

From [source distributed here][1]

Summary
-------

We present nrgrep (‘non-deterministic reverse grep’), a new pattern-matching tool designed for efficient
search of complex patterns. Unlike previous tools of the grep family, such as agrep and Gnu grep, nrgrep is
based on a single and uniform concept: the bit-parallel simulation of a non-deterministic suffix automaton.
As a result, nrgrep can find from simple patterns to regular expressions, exactly or allowing errors in the
matches, with an efficiency that degrades smoothly as the complexity of the searched pattern increases.
Another concept that is fully integrated into nrgrep and that contributes to this smoothness is the selection
of adequate subpatterns for fast scanning, which is also absent in many current tools. We show that the
efficiency of nrgrep is similar to that of the fastest existing string-matching tools for the simplest patterns,
and is by far unmatched for more complex patterns. Copyright © 2001 John Wiley & Sons, Ltd.

Source:

Navarro, G. (2001). NR‐grep: a fast and flexible pattern‐matching tool. Software: Practice and Experience, 31(13), 1265–1312.
doi:10.1002/spe.411

[Read full paper here][2]

[1]: https://users.dcc.uchile.cl/~gnavarro/pubcode/
[2]: https://sci-hub.se/10.1002/spe.411
