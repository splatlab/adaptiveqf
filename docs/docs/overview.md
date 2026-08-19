# Overview


Filters are probabilistic data structures that are used to check if an element
is present in a set. 

Filters answer membership queries with one-sided errors: they never return false
negatives, but may return false positives. If the filter returns NO (element not
in set), the filter is always correct. However, if the filter returns YES
(element in set), there is a small chance of it being wrong.

> **TODO**: Describe why filters are useful.

> **TODO:** Describe the problem of repeated false positives. 

> **TODO:** Describe adaptivity — see [The Quotient Filter](qf.md).

> **TODO:** Further reading: point to papers, add citations.

> **TODO:** If you just want to get started, go to getting started.

