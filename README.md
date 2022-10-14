# cqf
A General-Purpose Counting Filter: Counting Quotient Filter (CQF)

This work appeared at SIGMOD 2017. If you use this software please cite us:
```
@inproceedings{DBLP:conf/sigmod/PandeyBJP17,
  author    = {Prashant Pandey and
               Michael A. Bender and
               Rob Johnson and
               Robert Patro},
  title     = {A General-Purpose Counting Filter: Making Every Bit Count},
  booktitle = {Proceedings of the 2017 {ACM} International Conference on Management
               of Data, {SIGMOD} Conference 2017, Chicago, IL, USA, May 14-19, 2017},
  pages     = {775--787},
  year      = {2017},
  crossref  = {DBLP:conf/sigmod/2017},
  url       = {http://doi.acm.org/10.1145/3035918.3035963},
  doi       = {10.1145/3035918.3035963},
  timestamp = {Wed, 10 May 2017 22:12:12 +0200},
  biburl    = {http://dblp.org/rec/bib/conf/sigmod/PandeyBJP17},
  bibsource = {dblp computer science bibliography, http://dblp.org}
}
```

Overview
--------
 The CQF supports approximate membership testing and counting the occurrences of
 items in a data set. This general-purpose AMQ is small and fast, has good
 locality of reference, scales out of RAM to SSD, and supports deletions,
 counting (even on skewed data sets), resizing, merging, and highly concurrent
 access.

API
--------
* 'qf_insert(item, count)': insert an item to the filter
* 'qf_count_key_value(item)': return the count of the item. Note that this
  method may return false positive results like Bloom filters or an over count.
* 'qf_remove(item, count)': decrement the count of the item by count. If count
  is 0 then completely remove the item.

Build
-------
This library depends on libssl. 

The code uses two new instructions to implement select on machine words introduced 
in intel's Haswell line of CPUs. However, there is also an alternate implementation
of select on machine words to work on CPUs older than Haswell.

To build on a Haswell or newer hardware:
```bash
 $ make test
 $ ./test 8 7 100000000 1000 1000000 20
```

To build on an older hardare (older than Haswell):
```bash
 $ make NH=1 test
 $ ./test 8 7 100000000 1000 1000000 20
 ```
 
 The arguments for the test program are [log of filter size] [number of remainder bits] [universe size] [number of inserts] [number of queries] [number of trials]
 If the number of inserts exceeds the filter size, the filter will perform inserts until full
 For example, ./test 16 7 100000000 20000000 1000000 20 will make a filter of size 2^16 = 65536 using 7 bits per slot,
              then insert 20000000 random items from a universe of size 100000000 (or until full),
              then make 1000000 queries to test the false positive rate
              This will be done 20 times and the info from the trials will be averaged and reported
 For testing purposes, an optional additional numerical argument can be provided to function as a set seed

Contributing
------------
Contributions via GitHub pull requests are welcome.


Authors
-------
- Prashant Pandey <ppandey@cs.stonybrook.edu>
- Rob Johnson <rob@cs.stonybrook.edu>
