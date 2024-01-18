# AQF
Adaptive Quotient Filter (AQF)
```
```

Overview
--------
 The AQF supports approximate membership testing and counting the occurrences of
 items in a data set. Like other AMQs, the AQF has a chance for false positives
 during queries. However, the AQF has the ability to adapt to false positives
 after they have occurred so they are not repeated. At the same time, the AQF
 maintains the benefits of a quotient filter, as it is small and fast, has good
 locality of reference, scales out of RAM to SSD, and supports deletions,
 counting, resizing, merging, and highly concurrent access.

API
--------
* 'qf_malloc(nslots, pbits)': initialize an AQF
* 'qf_insert_ret(item, count)': insert an item to the filter. Returns 0 if the
  item to insert shares a fingerprint with an existing item in the filter.
* 'insert_and_extend(item, count, other)': insert an item and extend it. Should be called
  after qf_insert_ret returns 0, where other is the key of the previously inserted
  item which caused qf_insert_ret to return 0. If item and other are equal, this
  updates the count for the item in the filter.
* 'qf_query(item)': return the count of the item. Note that this
  method may return false positive results like Bloom filters or an over count.
* 'qf_remove(item)': remove the item from the filter
* 'qf_adapt(item, other)': extend the item in the filter. Afterwards, other is
  guaranteed to not be a false positive query unless the two items are identical.
Functions may need additional parameters for the AQF, flags, and/or return values.
See src/test_micro_throughput.c for an example of an in-memory reverse map.
See src/test_splinter_throughput.c for an example using SplinterDB as reverse map.

Build
-------
This library depends on libssl. 

The code uses two new instructions to implement select on machine words introduced 
in intel's Haswell line of CPUs. However, there is also an alternate implementation
of select on machine words to work on CPUs older than Haswell.

To build on a Haswell or newer hardware:
```bash
 $ make test_micro
 $ ./test_micro 26 9 200000000
```

To build on an older hardware (older than Haswell):
```bash
 $ make NH=1 test_micro
 $ ./test_micro 26 9 200000000
 ```

Contributing
------------
Contributions via GitHub pull requests are welcome.

