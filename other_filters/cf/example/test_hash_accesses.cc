#include "cuckoofilter.h"

#include <assert.h>
#include <math.h>
#include <time.h>

#include <iostream>
#include <vector>

using cuckoofilter::CuckooFilter;

double rand_zipfian(double s, double max) {
        double p = (double)rand() / RAND_MAX;

        double pD = p * (12 * (pow(max, -s + 1) - 1) / (1 - s) + 6 + 6 * pow(max, -s) + s - s * pow(max, -s + 1));
        double x = max / 2;
        while (true) {
                double m = pow(x, -s - 2);
                double mx = m * x;
                double mxx = mx * x;
                double mxxx = mxx * x;

                double b = 12 * (mxxx - 1) / (1 - s) + 6 + 6 * mxx + s - (s * mx) - pD;
                double c = 12 * mxx - (6 * s * mx) + (m * s * (s + 1));
                double newx = x - b / c > 1 ? x - b / c : 1;
                if (abs(newx - x) <= 0.01) { // this is the tolerance for approximation
                        return newx;
                }
                x = newx;
        }
}

int main(int argc, char **argv) {
  size_t total_items = 1 << 23;

  // Create a cuckoo filter where each item is of type size_t and
  // use 12 bits for each item:
  //    CuckooFilter<size_t, 12> filter(total_items);
  // To enable semi-sorting, define the storage of cuckoo filter to be
  // PackedTable, accepting keys of size_t type and making 13 bits
  // for each key:
  //   CuckooFilter<size_t, 13, cuckoofilter::PackedTable> filter(total_items);
  CuckooFilter<size_t, 12> filter((int)(total_items * 0.9));

  // Insert items to this cuckoo filter
  clock_t start_time = clock(), end_time;
  size_t num_inserted = 0;
  for (size_t i = 0; i < total_items; i++, num_inserted++) {
    if (filter.Add(rand()) != cuckoofilter::Ok) {
      break;
    }
  }
  end_time = clock();
  printf("time per insert: %f\n", (double)(end_time - start_time) / total_items);

  // Check if previously inserted items are in the filter, expected
  // true for all items
  /*for (size_t i = 0; i < num_inserted; i++) {
    assert(filter.Contain(i) == cuckoofilter::Ok);
  }*/

  /*uint64_t *query_set = (uint64_t*)calloc(total_items, sizeof(uint64_t));
  for (size_t i = 0; i < total_items; i++) {
	  //query_set[i] = rand_zipfian(1.5f, 1lu << 30);
	  query_set[i] = i + total_items;
  }*/

  FILE *fp = fopen("hash_accesses.txt", "w");
  fprintf(fp, "queries\taccesses\n");
  fclose(fp);
  fp = fopen("hash_accesses.txt", "a");

  // Check non-existing items, a few false positives expected
  start_time = clock();
  size_t total_queries = 0;
  size_t false_queries = 0;
  /*for (size_t i = total_items; i < 2 * total_items; i++) {
    if (filter.Contain(i) == cuckoofilter::Ok) {
      false_queries++;
    }
    total_queries++;
  }*/
  for (size_t i = 0; i < 10000000; i++) {
    if (filter.Contain(rand_zipfian(1.5f, 1ull << 30)) == cuckoofilter::Ok) {
      false_queries++;
    }
    total_queries++;
    if (i % 100000 == 0) {
	    fprintf(fp, "%lu\t%lu\n", i, false_queries);
    }
  }
  end_time = clock();
  printf("time per query: %f\n", (double)(end_time - start_time) / total_queries);
  fclose(fp);

  // Output the measured false positive rate
  std::cout << "false positive rate is "
            << 100.0 * false_queries / total_queries << "%\n";

  return 0;
}
