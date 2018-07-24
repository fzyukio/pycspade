#ifndef __EXTTPOSE_H
#define __EXTTPOSE_H

#include <fstream>

//using namespace std;

/**
 *   -s 0 -l -x
 * @param dbname
 * @param num_partitions
 * @param min_support 0
 * @param twoseq
 * @param use_diff
 * @param do_l2 = false
 * @param do_invert
 * @param use_newformat
 * @param maxmem
 * @param no_minus_off = true
 * @return
 */
extern std::vector<std::string> _exttpose(const std::string& dbname, int num_partitions = 1, double min_support = 0., bool twoseq = false,
                     bool use_diff = false, bool do_l2 = true, bool do_invert = true, bool use_newformat = true,
                     int maxmem = 32, bool no_minus_off = false);

#endif //__EXTTPOSE_H