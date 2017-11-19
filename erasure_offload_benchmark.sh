#!/bin/bash

core_num=$(nproc --all)

#
# Configuration
#
CEPH_BIN=${CEPH_BIN:-"/mnt/data/sashakot/ceph-install/bin"}
CEPH_ERASURE_BENCHMARK=${CEPH_BIN}/ceph_erasure_code_benchmark

#
# Runtime parameters
#

iter=${iter:-10000}
k=${k:-2}
m=${m:-1}
min_kb=${min_kb:-4}
max_kb=${max_kb:-16}
step_kb=${step_kb:-2}
core=${core:-$[$core_num/2]}

echo "Configuration:"
echo $CEPH_BIN
echo $CEPH_ERASURE_BENCHMARK

echo ""
echo "Runtime parameters"
echo "Iterations: $iter"
echo "k: $k"
echo "m: $m"
echo "min data size in: $min_kb (KB)"
echo "max data size in: $max_kb (KB)"
echo "core $core"

taskset="taskset -c $core "

echo ""
echo "Size (KB), Time (sec), CPU (%)"
for kb_bytes in $(seq $min_kb $step_kb $max_kb); do
	bytes=$[$kb_bytes*1024]
	params="-w encode  -P k=$k -P m=$m -P technique=reed_sol_van -s $bytes -i $iter"
	output=$(${taskset} ${CEPH_ERASURE_BENCHMARK} ${params})

	printf '%d,\t%f\n' $kb_bytes `echo $output | awk '{print $1}'`
done
