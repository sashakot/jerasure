#!/bin/bash

core_num=$(nproc --all)

#
# Configuration
#
CEPH_BIN=${CEPH_BIN:-"/mnt/data/sashakot/ceph-install/bin"}
CEPH_ERASURE_BENCHMARK=${CEPH_BIN}/ceph_erasure_code_benchmark
PIDSTAT=${PIDSTAT:-"/hpc/local/work/sashakot/sysstat/pidstat"}

#
# Runtime parameters
#

iter=${iter:-1000000}
k=${k:-2}
m=${m:-1}
min_kb=${min_kb:-102}
max_kb=${max_kb:-500}
step_kb=${step_kb:-2}
core=${core:-$[$core_num/2]}
tmpfile=$(mktemp /tmp/erasure-benchmark-script.XXXXXX)

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
pidstat="$PIDSTAT 1 -u -e"

function clean
{
	rm -f ${tmpfile}
}

trap clean EXIT

echo ""
printf "%s,\t%s,\t%s\n" "Size (KB)" "Time (sec)" "CPU (%)"
for kb_bytes in $(seq $min_kb $step_kb $max_kb); do
	rm -f ${tmpfile}

	bytes=$[$kb_bytes*1024]
	params="-w encode  -P k=$k -P m=$m -P technique=reed_sol_van -s $bytes -i $iter"
	${pidstat}  ${taskset} ${CEPH_ERASURE_BENCHMARK} ${params} > ${tmpfile}

	reported_time=$(cat ${tmpfile} | awk ' NF==2 {print $1} ')
	cpu=$(grep Average ${tmpfile} | awk '{print $8}')

	printf '%d,\t\t%.2f,\t\t%.0f\n' $kb_bytes $reported_time $cpu
done
