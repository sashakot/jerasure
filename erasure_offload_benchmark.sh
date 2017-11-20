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

iter=${iter:-10000}
k=${k:-2}
m=${m:-1}
min_kb=${min_kb:-512}
max_kb=${max_kb:-$[ 10*1024 ]}
step_kb=${step_kb:-512}
core=${core:-$[$core_num/2]}
tmpfile=$(mktemp /tmp/erasure-benchmark-script.XXXXXX)
ib_dev=${ib_dev:-"mlx5_2"}

echo "Configuration:"
echo $CEPH_BIN
echo $CEPH_ERASURE_BENCHMARK

echo ""
echo "Runtime parameters"
echo "Iterations: $iter"

echo "k: $k"

if [[ (( $k != 2 )) ]]; then
	echo "WARNING: k > 2 is not supported. Mememory registarion problem"
fi

echo "m: $m"

if [[ (( $m != 1 )) ]]; then
	echo "WARNING: m > 1 is not supported. Mememory registarion problem"
fi

echo "min data size in: $min_kb (KB)"
echo "max data size in: $max_kb (KB)"
echo "data size step: $step_kb (KB)"
echo "ib device: $ib_dev"
echo "core $core"

taskset="taskset -c $core "
pidstat="$PIDSTAT 1 -u -e"

function clean
{
	rm -f ${tmpfile}
}

trap clean EXIT

function run_test()
{
	local bytes=$1
	local iter=$2
	local offload=$3

	rm -f ${tmpfile}

	params="-w encode  -P k=$k -P m=$m -P technique=reed_sol_van -s $bytes -i $iter"
	MLNX_OFFLOAD=$offload MLNX_K=$k MLNX_M=$m MLNX_IB_DEV=$ib_dev ${pidstat} ${taskset} ${CEPH_ERASURE_BENCHMARK} ${params} > ${tmpfile} 2>/dev/null

	reported_time=$(cat ${tmpfile} | awk ' NF==2 {print $1} ')
	cpu=$(grep Average ${tmpfile} | awk '{print $8}')

#	cat ${tmpfile}

	printf '%.2f,\t\t\t%.0f\n' $reported_time $cpu
}

echo ""
printf "%s,\t%s,\t%s,\t%s,\t%s\n" "Size (KB)" "Baseline time (sec)" "CPU (%)" "Offload time (sec)" "CPU (%)"
for kb_bytes in $(seq $min_kb $step_kb $max_kb); do
	bytes=$[$kb_bytes*1024]

	printf '%d,\t\t%s,\t\t%s\n' $kb_bytes "$(run_test $bytes $iter OFF)" "$(run_test $bytes $iter ON)"
done
