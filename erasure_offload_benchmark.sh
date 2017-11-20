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

iter=${iter:-1}
k=${k:-2}
m=${m:-1}
min_kb=${min_kb:-102}
max_kb=${max_kb:-102}
step_kb=${step_kb:-2}
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
echo "ib device: $ib_dev"
echo "core $core"

taskset="taskset -c $core "
pidstat="$PIDSTAT 1 -u -e"

function clean
{
	rm -f ${tmpfile}
}

trap clean EXIT

export MLNX_K=$k
export MLNX_M=$m
export MLNX_IB_DEV=$ib_dev

function run_test()
{
	local bytes=$1
	local iter=$2
	local offload=$3

	params="-w encode  -P k=$k -P m=$m -P technique=reed_sol_van -s $bytes -i $iter"
	MLNX_OFFLOAD=$offload ${pidstat}  ${taskset} ${CEPH_ERASURE_BENCHMARK} ${params} > ${tmpfile}

	reported_time=$(cat ${tmpfile} | awk ' NF==2 {print $1} ')
	cpu=$(grep Average ${tmpfile} | awk '{print $8}')

	printf '%.2f,\t\t%.0f\n' $reported_time $cpu
}

echo ""
printf "%s,\t%s,\t%s\n" "Size (KB)" "Time (sec)" "CPU (%)"
for kb_bytes in $(seq $min_kb $step_kb $max_kb); do
	rm -f ${tmpfile}

	bytes=$[$kb_bytes*1024]

	printf '%d,\t%s,\t%s\n' $kb_bytes $(run_test $bytes $iter OFF) $(run_test $bytes $iter ON)
done
