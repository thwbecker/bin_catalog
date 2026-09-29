#!/bin/bash
#
# checks for binning and sampling conventions of bin_catalog, nsample_catalog
# and psample_catalog. run from anywhere:
#
#   test/run_tests.sh [bin_dir]
#
# prints PASS/FAIL per check and returns the number of failures
#
bdir=${1-$(dirname "$(readlink -f "$0")")/../bin}

tdir=$(mktemp -d)
trap "rm -rf $tdir" EXIT
cd $tdir || exit 1
nfail=0

check(){			# name condition(0/1)
    if [ "$2" = 1 ];then
	echo "PASS: $1"
    else
	echo "FAIL: $1"
	((nfail++))
    fi
}
# a cluster of events with a range of strike-slip mechanisms around lon lat
cluster(){			# lon lat n
    gawk -v x=$1 -v y=$2 -v n=$3 'BEGIN{srand(1);
	for(i=0;i<n;i++){
	    lon=x+(rand()-0.5)*0.004;lat=y+(rand()-0.5)*0.004;
	    printf("%.5f %.5f 8.0 %.1f %.1f %.1f 3.0 %.5f %.5f %i\n",
		   lon,lat,rand()*360,70+rand()*20,(rand()<0.5)?(rand()*20-10):(180-rand()*20),lon,lat,1e9+i*3600)}}'
}
common="-m 0 -M 10 -z 30 -Z -10"
#
# 1: bin_catalog reports the node that events are assigned to
#
cluster -117.5 34.0 30 > c1.aki
$bdir/bin_catalog --dx 0.5 $common -l -118 -r -116 -b 33 -t 35 -o b1 c1.aki 2> /dev/null
check "bin_catalog cluster at node" \
      $(gawk '{if($7==-117.5 && $8==34 && $9==30)ok=1}END{print(ok+0)}' b1.0.5.0.5.0.norm.dat)
#
# 2: nsample_catalog reports the location it searched around
#
$bdir/nsample_catalog --dx 0.5 $common -l -118 -r -116 -b 33 -t 35 -p 10 -D 5 -o n1 c1.aki 2> /dev/null
check "nsample_catalog cluster at search point" \
      $(gawk '{n++;if($7==-117.5 && $8==34 && $9==30)ok=1}END{print((ok && n==1)?1:0)}' n1.0.5.0.5.0.norm.dat)
check "nsample_catalog stress at search point" \
      $(gawk '{n++;if($7==-117.5 && $8==34)ok=1}END{print((ok && n==1)?1:0)}' n1.0.5.0.5.s.dat)
#
# 3: psample_catalog at the cluster, with and without distance weighting, finite scaled output
#
for w in 0 1;do
    echo -117.5 34.0 | timeout 120 $bdir/psample_catalog $common -p 10 -D 5 -w $w -o p$w c1.aki 2> /dev/null
    check "psample_catalog weighting $w" \
	  $(gawk '{if($7==-117.5 && $8==34 && $9==30)ok=1}END{print(ok+0)}' p$w.0.norm.dat)
    check "psample_catalog weighting $w finite scaled and stress" \
	  $(cat p$w.0.scaled.dat p$w.s.dat | gawk '{for(i=1;i<=6;i++)if(tolower($i)~/nan/)bad=1}END{print((NR==2 && !bad)?1:0)}')
done
#
# 4: an event 0.8 dy south of the southern boundary is not used
#
echo "-117.0 32.60 8 0 85 0 3 -117.0 32.60 1000000000" > e1.aki
$bdir/bin_catalog --dx 0.5 $common -l -118 -r -116 -b 33 -t 35 -o e1 e1.aki 2> e1.log
check "bin_catalog southern edge" $(grep -c "used 0 out of 1" e1.log)
#
# 5: number of bins for a region that is not exactly representable
#
$bdir/bin_catalog --dx 0.1 $common -l -117.5 -r -117.2 -b 34.0 -t 34.3 -o g1 c1.aki 2> g1.log
check "bin_catalog nx=3 for 0.3/0.1" $(grep -c "nx: 3 ny 3" g1.log)
#
# 6: k nearest with fewer events than k gives no bins
#
cluster -117.5 34.0 20 > c2.aki
$bdir/nsample_catalog --dx 0.5 $common -l -118 -r -116 -b 33 -t 35 -p -50 -D 50 -o k1 c2.aki 2> /dev/null
check "nsample_catalog fewer events than k" $( [ ! -s k1.0.5.0.5.0.norm.dat ] && echo 1 || echo 0)
#
# 7: undefined friction mode is rejected
#
$bdir/bin_catalog --dx 0.5 $common -l -118 -r -116 -b 33 -t 35 -F 7 -o f1 c1.aki > /dev/null 2>&1
check "bin_catalog rejects -F 7" $( [ $? -ne 0 ] && echo 1 || echo 0)

echo "$nfail failure(s)"
exit $nfail
