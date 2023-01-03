
initTPRFile=wt_ca_domain.tpr
initXTCFile=all_combined.xtc
initNDXFile=wt_ca_domains.ndx

prefix=both_tmd
TPRFile=backbone.tpr
XTCFile=tetramer_combined.xtc
output=cov-domain-tetramer
selText=both_tmd
## = = Extract alpha-Carbons for quicker processing.
#echo "$selText $selText" | gmx trjconv \
#-s $initTPRFile -f $initXTCFile -n $initNDXFile \
#-pbc mol -ur compact -center \
#-o $prefix 
#
#echo "$selText" | gmx convert-tpr \
#-s $initTPRFile -n $initNDXFile \
#-o $prefix

## = = Generate backbone PCA
echo "System System" | gmx covar \
-s $TPRFile -f $XTCFile \
-o ${output}-eigval.xvg \
-v ${output}-eigvec.trr \
-av ${output}-average.pdb \
-l ${output}.log

# = = 
echo "System System" | gmx anaeig \
-s $TPRFile -f $XTCFile \
-eig ${output}-eigval.xvg \
-v ${output}-eigvec.trr \
-nframes 81 -first 1 -last 9 \
-extr ${output}-extr.pdb \
-proj ${output}-proj1to9.xvg
