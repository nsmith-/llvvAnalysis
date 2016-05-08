# source me
inputdir=$1
shift

mkdir -p /nfs_scratch/nsmith/llvvAnalysis-$inputdir/output

for name in $inputdir/*; do
  sample=${name##*/}
  farmoutAnalysisJobs \
    --infer-cmssw-path \
    --fwklite \
    --submit-dir=/nfs_scratch/nsmith/llvvAnalysis-$inputdir/$sample \
    --input-dir=root://cmsxrootd.hep.wisc.edu//store/user/nsmith/$inputdir/$sample \
    --output-dir=srm://cmssrm2.hep.wisc.edu:8443/srm/v2/server?SFN=/hdfs/store/user/nsmith/llvvAnalysis-$inputdir \
    --job-generates-output-name \
    --extra-inputs=runAnalysis_template.py \
    $@ \
    llvvAnalysis \
    $CMSSW_BASE/bin/$SCRAM_ARCH/run2015_WIMPAnalysis \
    runAnalysis_template.py inputFiles='$inputFileNames' jobNumber='$jobNumber'
done
