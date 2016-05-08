# source me
inputdir=$1
shift

for name in $inputdir/MC13TeV_DYJetsToLL_M50Pt100_amcatnlo; do
  sample=${name##*/}
  farmoutAnalysisJobs \
    --infer-cmssw-path \
    --fwklite \
    --submit-dir=/nfs_scratch/nsmith/llvvAnalysis-$inputdir/$sample \
    --input-dir=root://cmsxrootd.hep.wisc.edu//store/user/nsmith/$inputdir/$sample \
    --output-dir=/nfs_scratch/nsmith/llvvAnalysis-$inputdir/output \
    --job-generates-output-name \
    --extra-inputs=runAnalysis_template.py \
    $@ \
    llvvAnalysis \
    $CMSSW_BASE/bin/$SCRAM_ARCH/run2015_WIMPAnalysis \
    runAnalysis_template.py inputFiles='$inputFileNames' jobNumber='$jobNumber'
done
