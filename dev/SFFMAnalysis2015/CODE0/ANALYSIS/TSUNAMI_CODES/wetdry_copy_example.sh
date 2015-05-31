#for i in $(ls EIGHT_FFI/*/*/*/MODEL_OUTPUTS/*/wet_dry*); do mkdir -p WDOUT/$(dirname $i); cp $i WDOUT/$i; done
for i in $(ls BEST_MODEL_REMAINING_FFI/*/*/*/MODEL_OUTPUTS/*/wet_dry*); do mkdir -p WDOUT/$(dirname $i); cp $i WDOUT/$i; done
