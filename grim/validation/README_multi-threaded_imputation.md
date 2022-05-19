# Multi-threaded imputation
  
```runfile_mt.py``` is a multi-threaded version of ```runfile.py```. 
```
cd validation/
python  runfile_mt.py -d num_threads -c  CONFIG_FILE
```
where ```d``` is number of threads (batches) and by default is 4

output of imputation for individual batches can be merged to get one file for imputation results. for example
```
cd validation/output  
cat don.a?_outdeltaFunc_muug > don_outdeltaFunc_muug
```
