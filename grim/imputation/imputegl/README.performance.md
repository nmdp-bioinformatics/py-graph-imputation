
# Performance enhancements

Timings on 100 donors: data/testwmda.100.txt
	testwmda.100.txt D001000 - D001100

```
$ time python runfile.py

real	1m49.707s
user	0m11.509s
sys	0m0.632s
```


## changing queries to use node id
```
$ time python runfile.py

real	0m20.121s
user	0m9.628s
sys	0m0.422s
```

## create indexes

```
$ python create_indexes.py
```

## confirm that the indexes exist
```
 (cypher) 
 $ CALL db.indexes
 (returns 31 rows)
```

```

$ time python runfile.py

real	0m11.900s
user	0m9.515s
sys	0m0.420s
```

So having the queries use node id resulted in a 5.45x speedup and using indexes added another 1.68x speedup so overall a 9.15x speedup.
