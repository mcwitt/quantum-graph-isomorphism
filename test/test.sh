GSL_RNG_SEED=137 ./cg -s 0 -S 0.95 -n 20 -m 11 graphs/VMB-2013/g13a.txt > test.out
diff -q test.{out,cmp} && rm test.out && echo "success!"
