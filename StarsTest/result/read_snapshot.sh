#mkdir ./result
#mkdir ./result/dat
for file in `ls out*`
do
  ./read_snapshot.exe ${file}
done
