mkdir empty
for i in results_*; do rsync -avu empty/ $i/ --delete; done
rmdir empty
rmdir results_*
