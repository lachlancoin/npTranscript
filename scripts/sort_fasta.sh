grep '>' $1  | cut -f 2  -d ' ' | cut -f 5 -d ';' | sort -k 2,2gr -t '=' | head -n $2 > tmpfile.txt
while read line; do
	#echo $line
	a=$(cat -n  $1 | grep $line | cut -f 1 | tr -d ' ')
	#echo $a
	b=$(tail -n +$a $1 | cat -n | grep '>' | head  -n 2 | tail -n +2 | cut -f 1 | tr -d ' ')
	#echo $a $b
	c=$(($b-1))

	tail -n +$a  $1 | head -n $c
done < tmpfile.txt
rm tmpfile.txt

