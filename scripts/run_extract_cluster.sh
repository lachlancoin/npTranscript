JSA_MEM1=8000m

npTranscript=$HOME/github/npTranscript
classp=$(ls ${npTranscript}/libs | xargs -I {} echo ${npTranscript}/libs/{} )
classpath=$(echo $classp | sed 's/ /:/g')
JSA_CP=${npTranscript}/npTranscript.jar:${classpath}

str="java -Xmx${JSA_MEM1} -ea -Djava.awt.headless=true -Dfile.encoding=UTF-8 -classpath ${JSA_CP} npTranscript.run.ExtractClusterCmd $@"

echo "running.."
echo $str
$str

echo "finished"
