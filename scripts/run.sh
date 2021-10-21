##general script for running analysis
##ASSUMES THE PROJECT IS STORED IN $HOME/github/npTranscript
if [ ! $JSA_MEM ] ; then
	JSA_MEM=8000m
fi
if [ ! $npTranscript ] ; then
 export npTranscript=$HOME/github/npTranscript
fi

if [ ! $mainclass ] ; then
 mainclass="npTranscript.run.ViralTranscriptAnalysisCmd2"
fi

#npTranscript=$HOME/github/npTranscript
#classp=$(ls ${npTranscript}/libs | xargs -I {} echo ${npTranscript}/libs/{} )
#classpath=$(echo $classp | sed 's/ /:/g')
JSA_CP=${npTranscript}/java/target/npTranscript-1.0.jar:${npTranscript}/java/target/classes/
#echo $JSA_CP



str="java -Xmx${JSA_MEM} -ea -Djava.awt.headless=true -Dfile.encoding=UTF-8 -classpath ${JSA_CP} ${mainclass} $@"
echo "running .."
echo $str
$str
echo "finished"


