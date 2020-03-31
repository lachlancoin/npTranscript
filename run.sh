JSA_MEM=8000m
npTranscript=$HOME/github/npTranscript
classp=$(ls ${npTranscript}/libs | xargs -I {} echo ${npTranscript}/libs/{} )
classpath=$(echo $classp | sed 's/ /:/g')
JSA_CP=${npTranscript}/npTranscript.jar:${classpath}
echo $JSA_CP



java -Xmx${JSA_MEM} -ea -Djava.awt.headless=true -Dfile.encoding=UTF-8 -classpath ${JSA_CP} npTranscript.run.ViralTranscriptAnalysisCmd2 "$@"


