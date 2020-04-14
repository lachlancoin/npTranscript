
#this makes the jar from the binary files but does not compile the binary files
#should only be used if remaking npTranscript.jar from compiled binary

rm npTranscript.jar
cd java/bin
jar -cf npTranscript.jar npTranscript
mv npTranscript.jar ..
cd ../..
git add npTranscript.jar
