
#this makes the jar from the binary files but does not compile the binary files
rm npTranscript.jar
cd bin
jar -cf npTranscript.jar npTranscript
mv npTranscript.jar ..
